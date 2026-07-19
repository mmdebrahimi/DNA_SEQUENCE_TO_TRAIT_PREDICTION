"""CDC & FDA AR Isolate Bank ingester — public isolate HTML -> BioSample-keyed MIC labels.

The AR Isolate Bank (wwwn.cdc.gov/ARIsolateBank) is PUBLIC (no login, no order, no biosafety
attestation for the DATA — only the physical isolates need that). It server-renders, per isolate:
  - the NCBI **BioSample** accession (SAMN...) linking the whole-genome sequence,
  - the **organism** (genus species),
  - reference **broth-microdilution MICs** + CLSI/FDA **S/I/R** interpretation per drug.

That makes it a rare thing for this project: measured **BMD-MIC ⋈ WGS**, per isolate, FREE — the
CRyPTIC-analogue shape the label wall needed (`wiki/cryptic_analogue_consortia_contacts_2026-07-18.md`).
AR Bank isolates are curated CDC outbreak/surveillance isolates -> **provenance-separable** from the
decoder's broad NCBI-PD tuning pull (the external-cohort BioSample-level preflight verifies
disjointness at run time).

Enumeration (all discovered live 2026-07-18, not guessed):
  - 34 panels at ``Panel/PanelDetail.aspx?ID=<panel_id>``; each panel ``<tbody>`` row carries
    AR# + organism + BioSample + a link ``IsolateDetail.aspx?IsolateID=<i>&PanelID=<p>``.
  - the per-isolate MIC table lives on ``IsolateDetail`` (``<tr><td>Drug<td>MIC<td>INT``).

This module bridges STRAIGHT into the FROZEN external-cohort arm: ``to_label_inputs(details, drug)``
yields the ``({biosample:[mic tokens]}, {biosample:{calls}})`` that
``dna_decode.data.external_mic_labels.build_drug_labels`` consumes DIRECTLY — no crosswalk
(the page gives the BioSample), unlike the Oxford ingester.

The PARSERS are pure + offline-tested (fixtures under ``tests/data/ar_isolate_bank/``). The FETCH
layer (``fetch`` / ``enumerate_panel`` / ``fetch_isolate_detail``) is lazy, cache-first, stdlib-only.
The frozen decoder surface is untouched — this is a READ-only ingestion adapter.
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from pathlib import Path

BASE = "https://wwwn.cdc.gov/ARIsolateBank"
PANEL_URL = BASE + "/Panel/PanelDetail.aspx?ID={panel_id}"
ISOLATE_URL = BASE + "/Panel/IsolateDetail.aspx?IsolateID={isolate_id}&PanelID={panel_id}"

# The 34 established panels (custom panels excluded — they only re-slice established-panel isolates,
# so every isolate is reachable via an established panel; discovered from Panel/AllIsolate 2026-07-18).
PANEL_IDS: tuple[int, ...] = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
    1034, 1035, 1036, 1125, 1126, 1127, 1128, 1130, 1132, 1133, 1140,
    1153, 1156, 1157, 1162, 1164, 1174, 1222, 1247, 1368,
)

_BIOSAMPLE_RE = re.compile(r"(SAM[NED][A-Z]?\d+)")
_TAG_RE = re.compile(r"<[^>]+>")
_ISOLATE_LINK_RE = re.compile(
    r'IsolateDetail\.aspx\?IsolateID=(\d+)&(?:amp;)?PanelID=(\d+)', re.I)


@dataclass(frozen=True)
class PanelRow:
    """One isolate as listed on a panel page (no MIC — that's on the detail page)."""
    isolate_id: int
    panel_id: int
    ar_number: str          # e.g. "0021"
    organism: str           # e.g. "Citrobacter freundii"
    biosample: str          # e.g. "SAMN04014862" ("" if none shown)
    mechanisms: str = ""     # free text, e.g. "CMY-108, QnrB13"


@dataclass
class IsolateDetail:
    """One isolate's parsed detail page: BioSample + organism + per-drug (MIC, INT)."""
    isolate_id: int
    panel_id: int
    ar_number: str
    organism: str
    biosample: str
    mics: dict[str, str] = field(default_factory=dict)   # drug (as shown) -> raw MIC token
    calls: dict[str, str] = field(default_factory=dict)  # drug (as shown) -> S/I/R


def _strip_tags(s: str) -> str:
    return re.sub(r"\s+", " ", _TAG_RE.sub("", s)).strip()


def _first_biosample(html_fragment: str) -> str:
    m = _BIOSAMPLE_RE.search(html_fragment)
    return m.group(1) if m else ""


# --------------------------------------------------------------------------- panel parsing

def parse_panel_isolates(html: str) -> list[PanelRow]:
    """Parse a PanelDetail page ``<tbody>`` into PanelRows.

    Each row: ``<a href="IsolateDetail.aspx?IsolateID=21&PanelID=1177">0021</a>``, then
    ``<td><i>organism</i>``, a mechanisms cell, and an NCBI BioSample link. Robust to missing
    BioSample / mechanisms cells; skips header/blank rows.
    """
    rows: list[PanelRow] = []
    for m in re.finditer(r"<tr>(.*?)</tr>", html, re.S | re.I):
        tr = m.group(1)
        link = _ISOLATE_LINK_RE.search(tr)
        if not link:
            continue
        isolate_id, panel_id = int(link.group(1)), int(link.group(2))
        cells = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S | re.I)
        ar_number = _strip_tags(cells[0]) if cells else str(isolate_id)
        organism = _strip_tags(cells[1]) if len(cells) > 1 else ""
        mechanisms = _strip_tags(cells[2]) if len(cells) > 2 else ""
        biosample = _first_biosample(tr)
        rows.append(PanelRow(isolate_id=isolate_id, panel_id=panel_id, ar_number=ar_number,
                             organism=organism, biosample=biosample, mechanisms=mechanisms))
    return rows


# --------------------------------------------------------------------------- detail parsing

def _parse_mic_table(html: str) -> tuple[dict[str, str], dict[str, str]]:
    """Extract {drug: raw MIC} + {drug: INT} from the MIC/INT table on a detail page.

    The table header is ``Drug | MIC (μg/ml) | INT``; each row is
    ``<tr><td>Drug <sup>footnotes</sup></td><td...>MICvalue</td><td...>INT</td></tr>``.
    Footnote superscripts are stripped from the drug name; MIC tokens keep their operator
    (``<=2``, ``>8``) so downstream operator-aware tiering works. Empty INT/MIC cells are skipped.
    """
    mics: dict[str, str] = {}
    calls: dict[str, str] = {}
    # anchor on the header cell so we don't pick up unrelated tables
    hdr = re.search(r"MIC\s*\(&mu;g/ml\)|MIC\s*\(.g/ml\)|MIC\s*\(&#181;g/ml\)", html, re.I)
    region = html[hdr.start():] if hdr else html
    body = re.search(r"<tbody>(.*?)</tbody>", region, re.S | re.I)
    region = body.group(1) if body else region
    for m in re.finditer(r"<tr>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*</tr>",
                         region, re.S | re.I):
        drug = _strip_tags(re.sub(r"<sup>.*?</sup>", "", m.group(1), flags=re.S | re.I))
        mic = _strip_tags(m.group(2)).replace("&lt;", "<").replace("&gt;", ">").replace("&le;", "<=").replace("&ge;", ">=")
        intr = _strip_tags(m.group(3)).upper()
        if not drug:
            continue
        if mic:
            mics[drug] = mic
        if intr in ("S", "I", "R", "SDD", "NS"):
            calls[drug] = intr
    return mics, calls


def parse_isolate_detail(html: str, *, isolate_id: int = 0, panel_id: int = 0) -> IsolateDetail | None:
    """Parse an IsolateDetail page. Returns None for a blank/template page (no BioSample AND no MIC).

    - BioSample: the ``MainContent_lblSeqAcc`` NCBI link (falls back to any SAMN in the page).
    - AR number + organism: the ``AR Bank #NNNN`` header + the following ``<i>genus species</i>``.
    """
    seq = re.search(r'lblSeqAcc"?>(.*?)</span>', html, re.S | re.I)
    biosample = _first_biosample(seq.group(1)) if seq else _first_biosample(html)
    hdr = re.search(r"AR\s*Bank\s*#\s*(\d+)\s*</b>.*?<i>(.*?)</i>", html, re.S | re.I)
    ar_number = hdr.group(1) if hdr else ""
    organism = _strip_tags(hdr.group(2)) if hdr else ""
    mics, calls = _parse_mic_table(html)
    if not biosample and not mics:
        return None   # template/blank render (no PanelID etc.)
    return IsolateDetail(isolate_id=isolate_id, panel_id=panel_id, ar_number=ar_number,
                         organism=organism, biosample=biosample, mics=mics, calls=calls)


# --------------------------------------------------------------------------- bridge to the frozen arm

def organism_matches(organism: str, needles: tuple[str, ...]) -> bool:
    """Case-insensitive substring match of an organism against any of `needles`."""
    o = organism.lower()
    return any(n.lower() in o for n in needles)


def to_label_inputs(details, drug: str) -> tuple[dict[str, list[str]], dict[str, dict]]:
    """Bridge parsed IsolateDetails -> the two build_drug_labels inputs for ONE drug.

    Returns ({biosample: [raw MIC token]}, {biosample: {S/I/R call}}), keyed by BioSample and
    filtered to isolates that (a) have a BioSample and (b) list that drug (canonicalized). The
    drug match is via `external_mic_labels.canonical_drug` so AR Bank names (Ciprofloxacin,
    Ceftriaxone, Gentamicin) map onto the pilot canon; non-pilot drugs yield empty inputs.
    """
    from dna_decode.data.external_mic_labels import canonical_drug
    canon = canonical_drug(drug)
    if canon is None:
        return {}, {}     # non-pilot drug is not scorable (every non-pilot name also canon->None)
    iso_mics: dict[str, list[str]] = {}
    iso_calls: dict[str, set] = {}
    for d in details:
        if not d.biosample:
            continue
        for raw_drug, mic in d.mics.items():
            if canonical_drug(raw_drug) != canon:
                continue
            iso_mics.setdefault(d.biosample, []).append(mic)
            call = d.calls.get(raw_drug)
            if call in ("R", "S"):   # I/SDD/NS are not decisive calls; MIC tiering decides
                iso_calls.setdefault(d.biosample, set()).add(call)
    return iso_mics, {bs: c for bs, c in iso_calls.items()}


# --------------------------------------------------------------------------- fetch layer (lazy)

def fetch(url: str, cache_dir: str | Path, *, name: str | None = None,
          sleep: float = 0.5, offline_ok: bool = False) -> str | None:
    """Cache-first GET. Reads ``<cache_dir>/<name>.html`` if present; else downloads + caches.

    Returns None if offline_ok and the fetch fails (so an offline caller degrades rather than
    raising). A one-request politeness ``sleep`` is applied only on a real network hit.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    fname = (name or re.sub(r"[^A-Za-z0-9]+", "_", url)[-120:]) + ".html"
    fp = cache_dir / fname
    if fp.exists():
        return fp.read_text(encoding="utf-8", errors="replace")
    import urllib.request
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 dna_decode-research"})
        html = urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")
    except Exception:
        if offline_ok:
            return None
        raise
    fp.write_text(html, encoding="utf-8")
    if sleep:
        time.sleep(sleep)
    return html


def enumerate_panel(panel_id: int, cache_dir: str | Path, **kw) -> list[PanelRow]:
    """Fetch + parse one panel's isolate rows."""
    html = fetch(PANEL_URL.format(panel_id=panel_id), cache_dir, name=f"panel_{panel_id}", **kw)
    return parse_panel_isolates(html) if html else []


def enumerate_all(cache_dir: str | Path, panel_ids: tuple[int, ...] = PANEL_IDS,
                  organism_needles: tuple[str, ...] | None = None, **kw) -> list[PanelRow]:
    """Enumerate every isolate across `panel_ids`, deduped by BioSample (an isolate recurs across
    panels). Optionally filter to organisms matching `organism_needles`. Rows without a BioSample
    are kept only if no filter is given (can't dedup/score them, but visible for auditing)."""
    seen: dict[str, PanelRow] = {}
    orphans: list[PanelRow] = []
    for pid in panel_ids:
        for row in enumerate_panel(pid, cache_dir, **kw):
            if organism_needles and not organism_matches(row.organism, organism_needles):
                continue
            if row.biosample:
                seen.setdefault(row.biosample, row)
            else:
                orphans.append(row)
    out = list(seen.values())
    if organism_needles is None:
        out += orphans
    return out


def fetch_isolate_detail(ref: PanelRow, cache_dir: str | Path, **kw) -> IsolateDetail | None:
    """Fetch + parse one isolate's MIC detail page (keyed by the panel row's IDs)."""
    html = fetch(ISOLATE_URL.format(isolate_id=ref.isolate_id, panel_id=ref.panel_id),
                 cache_dir, name=f"isolate_{ref.isolate_id}_{ref.panel_id}", **kw)
    if not html:
        return None
    return parse_isolate_detail(html, isolate_id=ref.isolate_id, panel_id=ref.panel_id)
