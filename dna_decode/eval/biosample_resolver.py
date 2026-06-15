"""Bidirectional BioSample resolver — the cross-accession-type identity bridge.

The leakage gate (`cohort_manifest`) matches accession STRINGS only, so a GCA
assembly and an ERR run that are the SAME physical isolate are invisible to it.
This module resolves both directions through the BioSample (SAMEA…/SAMN…) that
both an NCBI assembly and an ENA run carry:

  - `runs_for_project(project)`     -> [(run_accession, biosample)]   (ENA read_run)
  - `biosample_to_assemblies(bs)`   -> [gca,...]                       (Step 3 genome resolution)
  - `assembly_to_biosample(gca)`    -> biosample | None                (leakage cross-check)

PRIMARY resolver for the BioSample<->Assembly directions is NCBI Entrez
`esearch`+`esummary` (works for both SAMEA and SAMN). ENA portal
`filereport?accession=<biosample>&result=assembly&fields=assembly_accession,sample_accession`
is the FALLBACK (the filereport `accession=` form — UNVERIFIED against the live API
until the pre-scoring smoke runs). When BOTH sources
answer and DISAGREE, the result is None with source="disagreement" — the caller
counts a disagreement as UNRESOLVED for the leakage check, never as disjoint
(C4 from the plan brainstorm).

A BioSample with no linked GCA/GCF is a VALID cohort outcome (-> ASSEMBLY-REQUIRED),
NOT an error.

The HTTP calls are isolated behind an injectable `fetch(url) -> str` so the pure
parse/compose logic is unit-testable offline with no network. Results are cached
to a JSON file (idempotent; network only on cache miss), each entry recording its
source for provenance.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Callable

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENTREZ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

Fetch = Callable[[str], str]


def _default_fetch(url: str, *, timeout: int = 120) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


# --------------------------------------------------------------------------- #
# Pure parsers (no network) — unit-tested directly against fixture strings.
# --------------------------------------------------------------------------- #
def parse_ena_read_run(tsv: str) -> list[tuple[str, str]]:
    """Parse an ENA filereport `result=read_run` TSV into [(run, biosample)].

    Expects a header row with `run_accession` + `sample_accession`. Rows missing
    either field are skipped. NO assembly column exists in read_run.
    """
    lines = [ln for ln in tsv.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        ri = header.index("run_accession")
        si = header.index("sample_accession")
    except ValueError:
        return []
    out: list[tuple[str, str]] = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) <= max(ri, si):
            continue
        run, bs = cells[ri].strip(), cells[si].strip()
        if run and bs:
            out.append((run, bs))
    return out


def parse_ena_read_run_records(tsv: str, fields: tuple[str, ...]) -> list[dict]:
    """Parse an ENA filereport `result=read_run` TSV into a list of per-run dicts.

    Unlike `parse_ena_read_run` (which returns just (run, biosample) tuples), this
    keeps EVERY requested field present in the header — the crosswalk's candidate-key
    source (run_accession / sample_accession / sample_alias / secondary_sample_accession).
    Columns absent from the header are simply omitted from each record (tolerant).
    A record is emitted only if it has at least one non-empty requested field.
    """
    lines = [ln for ln in tsv.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    col_idx = {f: header.index(f) for f in fields if f in header}
    if not col_idx:
        return []
    out: list[dict] = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        rec = {}
        for f, i in col_idx.items():
            if i < len(cells) and cells[i].strip():
                rec[f] = cells[i].strip()
        if rec:
            out.append(rec)
    return out


def parse_ena_assembly(tsv: str) -> list[str]:
    """Parse an ENA filereport `result=assembly` TSV into [assembly_accession,...]."""
    lines = [ln for ln in tsv.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    try:
        ai = header.index("assembly_accession")
    except ValueError:
        return []
    out: list[str] = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) <= ai:
            continue
        acc = cells[ai].strip()
        if acc:
            out.append(acc)
    return out


def parse_entrez_esearch(json_text: str) -> list[str]:
    """Parse an Entrez esearch JSON response into the UID id-list."""
    try:
        d = json.loads(json_text)
    except (ValueError, TypeError):
        return []
    return list(d.get("esearchresult", {}).get("idlist", []) or [])


def parse_entrez_assembly_summary(json_text: str) -> list[dict]:
    """Parse an Entrez esummary (db=assembly) JSON response.

    Returns [{"assembly": <AssemblyAccession>, "biosample": <BioSampleAccn>}, ...]
    over every UID in the result set.
    """
    try:
        d = json.loads(json_text)
    except (ValueError, TypeError):
        return []
    res = d.get("result", {})
    uids = res.get("uids", []) or []
    out: list[dict] = []
    for uid in uids:
        rec = res.get(uid, {})
        out.append({
            "assembly": (rec.get("assemblyaccession") or rec.get("AssemblyAccession") or "").strip(),
            "biosample": (rec.get("biosampleaccn") or rec.get("BioSampleAccn") or "").strip(),
        })
    return out


# --------------------------------------------------------------------------- #
# Resolver — network composition over the pure parsers + a persisted cache.
# --------------------------------------------------------------------------- #
class BioSampleResolver:
    """Bidirectional resolver with a JSON cache. Network only on cache miss."""

    def __init__(self, cache_path: str | Path = "data/raw/_biosample_cache.json",
                 fetch: Fetch | None = None):
        self.cache_path = Path(cache_path)
        self.fetch = fetch or _default_fetch
        self._cache = self._load_cache()

    # -- cache ------------------------------------------------------------- #
    def _load_cache(self) -> dict:
        if self.cache_path.exists():
            try:
                return json.loads(self.cache_path.read_text(encoding="utf-8"))
            except (ValueError, OSError):
                pass
        return {"biosample_to_assemblies": {}, "assembly_to_biosample": {}}

    def save_cache(self) -> None:
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.cache_path.with_suffix(self.cache_path.suffix + ".tmp")
        tmp.write_text(json.dumps(self._cache, indent=2), encoding="utf-8")
        tmp.replace(self.cache_path)

    # -- ENA read_run (run -> biosample) ----------------------------------- #
    def runs_for_project(self, project: str) -> list[tuple[str, str]]:
        url = (f"{ENA_PORTAL}?accession={urllib.parse.quote(project)}"
               f"&result=read_run&fields=run_accession,sample_accession&format=tsv")
        return parse_ena_read_run(self.fetch(url))

    # -- ENA read_run multi-field records (crosswalk candidate-key source) -- #
    def read_run_records_for_project(
        self, project: str,
        fields: tuple[str, ...] = ("run_accession", "sample_accession",
                                   "sample_alias", "secondary_sample_accession"),
    ) -> list[dict]:
        """Fetch per-run records with extra ENA fields for crosswalk construction.

        Additive — does NOT change `runs_for_project`'s (run, sample) tuple contract.
        """
        url = (f"{ENA_PORTAL}?accession={urllib.parse.quote(project)}"
               f"&result=read_run&fields={','.join(fields)}&format=tsv")
        return parse_ena_read_run_records(self.fetch(url), fields)

    # -- biosample -> assemblies (Entrez primary, ENA fallback) ------------ #
    def _entrez_assemblies_for_biosample(self, biosample: str) -> list[str]:
        es = self.fetch(f"{ENTREZ}/esearch.fcgi?db=assembly"
                        f"&term={urllib.parse.quote(biosample)}[BioSample]&retmode=json")
        uids = parse_entrez_esearch(es)
        if not uids:
            return []
        summ = self.fetch(f"{ENTREZ}/esummary.fcgi?db=assembly"
                          f"&id={','.join(uids)}&retmode=json")
        return [r["assembly"] for r in parse_entrez_assembly_summary(summ) if r["assembly"]]

    def _ena_assemblies_for_biosample(self, biosample: str) -> list[str]:
        url = (f"{ENA_PORTAL}?accession={urllib.parse.quote(biosample)}"
               f"&result=assembly&fields=assembly_accession,sample_accession&format=tsv")
        return parse_ena_assembly(self.fetch(url))

    def biosample_to_assemblies(self, biosample: str) -> list[str]:
        cache = self._cache.setdefault("biosample_to_assemblies", {})
        if biosample in cache:
            return list(cache[biosample]["value"])
        try:
            gcas = self._entrez_assemblies_for_biosample(biosample)
            source = "entrez"
        except Exception:  # noqa: BLE001 — fall back to ENA on any Entrez failure
            gcas, source = [], "entrez_error"
        if not gcas:
            try:
                gcas = self._ena_assemblies_for_biosample(biosample)
                source = "ena" if gcas else source
            except Exception:  # noqa: BLE001
                pass
        cache[biosample] = {"value": sorted(set(gcas)), "source": source}
        return list(cache[biosample]["value"])

    # -- assembly -> biosample (Entrez primary, ENA fallback; disagree->None) #
    def _entrez_biosample_for_assembly(self, gca: str) -> str | None:
        es = self.fetch(f"{ENTREZ}/esearch.fcgi?db=assembly"
                        f"&term={urllib.parse.quote(gca)}&retmode=json")
        uids = parse_entrez_esearch(es)
        if not uids:
            return None
        summ = self.fetch(f"{ENTREZ}/esummary.fcgi?db=assembly&id={uids[0]}&retmode=json")
        recs = parse_entrez_assembly_summary(summ)
        return (recs[0]["biosample"] or None) if recs else None

    def _ena_biosample_for_assembly(self, gca: str) -> str | None:
        url = (f"{ENA_PORTAL}?accession={urllib.parse.quote(gca)}"
               f"&result=assembly&fields=assembly_accession,sample_accession&format=tsv")
        tsv = self.fetch(url)
        lines = [ln for ln in tsv.splitlines() if ln.strip()]
        if len(lines) < 2:
            return None
        header = lines[0].split("\t")
        if "sample_accession" not in header:
            return None
        si = header.index("sample_accession")
        cells = lines[1].split("\t")
        return cells[si].strip() if len(cells) > si and cells[si].strip() else None

    def assembly_to_biosample(self, gca: str) -> str | None:
        """Resolve GCA/GCF -> BioSample. Entrez primary; ENA fallback/cross-check.

        If BOTH sources answer and DISAGREE -> None (source="disagreement"); the
        caller treats this as UNRESOLVED for the leakage bound, never disjoint.
        """
        cache = self._cache.setdefault("assembly_to_biosample", {})
        if gca in cache:
            return cache[gca]["value"]
        entrez_bs = ena_bs = None
        try:
            entrez_bs = self._entrez_biosample_for_assembly(gca)
        except Exception:  # noqa: BLE001
            entrez_bs = None
        try:
            ena_bs = self._ena_biosample_for_assembly(gca)
        except Exception:  # noqa: BLE001
            ena_bs = None
        value, source = _reconcile_biosample(entrez_bs, ena_bs)
        cache[gca] = {"value": value, "source": source}
        return value


def _reconcile_biosample(entrez_bs: str | None, ena_bs: str | None) -> tuple[str | None, str]:
    """Reconcile a BioSample answer from the two sources.

    - both present + equal  -> (value, "entrez+ena")
    - both present + differ  -> (None, "disagreement")   [counts as UNRESOLVED]
    - only one present       -> (that one, "entrez" | "ena")
    - neither               -> (None, "unresolved")
    """
    if entrez_bs and ena_bs:
        if entrez_bs == ena_bs:
            return entrez_bs, "entrez+ena"
        return None, "disagreement"
    if entrez_bs:
        return entrez_bs, "entrez"
    if ena_bs:
        return ena_bs, "ena"
    return None, "unresolved"
