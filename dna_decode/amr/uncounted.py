"""Determinants that are RELEVANT to the drug's class but which the deployed rule did NOT count.

WHY THIS EXISTS
The frozen `call_resistance` rules match a specific AMRFinder **Subclass** string (gentamicin counts
`GENTAMICIN`, ceftriaxone counts `CARBAPENEM`/`CEPHALOSPORIN`, ...). Most of what that excludes is
excluded CORRECTLY -- `aph(3')-Ia` is a kanamycin gene and does not confer gentamicin resistance, so
counting it would over-call. But the same narrowness has a real blind spot, measured on the first
prospective cohort (2026-08-24, `wiki/prospective_lock_first_accrual_2026-08-24.md`):

    E. coli x gentamicin, prospective sens 0.429 -- 24 of the 28 false negatives carry a 16S rRNA
    methyltransferase (rmtE1 / armA), which AMRFinder files under the GENERIC `AMINOGLYCOSIDE`
    subclass. Those isolates decode as S with no hint that a pan-aminoglycoside mechanism is present.

This module makes that visible AT THE POINT OF USE. It is a DISCLOSURE, not a rule:

  * it NEVER changes a prediction -- the frozen decoder's calls stay byte-identical;
  * it lives OUTSIDE the sha256-pinned surface (`amr_rules.py` / `calibrated_amr_rules.json` /
    `mic_tiers.py` / `shipped_decoder_surface.py` / `cohort_manifest.py`), so the prospective lock and
    the reproducibility freeze are untouched. It only READS `mic_tiers.amrfinder_classes_for`.

Same idiom the genome-map layer already uses (`DETERMINANT_PRESENT` for determinants the deployed rule
excludes): surface the evidence, keep the call honest about what it counted.
"""
from __future__ import annotations

_SYMBOL_HEADERS = ("Element symbol", "Gene symbol")   # AMRFinder renamed this between versions


def parse_main_tsv_rows(text: str) -> list[dict]:
    """Parse an AMRFinder `main.tsv` into {symbol, cls, subclass} rows. PURE.

    Header-driven (not positional) because AMRFinder's column set has shifted across versions; a file
    whose header carries none of the known symbol names yields no rows rather than mis-indexed junk.
    """
    lines = [ln for ln in str(text).splitlines() if ln.strip()]
    if not lines:
        return []
    hdr = lines[0].split("\t")
    sym_i = next((hdr.index(h) for h in _SYMBOL_HEADERS if h in hdr), None)
    cls_i = hdr.index("Class") if "Class" in hdr else None
    sub_i = hdr.index("Subclass") if "Subclass" in hdr else None
    if sym_i is None or cls_i is None:
        return []
    out = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= max(sym_i, cls_i, sub_i if sub_i is not None else 0):
            continue
        out.append({"symbol": f[sym_i].strip(),
                    "cls": f[cls_i].strip().upper(),
                    "subclass": (f[sub_i].strip().upper() if sub_i is not None else "")})
    return out


def uncounted_class_determinants(rows: list[dict], drug: str, counted: list | None) -> list[dict]:
    """Rows whose AMRFinder Class is relevant to `drug` but whose symbol the rule did NOT count. PURE.

    `counted` is the deployed call's determinant list (dicts with a 'gene'/'symbol', or bare strings).
    Comparison is by SYMBOL, deduped, order-stable.
    """
    try:
        from dna_decode.data.mic_tiers import amrfinder_classes_for   # frozen; READ-only
        relevant = {c.upper() for c in amrfinder_classes_for(drug)}
    except Exception:       # noqa: BLE001 -- an unknown drug simply has no class list to disclose against
        return []
    if not relevant:
        return []

    counted_syms = set()
    for d in (counted or []):
        s = d.get("gene") or d.get("symbol") or d.get("element_symbol") if isinstance(d, dict) else d
        if s:
            counted_syms.add(str(s).strip())

    out, seen = [], set()
    for r in rows:
        sym = r.get("symbol", "")
        if not sym or sym in counted_syms or sym in seen:
            continue
        if r.get("cls", "") in relevant:
            seen.add(sym)
            out.append({"symbol": sym, "class": r.get("cls", ""), "subclass": r.get("subclass", "")})
    return out


def has_16s_methyltransferase(uncounted: list[dict]) -> bool:
    """Is a 16S rRNA methyltransferase among the uncounted determinants? PURE.

    rmtA-H / armA / npmA confer HIGH-LEVEL resistance to the 4,6-disubstituted deoxystreptamines
    (gentamicin / tobramycin / amikacin), but AMRFinder files them under the generic `AMINOGLYCOSIDE`
    subclass, so a Subclass=GENTAMICIN rule cannot see them. This is the ONE exclusion measured to be a
    genuine miss rather than a deliberate one.
    """
    for u in uncounted:
        s = str(u.get("symbol", "")).lower()
        if s.startswith("rmt") or s.startswith("arma") or s.startswith("npma"):
            return True
    return False


def render_note(uncounted: list[dict], drug: str) -> str:
    """Human-readable disclosure (empty when there is nothing to disclose).

    The specific 16S-methyltransferase warning is emitted ONLY when one is actually detected. An earlier
    version printed it unconditionally, so a ciprofloxacin call that correctly flagged `qnrB19` also
    carried an irrelevant paragraph about gentamicin -- noise that trains a reader to skip the note.
    """
    if not uncounted:
        return ""
    shown = ", ".join(f"{u['symbol']} [{u['subclass'] or u['class']}]" for u in uncounted[:6])
    more = f" (+{len(uncounted) - 6} more)" if len(uncounted) > 6 else ""
    head = (f"  note: {len(uncounted)} determinant(s) of a {drug}-relevant CLASS are present but were NOT "
            f"counted by this rule: {shown}{more}.\n"
            f"        Most such exclusions are DELIBERATE (a broader-class gene that does not confer this "
            f"drug's resistance would over-call). This note does NOT change the call.")
    if has_16s_methyltransferase(uncounted):
        head += ("\n        ** One of these is a 16S rRNA METHYLTRANSFERASE (rmt*/armA/npmA) -- a real "
                 "high-level aminoglycoside mechanism this rule MISSES, because AMRFinder files it under "
                 "the generic AMINOGLYCOSIDE subclass. Measured prospective sens 0.429 on gentamicin; "
                 "24/28 false negatives carried one. See "
                 "wiki/prospective_lock_first_accrual_2026-08-24.md. Treat an S call on this isolate as "
                 "UNRELIABLE. **")
    return head
