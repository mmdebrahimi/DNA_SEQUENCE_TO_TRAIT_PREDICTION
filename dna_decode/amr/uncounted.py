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


def primary_mechanism_misses(uncounted: list[dict], drug: str) -> list[dict]:
    """The uncounted determinants that are a GENUINE miss rather than a deliberate exclusion. PURE.

    Two tests, both drug-general -- no gene name is hardcoded:

    1. **The subclass must be GENERIC.** The frozen rules match a SPECIFIC drug Subclass. A determinant
       whose Subclass names a DIFFERENT specific drug (`aph(3')-Ia` -> KANAMYCIN, `aadA5` -> STREPTOMYCIN)
       is excluded CORRECTLY -- counting it would over-call. A determinant filed under the bare CLASS name
       (`rmtE1` -> Class AMINOGLYCOSIDE, Subclass AMINOGLYCOSIDE) names no drug at all, and is therefore
       invisible to ANY subclass rule by construction. That is the shape of a real miss.
    2. **The mechanism must be a PRIMARY one for this drug**, via the curated
       `mic_tiers.classify_gene_symbol` + `primary_mechanisms_for` catalogs.

    Test 1 is load-bearing on its own: `classify_gene_symbol('gentamicin', "aph(3')-Ia")` returns
    `aminoglycoside_modifying_enzymes`, which IS a primary gentamicin mechanism -- so the mechanism test
    ALONE would flag a kanamycin gene as a miss and re-introduce exactly the noise this replaces.

    Replaces an earlier `startswith('rmt'/'arma'/'npma')` prefix match, which hardcoded one gene family
    and would mis-fire on any future symbol merely starting with those letters.
    """
    try:
        from dna_decode.data.mic_tiers import classify_gene_symbol, primary_mechanisms_for  # frozen; READ-only
        primary = {m for m in primary_mechanisms_for(drug)}
    except Exception:  # noqa: BLE001 -- a drug with no curated catalog simply has no primary-miss concept
        return []
    if not primary:
        return []

    out = []
    for u in uncounted:
        cls, sub = str(u.get("class", "")), str(u.get("subclass", ""))
        if sub and sub != cls:          # names a specific, DIFFERENT drug -> deliberate exclusion
            continue
        try:
            mech = classify_gene_symbol(drug, str(u.get("symbol", "")))
        except Exception:  # noqa: BLE001
            mech = ""
        if mech and mech in primary:
            out.append({**u, "mechanism": mech})
    return out


def disclosure_provenance() -> dict:
    """Which catalog produced the warning. PURE-ish (hashes a frozen file).

    The disclosure READS a sha256-pinned file (`mic_tiers.py`) without altering the frozen decision. That
    boundary has to be auditable: without recording WHICH catalog spoke, two decoder generations could
    silently share or diverge on the mechanism catalog and no consumer could tell.
    """
    import hashlib
    from pathlib import Path
    p = Path(__file__).resolve().parent.parent / "data" / "mic_tiers.py"
    try:
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
    except OSError:
        sha = ""
    return {"disclosure_schema": "amr-uncounted-disclosure-v1",
            "disclosure_catalog": "dna_decode/data/mic_tiers.py",
            "disclosure_catalog_sha256": sha,
            "alters_frozen_decision": False}


# (drug, mechanism) pairs where the project has MEASURED that the frozen rule misses real resistance.
# Deliberately evidence-gated rather than inferred, and this is the third narrowing -- each earlier,
# broader trigger was rejected by measurement, not by taste:
#   * every uncounted class-relevant determinant       -> fired on 70% of gentamicin calls
#   * every uncounted PRIMARY-mechanism determinant    -> still fired on 48% of CEFTRIAXONE calls,
#     because a narrow-spectrum `blaTEM-1` also files under the generic BETA-LACTAM subclass and
#     excluding it from ceftriaxone is CORRECT. A warning that fires on half of all calls is not a
#     warning.
# So the bar is a measured gap with a citable artifact. Adding a row is a claim and needs evidence.
_MEASURED_GAPS: dict[str, dict[str, str]] = {
    "gentamicin": {
        "16S_rRNA_methyltransferase":
            "measured prospective sens 0.429 -- 24 of 28 false negatives carried one "
            "(wiki/prospective_lock_first_accrual_2026-08-24.md)",
    },
}


def measured_gap_misses(uncounted: list[dict], drug: str) -> list[dict]:
    """Primary-mechanism misses for which a MEASURED gap is on record. PURE."""
    gaps = _MEASURED_GAPS.get(str(drug).lower(), {})
    if not gaps:
        return []
    return [{**m, "evidence": gaps[m["mechanism"]]}
            for m in primary_mechanism_misses(uncounted, drug) if m.get("mechanism") in gaps]


def render_note(uncounted: list[dict], drug: str, prediction: str | None = None) -> str:
    """Human-readable disclosure. Prints ONLY measured gaps (empty otherwise).

    The full class-relevant list and the wider primary-mechanism list both remain in the JSON record as
    audit metadata; only the HUMAN line is evidence-gated, because that is the surface that goes numb
    when it cries wolf. Measured fire rates: 2% on gentamicin (the real gap) vs 70% / 48% for the two
    broader triggers this replaced.
    """
    misses = measured_gap_misses(uncounted, drug)
    if not misses:
        return ""
    shown = ", ".join(f"{m['symbol']} [{m['mechanism']}]" for m in misses[:5])
    more = f" (+{len(misses) - 5} more)" if len(misses) > 5 else ""
    # evidence is per-(drug, mechanism), so a note can never cite another drug's measurement
    evidence = "; ".join(sorted({m["evidence"] for m in misses}))
    lead = (f"  ** {len(misses)} determinant(s) implicating a mechanism this {drug} rule is MEASURED to "
            f"miss are present but were NOT counted: {shown}{more}.")
    why = (f"\n     AMRFinder files these under the bare CLASS name rather than a specific drug subclass, "
           f"so a subclass-matching rule cannot see them. Evidence: {evidence}.")
    if str(prediction).upper() == "S":
        return lead + why + ("\n     The call above is S. Treat that S as UNRELIABLE for this isolate. "
                             "The prediction is left unchanged because the decoder rule is frozen. **")
    return lead + why + "\n     This does NOT change the call. **"
