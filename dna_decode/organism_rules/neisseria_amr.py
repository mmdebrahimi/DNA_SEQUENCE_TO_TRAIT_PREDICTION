"""Neisseria gonorrhoeae ciprofloxacin — curated organism-specific determinant rule (Tier-3 cell, NON-frozen).

First AMR-Portal Tier-3/4 curated cell (plan `plans/AMR_Portal_Tier34_Curation_Plan_2026-07-01.md`).
Mechanism class A (target-site point mutation) — genuinely decodable. NG cipro resistance is driven by
**gyrA QRDR** substitutions at Ser91 and Asp95 (the canonical first-/second-step fluoroquinolone-resistance
mutations; Belland 1994, WHO GC surveillance): gyrA S91F is by far the dominant determinant, with S91Y and
D95N/A/G also resistance-conferring. parC (Ser87/Ser88/Glu91) mutations are ACCESSORY — they raise MIC on a
gyrA background but rarely confer resistance alone — so the primary binary rule is gyrA-QRDR-driven.

NON-FROZEN + scorer-local (like `organism_rules/tb_*` + `experimental_drug_rules`): this is NOT in the frozen
`amr_rules.py` / `calibrated_amr_rules.json` deployed surface (those stay byte-pinned). Scored on the EBI AMR
Portal via `scripts/neisseria_cipro_amr_portal_validate.py`; endorsed only if it clears the spec>=0.85
falsifier on the powered provenance-disjoint set.

Input = the isolate's AMRFinderPlus point-mutation symbols (the Portal reports them in `amr_element_symbol`,
e.g. `gyrA_S91F`, with `element_subtype=POINT`).
"""
from __future__ import annotations

import re

DRUG = "ciprofloxacin"
ORGANISM = "Neisseria gonorrhoeae"

# QRDR resistance codons (any substitution at these positions confers/contributes to FQ resistance).
GYRA_QRDR_CODONS = frozenset({91, 95})          # Ser91, Asp95 — primary
PARC_QRDR_CODONS = frozenset({87, 88, 91})      # Ser87, Ser88, Glu91 — accessory (recorded, not required)
_POINT_RE = re.compile(r"^(gyrA|parC)_([A-Z])(\d+)([A-Z*])$")


def _qrdr_hits(symbols: list[str]) -> dict[str, list[str]]:
    """Parse AMRFinder point-mutation symbols -> {'gyrA': [...], 'parC': [...]} of QRDR-codon substitutions."""
    hits = {"gyrA": [], "parC": []}
    for s in symbols:
        m = _POINT_RE.match((s or "").strip())
        if not m:
            continue
        gene, _wt, pos, _sub = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        codons = GYRA_QRDR_CODONS if gene == "gyrA" else PARC_QRDR_CODONS
        if pos in codons:
            hits[gene].append(s.strip())
    return hits


DRUG_TET = "tetracycline"


def call_ng_tetracycline(symbols: list[str]) -> dict:
    """Predict tetracycline R/S for N. gonorrhoeae.

    Two determinants: acquired **tet(M)** (plasmid ribosomal-protection -> HIGH-level R) and the **rpsJ V57M**
    ribosomal S10 point mutation (chromosomal -> LOW-level R, often at/above the CLSI R>=2 breakpoint). Either
    -> R. (mtrR/porB efflux variants contribute to multidrug tolerance but are not the primary tet
    determinant.) Non-frozen / scorer-local; endorsed only if it clears the spec>=0.85 falsifier on the AMR
    Portal (the rpsJ-only isolates are the spec risk if V57M sits below the breakpoint)."""
    tetm = [s.strip() for s in symbols if s and s.strip().startswith("tet(M)")]
    rpsj = [s.strip() for s in symbols if (s or "").strip() == "rpsJ_V57M"]
    # tet(M)-ONLY: rpsJ V57M is common + LOW-level (frequently below the CLSI R>=2 breakpoint), so including
    # it collapses specificity (empirically 0.35 on the AMR Portal). Demoted to accessory context.
    return {
        "prediction": "R" if tetm else "S",
        "matched_tetM": tetm, "accessory_rpsJ_V57M": rpsj,
        "rule": "tet(M) (high-level ribosomal protection) -> R (rpsJ V57M accessory-only: low-level, over-calls)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


def call_ng_ciprofloxacin(symbols: list[str]) -> dict:
    """Predict cipro R/S for N. gonorrhoeae from its determinant symbols.

    R iff >=1 gyrA QRDR (Ser91/Asp95) substitution is present (the deterministic primary rule). parC QRDR
    hits are surfaced as accessory context but do NOT change the binary call on their own."""
    hits = _qrdr_hits(symbols)
    resistant = bool(hits["gyrA"])
    return {
        "prediction": "R" if resistant else "S",
        "matched_gyrA_qrdr": hits["gyrA"],
        "accessory_parC_qrdr": hits["parC"],
        "rule": "gyrA QRDR (Ser91/Asp95) point mutation -> R (parC accessory-only)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


# --------------------------------------------------------------------------------------------------
# Extension 2026-07-20: the remaining AR-Bank gonococcal panel drugs (azithromycin / ceftriaxone /
# cefixime / penicillin) + an explicit gentamicin ABSTAIN. Determinants sourced from the Pathogenwatch
# N. gonorrhoeae AMR library (`data/gono_catalogue/485.toml`, sha256 00f44d16..., taxid 485) +
# `wiki/gonorrhoeae_amr_determinant_catalog_2026-07-19.md`. Same NON-FROZEN / scorer-local posture as
# the cipro/tet rules; primary high-level determinant -> R, low-level chromosomal determinants demoted
# to accessory (same over-call discipline that keeps rpsJ/parC out of the binary call). Endorsed only
# if each clears the spec>=0.85 falsifier on an independent measured-MIC cohort (AR Bank / AMR Portal).
# --------------------------------------------------------------------------------------------------

DRUG_AZM = "azithromycin"
DRUG_CRO = "ceftriaxone"
DRUG_CFM = "cefixime"
DRUG_PEN = "penicillin"
DRUG_GEN = "gentamicin"

# High-level macrolide determinant: 23S rRNA peptidyl-transferase mutations. The 4 rDNA copies collapse
# in short-read assemblies, so a call is consensus-base dependent (a known sensitivity caveat). Accept
# BOTH the gonococcal WHO-2016 coordinates (A2045G / C2597T) and the E. coli-equivalent coordinates
# (A2059G / C2611T) that some callers emit.
_23S_AZM_MUTS = ("A2045G", "C2597T", "A2059G", "C2611T")
# penA PBP2 determinant. AMRFinderPlus (-O Neisseria_gonorrhoeae) reports ONLY resistance-CURATED penA
# point mutations, each tagged Subclass=CEPHALOSPORIN (verified on real output: penA_I312M/V316T/F504L/
# A510V/N512Y/G545S). So we match ANY penA point mutation (gene-level style, like tet(M)) rather than a
# hard-coded codon list -- a partial codon set silently missed real mosaic positions 504/510 (R3 catch).
_PENA_POINT_RE = re.compile(r"^penA_[A-Z]\d+[A-Z*]$")


def _penA_esc_hits(symbols: list[str]) -> list[str]:
    """Any curated penA point mutation (AMRFinder emits only ESC/penicillin-relevant ones) or a mosaic tag."""
    out = []
    for s in symbols:
        t = (s or "").strip()
        if _PENA_POINT_RE.match(t) or (t.startswith("penA") and "mosaic" in t.lower()):
            out.append(t)
    return out


def call_ng_azithromycin(symbols: list[str]) -> dict:
    """Predict azithromycin R/S. Primary = 23S rRNA A2045G/C2597T (WHO coords; = A2059G/C2611T E. coli)
    -> high-level macrolide R. mtrR mosaic / promoter (`a-57del`) efflux variants raise MIC but are
    LOW-level and over-call, so they are accessory-only (same demotion as rpsJ for tet)."""
    hits_23s = [s.strip() for s in symbols if any(mut in (s or "") for mut in _23S_AZM_MUTS)]
    mtr = [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")]
    return {
        "prediction": "R" if hits_23s else "S",
        "matched_23S": hits_23s, "accessory_mtr": mtr,
        "rule": "23S rRNA A2045G/C2597T (macrolide target) -> R (mtrR efflux accessory-only: low-level, over-calls)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


def call_ng_penicillin(symbols: list[str]) -> dict:
    """Predict penicillin R/S. Primary = plasmid **blaTEM** (TEM-1/135 penicillinase) -> high-level R.
    Chromosomal penA/mtrR/ponA/porB confer LOW-to-intermediate resistance and over-call, so accessory-only."""
    tem = [s.strip() for s in symbols if (s or "").strip().lower().startswith("blatem")]
    chrom = ([s.strip() for s in symbols if (s or "").strip().startswith(("penA", "ponA", "porB"))]
             + [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")])
    return {
        "prediction": "R" if tem else "S",
        "matched_blaTEM": tem, "accessory_chromosomal": chrom,
        "rule": "blaTEM plasmid penicillinase -> R (chromosomal penA/mtrR/ponA/porB accessory-only: low-level)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


def _call_ng_esc(symbols: list[str], drug: str) -> dict:
    """Shared ceftriaxone/cefixime rule: primary = penA mosaic / ESC-associated PBP2 substitution ->
    decreased ESC susceptibility / R. ponA L421P + porB + mtrR modulate the level (accessory). ESC
    resistance in NG is subtle (mostly reduced-susceptibility, few fully-R) -> flagged for validation."""
    pena = _penA_esc_hits(symbols)
    accessory = ([s.strip() for s in symbols if (s or "").strip().startswith(("ponA", "porB"))]
                 + [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")])
    return {
        "prediction": "R" if pena else "S",
        "matched_penA_esc": pena, "accessory_ponA_porB_mtr": accessory,
        "rule": f"penA mosaic/ESC PBP2 substitution -> {drug} R (ponA/porB/mtrR accessory; ESC-R is "
                "subtle/reduced-susceptibility -> validate spec on measured MIC)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local",
    }


def call_ng_ceftriaxone(symbols: list[str]) -> dict:
    return _call_ng_esc(symbols, DRUG_CRO)


def call_ng_cefixime(symbols: list[str]) -> dict:
    return _call_ng_esc(symbols, DRUG_CFM)


def call_ng_gentamicin(symbols: list[str]) -> dict:
    """Gentamicin has NO validated gonococcal genetic determinant (an alternative/newer agent; absent from
    the Pathogenwatch 485.toml catalogue and the literature). A determinant rule cannot decode it, so this
    ABSTAINS (INDETERMINATE) rather than calling a vacuous S-by-absence. The external scorer excludes
    INDETERMINATE from sens/spec (honest non-coverage, mirroring the fungal/antiviral abstainers)."""
    return {
        "prediction": "INDETERMINATE",
        "rule": "no validated gonococcal gentamicin determinant -> ABSTAIN (not decodable)",
        "rule_status": "ABSTAIN_NO_DETERMINANT", "rule_scope": "scorer_local",
    }


_NG_DISPATCH = {
    DRUG: call_ng_ciprofloxacin, "ciprofloxacin": call_ng_ciprofloxacin,
    DRUG_TET: call_ng_tetracycline,
    DRUG_AZM: call_ng_azithromycin,
    DRUG_CRO: call_ng_ceftriaxone,
    DRUG_CFM: call_ng_cefixime,
    DRUG_PEN: call_ng_penicillin,
    DRUG_GEN: call_ng_gentamicin,
}


def call_ng_amr(drug: str, symbols: list[str]) -> dict:
    """Dispatch to the per-drug gonococcal rule for the AR-Bank panel (azithromycin / cefixime /
    ceftriaxone / ciprofloxacin / gentamicin / penicillin / tetracycline). An unsupported drug ABSTAINS
    (INDETERMINATE) rather than guessing."""
    fn = _NG_DISPATCH.get((drug or "").strip().lower())
    if fn is None:
        return {"prediction": "INDETERMINATE", "rule": f"no gonococcal rule for {drug!r} -> ABSTAIN",
                "rule_status": "ABSTAIN_UNSUPPORTED_DRUG", "rule_scope": "scorer_local"}
    return fn(symbols)
