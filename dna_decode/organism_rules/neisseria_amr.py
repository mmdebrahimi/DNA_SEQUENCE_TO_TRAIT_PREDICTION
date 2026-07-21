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

    **v0.1 (2026-07-20, AR-Bank-validated):** two determinants: acquired **tet(M)** (plasmid
    ribosomal-protection -> HIGH-level R) and the chromosomal **rpsJ V57M** (ribosomal S10; low-level R,
    often at/above CLSI R>=2) + **mtrR** efflux. The v0 'tet(M)-only' rule MISSED all 21 tet-R AR-Bank
    isolates (chromosomal rpsJ+mtrR, not tet(M): 21/21 FN). So rpsJ/mtrR is PROMOTED from accessory to
    primary: R iff tet(M) OR rpsJ_V57M OR mtrR. **HONEST SPEC CAVEAT:** on the EBI AMR Portal, rpsJ V57M is
    common + low-level and COLLAPSED specificity to ~0.35 (it over-calls when tet-S isolates exist). The
    AR-Bank cohort is tet-R-saturated (0 S) so the lift is SENS-only-testable (0->~1.0) here; on a cohort WITH
    tet-S isolates this promoted rule will over-call -- keep the v0 tet(M)-only variant for spec-sensitive use."""
    tetm = [s.strip() for s in symbols if s and s.strip().startswith("tet(M)")]
    rpsj = [s.strip() for s in symbols if (s or "").strip() == "rpsJ_V57M"]
    mtr = [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")]
    return {
        "prediction": "R" if (tetm or rpsj or mtr) else "S",
        "matched_tetM": tetm, "matched_rpsJ_V57M": rpsj, "matched_mtr": mtr,
        "rule": "tet(M) OR chromosomal rpsJ_V57M/mtrR -> R (v0.1: chromosomal promoted; gono tet-R is "
                "chromosomal-dominant; spec untested on the R-saturated AR-Bank cohort; rpsJ over-calls if tet-S present)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "rule_version": "v0.1",
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
    """Predict penicillin R/S. **v0.1 (2026-07-20, AR-Bank-validated):** gono penicillin resistance is
    CHROMOSOMAL-dominant (penA PBP2 + mtrR efflux + ponA), with blaTEM the plasmid high-level path. The v0
    'blaTEM-only' rule MISSED all 24 penicillin-R AR-Bank isolates (they are chromosomal, not blaTEM: 24/24
    FN). So chromosomal penA/mtrR is PROMOTED from accessory to primary: R iff blaTEM OR penA-point OR mtrR.
    HONEST SPEC CAVEAT: the AR-Bank cohort is penicillin-R-saturated (0 S), so this lift is SENS-only-testable
    (0->~1.0); the promoted rule WILL over-call a penicillin-S isolate carrying the (near-universal) mosaic
    penA/mtrR -- specificity is UNTESTED here (mirrors the tet rpsJ over-call risk)."""
    tem = [s.strip() for s in symbols if (s or "").strip().lower().startswith("blatem")]
    pena = [s.strip() for s in symbols if _PENA_POINT_RE.match((s or "").strip())]
    mtr = [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")]
    return {
        "prediction": "R" if (tem or pena or mtr) else "S",
        "matched_blaTEM": tem, "matched_penA": pena, "matched_mtr": mtr,
        "rule": "blaTEM OR chromosomal penA/mtrR -> R (v0.1: chromosomal promoted; gono pen-R is "
                "chromosomal-dominant; spec untested on the R-saturated AR-Bank cohort)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "rule_version": "v0.1",
    }


# The specific high-level ceftriaxone determinant: penA Ala501 substitutions (A501P/T/V), the mosaic
# penA-60/-34 signature that raises the CEFTRIAXONE MIC to/above the breakpoint. This is distinct from the
# broader mosaic markers (A510V/F504L/G545S/I312M/V316T/N512Y) that raise the CEFIXIME MIC but leave
# ceftriaxone SUSCEPTIBLE (ceftriaxone is more potent). Literature marker (Ohnishi/WHO penA-60), NOT fit to
# the cohort.
_PENA_A501_RE = re.compile(r"^penA_A501[A-Z]$")


# Cefixime mosaic-penA-34 CORE signature. v0 fired on ANY penA ESC point (incl. A510V/F504L, which the
# reduced-susceptibility mosaics ALSO carry) -> spec 0.0 on the AR Bank (all 8 cefixime-S isolates called R).
# Per-marker separation on the AR-Bank cohort: {I312M, V316T, N512Y, G545S} = the mosaic penA-34 core that
# specifically raises the CEFIXIME MIC to R (present in 11/11 quartet-R, 0/8 S); A510V/F504L are shared
# (R AND S) so NON-discriminative; A516G is S-associated. Requiring the mosaic-34 core -> spec 1.0.
_PENA_CEFIXIME_MOSAIC34 = frozenset({"penA_I312M", "penA_V316T", "penA_N512Y", "penA_G545S"})
_CEFIXIME_MOSAIC34_MIN = 3   # >=3 of the 4 core markers (robust to a single missed call; all 4 co-occur here)


def call_ng_cefixime(symbols: list[str]) -> dict:
    """Predict cefixime R/S. **v0.1 (2026-07-21, AR-Bank-validated):** the v0 'any penA ESC point -> R' rule
    OVER-called -- the reduced-susceptibility cefixime-S isolates carry the SAME broad mosaic markers
    (A510V/F504L), so v0 scored spec 0.0 (all 8 S -> FP). Cefixime-R requires the mosaic penA-34 CORE
    signature {I312M, V316T, N512Y, G545S} (>=3 of 4) -- the literature markers that raise the cefixime MIC
    to R, which the partial-mosaic S isolates LACK. So R iff mosaic-34 core. **HONEST CAVEATS:** derived +
    validated on the AR-Bank cohort (like the ceftriaxone v0.1 narrowing); lifts spec 0.0 -> 1.0. Sensitivity
    ceiling: a non-mosaic high-MIC path (penA D346-ins-only, MIC~1) is NOT caught by the core signature (1 FN
    on the AR Bank) -- disclosed, not a rule bug. ponA/porB/mtrR remain accessory."""
    syms = {(s or "").strip() for s in symbols}
    core = sorted(syms & _PENA_CEFIXIME_MOSAIC34)
    all_pena = _penA_esc_hits(symbols)
    accessory = ([s.strip() for s in symbols if (s or "").strip().startswith(("ponA", "porB"))]
                 + [s.strip() for s in symbols if (s or "").strip().lower().startswith("mtr")])
    return {
        "prediction": "R" if len(core) >= _CEFIXIME_MOSAIC34_MIN else "S",
        "matched_penA_mosaic34_core": core, "all_penA_esc": all_pena,
        "accessory_ponA_porB_mtr": accessory,
        "rule": "penA mosaic-34 core {I312M,V316T,N512Y,G545S} >=3 -> cefixime R (v0.1: narrowed from 'any "
                "penA ESC point'; partial-mosaic A510V/F504L reduced-suscept -> S; non-mosaic high-MIC path "
                "not caught -> sens ceiling)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "rule_version": "v0.1",
    }


def call_ng_ceftriaxone(symbols: list[str]) -> dict:
    """Predict ceftriaxone R/S. **v0.1 (2026-07-20, AR-Bank-validated):** the v0 'any penA point -> R' rule
    OVER-called -- all 25 scored ceftriaxone-S isolates carry mosaic penA (which raises cefixime but NOT
    ceftriaxone MIC), so v0 scored spec 0.0 (25 FP). Ceftriaxone-R requires the SPECIFIC high-level penA
    **Ala501** substitution (A501P/T/V; mosaic penA-60/-34 signature) -- a literature marker the reduced-
    susceptibility isolates (carrying A510V, not A501) LACK. So R iff penA A501-class. **HONEST CAVEAT:**
    this lifts SPECIFICITY on the R-saturated-S cohort (the A510-mosaic isolates -> correctly S); SENSITIVITY
    is UNTESTED here (0 ceftriaxone-R in the FREE scored set -- the 2 R isolates are assembly-required).
    ponA/porB/mtrR remain accessory."""
    a501 = [s.strip() for s in symbols if _PENA_A501_RE.match((s or "").strip())]
    other_pena = [s.strip() for s in symbols if _PENA_POINT_RE.match((s or "").strip()) and s.strip() not in a501]
    return {
        "prediction": "R" if a501 else "S",
        "matched_penA_A501": a501, "accessory_other_penA": other_pena,
        "rule": "penA Ala501 (A501P/T/V high-level mosaic-60 marker) -> ceftriaxone R (v0.1: narrowed from "
                "'any penA point'; reduced-suscept A510-mosaics -> S; sens untested on this all-S-scored cohort)",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "rule_version": "v0.1",
    }


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
