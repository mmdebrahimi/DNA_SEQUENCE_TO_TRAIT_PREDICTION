"""Measured catalog-completeness gaps in the target-site catalogs -- the L2 index, not a call surface.

WHAT THIS IS. The doubt layer's AMR arm ranks determinant FAMILIES the deployed rule cannot represent by
a purity signature (carriers labelled resistant, never susceptible). The 2026-09-02 probe measured whether
that same signature is well-formed on the target-site arm and found that it is -- one vocabulary, not two
(`wiki/doubt_target_site_denominator_2026-09-02.md`). This module is the committed result of that
measurement, so the finding reaches a CALL instead of living only in a probe artifact.

WHY IT IS SEPARATE FROM THE POSITION-NOVELTY FLAG, and this is load-bearing rather than tidy: the two
cover DIFFERENT blind spots and neither subsumes the other.

  position-novelty   fires on a novel substitution AT a catalogued position (K103R -> fires).
  completeness       fires on a substitution at a position the catalog does NOT carry (V179F -> fires).

Verified: `flag_for_cell(['V179F'], 'hiv-nnrti-rt').position_novel` is False, so without this index a
V179F carrier receiving a susceptible call gets NO doubt at all from L2.

WHAT IT IS NOT. Not a catalog, not a call, and not a curation recommendation. Data-derived NNRTI curation
was measured three ways on 2026-09-01 and DECLINED (every variant recovered less of the blind spot than
the free position-novelty flag, and the best-scoring one deleted canonical Y181C). L2 qualifies a call
without competing with it; entries here can only ever attach doubt to a susceptible call.

ENTRY BAR. A unit is listed only if it survived a FAMILY-WISE correction over the units actually tested
in its cell -- 638 for HIV NNRTI/EFV -- with at least 5 carriers and zero susceptible ones. Purity is the
signature: one susceptible carrier is positive evidence the exclusion is deliberate and ENDS the signal.
"""
from __future__ import annotations

# cell -> substitution -> the measurement that put it here. Every field is traceable to the probe
# artifact; nothing is asserted from memory or from literature.
TARGET_SITE_COMPLETENESS: dict[str, dict[str, dict]] = {
    "hiv-nnrti-rt": {
        "V179F": {
            "carriers_labelled_r": 15,
            "carriers_labelled_s": 0,
            "purity_surprise_p": 8.83e-06,
            "n_units_tested": 638,
            "base_s_rate": 0.5397,
            "scored_on": "efavirenz",
            "label": "Stanford HIVDB PhenoSense fold-change, DRMcv.R cutoff 3.0",
            "artifact": "wiki/doubt_target_site_denominator_probe.json",
            # Corroboration, NOT the reason it is listed: the 2026-09-01 curation measurement reached
            # position 179 independently by multivariate OLS, naming V179D (12 carriers) and V179E (3).
            "independent_corroboration": "V179D/V179E named by the 2026-09-01 OLS curation measurement",
        },
    },
}

# Cells where the screen COULD apply (mutant-level catalogs) but has not been measured. Reporting these
# as "no doubt" would be a false clean bill -- the same three-state discipline the position-novelty flag
# already holds. They are listed so the gap between "measured clean" and "never measured" stays visible.
UNMEASURED_CELLS: frozenset[str] = frozenset({
    "sarscov2-mpro",                  # CoV-RDB fold-change exists but is TN-starved (37R/5S)
    "fungal-fluconazole-erg11",       # no free isolate-level phenotype source
    "fungal-voriconazole-erg11",
})


def completeness_units_for(cell: str) -> dict[str, dict]:
    """Measured gap units for a cell. Empty dict means measured-and-empty OR never measured -- callers
    must distinguish those via `is_measured`, never treat empty as clean."""
    return TARGET_SITE_COMPLETENESS.get(cell, {})


def is_measured(cell: str) -> bool:
    return cell in TARGET_SITE_COMPLETENESS


def matching_units(observed_substitutions, cell: str) -> list[tuple[str, dict]]:
    """Observed substitutions that hit a measured completeness gap, in stable order."""
    units = completeness_units_for(cell)
    seen = {str(s).strip().upper() for s in (observed_substitutions or [])}
    return [(s, units[s]) for s in sorted(units) if s.upper() in seen]
