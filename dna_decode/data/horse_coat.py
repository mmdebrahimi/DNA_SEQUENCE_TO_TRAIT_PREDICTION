"""Horse base coat-colour decoder (off-pathogen, non-human, MEASURED/OBSERVED-label cell).

The deployed Mendelian two-locus rule (Rieder et al. 2001; UC Davis VGL): MC1R "Extension" (E/e) is
EPISTATIC to ASIP "Agouti" (A/a). Molecular basis: chestnut = MC1R c.901C>T (recessive e); black = ASIP
11-bp exon-2 deletion (recessive a).

    e/e            -> chestnut   (no eumelanin at all, regardless of agouti)
    E_ (E/E,E/e) + A_ (A/A,A/a) -> bay     (black restricted to points)
    E_           + a/a          -> black   (uniform black)

This is a DEPLOYED published rule (VGL sells these tests) -> a rule-INTEGRATION cell, not a novel finding.
Validation must use INDEPENDENTLY OBSERVED colour (not colour assigned FROM genotype = circular). See
`scripts/horse_coat_validate.py` for the published-observed-contingency validation + its honest caveats.
"""
from __future__ import annotations

INDETERMINATE = "INDETERMINATE"


def _locus(genotype: str, dominant: str, recessive: str) -> str | None:
    """Normalise a 2-allele locus string -> 'hom_dom' / 'het' / 'hom_rec' / None.

    Accepts e.g. 'EE'/'Ee'/'eE'/'ee' (or 'E/e', 'E E'); case-sensitive on the allele letters."""
    g = "".join(c for c in (genotype or "") if c in (dominant + recessive))
    if len(g) != 2:
        return None
    n_dom = g.count(dominant)
    return "hom_dom" if n_dom == 2 else ("het" if n_dom == 1 else "hom_rec")


def call_horse_base_colour(mc1r: str, asip: str) -> str:
    """Predict base coat colour (chestnut/bay/black) from MC1R (E/e) + ASIP (A/a) genotypes.

    INDETERMINATE if either locus is not two callable alleles. Epistasis: e/e overrides agouti."""
    e = _locus(mc1r, "E", "e")
    if e is None:
        return INDETERMINATE
    if e == "hom_rec":            # e/e -> chestnut, agouti irrelevant (so ASIP may be absent)
        return "chestnut"
    a = _locus(asip, "A", "a")    # has >=1 E -> agouti decides bay vs black
    if a is None:
        return INDETERMINATE
    return "bay" if a in ("hom_dom", "het") else "black"
