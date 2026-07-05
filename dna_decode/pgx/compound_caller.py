"""Compound-allele PGx caller — resolves stars DEFINED BY >1 co-occurring SNP on the same haplotype.

The core `caller.py:assemble_diplotype` maps ONE defining SNP -> one star. Some clinically dominant
star alleles are COMPOUND: they require two (or more) SNPs on the SAME haplotype, and each component SNP
ALONE defines a DIFFERENT allele. Examples (the genuine "multi-SNP calling" validation targets):
  * CYP2B6*6  = rs3745274 (516G>T) + rs2279343 (785A>G) in cis; 516T alone = *9, 785G alone = *4.
  * TPMT*3A   = rs1800460 (*3B) + rs1142345 (*3C) in cis; each alone = *3B / *3C respectively.
A single-SNP tag would MIS-call *4/*9 (or *3B/*3C) as the compound. So this needs true per-haplotype
resolution of a component-SNP SET -> star label. That resolution IS the added validation value.

This module is SEPARATE from the frozen `caller.py` path (CYP2C19/2C9/2C8/3A5 are untouched). It reuses
`caller.scan_vcf` to read the component sites (with phasing), then applies a per-haplotype compound rule.

Honest scope: phased input (1000G 30x is phased) -> exact cis/trans resolution. Unphased with >=2 het
component sites -> cis/trans ambiguous -> flagged (`phase_ambiguous`), the standard-call kept with low
confidence.
"""
from __future__ import annotations

from dataclasses import dataclass

from dna_decode.pgx.caller import DiplotypeResult, _vc_dict, scan_vcf
from dna_decode.pgx.cyp2c19_catalog import DefiningVariant


@dataclass(frozen=True)
class CompoundAllele:
    """A star allele defined by the SET of component-SNP tags that must co-occur on one haplotype."""
    star: str
    components: frozenset[str]   # e.g. frozenset({"516T", "785G"}) for CYP2B6*6


def _haplotype_star(present: frozenset[str], rules: list[CompoundAllele], reference: str) -> str:
    """Given the component tags present on ONE haplotype, pick the best-matching star.

    Rules are matched most-specific-first (largest component set that is a SUBSET of `present`). A
    haplotype with components not matching any rule beyond the reference falls back to `reference`.
    """
    best = None
    for r in sorted(rules, key=lambda r: -len(r.components)):
        if r.components <= present:
            best = r
            break
    return best.star if best else reference


def assemble_compound_diplotype(
    vcf,
    components: list[DefiningVariant],   # one per component SNP; `.star` holds the component TAG
    rules: list[CompoundAllele],
    *,
    reference_allele: str,
    phenotype_fn,
    gene: str,
    sample: str | None = None,
) -> DiplotypeResult:
    """Resolve a compound-allele diplotype from a VCF. `components[i].star` is the component tag used in
    the `rules`; `rules` map component-tag SETS -> star labels. Returns the same DiplotypeResult shape as
    the core caller so the runner/report wiring is uniform."""
    calls = scan_vcf(vcf, defining=components, sample=sample)
    present_calls = [c for c in calls.values() if c.found]
    flags: list[str] = []
    if not present_calls:
        return DiplotypeResult("no_input", None, None, None, None, None,
                               flags=["no_defining_site_in_vcf"],
                               variant_calls=[_vc_dict(c) for c in calls.values()],
                               reason=f"No {gene} component site found in the VCF (wrong region/assembly?).")

    absent = [c.star for c in calls.values() if not c.found]
    if absent:
        flags.append(f"assumed_reference_at_uncalled_sites:{','.join(absent)}")
    nocalls = [c.rsid for c in calls.values() if c.no_call]
    if nocalls:
        flags.append(f"no_call_at:{','.join(nocalls)}")

    # tag -> VariantCall (star field carries the component tag)
    by_tag = {c.star: c for c in calls.values()}
    alt_calls = [c for c in calls.values() if c.found and c.alt_count > 0]
    hets = [c for c in alt_calls if c.alt_count == 1]

    phased = bool(hets) and all(c.phased and c.hap is not None for c in hets)
    if not hets:
        phased = True  # only homs/refs -> phasing-agnostic

    def present_on(hap_idx: int) -> frozenset[str]:
        out = set()
        for tag, c in by_tag.items():
            if not c.found or c.alt_count == 0:
                continue
            if c.alt_count == 2:
                out.add(tag)                       # hom -> present on both haplotypes
            elif c.hap is not None and c.hap[hap_idx] == 1:
                out.add(tag)
        return frozenset(out)

    if phased:
        a1 = _haplotype_star(present_on(0), rules, reference_allele)
        a2 = _haplotype_star(present_on(1), rules, reference_allele)
        phasing = "phased"
    else:
        # unphased fallback: parsimonious trans (each het on a separate haplotype); homs on both.
        homset = frozenset(c.star for c in alt_calls if c.alt_count == 2)
        if len(hets) <= 1:
            hapA = homset | frozenset(c.star for c in hets)
            hapB = homset
        else:
            # >=2 unphased hets -> cis/trans ambiguous. Standard call = trans (all hets split one-each is
            # not meaningful for compounds); default to CIS-compound on hapA (the clinically-standard *6/*3A
            # call) and expose the trans alternate.
            hapA = homset | frozenset(c.star for c in hets)   # cis: all components together -> compound
            hapB = homset
            flags.append("unphased_compound_cis_assumption")
        a1 = _haplotype_star(hapA, rules, reference_allele)
        a2 = _haplotype_star(hapB, rules, reference_allele)
        phasing = "unphased"

    a1, a2 = sorted((a1, a2), key=_star_sort_key)
    pheno = phenotype_fn(a1, a2)
    proxy = f"{a1}/{a2}"
    res = DiplotypeResult("ok", proxy, a1, a2, pheno, phasing,
                          core_proxy_diplotype=proxy, flags=flags,
                          variant_calls=[_vc_dict(c) for c in calls.values()])

    # phase-ambiguity surfacing: >=2 unphased hets whose cis vs trans resolutions differ in phenotype.
    if phasing == "unphased" and len(hets) >= 2:
        # trans alternate: each het component on its own haplotype (singles), homs on both
        homset = frozenset(c.star for c in alt_calls if c.alt_count == 2)
        t1 = _haplotype_star(homset | {hets[0].star}, rules, reference_allele)
        t2 = _haplotype_star(homset | frozenset(c.star for c in hets[1:]), rules, reference_allele)
        t1, t2 = sorted((t1, t2), key=_star_sort_key)
        alt_pheno = phenotype_fn(t1, t2)
        if alt_pheno != pheno:
            res.phenotype_status = "phase_ambiguous"
            res.phenotype_confidence = "low"
            res.alternate_diplotype = f"{t1}/{t2}"
            res.alternate_phenotype = alt_pheno
            res.flags.append("phase_ambiguous_phenotype_differs")
    return res


def _star_sort_key(a: str):
    core = a.lstrip("*")
    num = ""
    for ch in core:
        if ch.isdigit():
            num += ch
        else:
            break
    return (int(num) if num else 999, a)
