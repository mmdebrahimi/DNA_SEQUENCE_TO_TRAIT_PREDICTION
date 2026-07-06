"""CYP2D6 diplotype caller — PRIORITY-ORDERED per-haplotype star resolution (v0, SNP surface).

Distinct from BOTH the single-SNP core caller (`caller.py`) and the subset-largest compound caller
(`compound_caller.py`). CYP2D6 star alleles share a SNP BACKGROUND (*2's 2851C>T / 4181G>C is carried by
*4/*17/*29/*35/*41; *4 carries *10's 100C>T), so neither a one-SNP-per-star map nor a subset-largest match
resolves them correctly — a shared background SNP would over-match. This resolver picks the MOST-SPECIFIC
defining SNP present on each haplotype via an explicit priority order (`cyp2d6_catalog.STAR_PRIORITY`), so a
background SNP never mis-calls the more-specific allele.

It reuses `caller.scan_vcf` to read the component sites (with phasing). Phased input (1000G 30x is phased)
-> exact per-haplotype resolution; unphased with >=2 het specific sites -> parsimonious trans, flagged.

HONESTY: structural alleles (*5/*13/*36/*61/*63/*68/*xN) are INVISIBLE to a SNP VCF and are NOT withheld
(the proxy cannot see them) -> they may be SILENTLY MIS-CALLED; the runner stamps `cnv_hybrid_unassessed`.
"""
from __future__ import annotations

from dna_decode.pgx.caller import DiplotypeResult, _vc_dict, scan_vcf
from dna_decode.pgx.cyp2c19_catalog import DefiningVariant
from dna_decode.pgx.cyp2d6_catalog import BACKGROUND_TAGS, SPECIFIC_TAGS


def _haplotype_star_priority(present: frozenset[str], priority: list[tuple[str, str]],
                             reference: str) -> str:
    """Given the component tags present on ONE haplotype, return the star of the FIRST (star, tag) in
    `priority` whose tag is present. Priority is ordered most-specific -> background, so a shared
    background SNP never wins over a more-specific allele-defining SNP. Falls back to `reference`."""
    for star, tag in priority:
        if tag in present:
            return star
    return reference


def _star_sort_key(a: str):
    core = a.lstrip("*")
    num = ""
    for ch in core:
        if ch.isdigit():
            num += ch
        else:
            break
    return (int(num) if num else 999, a)


def assemble_cyp2d6_diplotype(
    vcf,
    components: list[DefiningVariant],
    priority: list[tuple[str, str]],
    *,
    reference_allele: str,
    phenotype_fn,
    gene: str,
    sample: str | None = None,
) -> DiplotypeResult:
    """Resolve a CYP2D6 diplotype from a VCF via per-haplotype priority resolution. Returns the same
    DiplotypeResult shape as the core caller so the runner/report wiring is uniform."""
    calls = scan_vcf(vcf, defining=components, sample=sample)
    present_calls = [c for c in calls.values() if c.found]
    flags: list[str] = []
    if not present_calls:
        return DiplotypeResult("no_input", None, None, None, None, None,
                               flags=["no_defining_site_in_vcf"],
                               variant_calls=[_vc_dict(c) for c in calls.values()],
                               reason=f"No {gene} defining site found in the VCF (wrong region/assembly?).")

    absent = [c.star for c in calls.values() if not c.found]
    if absent:
        flags.append(f"assumed_reference_at_uncalled_sites:{','.join(absent)}")
    nocalls = [c.rsid for c in calls.values() if c.no_call]
    if nocalls:
        flags.append(f"no_call_at:{','.join(nocalls)}")

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
                out.add(tag)                      # hom -> present on both haplotypes
            elif c.hap is not None and c.hap[hap_idx] == 1:
                out.add(tag)
        return frozenset(out)

    if phased:
        h0, h1 = present_on(0), present_on(1)
        a1 = _haplotype_star_priority(h0, priority, reference_allele)
        a2 = _haplotype_star_priority(h1, priority, reference_allele)
        phasing = "phased"
        hap_tags = (h0, h1)
    else:
        # unphased fallback: parsimonious trans (each het on its own haplotype); homs on both.
        homset = frozenset(c.star for c in alt_calls if c.alt_count == 2)
        het_tags = [c.star for c in hets]
        if len(hets) <= 1:
            hA = homset | frozenset(het_tags)
            hB = homset
        else:
            hA = homset | {het_tags[0]}
            hB = homset | frozenset(het_tags[1:])
            flags.append("unphased_trans_assumption")
        a1 = _haplotype_star_priority(hA, priority, reference_allele)
        a2 = _haplotype_star_priority(hB, priority, reference_allele)
        phasing = "unphased"
        hap_tags = (hA, hB)

    # Honesty flag (brainstorm): >=2 allele-SPECIFIC tags on ONE haplotype is a data anomaly the priority
    # resolver silently collapses -> surface it rather than hide it.
    for h in hap_tags:
        if len(h & SPECIFIC_TAGS) >= 2:
            flags.append("multi_specific_haplotype:" + "+".join(sorted(h & SPECIFIC_TAGS)))

    a1, a2 = sorted((a1, a2), key=_star_sort_key)
    pheno = phenotype_fn(a1, a2)
    proxy = f"{a1}/{a2}"
    res = DiplotypeResult("ok", proxy, a1, a2, pheno, phasing,
                          core_proxy_diplotype=proxy, flags=flags,
                          variant_calls=[_vc_dict(c) for c in calls.values()])

    # phase-ambiguity surfacing: >=2 unphased hets whose trans vs cis resolutions differ in phenotype.
    if phasing == "unphased" and len(hets) >= 2:
        homset = frozenset(c.star for c in alt_calls if c.alt_count == 2)
        het_tags = [c.star for c in hets]
        c1 = _haplotype_star_priority(homset | frozenset(het_tags), priority, reference_allele)  # cis
        c2 = _haplotype_star_priority(homset, priority, reference_allele)
        c1, c2 = sorted((c1, c2), key=_star_sort_key)
        alt_pheno = phenotype_fn(c1, c2)
        if alt_pheno != pheno:
            res.phenotype_status = "phase_ambiguous"
            res.phenotype_confidence = "low"
            res.alternate_diplotype = f"{c1}/{c2}"
            res.alternate_phenotype = alt_pheno
            res.flags.append("phase_ambiguous_phenotype_differs")

    return res
