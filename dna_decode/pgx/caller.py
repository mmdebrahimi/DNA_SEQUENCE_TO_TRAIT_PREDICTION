"""CYP2C19 diplotype caller: VCF -> defining-SNP genotypes -> star alleles -> diplotype -> CPIC phenotype.

Pure-stdlib (no pysam): a minimal single-sample VCF reader matching defining sites by (chrom,pos) with an
rsID fallback. Handles phased (`|`) and unphased (`/`) GT, multiallelic ALT, no-call (`./.`), and -- the
documented hard case -- absence of a record for a defining site (population variant-only VCFs only list
ALT-bearing sites). Absence is treated as reference BUT surfaced explicitly (never silent ref-by-absence).

Diplotype assembly:
  * fully-phased het sites -> exact per-haplotype resolution.
  * unphased: parsimonious trans assembly (homs on both haplotypes; <=1 het -> unambiguous; 2 hets ->
    trans assumption, flagged; >2 hets -> flagged ambiguous). The cis/trans ambiguity at >=2 unphased het
    sites is the NA19122-class edge (a known GeT-RM phasing case) -- flagged, not hidden.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from dna_decode.pgx.cyp2c19_catalog import (
    ASSEMBLY,
    CORE_DEFINING,
    GENE,
    REFERENCE_ALLELE,
    DefiningVariant,
    diplotype_phenotype,
)


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


@dataclass
class VariantCall:
    star: str
    rsid: str
    pos: int
    found: bool                 # was a VCF record present for this defining site?
    gt: str | None             # raw GT string
    phased: bool
    alt_count: int             # number of ALT copies of this allele (0/1/2)
    no_call: bool              # GT contained a missing allele (".")
    hap: tuple[int, int] | None  # per-haplotype ALT presence (1=alt) when GT is diploid


@dataclass
class DiplotypeResult:
    status: str                # "ok" | "no_input"
    diplotype: str | None      # "*2/*17"
    allele1: str | None
    allele2: str | None
    phenotype: str | None      # CPIC metabolizer phenotype
    phasing: str | None        # "phased" | "unphased" | None
    flags: list[str] = field(default_factory=list)
    variant_calls: list[dict] = field(default_factory=list)
    reason: str | None = None


def _interpret(d: DefiningVariant, ref: str, alt_field: str, gt: str | None) -> VariantCall:
    """Build a VariantCall from one matched VCF record's REF/ALT/GT for defining variant `d`."""
    alts = alt_field.split(",")
    # allele index for d.alt within [REF, ALT1, ALT2, ...]; -1 if d.alt isn't an ALT here
    alt_index = alts.index(d.alt) + 1 if d.alt in alts else -1
    if gt is None or gt in (".", "./.", ".|."):
        no_call = gt is not None
        return VariantCall(d.star, d.rsid, d.pos, True, gt, "|" in (gt or ""),
                           0, no_call, None)
    phased = "|" in gt
    raw = gt.replace("|", "/").split("/")
    no_call = any(a == "." for a in raw)
    nums = [int(a) for a in raw if a.isdigit()]
    if alt_index < 0:
        # d.alt not represented at this record (different ALT) -> treat as reference for d
        return VariantCall(d.star, d.rsid, d.pos, True, gt, phased, 0, no_call, None)
    alt_count = sum(1 for n in nums if n == alt_index)
    hap = None
    if len(raw) == 2 and all(a.isdigit() for a in raw):
        hap = (1 if int(raw[0]) == alt_index else 0, 1 if int(raw[1]) == alt_index else 0)
    return VariantCall(d.star, d.rsid, d.pos, True, gt, phased, alt_count, no_call, hap)


def scan_vcf(path: str | Path, defining: list[DefiningVariant] = CORE_DEFINING,
             sample: str | None = None) -> dict[str, VariantCall]:
    """Parse a single-sample VCF; return {star_allele: VariantCall} for every defining site.

    Sites with no record in the VCF are returned with found=False (assumed reference downstream, flagged).
    `sample` selects a column by name; default = the first sample column.
    """
    by_pos = {(d.chrom, d.pos): d for d in defining}
    by_rsid = {d.rsid: d for d in defining}
    calls: dict[str, VariantCall] = {}
    sample_idx = 0
    text = Path(path).read_text(encoding="utf-8")
    for line in text.splitlines():
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.split("\t")
            samples = header[9:] if len(header) > 9 else []
            if sample and sample in samples:
                sample_idx = samples.index(sample)
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 8:
            continue
        chrom, pos_s, vid, ref, alt = _norm_chrom(cols[0]), cols[1], cols[2], cols[3], cols[4]
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        d = by_pos.get((chrom, pos))
        if d is None and vid and vid != ".":
            d = by_rsid.get(vid)
        if d is None:
            continue
        gt = None
        if len(cols) >= 10:
            fmt = cols[8].split(":")
            col = 9 + sample_idx
            if "GT" in fmt and col < len(cols):
                parts = cols[col].split(":")
                gi = fmt.index("GT")
                if gi < len(parts):
                    gt = parts[gi]
        calls[d.star] = _interpret(d, ref, alt, gt)
    # fill absent sites
    for d in defining:
        if d.star not in calls:
            calls[d.star] = VariantCall(d.star, d.rsid, d.pos, False, None, False, 0, False, None)
    return calls


def _combine_haplotype(stars: list[str], flags: list[str]) -> str:
    """Reduce a haplotype's ALT-bearing star list to a single allele label."""
    if not stars:
        return REFERENCE_ALLELE
    if len(stars) == 1:
        return stars[0]
    label = "+".join(sorted(set(stars)))
    if "compound_haplotype_outside_core" not in flags:
        flags.append("compound_haplotype_outside_core")
    return label


def assemble_diplotype(calls: dict[str, VariantCall]) -> DiplotypeResult:
    """Combine per-site VariantCalls into a diplotype + CPIC phenotype + honesty flags."""
    flags: list[str] = []
    present = [c for c in calls.values() if c.found]
    if not present:
        return DiplotypeResult("no_input", None, None, None, None, None,
                               flags=["no_defining_site_in_vcf"],
                               variant_calls=[_vc_dict(c) for c in calls.values()],
                               reason="No CYP2C19 defining site found in the VCF (wrong region/assembly?).")

    absent = [c.star for c in calls.values() if not c.found]
    if absent:
        flags.append(f"assumed_reference_at_uncalled_sites:{','.join(absent)}")
    nocalls = [c.rsid for c in calls.values() if c.no_call]
    if nocalls:
        flags.append(f"no_call_at:{','.join(nocalls)}")

    alts = [c for c in calls.values() if c.found and c.alt_count > 0]
    homs = [c for c in alts if c.alt_count == 2]
    hets = [c for c in alts if c.alt_count == 1]

    phasing: str
    if hets and all(c.phased and c.hap is not None for c in hets):
        # exact phased resolution (homs are phasing-agnostic; their phased GT 1|1 -> hap (1,1))
        hap0 = [c.star for c in alts if c.hap and c.hap[0] == 1]
        hap1 = [c.star for c in alts if c.hap and c.hap[1] == 1]
        a1, a2 = _combine_haplotype(hap0, flags), _combine_haplotype(hap1, flags)
        phasing = "phased"
    else:
        hap_a = [c.star for c in homs]
        hap_b = [c.star for c in homs]
        if len(hets) == 1:
            hap_a.append(hets[0].star)
        elif len(hets) == 2:
            hap_a.append(hets[0].star)
            hap_b.append(hets[1].star)
            flags.append("unphased_trans_assumption")
        elif len(hets) > 2:
            flags.append("ambiguous_multi_het_unphased")
            for i, c in enumerate(hets):
                (hap_a if i % 2 == 0 else hap_b).append(c.star)
        a1, a2 = _combine_haplotype(hap_a, flags), _combine_haplotype(hap_b, flags)
        phasing = "unphased"

    a1, a2 = sorted((a1, a2), key=_allele_sort_key)
    pheno = diplotype_phenotype(a1, a2)
    return DiplotypeResult("ok", f"{a1}/{a2}", a1, a2, pheno, phasing, flags=flags,
                           variant_calls=[_vc_dict(c) for c in calls.values()])


def _allele_sort_key(a: str):
    """Sort *1 first, then by numeric star where possible, else lexicographic."""
    core = a.lstrip("*").split("+")[0]
    try:
        return (0, int(core), a)
    except ValueError:
        return (1, 0, a)


def _vc_dict(c: VariantCall) -> dict:
    return {"star": c.star, "rsid": c.rsid, "pos": c.pos, "found": c.found,
            "gt": c.gt, "phased": c.phased, "alt_count": c.alt_count, "no_call": c.no_call}


def call_diplotype(path: str | Path, sample: str | None = None,
                   defining: list[DefiningVariant] = CORE_DEFINING) -> DiplotypeResult:
    """End-to-end: VCF -> CYP2C19 diplotype + CPIC phenotype. Convenience wrapper over scan + assemble."""
    return assemble_diplotype(scan_vcf(path, defining=defining, sample=sample))
