"""HLA pharmacogenomic tag-SNP catalog — GRCh38, NCBI-verified + AF-confirmed on 1000G.

Each entry: a clinically-actionable HLA allele + its near-perfect LD TAG SNP + the drug / adverse reaction /
CPIC action. The tag SNP's ALT allele marks CARRIAGE of the HLA allele (a PROXY via linkage disequilibrium;
concordance vs real HLA truth is the deliverable number, see scripts/hla_concordance.py).

PROVENANCE (grounded, NO fabrication):
  * B*57:01 / rs2395029 (HCP5 g.31464003 T>G, GRCh38 verified via NCBI): the G allele tags HLA-B*57:01 with
    ~100% sensitivity + high specificity (Colombo 2008; the PREDICT-1 / Mallal 2008 abacavir screen; the
    clinical rs2395029 abacavir test). CPIC abacavir guideline (Martin 2014): B*57:01 carrier -> AVOID
    abacavir (hypersensitivity reaction, HSR). AF(G) on 1000G = 0.031.
  * B*58:01 / rs9263726 (g.31138722 G>A): the A allele tags HLA-B*58:01 (Asian-common); allopurinol
    SJS/TEN. CPIC allopurinol (Saito 2016): B*58:01 carrier -> avoid allopurinol. Proxy quality is
    population-dependent (strong in East Asians) -> flagged PROVISIONAL.
  * A*31:01 / rs1061235 (g.29945521 A>T): tags HLA-A*31:01 region; carbamazepine HSR/DRESS. CPIC
    carbamazepine (Phillips 2018). Proxy quality less established than B*57:01 -> flagged PROVISIONAL.

NOT a clinical tool. v0 anchor = B*57:01/abacavir (the gold-standard tag); B*58:01 + A*31:01 ship as
PROVISIONAL proxies (validate the LD in the target population before any use).
"""
from __future__ import annotations

from dataclasses import dataclass

ASSEMBLY = "GRCh38"


@dataclass(frozen=True)
class HLATagAllele:
    key: str            # CLI key, e.g. "b5701"
    allele: str         # "HLA-B*57:01"
    rsid: str
    chrom: str          # GRCh38, normalized (no "chr")
    pos: int
    ref: str
    tag_alt: str        # the ALT allele that marks CARRIAGE of the HLA allele
    drug: str
    reaction: str
    cpic_action: str
    proxy_tier: str     # "gold_standard" | "provisional"
    proxy_note: str
    source: str


CATALOG: dict[str, HLATagAllele] = {
    "b5701": HLATagAllele(
        key="b5701", allele="HLA-B*57:01", rsid="rs2395029", chrom="6", pos=31464003, ref="T", tag_alt="G",
        drug="abacavir", reaction="hypersensitivity reaction (HSR)",
        cpic_action="AVOID abacavir in a B*57:01 carrier (CPIC Martin 2014)",
        proxy_tier="gold_standard",
        proxy_note="rs2395029(G) tags B*57:01 with ~100% sensitivity + high specificity (Colombo 2008; "
                   "PREDICT-1/Mallal 2008); the deployed clinical abacavir screen.",
        source="NCBI dbSNP (GRCh38) + CPIC abacavir guideline (Martin 2014); AF-confirmed on 1000G"),
    "b5801": HLATagAllele(
        key="b5801", allele="HLA-B*58:01", rsid="rs9263726", chrom="6", pos=31138722, ref="G", tag_alt="A",
        drug="allopurinol", reaction="SJS/TEN (severe cutaneous)",
        cpic_action="AVOID allopurinol in a B*58:01 carrier (CPIC Saito 2016)",
        proxy_tier="provisional",
        proxy_note="rs9263726(A) tags B*58:01; LD is strong in East Asians, weaker elsewhere -> PROVISIONAL, "
                   "validate in the target population.",
        source="NCBI dbSNP (GRCh38) + CPIC allopurinol guideline (Saito 2016)"),
    "a3101": HLATagAllele(
        key="a3101", allele="HLA-A*31:01", rsid="rs1061235", chrom="6", pos=29945521, ref="A", tag_alt="T",
        drug="carbamazepine", reaction="HSR / DRESS",
        cpic_action="consider avoiding carbamazepine in an A*31:01 carrier (CPIC Phillips 2018)",
        proxy_tier="provisional",
        proxy_note="rs1061235(T) tags the HLA-A*31:01 region; proxy quality less established than B*57:01 -> "
                   "PROVISIONAL.",
        source="NCBI dbSNP (GRCh38) + CPIC carbamazepine guideline (Phillips 2018)"),
}


def get(key: str) -> HLATagAllele:
    if key not in CATALOG:
        raise KeyError(f"unknown HLA allele key {key!r} (known: {', '.join(CATALOG)})")
    return CATALOG[key]
