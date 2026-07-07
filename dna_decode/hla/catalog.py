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

NOT a clinical tool. VALIDATION OUTCOME (2026-07-06, vs the free 1000G HLA truth `20140702_hla_diversity`):
only **B*57:01/abacavir (rs2395029)** cleared deployment — sens 0.979 / spec 0.992 / PPV 0.855 (n=1103).
The provisional **B*58:01 (rs9263726)** measured a WEAK sens 0.61 / PPV 0.18 (mixed-population LD) and
**A*31:01 (rs1061235)** is NOT paneled on 1000G (sens 0.0) — both DEMOTED to a documented negative
(`_UNVALIDATED_TAGS`), NOT shipped as routable cells. The single-SNP LD-proxy approach that works for
B*57:01 does NOT generalize (`wiki/hla_validation_2026-07-06.md`).
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
                   "PREDICT-1/Mallal 2008); the deployed clinical abacavir screen. VALIDATED vs 1000G HLA "
                   "truth 2026-07-06: sens 0.979 / spec 0.992 / PPV 0.855 (n=1103, 47TP/1FN/8FP).",
        source="NCBI dbSNP (GRCh38) + CPIC abacavir guideline (Martin 2014); AF-confirmed on 1000G"),
}

# DEMOTED — provisional tags that FAILED real 1000G-HLA-truth validation (2026-07-06); NOT shipped as
# routable cells (would mislead a clinical screen). Kept as a documented negative + the measured numbers.
# The B*58:01/A*31:01 lesson: a single-SNP LD proxy that works for B*57:01 does NOT generalize — these need
# either a population-specific tag (B*58:01: rs9263726 LD is strong only in East Asians) or sequence-based
# typing (A*31:01: no clean 1000G-paneled single-SNP tag). See wiki/hla_validation_2026-07-06.md.
_UNVALIDATED_TAGS: dict[str, dict] = {
    "b5801": {"allele": "HLA-B*58:01", "rsid": "rs9263726", "drug": "allopurinol",
              "measured": "sens 0.609 / spec 0.824 / PPV 0.176 (n=1103) — WEAK proxy (mixed-population LD); "
                          "misses 39% of carriers -> unsafe for an SJS/TEN screen",
              "verdict": "weak_proxy_measured_not_deployable"},
    "a3101": {"allele": "HLA-A*31:01", "rsid": "rs1061235", "drug": "carbamazepine",
              "measured": "sens 0.0 (0TP/74FN) — rs1061235 is NOT paneled on 1000G (no variant record at "
                          "chr6:29945521) -> the tag cannot call carriers",
              "verdict": "no_valid_tag_on_1000g"},
}


def get(key: str) -> HLATagAllele:
    if key not in CATALOG:
        raise KeyError(f"unknown HLA allele key {key!r} (known: {', '.join(CATALOG)})")
    return CATALOG[key]
