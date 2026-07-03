"""Deterministic single-locus human-trait decoders (openSNP-scoreable, self-report tier).

Extends the visible/organismal Mendelian decoder set (eye colour rs12913832, ABO rs8176719) to more
single-SNP traits. Each rule is a SOURCED, textbook genotype->phenotype call. Label tier = self-reported
(non-circular, noisy) -> PILOT/DEMO, exactly like the eye-colour openSNP cell.

STRAND-AGNOSTIC design (openSNP mixes 23andMe forward + AncestryDNA strands): every SNP here is a
NON-complementary pair (C/T, A/G, A/C — never A/T or C/G), so the forward-strand allele set and its
reverse complement are DISJOINT and the strand is recoverable from the genotype itself (like the
strand-agnostic eye-colour cell). A genotype whose alleles fall outside both strand representations ->
INDETERMINATE (never guessed).

Sources (verbatim, per rule):
  * Earwax  — ABCC11 rs17822931, Yoshiura 2006 Nat Genet 38:324. c.538G>A; the derived A (=T on the
    23andMe C/T strand) allele is non-functional -> DRY earwax is HOMOZYGOUS derived (AA/TT); any ancestral
    G/C allele -> WET (dominant). Near-Mendelian (the strongest single-SNP human trait).
  * Lactase — MCM6/LCT rs4988235 (-13910C>T), Enattah 2002 Nat Genet 30:233. The derived T (=A on the
    23andMe A/G strand) allele confers lactase PERSISTENCE (dominant) -> a T/A carrier is lactose-TOLERANT;
    homozygous ancestral (CC/GG) -> non-persistent (lactose-INTOLERANT). European-calibrated (ancestry
    caveat, like eye colour); other persistence alleles exist in African/Arabian populations (v0 misses them).
  * Cilantro — near OR6A2 rs72921001, Eriksson 2012 (23andMe, Flavour 1:22). The C (=G on reverse) allele
    is associated with perceiving cilantro as SOAPY. DELIBERATELY a WEAK-EFFECT association (small OR), NOT a
    Mendelian rule -> included as an HONEST calibration contrast: the deterministic decoder should score
    NEAR-CHANCE here, distinguishing a real single-locus trait from a weak GWAS hit.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

INDETERMINATE = "INDETERMINATE"


def _alleles(genotype: str) -> set[str]:
    return {c for c in (genotype or "").strip().upper() if c in "ACGT"}


# --- Deterministic calls (strand-agnostic; homozygous-derived vs dominant-ancestral per source) ---

def call_earwax(rs17822931: str) -> str:
    """DRY iff homozygous derived (AA or TT); WET iff any ancestral G/C allele. Yoshiura 2006."""
    a = _alleles(rs17822931)
    if not a:
        return INDETERMINATE
    if a <= {"A"} or a <= {"T"}:          # homozygous derived (dry allele) on either strand
        return "dry"
    if a & {"G", "C"}:                    # >=1 ancestral (wet) allele -> dominant WET
        return "wet"
    return INDETERMINATE


def call_lactase(rs4988235: str) -> str:
    """PERSISTENT (lactose-tolerant) iff a derived A/T allele present; non-persistent (INTOLERANT) iff
    homozygous ancestral GG/CC. Enattah 2002; European -13910*T. Strand-agnostic."""
    a = _alleles(rs4988235)
    if not a:
        return INDETERMINATE
    if a & {"A", "T"}:                    # dominant persistence allele
        return "tolerant"
    if a <= {"G"} or a <= {"C"}:          # homozygous ancestral -> non-persistent
        return "intolerant"
    return INDETERMINATE


def call_cilantro(rs72921001: str) -> str:
    """SOAPY iff a C/G risk allele present (weak association); not-soapy iff homozygous A/T. Eriksson 2012.
    WEAK EFFECT by design -> expected near-chance (the honest calibration contrast)."""
    a = _alleles(rs72921001)
    if not a:
        return INDETERMINATE
    if a & {"C", "G"}:
        return "soapy"
    if a <= {"A"} or a <= {"T"}:
        return "not-soapy"
    return INDETERMINATE


# --- Self-report binners (raw openSNP free-text -> the binary label, or None = unscoreable) ---

def bin_earwax(raw: str) -> str | None:
    s = (raw or "").strip().lower()
    if not s or s in ("-", "unknown", "n/a", "rather not say"):
        return None
    if "dry" in s or "flaky" in s or "crumbl" in s or "brittle" in s:
        return "dry"
    if "wet" in s or "sticky" in s or "moist" in s or "soft" in s:
        return "wet"
    return None


def bin_lactase(raw: str) -> str | None:
    """'Lactose intolerance' free text -> tolerant / intolerant. Note the label POLARITY: this column is
    'intolerance', so 'yes'/'intolerant' -> intolerant; 'no'/'tolerant'/'none' -> tolerant."""
    s = (raw or "").strip().lower()
    if not s or s in ("-", "unknown", "n/a", "rather not say"):
        return None
    if "not intoler" in s or "no intoler" in s:      # 'not intolerant' -> tolerant (check before 'intoler')
        return "tolerant"
    if "intoler" in s or s in ("yes", "true", "y"):
        return "intolerant"
    if "toler" in s or "persist" in s or s in ("no", "false", "none", "n"):
        return "tolerant"
    return None


def bin_cilantro(raw: str) -> str | None:
    s = (raw or "").strip().lower()
    if not s or s in ("-", "unknown", "n/a", "rather not say"):
        return None
    if "soap" in s or s in ("yes", "true", "y"):
        return "soapy"
    if s in ("no", "false", "n") or "not" in s or "normal" in s or "fine" in s or "good" in s:
        return "not-soapy"
    return None


@dataclass(frozen=True)
class SingleSnpTrait:
    key: str
    rsid: str
    gene: str
    phenotype_keywords: tuple[str, ...]   # match the openSNP phenotype column header (lowercased contains)
    positive_label: str                   # the "R-like" positive class for the confusion matrix
    negative_label: str
    call: Callable[[str], str]
    binner: Callable[[str], "str | None"]
    tier: str                             # STRONG_MENDELIAN / WEAK_ASSOCIATION_CONTRAST
    source: str


TRAITS: dict[str, SingleSnpTrait] = {
    "earwax": SingleSnpTrait(
        key="earwax", rsid="rs17822931", gene="ABCC11", phenotype_keywords=("earwax",),
        positive_label="dry", negative_label="wet", call=call_earwax, binner=bin_earwax,
        tier="STRONG_MENDELIAN", source="Yoshiura 2006 Nat Genet 38:324"),
    "lactase": SingleSnpTrait(
        key="lactase", rsid="rs4988235", gene="MCM6/LCT", phenotype_keywords=("lactose", "lactase"),
        positive_label="intolerant", negative_label="tolerant", call=call_lactase, binner=bin_lactase,
        tier="STRONG_MENDELIAN", source="Enattah 2002 Nat Genet 30:233 (European -13910*T)"),
    "cilantro": SingleSnpTrait(
        key="cilantro", rsid="rs72921001", gene="OR6A2-region", phenotype_keywords=("cilantro", "coriander"),
        positive_label="soapy", negative_label="not-soapy", call=call_cilantro, binner=bin_cilantro,
        tier="WEAK_ASSOCIATION_CONTRAST", source="Eriksson 2012 Flavour 1:22 (23andMe GWAS, weak OR)"),
}
