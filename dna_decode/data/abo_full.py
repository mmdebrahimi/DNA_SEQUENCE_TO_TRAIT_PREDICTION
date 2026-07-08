"""Full ABO blood-type decoder (O/A/B/AB) — the 3-variant deterministic cell.

Completes the named-deferred A-vs-B extension of the O-vs-non-O ABO cell (wiki/abo_pgp_result_2026-06-30
flagged "A-vs-B ... NOT fabricated -> deferred"). Uses the STANDARD 3-variant method (the same one large
genomic cohorts use, e.g. UK Biobank ABO genotyping): the O-deletion + the two A/B-distinguishing exon-7 SNPs.

PROVENANCE (sourced, NOT memory — ClinVar RCV000019310 + the ABO 3-variant literature, 2026-07-07):
  * rs8176719  261delG (c.261delG) — the O frameshift deletion. 23andMe reports D (deletion) / I (insertion).
  * rs8176746  c.796C>A  p.Leu266Met — A vs B residue 266.
  * rs8176747  c.803G>C  p.Gly268Ala — A vs B residue 268.
  ABO is transcribed on the MINUS strand, so on a plus-strand array (23andMe) the cDNA alleles map to:
    rs8176746: A-allele = genomic G (cDNA 796C), B-allele = genomic T (cDNA 796A)
    rs8176747: A-allele = genomic C (cDNA 803G), B-allele = genomic G (cDNA 803C)
  Empirically cross-checked against the OpenSNP substrate (O+ -> GG/CC A-background tags; AB+ -> GT/CG het).

THE O-BACKGROUND SUBTLETY (load-bearing, why the deletion count is needed): the common O allele (O01) sits on
an A-type tag background (G/C), so tag-only A/B calling is ambiguous. A DI genotype (one O + one functional B)
shows HET tags but is phenotype B, NOT AB. The deletion count disambiguates:
  * 2 deletions (DD)            -> O
  * tags homozygous A (GG,CC)   -> A   (AA or AO)
  * tags homozygous B (TT,GG)   -> B   (BB or BO)
  * tags heterozygous, II       -> AB  (one A + one B functional allele)
  * tags heterozygous, DI       -> B   (O carries A-background tags; the functional allele is B)

HONESTY: near-Mendelian at the common alleles; rare/weak subgroup alleles (A2, cis-AB, Bw, non-deletional O)
+ self-report noise (~15%) are NOT captured -> Indeterminate on ambiguous input. NOT a clinical tool.
"""
from __future__ import annotations

# plus-strand genomic (23andMe) A-type vs B-type allele per tag SNP.
_A746, _B746 = "G", "T"   # rs8176746
_A747, _B747 = "C", "G"   # rs8176747


def _alleles(gt: str | None) -> str | None:
    if not gt:
        return None
    g = gt.upper().replace("/", "").replace("|", "").replace(" ", "").replace("-", "")
    return g if len(g) == 2 and g not in ("--", "..") else None


def _deletion_count(rs8176719: str | None) -> int | None:
    """rs8176719 D/I -> count of deletion (O) alleles; None if uncallable."""
    g = _alleles(rs8176719)
    if g is None:
        return None
    if set(g) <= {"D", "I"}:
        return g.count("D")
    return None


def _tag_state(gt: str | None, a_allele: str, b_allele: str) -> str | None:
    """A-hom / B-hom / het / None for a tag SNP genotype."""
    g = _alleles(gt)
    if g is None:
        return None
    s = set(g)
    if s == {a_allele}:
        return "A"
    if s == {b_allele}:
        return "B"
    if s == {a_allele, b_allele}:
        return "het"
    return None


def call_abo_full(rs8176719: str | None, rs8176746: str | None, rs8176747: str | None) -> str:
    """Deterministic O/A/B/AB call from the 3 ABO variants (plus-strand genotypes). Indeterminate on
    uncallable / inconsistent input (never guesses)."""
    ndel = _deletion_count(rs8176719)
    if ndel is None:
        return "Indeterminate"
    if ndel == 2:
        return "O"
    t746 = _tag_state(rs8176746, _A746, _B746)
    t747 = _tag_state(rs8176747, _A747, _B747)
    if t746 is None or t747 is None or t746 != t747:
        return "Indeterminate"          # the two tags must agree (same A/B/het state)
    state = t746
    if state == "A":
        return "A"
    if state == "B":
        return "B"
    # het tags: II -> AB (two functional, one A one B); DI -> B (O carries A-background)
    if ndel == 0:
        return "AB"
    if ndel == 1:
        return "B"
    return "Indeterminate"
