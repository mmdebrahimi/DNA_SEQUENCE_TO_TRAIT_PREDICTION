"""TB AMR decoder cell — RIF (rpoB) + INH (katG + inhA/...) determinant scoring (NON-FROZEN).

Organism-routed, deterministic, mirrors the TMP-SMX overlay branding. Ties together:
  - Step 1 (`tb_vcf`): masked-VCF determinant CALLS + regeno-VCF CALLABILITY, and
  - Step 2 (`tb_who_catalogue`): the WHO grade-1/2 determinant table (ratified A = all grade-1/2 loci).

Rule per drug:
  - R  iff the isolate's masked calls match >=1 grade-1/2 determinant (pos + ref + alt).
  - else, if a regeno VCF is supplied, S iff EVERY in-scope determinant position is callable; ABSTAIN if
    any is uncallable (brainstorm C3 — never susceptible-by-absence over a masked/uncallable window).
  - else (no regeno supplied — fixture/plumbing mode) S with `callability_assessed=False`.

This is the v1a plumbing cell: it makes per-isolate calls. It does NOT compute cohort sens/spec
(that is Step 5, lineage-collapsed + cohort/callability-gated).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from dna_decode.data.tb_who_catalogue import Determinant
from dna_decode.organism_rules import tb_vcf

# Branding — every emitted call carries these so a TB cell is never confused with a frozen claim.
RULE_STATUS = "KNOWLEDGE_BASELINE"
RULE_SCOPE = "organism_routed"
INPUT_TYPE = "vcf_h37rv"
CATALOGUE_COMMIT = "0bb3914348c5a4c981859601447834c08f03ee3d"

R, S, ABSTAIN = "R", "S", "ABSTAIN"


@dataclass(frozen=True)
class DrugCall:
    drug: str
    prediction: str                       # R / S / ABSTAIN
    matched: tuple[Determinant, ...] = ()
    coverage_scope: tuple[str, ...] = ()   # genes whose grade-1/2 determinants are in scope
    n_determinant_positions: int = 0
    n_uncallable_positions: int = 0
    callability_assessed: bool = True
    rule_status: str = RULE_STATUS
    rule_scope: str = RULE_SCOPE
    input_type: str = INPUT_TYPE
    catalogue_commit: str = CATALOGUE_COMMIT


def _isolate_snvs(masked_calls: dict[int, tb_vcf.VariantCall]) -> set[tuple[int, str, str]]:
    """All per-base SNV components carried by the isolate (MNV records decomposed)."""
    out: set[tuple[int, str, str]] = set()
    for c in masked_calls.values():
        out |= tb_vcf.snv_components(c.pos, c.ref, c.alt)
    return out


def _matched_determinants(masked_calls, determinants) -> tuple[Determinant, ...]:
    """A determinant matches iff ALL its SNV components are present in the isolate's decomposed SNVs.

    Single-SNV determinants (the common case) reduce to "isolate carries that exact SNV" — but now the
    isolate side is MNV-decomposed, so a determinant encoded inside a larger multi-base isolate record
    still fires. Multi-base determinant encodings require all their components present (no partial match).
    """
    iso = _isolate_snvs(masked_calls)
    return tuple(d for d in determinants if tb_vcf.snv_components(d.pos, d.ref, d.alt) <= iso)


def score_drug(
    drug: str,
    masked_calls: dict[int, tb_vcf.VariantCall],
    determinants: list[Determinant],
    regeno_text: str | None = None,
    *,
    absent_is_uncallable: bool = True,
) -> DrugCall:
    """Score one drug. A matched determinant short-circuits to R, so callability only affects S calls.

    `absent_is_uncallable` (default True, behaviour UNCHANGED) forwards to
    `tb_vcf.callable_positions`. Pass False for minos panel-VCFs, where a determinant position ABSENT from
    the regeno VCF means "not a genotyping target", not "could not be called" -- under the default rule
    100% of S-by-absence calls ABSTAIN (measured; see `tb_vcf.position_states`).
    """
    matched = _matched_determinants(masked_calls, determinants)
    genes = tuple(sorted({d.gene for d in determinants}))
    positions = sorted({d.pos for d in determinants})

    if matched:
        return DrugCall(drug, R, matched=matched, coverage_scope=genes,
                        n_determinant_positions=len(positions))

    if regeno_text is None:
        return DrugCall(drug, S, coverage_scope=genes, n_determinant_positions=len(positions),
                        callability_assessed=False)

    flags = tb_vcf.callable_positions(regeno_text, positions,
                                      absent_is_uncallable=absent_is_uncallable)
    uncallable = [p for p, ok in flags.items() if not ok]
    if uncallable:
        return DrugCall(drug, ABSTAIN, coverage_scope=genes, n_determinant_positions=len(positions),
                        n_uncallable_positions=len(uncallable))
    return DrugCall(drug, S, coverage_scope=genes, n_determinant_positions=len(positions))


def score_rif(masked_calls, determinants, regeno_text=None, *, absent_is_uncallable=True) -> DrugCall:
    return score_drug("rifampicin", masked_calls, determinants, regeno_text,
                      absent_is_uncallable=absent_is_uncallable)


def score_inh(masked_calls, determinants, regeno_text=None, *, absent_is_uncallable=True) -> DrugCall:
    return score_drug("isoniazid", masked_calls, determinants, regeno_text,
                      absent_is_uncallable=absent_is_uncallable)
