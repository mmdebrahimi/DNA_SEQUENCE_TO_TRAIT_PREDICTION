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


def _matches(call: tb_vcf.VariantCall, det: Determinant) -> bool:
    return call.ref == det.ref and call.alt == det.alt


def score_drug(
    drug: str,
    masked_calls: dict[int, tb_vcf.VariantCall],
    determinants: list[Determinant],
    regeno_text: str | None = None,
) -> DrugCall:
    matched = tuple(
        d for d in determinants
        if d.pos in masked_calls and _matches(masked_calls[d.pos], d)
    )
    genes = tuple(sorted({d.gene for d in determinants}))
    positions = sorted({d.pos for d in determinants})

    if matched:
        return DrugCall(drug, R, matched=matched, coverage_scope=genes,
                        n_determinant_positions=len(positions))

    if regeno_text is None:
        return DrugCall(drug, S, coverage_scope=genes, n_determinant_positions=len(positions),
                        callability_assessed=False)

    flags = tb_vcf.callable_positions(regeno_text, positions)
    uncallable = [p for p, ok in flags.items() if not ok]
    if uncallable:
        return DrugCall(drug, ABSTAIN, coverage_scope=genes, n_determinant_positions=len(positions),
                        n_uncallable_positions=len(uncallable))
    return DrugCall(drug, S, coverage_scope=genes, n_determinant_positions=len(positions))


def score_rif(masked_calls, determinants, regeno_text=None) -> DrugCall:
    return score_drug("rifampicin", masked_calls, determinants, regeno_text)


def score_inh(masked_calls, determinants, regeno_text=None) -> DrugCall:
    return score_drug("isoniazid", masked_calls, determinants, regeno_text)
