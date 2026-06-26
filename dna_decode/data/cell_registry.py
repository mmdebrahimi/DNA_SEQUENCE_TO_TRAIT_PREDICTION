"""Evidence-Contract Registry (v0, thin cut) — one checked-in, test-enforced contract per shipped cell.

A `CellContract` declares, for each deployed decoder cell, WHAT it claims, at what HONEST evidence tier,
on what validation SLICE, with what label provenance, what abstention vocabulary it speaks, and (declared,
not executed) its falsifier / incoming-data gate / demotion rule. The validation report card can later read
its grid + tiers from here so a shipped decoder cannot ship invisibly and abstention has ONE vocabulary.

v0 SCOPE (descoped per the 2026-06-26 /brainstorm + user direction): the VALIDATED/DEPLOYED CORE only —
the AMR surface (projected verbatim from the frozen `shipped_decoder_surface`) + the 3 human-PGx genes.
DEFERRED to v0.1: a `finder`/`typing`/`viral` track + a `route` field + report-card rebuild.

INTEGRITY RAILS (load-bearing):
- `cell_id` is a DISPLAY string ONLY (brainstorm C2). The AMR join key is `cell_key.canonical_cell_key`
  (organism, drug) — `amr_projection_keys()` returns exactly that set, and the consistency test asserts it
  EQUALS the frozen surface's keys. Never join AMR cells by raw `cell_id` string.
- AMR contracts are a PROJECTION of `shipped_decoder_surface.shipped_decoder_rows()` built programmatically,
  so the projection == surface BY CONSTRUCTION; the test pins it against silent drift.
- NO numeric confidence field exists anywhere (anti-"trust-layer-theater" guardrail). `evidence_tier` is a
  categorical honesty label; `claim_status` carries the structural status separately (brainstorm M1).
- This module imports `shipped_decoder_surface` READ-ONLY and touches NO frozen file (`amr_rules.py` /
  `calibrated_amr_rules.json` unchanged) -> `tests/test_tb_leak_guard.py` stays green.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from dna_decode.data.cell_key import canonical_cell_key
from dna_decode.data.cell_registry_vocab import AbstentionVocab
from dna_decode.data.shipped_decoder_surface import shipped_decoder_rows


class EvidenceTier(str, Enum):
    """Honest evidence tier per cell. Categorical, NOT numeric (guardrail). Ordered conceptually
    strong->weak for reading, but the enum carries no arithmetic and tests assert no numeric scale."""

    INDEPENDENT_MEASURED = "independent_measured"   # free INDEPENDENT isolate-level wet-lab label (e.g. HIV PhenoSense)
    NEAR_INDEPENDENT = "near_independent"           # provenance-disjoint stress test / consensus panel (NCBI-PD, GeT-RM)
    FAITHFUL_TO_TOOL = "faithful_to_tool"           # faithful to a reference tool/guideline, not an independent label
    KNOWLEDGE_BASELINE = "knowledge_baseline"       # literature/catalogue assignment, in-distribution
    NO_FREE_SOURCE = "no_free_source"               # no free isolate-level phenotype source exists


@dataclass(frozen=True)
class CellContract:
    """One shipped decoder cell's evidence contract. Frozen; NO numeric confidence field by design."""

    cell_id: str               # DISPLAY string only ("track:organism:target"); NOT the join key
    track: str                 # "amr" | "pgx"  (v0; finder/typing/viral deferred to v0.1)
    organism: str
    target: str                # drug (amr) | gene (pgx)
    claim: str                 # one-line plain claim the cell makes
    evidence_tier: EvidenceTier
    claim_status: str          # structural status (phenotype_source_status for amr; calling-status for pgx) — M1 split
    validation_slice: str      # the slice the tier was earned on
    label_provenance: str      # where the labels came from
    abstention_vocab: AbstentionVocab  # this cell's abstention KIND, collapsed to the controlled vocab
    native_abstention: str     # the cell's own raw in-tree abstention term
    falsifier_ref: str         # path to a falsifier script, or "none" (DECLARED, not executed)
    incoming_data_gate: str    # which of the 8 rejection gates apply, or "n/a" (DECLARED)
    demotion_rule: str         # free-text v0: the trigger that would demote this cell's tier


# --- AMR phenotype_source_status -> (evidence_tier, abstention_vocab) mapping (the only judgment in the AMR
#     projection; everything else is verbatim from the frozen surface). ---
_AMR_STATUS_MAP: dict[str, tuple[EvidenceTier, AbstentionVocab]] = {
    "ncbi_pd":          (EvidenceTier.NEAR_INDEPENDENT, AbstentionVocab.SCORED),
    "label_confounded": (EvidenceTier.FAITHFUL_TO_TOOL, AbstentionVocab.LABEL_CONFOUNDED),
    "no_free_source":   (EvidenceTier.NO_FREE_SOURCE,   AbstentionVocab.NO_FREE_SOURCE),
}


def _amr_contracts() -> list[CellContract]:
    """Project every frozen `shipped_decoder_surface` row to an AMR CellContract (== surface by construction)."""
    out: list[CellContract] = []
    for r in shipped_decoder_rows():
        org, drug = r["organism"], r["drug"]
        status = r["phenotype_source_status"]
        tier, vocab = _AMR_STATUS_MAP[status]
        scoreable = status == "ncbi_pd"
        out.append(CellContract(
            cell_id=f"amr:{org}:{drug}",
            track="amr",
            organism=org,
            target=drug,
            claim=f"{r['engine']} R/S call for {org} x {drug}",
            evidence_tier=tier,
            claim_status=status,
            validation_slice=("NCBI-PD provenance-disjoint stress test (lineage-disclosed)" if scoreable
                              else "label-confounded surrogate (cefoxitin is the CLSI surrogate)"
                              if status == "label_confounded"
                              else "no free isolate-level phenotype source"),
            label_provenance=("NCBI Pathogen Detection AST_phenotypes" if scoreable else "none (structural non-cell)"),
            abstention_vocab=vocab,
            native_abstention=("SCORED" if scoreable
                               else "LABEL_CONFOUNDED" if status == "label_confounded"
                               else "NO_FREE_PHENOTYPE"),
            falsifier_ref="scripts/provenance_disjoint_validate.py" if scoreable else "none",
            incoming_data_gate=("G1,G7,G8" if scoreable else "n/a"),
            demotion_rule=("SCORED -> UNDERPOWERED if either class falls below the powering floor; "
                           "lineage-collapse can demote the disclosed metric" if scoreable
                           else "n/a (no free label to demote against)"),
        ))
    return out


# --- PGx cells: the 3 deployed human-PGx genes (dna-pgx --gene). Calling is independently validatable vs
#     GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC. ---
_PGX_CONTRACTS: list[CellContract] = [
    CellContract(
        cell_id="pgx:human:cyp2c19", track="pgx", organism="human", target="cyp2c19",
        claim="CYP2C19 star-allele diplotype + CPIC metabolizer phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT,
        claim_status="calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM consensus core-diplotype concordance on real 1000G + trio Mendelian QC",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) ⋈ 1000G",
        abstention_vocab=AbstentionVocab.WITHHELD_NONCORE, native_abstention="phenotype_withheld",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="a non-core *4/*35 sentinel hit -> phenotype withheld rather than mis-called",
    ),
    CellContract(
        cell_id="pgx:human:cyp2c9", track="pgx", organism="human", target="cyp2c9",
        claim="CYP2C9 star-allele diplotype + CPIC activity-score phenotype from a phased VCF",
        evidence_tier=EvidenceTier.NEAR_INDEPENDENT,
        claim_status="calling_validated_phenotype_faithful_to_cpic",
        validation_slice="GeT-RM consensus core-diplotype concordance 73/73 on real 1000G + trio Mendelian QC",
        label_provenance="GeT-RM consensus (Astrolabe+Stargazer+Aldy; Gaedigk 2022) ⋈ 1000G",
        abstention_vocab=AbstentionVocab.WITHHELD_NONCORE, native_abstention="phenotype_withheld",
        falsifier_ref="scripts/pgx_getrm_concordance.py", incoming_data_gate="n/a",
        demotion_rule="a non-core *5/*8/*9/*11 sentinel hit -> phenotype withheld",
    ),
    CellContract(
        cell_id="pgx:human:vkorc1", track="pgx", organism="human", target="vkorc1",
        claim="VKORC1 -1639G>A (rs9923231) warfarin-sensitivity genotype from a phased VCF",
        evidence_tier=EvidenceTier.KNOWLEDGE_BASELINE,
        claim_status="single_snp_genotype_to_sensitivity",
        validation_slice="direct genotype readout (minus-strand encoded); not a star/diplotype system",
        label_provenance="literature sensitivity assignment (no independent panel validation in-repo)",
        abstention_vocab=AbstentionVocab.ABSTAIN_BY_DESIGN, native_abstention="ABSTAIN",
        falsifier_ref="none", incoming_data_gate="n/a",
        demotion_rule="n/a (deterministic single-SNP readout)",
    ),
]


def cells() -> list[CellContract]:
    """Every v0 cell contract (AMR projection + PGx)."""
    return _amr_contracts() + list(_PGX_CONTRACTS)


def by_cell_id() -> dict[str, CellContract]:
    """cell_id (display) -> contract. cell_id is unique within v0 by construction."""
    return {c.cell_id: c for c in cells()}


def amr_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "amr"]


def pgx_cells() -> list[CellContract]:
    return [c for c in cells() if c.track == "pgx"]


def amr_projection_keys() -> set[tuple[str, str]]:
    """The AMR cells' canonical (organism, drug) join keys — for the surface-consistency test (NOT cell_id)."""
    return {canonical_cell_key(c.organism, c.target) for c in amr_cells()}


def cli_routable_cell_ids() -> set[str]:
    """The v0 CLI-routable cell-id set: AMR surface cells + the 3 PGx genes."""
    return set(by_cell_id())
