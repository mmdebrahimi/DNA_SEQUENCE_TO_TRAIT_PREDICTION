"""Evidence-Contract Registry v0 — the coverage + consistency contract IS the deliverable.

Pins (thin v0 scope = AMR surface + 3 PGx genes):
  - COVERAGE: every v0 CLI-routable cell (each frozen-surface drug + each dna-pgx --gene) has exactly one contract.
  - CONSISTENCY: the AMR projection EQUALS the frozen shipped_decoder_surface, joined via canonical_cell_key
    (NOT raw cell_id string — brainstorm C2).
  - HONESTY GUARDRAILS: no numeric confidence field anywhere; tiers + abstention are categorical enums.
  - LEAK-GUARD: the registry imports the frozen surface read-only (asserted indirectly via the separate
    tests/test_tb_leak_guard.py, run alongside in the MVP bar).
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

from dna_decode.data.cell_key import canonical_cell_key
from dna_decode.data.cell_registry import (
    CellContract,
    EvidenceTier,
    amr_cells,
    amr_projection_keys,
    by_cell_id,
    cells,
    cli_routable_cell_ids,
    pgx_cells,
)
from dna_decode.data.cell_registry_vocab import (
    NATIVE_TO_VOCAB,
    AbstentionVocab,
    to_vocab,
)
from dna_decode.data.shipped_decoder_surface import all_surface_drugs, surface_index

REPO = Path(__file__).resolve().parent.parent

# --- PGx CLI-routable genes (dna-pgx --gene choices) — the authoritative v0 PGx routable set. ---
_PGX_GENES = {"cyp2c19", "cyp2c9", "vkorc1"}


# ------------------------- coverage -------------------------

def test_cell_id_unique():
    ids = [c.cell_id for c in cells()]
    assert len(ids) == len(set(ids)), "duplicate cell_id in the registry"


def test_every_surface_drug_has_an_amr_contract():
    """COVERAGE: no shipped AMR drug is missing from the registry."""
    covered = {c.target.lower() for c in amr_cells()}
    missing = all_surface_drugs() - covered
    assert not missing, f"AMR drugs in the frozen surface but missing a contract: {sorted(missing)}"


def test_every_pgx_gene_has_exactly_one_contract():
    """COVERAGE: each dna-pgx --gene has exactly one contract."""
    pgx_targets = [c.target.lower() for c in pgx_cells()]
    assert set(pgx_targets) == _PGX_GENES, f"PGx gene coverage mismatch: {sorted(pgx_targets)}"
    assert len(pgx_targets) == len(set(pgx_targets)), "duplicate PGx gene contract"


def test_cli_routable_set_is_amr_plus_pgx():
    """The v0 routable cell-id set is exactly the AMR surface cells + the 3 PGx genes."""
    routable = cli_routable_cell_ids()
    assert routable == set(by_cell_id())
    assert len([c for c in cells() if c.track == "pgx"]) == len(_PGX_GENES)


# ------------------------- consistency (the load-bearing join) -------------------------

def test_amr_projection_equals_frozen_surface_via_canonical_key():
    """CONSISTENCY: registry AMR projection == frozen surface, joined via canonical_cell_key (NOT cell_id)."""
    assert amr_projection_keys() == set(surface_index().keys())


def test_amr_projection_is_one_to_one_with_surface_rows():
    """No surface row is dropped or duplicated in the projection."""
    surface_keys = sorted(surface_index().keys())
    proj_keys = sorted(canonical_cell_key(c.organism, c.target) for c in amr_cells())
    assert proj_keys == surface_keys


def test_amr_claim_status_is_verbatim_from_surface():
    """The structural status (M1 separate field) is carried verbatim from the frozen surface."""
    surf = surface_index()
    for c in amr_cells():
        row = surf[canonical_cell_key(c.organism, c.target)]
        assert c.claim_status == row["phenotype_source_status"]


# ------------------------- honesty guardrails -------------------------

def test_no_numeric_confidence_field():
    """GUARDRAIL: no float/int 'confidence' field exists on the contract — categorical only."""
    for f in dataclasses.fields(CellContract):
        assert f.type not in ("float", "int", float, int), f"numeric field on CellContract: {f.name}"
    # and no live value is a float/bool numeric confidence
    for c in cells():
        for f in dataclasses.fields(CellContract):
            v = getattr(c, f.name)
            assert not isinstance(v, (float, bool)), f"{c.cell_id}.{f.name} is numeric/bool: {v!r}"


def test_tiers_and_vocab_are_categorical_enums_no_scale():
    """GUARDRAIL: evidence_tier + abstention_vocab are string enums, carrying no numeric scale."""
    for c in cells():
        assert isinstance(c.evidence_tier, EvidenceTier)
        assert isinstance(c.abstention_vocab, AbstentionVocab)
        assert isinstance(c.evidence_tier.value, str)
        assert isinstance(c.abstention_vocab.value, str)


def test_native_abstention_maps_to_declared_vocab():
    """Each cell's native abstention term is a known native term AND collapses to its declared vocab."""
    for c in cells():
        assert c.native_abstention in NATIVE_TO_VOCAB, f"unknown native term: {c.native_abstention!r}"
        assert to_vocab(c.native_abstention) == c.abstention_vocab, (
            f"{c.cell_id}: native {c.native_abstention!r} -> {to_vocab(c.native_abstention)} "
            f"!= declared {c.abstention_vocab}"
        )


def test_all_native_terms_map_to_real_vocab_members():
    for native, vocab in NATIVE_TO_VOCAB.items():
        assert isinstance(vocab, AbstentionVocab), f"{native} -> non-vocab {vocab!r}"


# ------------------------- declared falsifier refs -------------------------

def test_falsifier_refs_resolve_or_none():
    """Every non-'none' falsifier_ref resolves to an existing file (declared, not executed)."""
    for c in cells():
        if c.falsifier_ref != "none":
            assert (REPO / c.falsifier_ref).exists(), f"{c.cell_id}: falsifier_ref missing: {c.falsifier_ref}"


def test_incoming_gate_subset_of_known_gates():
    """incoming_data_gate is 'n/a' or a subset of the 8 named rejection gates G1-G8."""
    known = {f"G{i}" for i in range(1, 9)}
    for c in cells():
        if c.incoming_data_gate != "n/a":
            tokens = {t.strip() for t in c.incoming_data_gate.split(",")}
            assert tokens <= known, f"{c.cell_id}: unknown gate token in {c.incoming_data_gate!r}"
