"""Evidence-Contract Registry v0.1 — the coverage + consistency contract IS the deliverable.

Pins (v0.1 scope = the FULL CLI-routable surface):
  - COVERAGE: every CLI-routable cell has >=1 contract — every `dna-amr --drug` (bacterial + fungal +
    antimalarial + influenza + HIV + SARS), every `dna-pgx --gene`, and every typing/finder trait. The
    manifest is derived LIVE from the CLI catalogs (brainstorm C1 per-route manifest), so it cannot drift.
  - CONSISTENCY: the AMR projection EQUALS the frozen shipped_decoder_surface, joined via canonical_cell_key
    (NOT raw cell_id string — brainstorm C2); and registry.surface_index() == the frozen surface_index (the
    report card now reads its grid from the registry).
  - HONESTY GUARDRAILS: no numeric confidence field anywhere; tiers + abstention are categorical enums.
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

from dna_decode.data import shipped_decoder_surface as sds
from dna_decode.data.cell_key import canonical_cell_key
from dna_decode.data.cell_registry import (
    TRACKS,
    CellContract,
    EvidenceTier,
    amr_cells,
    amr_projection_keys,
    by_cell_id,
    cells,
    cli_routable_cell_ids,
    cli_routable_manifest,
    pgx_cells,
    surface_index,
)
from dna_decode.data.cell_registry_vocab import NATIVE_TO_VOCAB, AbstentionVocab, to_vocab

REPO = Path(__file__).resolve().parent.parent
_PGX_GENES = {"cyp2c19", "cyp2c9", "vkorc1"}


# ------------------------- coverage (the per-route manifest) -------------------------

def test_cell_id_unique():
    ids = [c.cell_id for c in cells()]
    assert len(ids) == len(set(ids)), "duplicate cell_id in the registry"


def test_every_cli_amr_drug_has_a_contract():
    """COVERAGE: every dna-amr --drug (incl HIV/SARS — the invisibility gap) has >=1 amr/viral cell."""
    routable = cli_routable_manifest()["dna-amr"]
    covered = {c.target.lower() for c in cells() if c.route == "dna-amr"}
    missing = routable - covered
    assert not missing, f"CLI-routable dna-amr drugs missing a contract: {sorted(missing)}"


def test_every_pgx_gene_has_exactly_one_contract():
    pgx_targets = [c.target.lower() for c in pgx_cells()]
    assert set(pgx_targets) == _PGX_GENES
    assert len(pgx_targets) == len(set(pgx_targets))


def test_every_typing_finder_trait_has_exactly_one_contract():
    """COVERAGE: every dna_decode.cli TRAIT except amr/pgx is a typing/finder whole-tool cell."""
    traits = cli_routable_manifest()["traits"]
    tf = [c for c in cells() if c.track in ("typing", "finder")]
    covered = {c.target for c in tf}
    assert covered == traits, f"typing/finder coverage mismatch: extra={covered - traits} missing={traits - covered}"
    assert len(tf) == len(covered), "duplicate typing/finder cell"


def test_tracks_are_all_declared():
    for c in cells():
        assert c.track in TRACKS, f"{c.cell_id}: undeclared track {c.track!r}"


def test_routable_cell_id_set():
    assert cli_routable_cell_ids() == set(by_cell_id())


# ------------------------- consistency (the load-bearing join) -------------------------

def test_amr_projection_equals_frozen_surface_via_canonical_key():
    """CONSISTENCY: registry AMR projection == frozen surface, joined via canonical_cell_key (NOT cell_id)."""
    assert amr_projection_keys() == set(sds.surface_index().keys())


def test_registry_surface_index_equals_frozen_surface_index():
    """The report card reads its grid from registry.surface_index() — it must == the frozen surface_index."""
    assert surface_index() == sds.surface_index()


def test_amr_claim_status_is_verbatim_from_surface():
    surf = sds.surface_index()
    for c in amr_cells():
        assert c.claim_status == surf[canonical_cell_key(c.organism, c.target)]["phenotype_source_status"]


# ------------------------- honesty guardrails -------------------------

def test_no_numeric_confidence_field():
    for f in dataclasses.fields(CellContract):
        assert f.type not in ("float", "int", float, int), f"numeric field on CellContract: {f.name}"
    for c in cells():
        for f in dataclasses.fields(CellContract):
            v = getattr(c, f.name)
            assert not isinstance(v, (float, bool)), f"{c.cell_id}.{f.name} is numeric/bool: {v!r}"


def test_tiers_and_vocab_are_categorical_enums_no_scale():
    for c in cells():
        assert isinstance(c.evidence_tier, EvidenceTier)
        assert isinstance(c.abstention_vocab, AbstentionVocab)
        assert isinstance(c.evidence_tier.value, str)
        assert isinstance(c.abstention_vocab.value, str)


def test_native_abstention_maps_to_declared_vocab():
    for c in cells():
        assert c.native_abstention in NATIVE_TO_VOCAB, f"unknown native term: {c.native_abstention!r}"
        assert to_vocab(c.native_abstention) == c.abstention_vocab, f"{c.cell_id}: vocab mismatch"


def test_all_native_terms_map_to_real_vocab_members():
    for native, vocab in NATIVE_TO_VOCAB.items():
        assert isinstance(vocab, AbstentionVocab), f"{native} -> non-vocab {vocab!r}"


# ------------------------- declared falsifier refs + gates -------------------------

def test_falsifier_refs_resolve_or_none():
    for c in cells():
        if c.falsifier_ref != "none":
            assert (REPO / c.falsifier_ref).exists(), f"{c.cell_id}: falsifier_ref missing: {c.falsifier_ref}"


def test_incoming_gate_subset_of_known_gates():
    known = {f"G{i}" for i in range(1, 9)}
    for c in cells():
        if c.incoming_data_gate != "n/a":
            assert {t.strip() for t in c.incoming_data_gate.split(",")} <= known, f"{c.cell_id}: unknown gate"
