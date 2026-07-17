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
from dna_decode.pgx import PGX_GENES  # M1: single source -> a 4th gene can't pass coverage vacuously
_PGX_GENES = set(PGX_GENES)


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


def test_every_mendelian_target_has_a_contract():
    """COVERAGE: the dna-clinvar (Mendelian) route has exactly one contract per target -> a shipped
    Mendelian decoder cannot ship invisibly to the trust surface."""
    from dna_decode.data.cell_registry import cli_routable_manifest, mendelian_cells
    targets = [c.target for c in mendelian_cells()]
    assert set(targets) == cli_routable_manifest()["dna-clinvar"]
    assert len(targets) == len(set(targets))


def test_every_hla_allele_has_a_contract():
    """COVERAGE: the dna-hla route has exactly one contract per allele -> a shipped HLA decoder cannot ship
    invisibly to the trust surface."""
    from dna_decode.data.cell_registry import cli_routable_manifest, hla_cells
    targets = [c.target for c in hla_cells()]
    assert set(targets) == cli_routable_manifest()["dna-hla"]
    assert len(targets) == len(set(targets))


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


# ------------------------- C1: HIV tiers data-driven from the report card (no overclaim) -------------------------

def _hiv_card_drugs():
    from dna_decode.data import trust_surface
    card = trust_surface._load("hiv_decoder_report_card.json")
    return {c["drug"] for c in (card or {}).get("cells", []) if c.get("drug")}


def test_hiv_not_censused_set_equals_routable_minus_card():
    """The HIV cells tiered NOT_CENSUSED must be EXACTLY the CLI-routable HIV drugs absent from the card.
    A drug with no report-card row may NEVER be independent_measured (the delavirdine overclaim, C1)."""
    from dna_decode.data.hiv_amr import all_supported_hiv_drugs
    not_censused = {c.target for c in cells()
                    if c.organism == "HIV-1" and c.evidence_tier == EvidenceTier.NOT_CENSUSED}
    expected = set(all_supported_hiv_drugs()) - _hiv_card_drugs()
    assert not_censused == expected
    # no card-row HIV drug is NOT_CENSUSED; no uncensused HIV drug claims independent_measured
    for c in cells():
        if c.organism != "HIV-1":
            continue
        in_card = c.target in _hiv_card_drugs()
        if not in_card:
            assert c.evidence_tier == EvidenceTier.NOT_CENSUSED, f"{c.target}: uncensused but tier {c.evidence_tier}"


def test_delavirdine_is_not_censused_regression_pin():
    """Regression pin for the exact shipped overclaim: delavirdine (CLI-routable, no card row) -> NOT_CENSUSED."""
    d = by_cell_id().get("viral:HIV-1:delavirdine")
    assert d is not None and d.evidence_tier == EvidenceTier.NOT_CENSUSED
    assert d.claim_status == "cli_routable_not_validated"


def test_pgx_genes_single_source():
    """M1: the manifest's PGx set derives from the importable PGX_GENES constant, not a copied literal."""
    from dna_decode.data.cell_registry import cli_routable_manifest
    assert cli_routable_manifest()["dna-pgx"] == set(PGX_GENES)


# ------------------------- contract-vs-evidence (added 2026-07-16 after a real mis-registration) ----------

def test_forward_cell_is_registered_as_dms_measured_not_unvalidated():
    """REGRESSION PIN. The forward cell shipped for a day registered as NOT_CENSUSED / 'prior_only_no_
    phenotype_claim' / 'validation_slice: NONE'. That was FALSE: the cell is validated per-variant against
    measured ProteinGym DMS (blaTEM genome-edit Spearman 0.7611 over 1,715 real single-nt variants). The
    contract was written from memory instead of from the cell -- under-claiming a validated cell is as much
    a trust-surface falsehood as over-claiming an unvalidated one."""
    c = by_cell_id()["finder:any:forward"]
    assert c.evidence_tier is EvidenceTier.INDEPENDENT_MEASURED
    assert "0.7611" in c.validation_slice          # the real blaTEM genome-edit number
    assert "proteingym" in c.label_provenance.lower()
    # ...and the scope rails that make the number honest must survive any future edit:
    assert "RANK" in c.demotion_rule                    # rank != magnitude
    assert "REGIME B ONLY" in c.demotion_rule           # never a resistance predictor
    assert "blosum62" in c.demotion_rule                # the CLI default is not the headline


def test_every_wiki_artifact_a_contract_cites_actually_exists():
    """ANTI-DRIFT. A contract that cites its evidence by filename must cite a file that EXISTS -- otherwise
    the trust surface points at nothing and no one notices. Generalizes the mis-registration above: the fix
    for 'written from memory' is 'the evidence must resolve'."""
    import re

    missing = []
    for c in cells():
        for field in (c.validation_slice, c.label_provenance, c.demotion_rule):
            for ref in re.findall(r"wiki/[\w./-]+\.(?:json|md)", field or ""):
                if not (REPO / ref).exists():
                    missing.append(f"{c.cell_id} -> {ref}")
    assert not missing, f"contracts cite non-existent evidence artifacts: {missing}"
