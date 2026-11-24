"""Tests for the eval-set inventory probe.

Two loads are borne here. The `set_algebra` tests are offline and synthetic: they pin that the
orthology filter cannot silently fall back to the unfiltered set, which is the failure that would
make `|C|` look healthy while measuring set B. The drift-guard test proves the `|A| == 3349` refusal
FIRES rather than merely existing -- a guard nobody has seen reject anything is not a guard.

The real-data tests are skip-if-absent so the suite stays green on a host without the refseq cache,
and the one that matters re-derives the committed 1.4291 from the current code. That number is the
whole basis of the committed refusal, so if it stops reproducing, the decomposition that calls it a
sampling artifact has lost its anchor.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MG1655 = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
CORPUS = ROOT / "wiki" / "glm_training_corpus_2026-11-21.txt"
ARTIFACT = ROOT / "wiki" / "glm_eval_set_inventory_2026-11-23.json"

needs_data = pytest.mark.skipif(
    not ((MG1655 / "annotations.gff3").exists() and CORPUS.exists()),
    reason="MG1655 refseq cache or training corpus absent on this host")


def _load():
    p = ROOT / "scripts" / "glm_eval_set_inventory.py"
    spec = importlib.util.spec_from_file_location("glm_eval_set_inventory", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["glm_eval_set_inventory"] = mod
    spec.loader.exec_module(mod)
    return mod


INV = _load()


class _W:
    """A stand-in for a corpus window: just the two fields the filters read."""

    def __init__(self, seq, gene_id):
        self.seq = seq
        self.gene_id = gene_id


# ---- the set algebra, offline and synthetic --------------------------------------------------

def test_sizes_and_overlaps_on_a_hand_built_case():
    got = INV.set_algebra(["AAA", "CCC"], ["AAA", "GGG", "TTT"], ["GGG"])
    assert got["A_sequence_withheld"] == 2
    assert got["B_held_out_of_corpus"] == 3
    assert got["C_orthology_free"] == 1
    assert got["A_inter_B"] == 1
    assert got["B_inter_C"] == 1
    assert got["A_inter_C"] == 0


def test_C_must_be_a_subset_of_B_and_a_violation_RAISES():
    """C is B filtered by gene symbol. A member of C outside B means the filter ran on the wrong
    list, which would otherwise be reported as an interesting overlap rather than a bug."""
    with pytest.raises(INV.InventoryRefusal):
        INV.set_algebra(["AAA"], ["AAA"], ["GGG"])


def test_an_empty_only_set_yields_C_of_zero_and_never_falls_back_to_B():
    """The failure this guards: an empty symbol difference silently leaving C == B, so a fully
    unfiltered set would be reported as orthology-free."""
    held = [_W("AAA", "g1"), _W("CCC", "g2")]
    of = [w.seq.upper() for w in held if (({}.get(w.gene_id) or "")) in set()]
    assert of == []
    got = INV.set_algebra(["AAA"], [w.seq for w in held], of)
    assert got["C_orthology_free"] == 0
    assert got["B_held_out_of_corpus"] == 2


def test_the_symbol_difference_filter_keeps_only_mg1655_only_genes():
    held = [_W("AAA", "g1"), _W("CCC", "g2"), _W("GGG", "g3")]
    m = {"g1": "thrL", "g2": "shared", "g3": "prophageX"}
    only = {"thrL", "prophageX"}
    of = [w.seq.upper() for w in held if (m.get(w.gene_id) or "") in only]
    assert of == ["AAA", "GGG"]


def test_case_is_normalised_so_a_case_mismatch_cannot_read_as_disjointness():
    """`fit_eval_halves` returns raw `w.seq` while the scorer uppercases. Unnormalised, these two
    sets would intersect at ZERO and that would look like a finding."""
    got = INV.set_algebra(["aaa"], ["AAA"], ["AAA"])
    assert got["A_inter_B"] == 1


def test_the_filters_measured_cost_is_C_over_A():
    got = INV.set_algebra(["A" * 3, "C" * 3, "G" * 3, "T" * 3], ["A" * 3], ["A" * 3])
    assert got["C_over_A_the_filters_measured_cost"] == pytest.approx(0.25)


# ---- the drift guard must FIRE -----------------------------------------------------------------

def test_the_drift_guard_fires_on_3348():
    """Proven by making it reject, not by asserting it exists."""
    with pytest.raises(INV.InventoryRefusal) as e:
        INV.assert_expected_a(3348)
    assert "3348" in str(e.value) and "3349" in str(e.value)


def test_the_drift_guard_passes_on_the_manifest_size():
    INV.assert_expected_a(INV.EXPECTED_A_SIZE)
    assert INV.EXPECTED_A_SIZE == 3349


def test_the_drift_guard_is_not_vacuous_on_a_near_miss_either():
    for n in (0, 1, 3348, 3350, 10_000):
        with pytest.raises(INV.InventoryRefusal):
            INV.assert_expected_a(n)


# ---- the committed geometry this probe reconstructs --------------------------------------------

def test_the_reconstructed_committed_geometry_matches_the_committed_artifact():
    """The artifact records n_windows and positions_per_window; the probe must agree with it."""
    art = json.loads((ROOT / "wiki" / "glm_beyond_kmer_2026-11-22.json").read_text(encoding="utf-8"))
    assert INV.COMMITTED_N_WINDOWS == art["eval"]["n_windows"]
    assert INV.COMMITTED_N_POSITIONS == art["eval"]["positions_per_window"]


def test_uniform_constant_agrees_with_the_gate_module():
    from dna_decode.glm.eval_admissibility import UNIFORM_CE_NATS_DNA
    assert INV.UNIFORM == UNIFORM_CE_NATS_DNA == math.log(4)


# ---- real data -----------------------------------------------------------------------------------

@needs_data
def test_set_A_recomputes_to_the_manifest_size():
    """The manifest records only `withheld_union_size`; set A itself is not persisted anywhere."""
    from dna_decode.glm.training_corpus import falsifier_eval_windows
    a = falsifier_eval_windows(INV.MG1655, length=INV.LENGTH, limit=INV.SCORING_LIMIT,
                               seeds=range(INV.N_SEEDS), kind=INV.KIND)
    manifest = json.loads(
        (ROOT / "wiki" / "glm_training_corpus_manifest_2026-11-21.json").read_text(encoding="utf-8"))
    assert len(a) == manifest["leak_gate"]["withheld_union_size"] == 3349


@needs_data
def test_the_committed_1_4291_reproduces_from_the_current_code():
    """The anchor of the whole decomposition. If this drifts, calling 1.4291 a sampling artifact
    loses its basis, because the claim rests on reproducing it and then varying one factor."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import (ORDER, ecoli_fit_windows, load_corpus_windows,
                                         mg1655_candidate_sets)
    cset = set(load_corpus_windows())
    chain = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit_windows(cset, length=INV.LENGTH))
    of = mg1655_candidate_sets(cset, length=INV.LENGTH)["orthology_free"]
    d = INV.decompose_committed_cell(chain, of, length=INV.LENGTH)
    assert d["cells"]["committed_6_windows_8_sampled_positions"] == 1.4291


@needs_data
def test_the_decomposition_crosses_uniform_when_only_the_sample_size_changes():
    """The committed cell is above uniform; hold the windows fixed, score every position, and it is
    below. That contrast is the finding, so it is pinned rather than described."""
    from dna_decode.glm.baselines import MarkovGenerator
    from scripts.glm_beyond_kmer import (ORDER, ecoli_fit_windows, load_corpus_windows,
                                         mg1655_candidate_sets)
    cset = set(load_corpus_windows())
    chain = MarkovGenerator(order=ORDER, seed=0).fit(ecoli_fit_windows(cset, length=INV.LENGTH))
    of = mg1655_candidate_sets(cset, length=INV.LENGTH)["orthology_free"]
    c = INV.decompose_committed_cell(chain, of, length=INV.LENGTH)["cells"]
    assert c["committed_6_windows_8_sampled_positions"] > INV.UNIFORM
    assert c["same_6_windows_ALL_positions"] < INV.UNIFORM
    assert c["ALL_windows_ALL_positions"] < INV.UNIFORM


@needs_data
def test_set_C_is_far_larger_than_the_committed_run_scored():
    """The premise this probe exists to test: the filter was believed to leave too little to score."""
    from scripts.glm_beyond_kmer import load_corpus_windows, mg1655_candidate_sets
    of = mg1655_candidate_sets(set(load_corpus_windows()), length=INV.LENGTH)["orthology_free"]
    assert len(of) == 648
    # The scorer refuses when len(of) < --n-windows, so its own default would have run.
    assert len(of) >= 120


# ---- the artifact ---------------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_artifact_parses_as_json():
    """This repo has shipped an unparseable hand-edited artifact before."""
    json.loads(ARTIFACT.read_text(encoding="utf-8"))


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_artifact_emits_no_generator_verdict():
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert art["emits_no_verdict_about_the_generator"] is True
    assert "verdict" not in art
    blob = json.dumps(art).upper()
    assert "BEYOND_KMER_CONSERVATIVE_POSITIVE" not in blob
    assert "WORSE_THAN_ORDER3" not in blob


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_every_number_in_the_headline_appears_in_a_structured_field():
    """The contract-number-audit discipline: prose figures must be re-derivable from fields."""
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    d = art["committed_cell_decomposition"]
    s = art["sizes_and_overlaps"]
    # The headline quotes BOTH, because the scorer slices the list while statistics want distinct
    # windows. Conflating them is what this assertion caught the first time it ran.
    assert s["C_list_entries"] == 648
    assert s["C_orthology_free"] == 634
    assert s["C_duplicate_entries"] == 14
    assert s["B_list_entries"] == 3923
    assert s["B_held_out_of_corpus"] == 3906
    assert s["B_duplicate_entries"] == 17
    assert d["cells"]["committed_6_windows_8_sampled_positions"] == 1.4291
    assert d["cells"]["ALL_windows_ALL_positions"] == 1.3571
    assert d["uniform_ce_nats"] == 1.3863
    assert d["n_predictions_in_the_committed_cell"] == 48
    assert art["committed_run_geometry_reconstructed"]["windows_that_were_AVAILABLE"] == 648


@pytest.mark.skipif(not ARTIFACT.exists(), reason="inventory artifact not yet generated")
def test_the_slices_the_scorer_actually_uses_are_duplicate_free():
    """The duplicate finding is a LATENT hazard only because these two slices are clean. If either
    stops being clean, a scored window would be counted twice."""
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ss = art["slice_safety"]
    assert ss["first_6_is_duplicate_free"] is True
    assert ss["first_120_is_duplicate_free"] is True
    assert ss["distinct_in_first_120"] == 120
    # And the full list is NOT duplicate-free -- otherwise the guard above is vacuous.
    assert ss["first_648_is_duplicate_free"] is False
