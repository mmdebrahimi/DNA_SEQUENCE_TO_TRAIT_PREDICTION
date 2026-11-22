"""The reconciled memo's numbers must be the artifacts', and its corrections must survive as corrections.

Same two jobs as the memo this supersedes, with one addition that matters more here than there.

1. **The prose cannot drift from the computation**: every headline figure is re-read from the artifact and
   compared.
2. **The corrections must keep correcting.** This memo retracts figures from an earlier memo that was
   itself a set of retractions, so the record is two layers deep. A tidy-up that deleted the withdrawn
   −0.0121 or the withdrawn 1.515 would leave the arm looking as though it had always known the answer,
   and the sequence of corrections is the most instructive thing in the file.
3. **The four open decisions must stay open and must stay priced.** All four turn on
   `COMPARISON_GEOMETRY_FIELDS`, so a reader who takes any of them as a free fix would break the anchor
   chain. The price has to be in the text, not inferred.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MEMO = ROOT / "wiki" / "glm_instrument_reconciled_2026-11-21.md"
RECON = ROOT / "wiki" / "glm_grid_reconcile_2026-11-21.json"
SWEEP = ROOT / "wiki" / "glm_margin_sample_size_sweep_2026-11-21.json"
NOISE = ROOT / "wiki" / "glm_spread_vs_sampling_noise_2026-11-21.json"
LEVERS = ROOT / "wiki" / "glm_instrument_levers_2026-11-21.json"
ORGS = ROOT / "wiki" / "glm_margin_across_organisms_2026-11-21.json"

pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in (MEMO, RECON, SWEEP, NOISE, LEVERS, ORGS)),
    reason="the reconciled instrument arm has not been run on this host")


def _ascii_minus(text: str) -> str:
    """Normalize the typographic minus to an ASCII hyphen before comparing against artifact values.

    The memo is layman-readable prose and uses U+2212 MINUS SIGN, which is correct typography; the JSON
    artifacts serialise an ASCII hyphen. Without this, every negative figure would appear to be missing
    from the memo and the drift check would fail on a character rather than on a number. Normalizing in
    one place beats or-ing both spellings into each assertion, which is what the superseded memo's test
    file had to do.
    """
    return text.replace("−", "-")


@pytest.fixture(scope="module")
def memo():
    return _ascii_minus(MEMO.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def flat(memo):
    """Whitespace-collapsed memo, so a multi-word phrase does not break on a line wrap."""
    return re.sub(r"\s+", " ", memo)


def test_the_minus_sign_normalisation_is_not_vacuous():
    """If the memo ever stopped using a typographic minus, this helper would be dead weight.

    Proving it is live keeps the fixture honest: a normalisation that never fires looks like a safeguard
    and is really just noise, which is the vacuous-control pattern this repo keeps running into.
    """
    raw = MEMO.read_text(encoding="utf-8")
    assert "−" in raw, "the memo no longer uses a typographic minus; drop _ascii_minus"
    assert "-" in _ascii_minus(raw)


def _load(p):
    return json.loads(p.read_text(encoding="utf-8"))


# ------------------------------------------------------------------ ordering (reader action first)

def test_decisions_needed_comes_before_corrections_and_before_the_results(memo):
    d = memo.index("## DECISIONS NEEDED")
    c = memo.index("## CORRECTIONS")
    r = memo.index("## What the instrument can and cannot do")
    assert d < c < r, (d, c, r)


def test_the_status_header_is_in_the_first_lines(memo):
    """`doc_validity` scans the first 30 lines; a status below that is absent, not merely late."""
    first = "\n".join(memo.splitlines()[:30])
    assert "**Status:**" in first
    assert "Supersedes" in first


# ------------------------------------------------------------------ the numbers are the artifacts'

def test_the_like_for_like_margins_match_the_reconcile_artifact(memo):
    rec = _load(RECON)
    for mode in ("kmer", "both"):
        s = rec["summary"][mode]
        assert f"{s['like_for_like_margin_mean']:+.4f}" in memo, (mode, s)
        assert f"{s['perfect_floor_mean']:.4f}" in memo, (mode, s)


def test_the_kmer_interval_and_its_unresolvability_are_stated(memo):
    rec = _load(RECON)
    sig = rec["margin_significance"]["kmer"]
    lo, hi = sig["ci95"]
    assert f"{lo}" in memo and f"{hi}" in memo
    assert sig["resolvable_above_zero"] is False
    assert "unresolved" in memo.lower()


def test_the_sampling_benchmark_and_the_240_replicate_arm_match(memo):
    rec = _load(NOISE)
    assert f"{rec['benchmark']['committed_geometry_sd']}" in memo
    assert f"{rec['benchmark']['grid_observed_sd_both']}" in memo
    for mode in ("kmer", "both"):
        r = rec["replication"][mode]
        assert f"{r['sd_A']}" in memo, (mode, r["sd_A"])
        assert r["A_consistent_with_benchmark"] is True


def test_the_feature_block_comparison_matches_the_levers_artifact(memo):
    rec = _load(LEVERS)
    q1 = rec["q1_feature_blocks"]
    assert f"{q1['both_minus_positional']['mean_diff']:+.4f}" in memo
    assert f"{q1['positional_minus_kmer']['mean_diff']:+.4f}" in memo
    assert q1["both_minus_positional"]["resolvably_nonzero"] is False
    assert q1["positional_minus_kmer"]["resolvably_nonzero"] is True


def test_the_order_sweep_ranges_match(memo):
    rec = _load(LEVERS)
    for mode in ("kmer", "positional"):
        assert f"{rec['q2_null_order'][mode]['range']:.4f}" in memo, mode


def test_the_cross_organism_table_matches_its_artifact(memo):
    rec = _load(ORGS)
    ok = {k: v for k, v in rec["by_substrate"].items() if v.get("status") == "ok"}
    assert len(ok) >= 4
    for name, v in ok.items():
        for mode in ("kmer", "both"):
            assert f"{v[mode]['margin']['mean']:+.4f}" in memo, (name, mode)
        assert f"{v['gc_natural_mean']:.3f}" in memo, name


def test_the_sweep_figure_that_makes_kmer_resolve_is_quoted(memo):
    rec = _load(SWEEP)
    cfg = rec["by_config"]["limit1200_n400"]["kmer"]["margin"]
    assert f"{cfg['mean']:+.4f}" in memo
    assert f"{cfg['ci95'][0]}" in memo and f"{cfg['ci95'][1]}" in memo


# ------------------------------------------------------------------ the retractions must persist

def test_the_withdrawn_kmer_margin_is_still_named_as_withdrawn(flat):
    """−0.0121 must appear AS AN ERROR. Deleting it hides that the arm once believed it."""
    assert "0.0121" in flat
    assert "withdrawn" in flat.lower()
    assert "opposite sign" in flat


def test_the_withdrawn_geometry_is_named_with_both_values(flat):
    assert "400" in flat and "120" in flat
    assert "glm_paired_sweep_kernel.py:45" in flat
    assert "0.0140" in flat and "0.0410" in flat


def test_the_refuted_fit_set_mechanism_is_recorded_with_its_interval(flat):
    """The arm built to test my own stated mechanism returned an interval containing 1."""
    assert "1.129" in flat
    assert "0.891" in flat and "1.445" in flat
    assert "contains 1" in flat


def test_my_replacement_mechanism_is_ALSO_recorded_as_not_replicating(flat):
    """Two mechanisms asserted and both withdrawn -- the second must not quietly vanish."""
    assert "1.515" in flat
    assert "one of five" in flat or "1 of 5" in flat
    assert "0.0352" in flat and "0.0448" in flat


def test_the_overstated_pairing_claim_is_corrected_with_counts(flat):
    """'Pairing helps model-vs-null' became 2 of 4 pairs in one mode, which is a different claim."""
    assert "2 of 4" in flat and "1 of 6" in flat
    assert "0.684" in flat
    assert "not" in flat and "clean rule" in flat


def test_the_retired_generator_distance_figure_is_named(flat):
    """+0.2150 was a distance above the inflated ceiling and is replaced by +0.2561."""
    assert "0.2150" in flat and "0.2561" in flat


# ------------------------------------------------------------------ claims it must NOT make

def test_the_gc_correlation_is_reported_with_its_lineage_confound(flat):
    """r = -0.941 over five points is two lineage groups. This repo has lost four arms to that."""
    assert "0.941" in flat
    assert "effective n of 2" in flat
    assert "four arms" in flat


def test_the_memo_does_not_claim_anti_correlation_between_checkpoints(flat):
    """The negative correlations have intervals spanning zero and must be labelled as such."""
    assert "spanning zero" in flat or "span zero" in flat
    assert "not" in flat


def test_the_memo_keeps_the_necessary_not_sufficient_rail(flat):
    assert "necessary, not sufficient" in flat.lower()
    assert "LOCAL composition" in flat or "local composition" in flat.lower()
    assert "long-range" in flat


def test_the_memo_scopes_itself_to_the_markov_null(flat):
    assert "Markov null" in flat
    assert "not the HF generator" in flat.lower() or "not the hf generator" in flat.lower()


def test_the_benchmark_is_declared_a_lower_bound(flat):
    """It assumes a fixed scoring function; the classifier is refitted per split."""
    assert "lower bound" in flat.lower()
    assert "refitted" in flat


def test_the_source_grid_s_own_invalid_verdict_is_disclosed(flat):
    """I built on a grid whose own verdict is INSTRUMENT_INVALID and had not said so."""
    assert "INSTRUMENT_INVALID" in flat


# ------------------------------------------------------------------ the decisions stay open and priced

def test_all_four_decisions_are_priced_against_the_geometry_contract(flat):
    """A reader who takes any of these as free would break the anchor chain."""
    assert "COMPARISON_GEOMETRY_FIELDS" in flat
    assert "knob_contract.py:243" in flat
    for field in ("feature_modes", "markov_order", "n_generated_per_candidate"):
        assert field in flat, field
    assert "fresh slug" in flat or "fresh Kaggle" in flat


def test_the_geometry_contract_really_locks_those_three_fields():
    """Non-vacuity: the memo's central constraint is checked against the code, not asserted."""
    from dna_decode.glm.knob_contract import COMPARISON_GEOMETRY_FIELDS

    for field in ("feature_modes", "markov_order", "n_generated_per_candidate"):
        assert field in COMPARISON_GEOMETRY_FIELDS, field


def test_the_memo_does_not_self_execute_its_own_recommendations(flat):
    assert "I have **not** edited the plan" in flat or "have not edited the plan" in flat
    assert "authority" in flat.lower()
    assert "arms ledger" in flat.lower() or "arms-ledger" in flat.lower()


def test_the_runner_disjointness_caveat_is_recorded_in_proportion(flat):
    """A real caveat on a claim stamped into every artifact, and a small one. Both halves matter."""
    assert "4323" in flat and "4340" in flat
    assert "0 to 3 of 600" in flat
    assert "too small to move a score" in flat


# ------------------------------------------------------------------ the structural cross-check section

def test_the_memo_reports_the_structural_cross_check_with_the_artifact_s_numbers(memo):
    """Every figure in that section is read from the artifact by script, so drift is the only risk."""
    rec = _load(LEVERS)
    cc = rec["structural_cross_check_second_substrate"]
    assert cc is not None
    assert cc["verdict"] in memo
    assert f"{cc['kmer_over_positional_range_ratio']}" in memo
    for mode in ("kmer", "positional", "both"):
        assert f"{cc['profile'][mode]['range']:.4f}" in memo, mode
        for v in cc["profile"][mode]["margin_by_order"].values():
            assert f"{v:+.4f}" in memo, (mode, v)


def test_the_memo_states_the_falsifying_outcome_for_that_cross_check(flat):
    """A cross-check whose failure mode is unstated can only confirm."""
    assert "falsifying outcome was named" in flat or "falsifying outcome" in flat
    assert "quirk" in flat


def test_the_memo_does_not_claim_the_two_substrates_agree_numerically(flat):
    """The claim is about SHAPE. Campylobacter carries a larger margin and the memo must not blur that."""
    assert "same shape" in flat
