"""Schema tests for the Enformer training-budget sourcing record.

These enforce the refusal discipline the artifact claims, which is the same one shipped in
`dna_decode/data/ecoff_catalog.py`: a value without a locator and a quote is REJECTED, never
approximated. A wrong training budget would produce a fully self-consistent feasibility verdict
that is entirely wrong, with nothing downstream able to catch it -- so the schema is the catch.

The tests check STRUCTURE, not wording. The quotes themselves are guaranteed by construction:
the generator pulls each one out of the fetched source by regex and refuses to emit if one fails
to match, so a quote in this artifact cannot be a transcription error.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "wiki" / "encode_roadmap_training_budget_sources.json"

pytestmark = pytest.mark.skipif(not ARTIFACT.exists(), reason="budget-sources artifact absent")


@pytest.fixture(scope="module")
def rec() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def enformer(rec) -> dict:
    return rec["reference_architectures"]["enformer"]


# ------------------------------------------------------------------ the sourcing bar
def _numeric_keys(d: dict) -> list[str]:
    return [k for k, v in d.items() if isinstance(v, (int, float)) and not isinstance(v, bool)]


def test_every_figure_carrying_a_number_has_BOTH_a_url_and_a_quote(enformer):
    """The load-bearing test: a number without a locator AND a quote is not sourced."""
    for name, fig in enformer["figures"].items():
        assert _numeric_keys(fig), f"{name}: a sourced figure with no number is not a figure"
        assert fig.get("source_url"), f"{name}: carries a number but no source_url"
        assert fig.get("verbatim_quote"), f"{name}: carries a number but no verbatim_quote"
        assert str(fig["source_url"]).startswith("https://"), f"{name}: source_url not resolvable"


def test_an_unsourced_entry_carries_NO_number(enformer):
    """An unsourced figure must be a refusal, not a guess wearing a caveat."""
    for name, fig in enformer["unsourced_figures"].items():
        assert fig.get("unsourced") is True, f"{name}: not marked unsourced"
        assert not _numeric_keys(fig), f"{name}: marked unsourced yet carries a number"
        assert fig.get("reason"), f"{name}: unsourced without a stated reason"


def test_the_three_figures_that_block_a_cost_estimate_are_all_refused(enformer):
    """Params, main-run step count and wall-clock: each absent, each recorded as such."""
    uns = enformer["unsourced_figures"]
    for key in ("main_model_parameter_count", "main_model_training_steps",
                "wall_clock_training_duration"):
        assert key in uns, f"{key} must be present AS a refusal, not merely omitted"


def test_the_absent_wall_clock_is_recorded_as_a_MEASURED_absence(enformer):
    """'Not found' must be distinguishable from 'not looked for'."""
    reason = enformer["unsourced_figures"]["wall_clock_training_duration"]["reason"]
    assert "ZERO" in reason and "78,591" in reason, (
        "the wall-clock refusal must cite the scan that found none, so a later reader can tell a "
        "measured absence from an unchecked one")


# ------------------------------------------------------------------ the derived ratio
def test_the_requested_training_time_ratio_is_REFUSED_not_estimated(rec):
    d = rec["derived"]
    assert d["computable"] is False
    assert d["unsourced"] is True
    assert not _numeric_keys(d), "the refused ratio must carry no number at all"


def test_the_derived_substitute_recomputes_from_its_two_NAMED_inputs(rec):
    """So the ratio cannot drift from the inputs it claims to come from."""
    sub = rec["derived"]["what_IS_computable_instead"]
    a = sub["input_a"]["value"]
    b = sub["input_b"]["value"]
    assert b != 0, "denominator must not be zero"
    assert sub["value"] == pytest.approx(a / b), (
        f"stated {sub['value']} != {a}/{b} = {a / b} -- the ratio has drifted from its inputs")
    assert sub["arithmetic"] == f"{a} / {b}", "the stated arithmetic must match the named inputs"


def test_each_derived_input_declares_whether_it_is_sourced_or_measured(rec):
    sub = rec["derived"]["what_IS_computable_instead"]
    for side in ("input_a", "input_b"):
        prov = sub[side]["provenance"]
        assert ("sourced" in prov) or ("measured" in prov), (
            f"{side}: provenance must say whether the number was sourced or measured")


def test_the_substitute_ratio_states_what_it_is_NOT(rec):
    """The whole hazard here is this being read as a training-time ratio."""
    sub = rec["derived"]["what_IS_computable_instead"]
    assert "training-time ratio" in sub["is_NOT"]
    assert "batch-size normaliser" in sub["reading"] or "NOT a step" in sub["reading"]


def test_the_batch_normaliser_matches_the_sourced_hardware_figure(rec, enformer):
    """input_a is not free-floating: it must equal the sourced batch total."""
    hw = enformer["figures"]["training_hardware"]
    sub = rec["derived"]["what_IS_computable_instead"]
    assert sub["input_a"]["value"] == hw["batch_total"] == hw["n_devices"] * hw["batch_per_device"]


# ------------------------------------------------------------------ provenance + honesty
def test_the_extraction_method_is_recorded_as_programmatic(rec):
    m = rec["extraction_method"]
    assert "REFUSES" in m, "the generator's refusal behaviour is part of the quotes' warrant"


def test_the_ablation_scope_note_is_present_on_the_500k_figure(enformer):
    """The one figure a reader is most likely to misattribute to the main run."""
    note = enformer["figures"]["ablation_budget"]["note"]
    assert "MAIN" in note and "unsupported" in note


def test_honest_limits_name_the_single_architecture_scope(rec):
    joined = " ".join(rec["honest_limits"])
    assert "ONE reference architecture" in joined
    assert "hardware CLASS" in joined, "a hardware figure must not read as a compute bound"
