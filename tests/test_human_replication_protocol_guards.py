"""Pin the preconditions added 2026-11-22 after P2's statistics re-test.

WHY THESE EXIST. Each one closes a defect that produced, or would have produced, a number that LOOKS
like a measurement:

  * `spearman` returned 0.0 on zero rank variance. 0.0 is an ordinary Spearman value, so a dead
    conditioning knob, a harness defect where no target varied, and total plumbing failure all
    arrived as one indistinguishable reading -- and the middle one would have read as the frozen
    bar's NO_TARGET_CONTROL outcome. Five other copies of the same guard in this repo return NaN or
    None; this arm was the lone copy returning a plausible number.
  * `eval_loss` reseeded the GLOBAL torch RNG every --log-every steps, making the training
    corruption draws periodic with a LOGGING parameter and giving every seed the same noise between
    evals. The bar's tolerance IS the seed spread, so that attacked the decision quantity itself.
  * Nothing checked the conditioning signal reaches the logits. `--shuffle-labels` tests whether the
    LABEL carries information, never whether it ARRIVES; a dead embedding would make both arms read
    zero and look like a clean negative.
  * The no-information arm drew ONE permutation, so every contrast from it was conditional on that
    single draw.

Structural pins are used only where this host cannot exercise the path. They are the ones that fail
if someone reintroduces a defect on a machine that cannot run the real job.
"""
from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ORACLE = REPO / "scripts" / "glm_human_activity_oracle.py"
TRAINER = REPO / "scripts" / "server1" / "glm_human_diffusion_server1.py"
PROTOCOL = REPO / "wiki" / "human_replication_and_multiplicity_protocol_2026-11-22.json"
RESERVATION = REPO / "wiki" / "human_final_shard_reservation_2026-11-22.json"
EXPORTER = REPO / "scripts" / "glm_export_server1_shard.py"

torch = pytest.importorskip("torch")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------- the silent-zero, behaviourally

@pytest.fixture(scope="module")
def oracle():
    return _load("_oracle_under_test", ORACLE)


@pytest.mark.parametrize("x,y,side", [
    ([1, 2, 3, 4, 5], [7, 7, 7, 7, 7], "y"),       # constant predictions -> dead knob
    ([3, 3, 3, 3, 3], [1, 5, 2, 4, 3], "x"),       # constant targets -> harness defect
    ([3, 3, 3], [7, 7, 7], "both"),                # total plumbing failure
])
def test_spearman_refuses_zero_variance_and_names_the_side(oracle, x, y, side):
    with pytest.raises(oracle.DegenerateSpearman) as got:
        oracle.spearman(x, y)
    assert got.value.side == side, (
        "the degenerate SIDE must be named: 'y' is a dead conditioning knob (a reportable state), "
        "'x' is a broken harness (not a result). Collapsing them is what NaN would do.")


def test_spearman_never_returns_a_plausible_zero_on_degenerate_input(oracle):
    """The specific regression: returning 0.0, which is indistinguishable from a measured null."""
    for x, y in [([1, 2, 3], [5, 5, 5]), ([2, 2, 2], [1, 2, 3])]:
        try:
            got = oracle.spearman(x, y)
        except oracle.DegenerateSpearman:
            continue
        pytest.fail("spearman returned %r for zero-variance input instead of refusing; a reader "
                    "cannot tell that from a real null" % got)


def test_spearman_still_measures_the_non_degenerate_case(oracle):
    """Non-vacuity in the other direction: the fix must not have broken measurement."""
    assert oracle.spearman([1, 2, 3, 4, 5], [1, 2, 3, 4, 5]) == pytest.approx(1.0)
    assert oracle.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]) == pytest.approx(-1.0)
    assert oracle.spearman([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]) == pytest.approx(0.8)
    # heavy ties must still average ranks rather than let sort order become a finding
    assert oracle.spearman([0, 0, 1, 1], [0, 1, 0, 1]) == pytest.approx(0.0)


# --------------------------------------------------------- the knob precondition, behaviourally

def test_the_conditioning_knob_moves_the_logits_and_the_probe_can_detect_a_dead_one():
    """The bar: 'an invariance from an unverified knob is plumbing, not a null result.'"""
    tr = _load("_trainer_under_test", TRAINER)
    cfg = dict(tr.SIZES["s"], dropout=0.1)
    model = tr.build_model(torch, cfg, seed=0)
    model.eval()

    def probe(m) -> float:
        with torch.no_grad():
            x = torch.full((1, tr.ELEMENT_LEN), tr.MASK_ID, dtype=torch.long)
            lo = m(x, torch.tensor([0], dtype=torch.long))
            hi = m(x, torch.tensor([tr.N_BINS - 1], dtype=torch.long))
            return float((lo - hi).abs().max())

    assert probe(model) > 0.0, "the quintile embedding does not reach the logits as shipped"
    with torch.no_grad():                       # what a dead knob looks like at the output
        model.cond.weight.zero_()
    assert not probe(model) > 0.0, (
        "GUARD IS VACUOUS: a zeroed conditioning embedding still moved the logits, so the preflight "
        "cannot detect the state it exists to catch")


# --------------------------------------------------------- structural pins (no GPU, no shard)

def test_eval_does_not_leave_the_global_rng_reseeded():
    src = TRAINER.read_text(encoding="utf-8")
    i = src.index("def eval_loss()")
    body = src[i:i + 2200]
    # Strip the docstring FIRST. It explains the defect and therefore quotes the very call this
    # test pins, so an ordering check over the raw text matches the PROSE and not the CODE -- which
    # is exactly how this test failed on its first run.
    q = body.find('"""')
    code = body[body.find('"""', q + 3) + 3:] if q != -1 else body
    assert "torch.manual_seed(999)" in code, "the fixed eval corruption draw went missing"
    assert "get_rng_state()" in code and "set_rng_state(" in code, (
        "eval_loss must SAVE and RESTORE the global torch RNG around its manual_seed. Without the "
        "bracket the training noise is periodic with --log-every and every seed shares the "
        "999-derived stream, which corrupts the seed spread the bar uses as its tolerance")
    assert code.index("get_rng_state()") < code.index("torch.manual_seed(999)"), (
        "the state must be captured BEFORE the reseed, or there is nothing to restore")
    assert code.index("torch.manual_seed(999)") < code.index("set_rng_state("), (
        "the restore must come AFTER the reseed, or the bracket does nothing")
    assert "cuda.set_rng_state_all" in code, (
        "torch.manual_seed seeds CPU *and* every CUDA device, so both must be restored")


def test_checkpoints_carry_rng_state_so_a_resumed_run_is_the_same_experiment():
    src = TRAINER.read_text(encoding="utf-8")
    assert '"rng"' in src, (
        "the checkpoint must persist RNG state. Without it a resumed run replays batch order and "
        "masking from step 0 while the step counter carries on, so an interrupted seed is a "
        "different experiment -- and this job runs on a box that is EXPECTED to be reclaimed")
    for needed in ("g_cpu.getstate()", "torch.get_rng_state()", "g_cpu.setstate(", "set_rng_state("):
        assert needed in src, "missing RNG %s" % needed


def test_the_no_information_arm_can_be_run_at_several_permutation_draws():
    src = TRAINER.read_text(encoding="utf-8")
    assert "--shuffle-seed" in src, (
        "the permutation is a random draw; a single draw leaves every no-information number "
        "conditional on it, so the draw must be settable independently of --seed")
    assert "a.shuffle_seed" in src
    i = src.index("if a.shuffle_labels:")
    assert "ONE permutation draw" in src[i:i + 1200], (
        "when the draw is derived from --seed the run must SAY it is a single draw; silence there is "
        "how a one-draw control gets read as a measured contrast")


# --------------------------------------------------------- the artifact corrections

def test_the_exporter_no_longer_contradicts_the_frozen_bar_about_the_nulls_strata():
    src = EXPORTER.read_text(encoding="utf-8")
    assert "are NOT the scoring strata" not in src, (
        "the manifest claimed the scorer re-derives its strata and that reusing folds-1-6 "
        "boundaries would fit the readout to the training set. The frozen bar fits the "
        "activity-stratified null per activity quintile of folds 1-6 DELIBERATELY, so a reader "
        "following the old text would build a null the bar does not specify")
    assert "deliberate, not leakage" in src, "the corrected statement went missing"
    assert "CORRECTION 2026-11-22" in src, (
        "the retraction must stay in the generated manifest. The phrase 're-derives its strata' "
        "still appears in this file ON PURPOSE -- inside the correction, quoting what was wrong -- "
        "so this test pins the old ASSERTION ('are NOT the scoring strata') rather than a substring "
        "the retraction legitimately contains.")


def test_the_reservation_records_that_its_mechanism_is_retracted():
    d = json.loads(RESERVATION.read_text(encoding="utf-8"))
    block = d.get("RETRACTED_MECHANISM_READ_THIS_FIRST")
    assert block, "the reservation must carry its retraction where a reader meets it first"
    assert "Repeated interrogation" in block["retracted_text_verbatim"], (
        "the retracted text must be preserved verbatim; a retraction that deletes what it retracts "
        "cannot be checked")
    assert "SELECTION" in block["why_it_is_wrong"]
    assert "RETRACTED" in d["purpose"], (
        "the purpose field itself stated the retracted mechanism and must point at the retraction")


def test_the_protocol_resolves_the_bars_ambiguity_in_the_conservative_direction():
    d = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    r1 = d["RULE_1_THE_UNIT_OF_REPLICATION"]
    assert "ACROSS-TRAINED-RUN" in r1["which_variance_gates_the_verdict"]["resolved"], (
        "the tolerance must be the OUTER variance. Choosing the inner sampling-only spread is the "
        "significance-manufacturing direction -- the same error the peer reported, one level up")
    assert len(r1["the_nesting_is_THREE_deep"]) == 3
    assert d["honest_limits"], "a protocol with no stated limits is overclaiming"
