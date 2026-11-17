"""Guards for `scripts/glm_active_frame_gate.py`.

Three things here matter more than the rest.

1. The verdict is ORDERED and the two indeterminates outrank every graded outcome, so a run
   with a dirty null or a failed pooled reconcile cannot report a margin.
2. `scripts/glm_encoder_gate.py` must be byte-identical to HEAD after this module is imported.
   The active frame breaks that file's reconcile anchor by construction, so the only safe
   relationship is read-only; an edit there would rewrite a frozen pre-registration.
3. The Phase 3 binned table must NOT be able to move the verdict. It is a shape, not a second
   bar, and a table that could grade would be a threshold chosen after seeing the data.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from scripts.glm_active_frame_gate import (
    PREREGISTERED_ACTIVE_FRAME,
    VERDICTS,
    compute_verdict,
    framing_recommendation,
)


def _arms(*, reconcile=True, null=0.01, margins=None, ci=(0.25, 0.35)):
    return {
        "pooled_reconcile_ok": reconcile,
        "null_max_abs": null,
        "active_margin_per_seed": list(margins if margins is not None else [0.30] * 10),
        "active_margin_ci": ci,
    }


# ----------------------------------------------------------------------------------------
# every outcome is reachable, and the ordering holds
# ----------------------------------------------------------------------------------------

def test_frame_explains_the_gap_is_reachable():
    v, _d = compute_verdict(_arms(margins=[0.30] * 10))
    assert v == "FRAME_EXPLAINS_THE_GAP"


def test_frame_contributes_is_reachable():
    """Above the contributes bar, below the explains bar, CI excluding zero."""
    v, _d = compute_verdict(_arms(margins=[0.09] * 10, ci=(0.06, 0.12)))
    assert v == "FRAME_CONTRIBUTES"


def test_frame_does_not_explain_is_reachable_on_a_small_margin():
    v, _d = compute_verdict(_arms(margins=[0.02] * 10, ci=(0.01, 0.03)))
    assert v == "FRAME_DOES_NOT_EXPLAIN"


def test_frame_does_not_explain_is_reachable_on_a_straddling_ci():
    """A large median with a CI through zero must NOT grade as a positive."""
    v, _d = compute_verdict(_arms(margins=[0.09] * 10, ci=(-0.02, 0.20)))
    assert v == "FRAME_DOES_NOT_EXPLAIN"


def test_indeterminate_null_not_clean_is_reachable():
    v, _d = compute_verdict(_arms(null=0.20))
    assert v == "INDETERMINATE_NULL_NOT_CLEAN"


def test_indeterminate_pooled_reconcile_failed_is_reachable():
    v, _d = compute_verdict(_arms(reconcile=False))
    assert v == "INDETERMINATE_POOLED_RECONCILE_FAILED"


def test_a_dirty_null_forces_indeterminate_whatever_the_margin():
    """The plan's explicit case: a null at 0.2 outranks any candidate margin."""
    for margins in ([0.60] * 10, [0.30] * 10, [0.01] * 10):
        v, _d = compute_verdict(_arms(null=0.20, margins=margins))
        assert v == "INDETERMINATE_NULL_NOT_CLEAN", f"margins {margins[0]} escaped the null gate"


def test_a_failed_pooled_reconcile_outranks_even_a_dirty_null():
    """Ordering: if the read path is wrong, nothing downstream is interpretable."""
    v, _d = compute_verdict(_arms(reconcile=False, null=0.99, margins=[0.60] * 10))
    assert v == "INDETERMINATE_POOLED_RECONCILE_FAILED"


def test_the_explains_bar_needs_the_seed_count_not_just_the_median():
    """A median over the bar carried by a few seeds must not grade as EXPLAINS."""
    margins = [0.40, 0.40, 0.40, 0.40, 0.40, 0.16, 0.02, 0.02, 0.02, 0.02]
    v, _d = compute_verdict(_arms(margins=margins, ci=(0.10, 0.30)))
    assert v != "FRAME_EXPLAINS_THE_GAP"


def test_every_declared_verdict_is_reachable_from_some_input():
    """No dead verdicts: a verdict nothing can return is decoration."""
    reached = {
        compute_verdict(_arms(margins=[0.30] * 10))[0],
        compute_verdict(_arms(margins=[0.09] * 10, ci=(0.06, 0.12)))[0],
        compute_verdict(_arms(margins=[0.02] * 10, ci=(0.01, 0.03)))[0],
        compute_verdict(_arms(null=0.20))[0],
        compute_verdict(_arms(reconcile=False))[0],
    }
    assert reached == set(VERDICTS)


# ----------------------------------------------------------------------------------------
# the framing recommendation map is TOTAL
# ----------------------------------------------------------------------------------------

def test_the_framing_recommendation_map_is_total_over_the_verdict_set():
    """Asserted by key-set equality, so adding a verdict without a recommendation fails."""
    assert set(PREREGISTERED_ACTIVE_FRAME["framing_recommendation"].keys()) == set(VERDICTS)


def test_every_verdict_resolves_to_a_recommendation():
    for v in VERDICTS:
        assert isinstance(framing_recommendation(v), str)
        assert framing_recommendation(v).strip()


def test_an_unknown_verdict_refuses_rather_than_returning_a_default():
    with pytest.raises(KeyError, match="TOTAL over the verdict set"):
        framing_recommendation("FRAME_PROBABLY_SOMETHING")


def test_the_pessimistic_entry_recommends_no_reframing():
    """The load-bearing entry: no measured headroom means do not re-base future numbers."""
    rec = framing_recommendation("FRAME_DOES_NOT_EXPLAIN")
    assert "KEEP THE POOLED FRAME" in rec.upper()


def test_neither_indeterminate_recommends_a_framing_change():
    for v in ("INDETERMINATE_NULL_NOT_CLEAN", "INDETERMINATE_POOLED_RECONCILE_FAILED"):
        assert "NO FRAMING CHANGE" in framing_recommendation(v).upper()


# ----------------------------------------------------------------------------------------
# the binned table is a shape, not a bar
# ----------------------------------------------------------------------------------------

def test_the_binned_table_cannot_change_the_verdict():
    """Called with two very different binned tables and identical medians -> same verdict.

    `compute_verdict` takes no binned argument at all, which is the structural guarantee; this
    test pins that the signature stays that way by passing tables through the arms dict and
    confirming they are ignored.
    """
    base = _arms(margins=[0.30] * 10)
    flat = dict(base, binned={"bins": [{"bin": i, "margin": 0.0} for i in range(10)]})
    steep = dict(base, binned={"bins": [{"bin": i, "margin": 0.9} for i in range(10)]})
    assert compute_verdict(flat) == compute_verdict(steep) == compute_verdict(base)


# ----------------------------------------------------------------------------------------
# the frozen pooled pre-registration is read-only
# ----------------------------------------------------------------------------------------

def test_glm_encoder_gate_is_byte_identical_to_HEAD_after_importing_this_module():
    import scripts.glm_active_frame_gate  # noqa: F401  (import is the point)

    r = subprocess.run(
        ["git", "diff", "--quiet", "HEAD", "--", "scripts/glm_encoder_gate.py"],
        cwd=Path(__file__).resolve().parents[1], capture_output=True,
    )
    assert r.returncode == 0, (
        "scripts/glm_encoder_gate.py differs from HEAD. Its reconcile anchor is broken by the "
        "active frame by construction, so it must be imported read-only and never edited"
    )


def test_this_gate_defines_its_OWN_pre_registration_distinct_from_the_pooled_one():
    from scripts.glm_encoder_gate import PREREGISTERED as pooled

    assert PREREGISTERED_ACTIVE_FRAME is not pooled
    assert PREREGISTERED_ACTIVE_FRAME != pooled
    # the pooled anchor survives unchanged inside the pooled dict
    assert pooled["reconcile_anchor"]["expected"] == 0.3000
    # and this bar keeps it only as a POOLED-phase check, never as an active-frame one
    assert PREREGISTERED_ACTIVE_FRAME["pooled_reconcile_anchor"]["expected"] == 0.3000
    assert "why_pooled" in PREREGISTERED_ACTIVE_FRAME["pooled_reconcile_anchor"]


def test_the_bar_says_twice_that_it_is_not_blind():
    assert PREREGISTERED_ACTIVE_FRAME["not_a_blind_pre_registration"] is True
    assert PREREGISTERED_ACTIVE_FRAME["frozen_after_the_mechanism_analysis"] is True
    assert PREREGISTERED_ACTIVE_FRAME["why_not_blind"].strip()
    assert PREREGISTERED_ACTIVE_FRAME["why_frozen_after"].strip()


def test_the_bar_pins_the_published_parameter_count_and_protocol():
    assert PREREGISTERED_ACTIVE_FRAME["expected_n_parameters"] == 13217
    assert PREREGISTERED_ACTIVE_FRAME["registered_seeds"] == list(range(10))
    assert "REGISTERED_PROTOCOL" in PREREGISTERED_ACTIVE_FRAME["protocol"]
    assert PREREGISTERED_ACTIVE_FRAME["device"] == "cpu"


# ----------------------------------------------------------------------------------------
# thread pinning, asserted from the environment rather than from a timing
# ----------------------------------------------------------------------------------------

def test_importing_the_module_guarantees_the_four_thread_variables_are_set():
    """Asserted by reading the environment, NOT by timing anything.

    `os.environ.setdefault` means the module guarantees PRESENCE and leaves an operator's own
    value alone, so the claim under test is presence rather than a specific value.
    """
    import scripts.glm_active_frame_gate  # noqa: F401

    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        assert v in os.environ, f"{v} was not set by importing the gate module"


def test_the_pinning_happens_before_numpy_is_imported_in_the_source():
    """Order matters: BLAS reads these at load time, so a late setdefault is a no-op."""
    src = (Path(__file__).resolve().parents[1]
           / "scripts" / "glm_active_frame_gate.py").read_text(encoding="utf-8")
    pin = src.index("OMP_NUM_THREADS")
    imp = src.index("import numpy as np")
    assert pin < imp, "thread pinning must precede the numpy import or it does nothing"


# ----------------------------------------------------------------------------------------
# end-to-end, opt-in because it trains a conv on CPU
# ----------------------------------------------------------------------------------------

@pytest.mark.skipif(
    os.environ.get("DNA_DECODE_RUN_SLOW") != "1",
    reason="trains a conv on CPU; set DNA_DECODE_RUN_SLOW=1 to run. A manual --seeds 1 "
           "--subsample smoke is recorded in the plan's execution log instead",
)
def test_end_to_end_smoke_on_a_tile_subsample(tmp_path):
    from scripts.glm_active_frame_gate import main

    rc = main(["--seeds", "1", "--subsample", "1200", "--n-bins", "4",
               "--out", str(tmp_path / "smoke.json")])
    assert rc == 0
    import json
    out = json.loads((tmp_path / "smoke.json").read_text(encoding="utf-8"))
    assert out["smoke"] is True
    assert out["verdict"] == "SMOKE_NOT_GRADED"
    assert out["phase1_pooled"]["reconcile_ok"] is True
