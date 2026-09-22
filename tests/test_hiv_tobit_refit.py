"""Guards for the censoring-aware (Tobit) refit of the HIV v0.1 catalog design.

The load-bearing guard is NOT the drift heuristic the acceptance bar froze -- that heuristic was
mis-specified (see the memo). It is the SYNTHETIC RECOVERY test: generate data from known
coefficients, censor it, and check Tobit recovers the truth better than OLS while both agree exactly
when nothing is censored. That answers "is the implementation correct" directly, where a drift
threshold only proxies it.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.linear_model import LinearRegression

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ARTIFACT = ROOT / "wiki" / "hiv_tobit_refit_2026-09-22.json"
BAR = ROOT / "wiki" / "hiv_tobit_refit_acceptance_bar.json"


def _mod():
    spec = importlib.util.spec_from_file_location(
        "hiv_tobit_refit", ROOT / "scripts" / "hiv_tobit_refit.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _synthetic(censor, seed=0, n=1500):
    rng = np.random.default_rng(seed)
    true = np.array([0.9, 0.5, 0.25, 0.1, 0.0, -0.3])
    X = (rng.random((n, len(true))) < 0.25).astype(float)
    y = np.minimum(0.2 + X @ true + rng.normal(0, 0.5, n), censor)
    return X, y, true


# --- THE implementation control: recovery of known coefficients -----------------------------------

def test_tobit_recovers_known_coefficients_better_than_ols_under_censoring():
    """THE guard. At ~27% censoring OLS attenuates a true 0.9 to ~0.56; Tobit recovers ~0.89."""
    m = _mod()
    X, y, true = _synthetic(1.0)
    ols = LinearRegression().fit(X, y).coef_
    tob = m.fit_tobit(X, y, censor=1.0)["coef"]
    assert np.abs(tob - true).mean() < np.abs(ols - true).mean() / 3


def test_attenuation_is_worst_for_the_largest_coefficient():
    """The audit's stated MECHANISM, verified where the truth is known: right-censoring biases the
    STRONGEST effects most -- which is why a coefficient THRESHOLD can exclude true majors."""
    X, y, true = _synthetic(1.0)
    ols = LinearRegression().fit(X, y).coef_
    big = abs(ols[0] - true[0]) / abs(true[0])      # true 0.9
    small = abs(ols[3] - true[3]) / abs(true[3])    # true 0.1
    assert big > small * 2


def test_tobit_reduces_to_ols_EXACTLY_when_nothing_is_censored():
    """Tobit has nothing to correct at zero censoring. Observed on the real data too: tenofovir and
    bictegravir, the two drugs with literally 0% censoring, both drift exactly 0.0000."""
    m = _mod()
    X, y, _ = _synthetic(99.0)                       # a cap nothing reaches
    ols = LinearRegression().fit(X, y).coef_
    tob = m.fit_tobit(X, y, censor=99.0)["coef"]
    assert np.abs(tob - ols).max() < 1e-4


def test_the_correction_grows_with_the_censored_fraction():
    """Monotonicity is why an ABSOLUTE drift threshold was the wrong control: the quantity it polices
    is continuous in the censored fraction, so any fixed cut slices through a smooth gradient."""
    m = _mod()
    drifts = []
    for c in (99.0, 2.0, 1.0):
        X, y, _ = _synthetic(c)
        ols = LinearRegression().fit(X, y).coef_
        drifts.append(float(np.abs(m.fit_tobit(X, y, censor=c)["coef"] - ols).max()))
    assert drifts[0] < drifts[1] < drifts[2]


# --- design + constants ---------------------------------------------------------------------------

def test_the_threshold_and_min_carriers_are_imported_not_restated():
    """OLS-vs-Tobit must be compared on the design the shipped catalogs actually use."""
    m = _mod()
    from scripts.hiv_nrti_mutant_catalog import MIN_CARRIERS, RESIST_COEF_MIN
    assert m.RESIST_COEF_MIN is RESIST_COEF_MIN and m.MIN_CARRIERS is MIN_CARRIERS
    assert m.CENSOR_LOG10 == pytest.approx(math.log10(100.0))


def test_a_mutant_below_min_carriers_never_enters_the_design():
    m = _mod()
    records = [(["RARE"], 10.0)] * 2 + [([], 1.0)] * 60
    _, _, cands = m.build_design(records)
    assert cands is None or "RARE" not in cands


# --- the frozen verdict rule ----------------------------------------------------------------------

def _d(drug, cens, drift, changed=False, conv=True):
    return {"status": "ok", "drug": drug, "censored_fraction": cens, "max_abs_coef_delta": drift,
            "membership_changed": changed, "tobit_converged": conv, "mean_coef_delta": 0.0,
            "gained_under_tobit": [], "lost_under_tobit": []}


def test_the_control_vetoes_before_anything_is_interpreted():
    m = _mod()
    out = m.verdict([_d("tenofovir", 0.0, 0.5), _d("lamivudine", 0.45, 1.0, changed=True)])
    assert out["verdict"] == m.SUSPECT


def test_a_clean_control_with_a_membership_change_is_material():
    m = _mod()
    out = m.verdict([_d("tenofovir", 0.0, 0.0), _d("lamivudine", 0.45, 1.0, changed=True)])
    assert out["verdict"] == m.MATERIAL


def test_a_clean_control_with_no_membership_change_is_negligible():
    m = _mod()
    out = m.verdict([_d("tenofovir", 0.0, 0.0), _d("lamivudine", 0.45, 1.0)])
    assert out["verdict"] == m.NEGLIGIBLE


def test_nonconvergence_on_a_high_censoring_drug_is_indeterminate_not_negligible():
    """A failed fit is a different fact from 'no effect' and must not be reported as one."""
    m = _mod()
    out = m.verdict([_d("tenofovir", 0.0, 0.0), _d("lamivudine", 0.45, 1.0, conv=False)])
    assert out["verdict"] == m.INDETERMINATE


# --- the committed run ----------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated (HIV data gitignored)")
def test_the_committed_verdict_is_re_derivable_and_is_the_FROZEN_one():
    """The run is recorded as IMPLEMENTATION_SUSPECT because that is what the frozen rule returns.
    It was NOT rewritten after the control was diagnosed as mis-specified -- re-sizing a bar after
    seeing the result is an authority call, and this session already overrode two."""
    m = _mod()
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["verdict"] == m.verdict(d["results"])["verdict"] == m.SUSPECT


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_two_zero_censoring_drugs_drift_exactly_zero():
    """The evidence that the implementation is sound even though the control fired."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    zero = [r for r in d["results"] if r.get("status") == "ok" and r["censored_fraction"] == 0.0]
    assert {r["drug"] for r in zero} == {"tenofovir", "bictegravir"}
    assert all(r["max_abs_coef_delta"] == 0.0 for r in zero)


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_knife_edge_flips_are_separated_from_decisive_ones():
    """didanosine K70N moves 0.174 -> 0.177 at 0.1% censoring: a coin toss, not a censoring effect.
    Unseparated it would read exactly like lamivudine K70T going -0.202 -> +1.312."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ddi = next(r for r in d["results"] if r["drug"] == "didanosine")
    assert ddi["n_flips_knife_edge"] == 1 and ddi["n_flips_decisive"] == 0
    lam = next(r for r in d["results"] if r["drug"] == "lamivudine")
    assert lam["n_flips_decisive"] >= 3


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_the_bar_declares_it_produces_no_new_gain_figures():
    b = json.loads(BAR.read_text(encoding="utf-8"))
    assert any("does NOT re-run the 5-fold CV" in s for s in b["declared_in_advance"])
