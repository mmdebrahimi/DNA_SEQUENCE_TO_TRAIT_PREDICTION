"""The ElasticNetCV alpha-count kwarg must track the installed scikit-learn.

The drift this guards
---------------------
scikit-learn folded `n_alphas` into `alphas` (which gained int support in 1.7) and the CV
*estimators* dropped `n_alphas` outright, so `ElasticNetCV(n_alphas=20)` became a
`TypeError` under 1.9.1. It broke `tests/test_tb_mic_calibration.py` through
`scripts/tb_mic_calibration.py` -> `scripts/hiv_epistasis.py`.

Why a guard and not just the fix
--------------------------------
The second call site, `scripts/geuvadis_stage1b_elasticnet_cis.py`, builds its estimator
*inside* a bare `try: ... except Exception:` per-fold. There the same `TypeError` is
SWALLOWED -- every fold fails, nothing is raised, and the script still finishes. That is
the silent-variant-drop shape this repo already documents for its smoke runner, so the
regression would come back invisible rather than red. These tests are what make it loud.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

import pytest
from sklearn.linear_model import ElasticNetCV

from dna_decode.sklearn_compat import _accepts, alpha_count_kwargs

ROOT = Path(__file__).resolve().parents[1]

#: Every project call site that asks for an auto-generated alpha grid.
CALL_SITES = (
    "scripts/hiv_epistasis.py",
    "scripts/geuvadis_stage1b_elasticnet_cis.py",
)


def test_the_kwarg_is_the_one_this_sklearn_actually_accepts() -> None:
    kw = alpha_count_kwargs(20)
    assert list(kw) in (["n_alphas"], ["alphas"]), kw
    assert next(iter(kw.values())) == 20
    accepted = set(inspect.signature(ElasticNetCV.__init__).parameters)
    assert set(kw) <= accepted, (
        f"alpha_count_kwargs returned {set(kw)} but this ElasticNetCV accepts "
        f"{sorted(p for p in accepted if 'alpha' in p)}"
    )


def test_the_estimator_really_constructs_and_fits() -> None:
    """A kwarg that typechecks is not evidence; fit it on real arrays."""
    import numpy as np

    rng = np.random.default_rng(0)
    X = rng.normal(size=(40, 5))
    y = X @ np.array([2.0, 0.0, -1.0, 0.0, 0.5]) + rng.normal(scale=0.1, size=40)
    model = ElasticNetCV(l1_ratio=0.5, cv=3, max_iter=2000, **alpha_count_kwargs(8))
    model.fit(X, y)
    assert len(model.alphas_) == 8, (
        f"asked for 8 alphas, path has {len(model.alphas_)} -- the kwarg was accepted "
        "but did not mean what we intended"
    )


def test_no_call_site_passes_the_dropped_kwarg_literally() -> None:
    """`n_alphas=` as a literal keyword is the regression; find it in the AST."""
    offenders = []
    for rel in CALL_SITES:
        tree = ast.parse((ROOT / rel).read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                for kw in node.keywords:
                    if kw.arg == "n_alphas":
                        offenders.append(f"{rel}:{node.lineno}")
    assert not offenders, (
        "these call sites pass n_alphas= literally, which this scikit-learn rejects: "
        + ", ".join(offenders)
        + ". Use dna_decode.sklearn_compat.alpha_count_kwargs()."
    )


@pytest.mark.parametrize("rel", CALL_SITES)
def test_every_call_site_routes_through_the_shim(rel: str) -> None:
    src = (ROOT / rel).read_text(encoding="utf-8")
    assert "alpha_count_kwargs" in src, (
        f"{rel} constructs an ElasticNetCV but does not use the compat shim; it will "
        "break again on the next scikit-learn that moves this parameter"
    )


def test_the_shim_detects_by_signature_not_by_version_string() -> None:
    """Version parsing encodes a belief; the signature is the fact."""
    src = (ROOT / "dna_decode" / "sklearn_compat.py").read_text(encoding="utf-8")
    assert "inspect.signature" in src, (
        "the shim should ask the class what it accepts, not compare version strings -- "
        "backports and vendored builds make the version a worse oracle than the signature"
    )
    # A matching `"__version__" not in src` assertion was REMOVED: it fails on a legitimate
    # docstring mention of the version, and the behavioural guard below
    # (test_the_detector_is_not_vacuous) is what actually pins signature-based detection.


def test_the_detector_is_not_vacuous() -> None:
    """A detector that answers the same for every parameter proves nothing."""
    assert _accepts("ElasticNetCV", "cv") is True, "cv is a real ElasticNetCV parameter"
    assert _accepts("ElasticNetCV", "definitely_not_a_parameter") is False


def test_a_bad_alpha_count_is_refused_rather_than_passed_through() -> None:
    for bad in (0, -1, 2.5, True, "20"):
        with pytest.raises(ValueError):
            alpha_count_kwargs(bad)  # type: ignore[arg-type]
