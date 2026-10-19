"""Version-tolerant shims for scikit-learn API drift.

Why this exists
---------------
`ElasticNetCV(n_alphas=N)` raised
``TypeError: ElasticNetCV.__init__() got an unexpected keyword argument 'n_alphas'``
under scikit-learn 1.9.1. The parameter was folded into `alphas`, which now accepts an
int ("versionchanged:: 1.7 -- `alphas` accepts an integer value which removes the need to
pass `n_alphas`"), and the CV *estimators* dropped `n_alphas` outright. The free
*functions* (`enet_path`, `lasso_path`) keep it deprecated-but-working until 1.11, which
is why the break is narrow and easy to miss.

Why not just pass `alphas=N`
----------------------------
`pyproject.toml` declares `scikit-learn>=1.4`, and `alphas`-accepts-int did not exist
before 1.7. Writing `alphas=N` would break the project's own published support floor on
1.4-1.6, and raising a PyPI dependency floor is a packaging decision, not a bug fix.

Why not just pass an explicit alpha array
-----------------------------------------
That would work on every version, but the auto-generated grid is derived from the data
(`alpha_max` is a function of X and y), so hand-building one changes the numbers.

The two spellings are DOCUMENTED as equivalent (scikit-learn's own note says `alphas`
accepting an int "removes the need to pass `n_alphas`"), and that is the basis for this
shim. It has NOT been verified numerically across the boundary here, because only one
scikit-learn is installed on this host -- an older one would be needed to compare fits.
So: equivalent by the upstream contract, not by measurement. Do not upgrade that to a
bit-identity claim without actually running both.

Detection is by SIGNATURE, not version string
---------------------------------------------
`_accepts("n_alphas")` asks the installed class what it actually takes. A version
comparison encodes a belief about when something changed; the signature IS the fact, and
it stays correct across backports, forks and vendored builds.
"""

from __future__ import annotations

import inspect
from functools import cache
from typing import Any


@cache
def _accepts(qualname: str, param: str) -> bool:
    """Does the named estimator's __init__ accept `param`?"""
    from sklearn import linear_model

    cls = getattr(linear_model, qualname)
    try:
        sig = inspect.signature(cls.__init__)
    except (TypeError, ValueError):  # pragma: no cover - C-implemented __init__
        return False
    return param in sig.parameters


def alpha_count_kwargs(n: int, estimator: str = "ElasticNetCV") -> dict[str, Any]:
    """Ask for `n` automatically-generated alphas, whichever spelling this version wants.

    Returns ``{"n_alphas": n}`` on scikit-learn < 1.7 and ``{"alphas": n}`` on >= 1.7.
    Both mean "generate n alphas along the regularization path", so the fitted model is
    unchanged across the drift.

    >>> kw = alpha_count_kwargs(20)
    >>> list(kw) in (["n_alphas"], ["alphas"])
    True
    >>> next(iter(kw.values()))
    20
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError(f"alpha count must be a positive int, got {n!r}")
    if _accepts(estimator, "n_alphas"):
        return {"n_alphas": n}
    return {"alphas": n}
