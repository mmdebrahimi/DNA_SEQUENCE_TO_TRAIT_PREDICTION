"""HIV supervised blind-spot COMPLEMENT scorer — NNRTI + PI + INSTI (2026-07-12).

Ships the learned per-mutation weights (built by `scripts/build_hiv_complement_model.py` from the free
Stanford PhenoSense fold-change label) as OFFLINE scorers — no sklearn, no training data at inference.

WHAT IT IS: a per-class RANKING complement to the deployed HIV catalog. For an isolate the catalog calls
susceptible (its blind spot), `blind_spot_risk(observed, drug_class=...)` returns a probability the genotype
is nonetheless resistant. Deployability is proven per class (leave-one-study-out blind-spot AUROC: NNRTI
0.81, PI 0.89 on LPV, INSTI 0.89 on RAL; see `wiki/hiv_supervised_deployability_2026-07-12.json` +
`wiki/hiv_supervised_targetsite_panel_2026-07-12.json`).

WHAT IT IS NOT: a hard R/S rule and NOT a catalog replacement. The catalog fold-in (turning these weights
into binary rules) was tested and REJECTED (net -0.006 balanced accuracy); the VALUE is the continuous
weighting. Supervised (needed the free label to train) and in-distribution to the Stanford knowledge base.
The frozen decoder surface + `hiv_amr.py` catalog are untouched.

Usage:
    from dna_decode.data import hiv_supervised_complement as C
    C.blind_spot_risk({"103N", "179D"})                    # NNRTI (default)
    C.blind_spot_risk({"82F", "54V"}, drug_class="PI")     # protease
    C.blind_spot_risk({"148H"}, drug_class="INSTI")        # integrase
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

_REF = Path(__file__).resolve().parent.parent.parent / "data" / "hiv_ref"
_MODELS = {
    "NNRTI": _REF / "hiv_nnrti_supervised_complement.json",
    "PI": _REF / "hiv_pi_supervised_complement.json",
    "INSTI": _REF / "hiv_insti_supervised_complement.json",
}
SUPPORTED_CLASSES = tuple(_MODELS)
DEFAULT_THRESHOLD = 0.5
_SCHEMA = "hiv-supervised-complement-v1"


class ComplementUnavailable(RuntimeError):
    """Raised when the serialized complement model JSON is absent or unknown class."""


@lru_cache(maxsize=8)
def _model(drug_class: str = "NNRTI", path: str | None = None):
    if path:
        p = Path(path)
    else:
        if drug_class not in _MODELS:
            raise ComplementUnavailable(f"unknown drug_class {drug_class!r}; expected one of {SUPPORTED_CLASSES}")
        p = _MODELS[drug_class]
    if not p.exists():
        raise ComplementUnavailable(f"complement model not found at {p}; run build_hiv_complement_model.py")
    m = json.loads(p.read_text(encoding="utf-8"))
    if m.get("schema") != _SCHEMA:
        raise ComplementUnavailable(f"unexpected complement schema: {m.get('schema')}")
    return m


def _norm(observed) -> set[str]:
    """Accept {'103N'} or {'K103N'} or (pos, aa) tuples -> canonical '<pos><aa>' tokens (matching feature_key)."""
    out = set()
    for o in observed:
        if isinstance(o, (tuple, list)) and len(o) == 2:
            out.add(f"{int(o[0])}{o[1]}")
            continue
        s = str(o).strip()
        if s and s[0].isalpha() and len(s) > 1 and s[1].isdigit():   # strip a leading WT letter (K103N -> 103N)
            s = s[1:]
        out.add(s)
    return out


def blind_spot_risk(observed, drug_class: str = "NNRTI", model_path: str | None = None) -> float:
    """P(resistant) for an HIV genotype from the class's supervised complement. `observed` = the isolate's
    non-WT residues as '<pos><aa>' (e.g. {'103N'} for NNRTI, {'82F'} for PI); unknown tokens contribute 0."""
    m = _model(drug_class, model_path)
    w = m["weights"]
    z = float(m["intercept"]) + sum(w.get(tok, 0.0) for tok in _norm(observed))
    return 1.0 / (1.0 + math.exp(-z))


def is_flagged(observed, drug_class: str = "NNRTI", threshold: float = DEFAULT_THRESHOLD,
               model_path: str | None = None) -> bool:
    """True iff the complement flags this (catalog-negative) genotype as likely-resistant at `threshold`."""
    return blind_spot_risk(observed, drug_class, model_path) >= threshold


def model_info(drug_class: str = "NNRTI", model_path: str | None = None) -> dict:
    m = _model(drug_class, model_path)
    return {k: m[k] for k in ("schema", "drug_class", "gene", "date", "drug_trained", "cutoff_fold",
                              "n_train", "n_features", "deployability", "honest_scope") if k in m}
