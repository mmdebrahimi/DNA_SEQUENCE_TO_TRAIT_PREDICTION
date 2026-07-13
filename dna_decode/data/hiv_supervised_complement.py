"""HIV NNRTI supervised blind-spot COMPLEMENT scorer (2026-07-12).

Ships the learned per-mutation weights (built by `scripts/build_hiv_complement_model.py` from the free
Stanford PhenoSense fold-change label) as an OFFLINE scorer — no sklearn, no training data at inference.

WHAT IT IS: a RANKING complement to the deployed NNRTI major-DRM catalog. For an isolate the catalog calls
susceptible (its blind spot), `blind_spot_risk` returns a probability that the genotype is nonetheless
resistant. Deployability is proven (leave-one-study-out blind-spot AUROC 0.81, patient-grouped 0.824 —
`wiki/hiv_supervised_deployability_2026-07-12.json`).

WHAT IT IS NOT: a hard R/S rule and NOT a catalog replacement. The catalog fold-in (turning these weights
into binary rules) was tested and REJECTED — it trades sensitivity for specificity (net -0.006 balanced
accuracy, `wiki/hiv_catalog_accessory_extension_2026-07-12.json`); the VALUE is the continuous weighting,
which is why it ships as a scorer. It is supervised (needed the free label to train) and in-distribution to
the Stanford knowledge base. The frozen decoder surface + `hiv_amr.py` catalog are untouched.

Usage:
    from dna_decode.data import hiv_supervised_complement as C
    C.blind_spot_risk({"103N", "179D"})       # -> probability in [0, 1]
    C.is_flagged({"103N"})                     # -> bool at the default 0.5 threshold
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

_MODEL_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "hiv_ref" / \
    "hiv_nnrti_supervised_complement.json"
DEFAULT_THRESHOLD = 0.5


class ComplementUnavailable(RuntimeError):
    """Raised when the serialized complement model JSON is absent."""


@lru_cache(maxsize=1)
def _model(path: str | None = None):
    p = Path(path) if path else _MODEL_PATH
    if not p.exists():
        raise ComplementUnavailable(f"complement model not found at {p}; run build_hiv_complement_model.py")
    m = json.loads(p.read_text(encoding="utf-8"))
    if m.get("schema") != "hiv-nnrti-supervised-complement-v1":
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
        # strip a leading WT letter if present (K103N -> 103N); keep pure '103N'
        if s and s[0].isalpha() and len(s) > 1 and s[1].isdigit():
            s = s[1:]
        out.add(s)
    return out


def blind_spot_risk(observed, model_path: str | None = None) -> float:
    """P(resistant) for an NNRTI genotype from the supervised complement. `observed` = the isolate's non-WT
    RT residues as '<pos><aa>' (e.g. {'103N','179D'}); unknown tokens contribute 0 (their weight was 0)."""
    m = _model(model_path)
    w = m["weights"]
    z = float(m["intercept"]) + sum(w.get(tok, 0.0) for tok in _norm(observed))
    return 1.0 / (1.0 + math.exp(-z))


def is_flagged(observed, threshold: float = DEFAULT_THRESHOLD, model_path: str | None = None) -> bool:
    """True iff the complement flags this (catalog-negative) genotype as likely-resistant at `threshold`."""
    return blind_spot_risk(observed, model_path) >= threshold


def model_info(model_path: str | None = None) -> dict:
    m = _model(model_path)
    return {k: m[k] for k in ("schema", "date", "drug_trained", "cutoff_fold", "n_train", "n_features",
                              "deployability", "honest_scope") if k in m}
