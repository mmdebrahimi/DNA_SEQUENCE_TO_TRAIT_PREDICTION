"""Validate the horse base coat-colour rule on an INDEPENDENTLY-OBSERVED dataset (non-circular).

THE anti-circular gate (load-bearing): the colour label MUST be independently observed, NOT assigned from
the MC1R/ASIP genotype. The one clean open CSV (Dryad 10.5061/dryad.3q111) FAILS this twice over -- its
colour was "determined based on the genotypes" (circular) AND the file is auth-gated (401 headless). Published
observed-colour x genotype contingencies (Rieder 2001; the 709-horse Synergy study; Noma CIE-Lab) are the
valid non-circular sources but live in paywalled/PDF supplementary tables (not headlessly fetchable). So this
harness scores a USER-PROVIDED observed-colour dataset and otherwise honestly reports the data wall -- never a
synthetic/circular number.

Input TSV (--data): columns `mc1r`, `asip`, `observed_colour` (chestnut/bay/black) per horse. The colour
column MUST be independently observed (visually classified), not genotype-derived.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.horse_coat import call_horse_base_colour  # noqa: E402

_BASE = {"chestnut", "bay", "black"}


def _norm_colour(s: str) -> str | None:
    t = (s or "").strip().lower()
    if t in _BASE:
        return t
    if t in ("sorrel", "red"):        # sorrel = US term for chestnut
        return "chestnut"
    return None                        # dark-bay/brown/grey/etc excluded from the base-3 test


def score_rows(rows: list[dict]) -> dict:
    """rows: dicts with mc1r/asip/observed_colour. Returns concordance of the rule vs OBSERVED colour."""
    n = correct = n_indet = n_excluded = 0
    confusion: dict[str, dict[str, int]] = {}
    for r in rows:
        obs = _norm_colour(r.get("observed_colour", ""))
        if obs is None:
            n_excluded += 1
            continue
        pred = call_horse_base_colour(r.get("mc1r", ""), r.get("asip", ""))
        if pred not in _BASE:
            n_indet += 1
            continue
        n += 1
        confusion.setdefault(obs, {}).setdefault(pred, 0)
        confusion[obs][pred] += 1
        if pred == obs:
            correct += 1
    return {"n_scored": n, "n_correct": correct,
            "concordance": round(correct / n, 3) if n else None,
            "n_indeterminate": n_indet, "n_excluded_nonbase": n_excluded,
            "confusion_observed_to_predicted": confusion}


def run(data_path: Path | None) -> dict:
    base = {
        "schema": "horse-coat-validation-v1", "date": _date.today().isoformat(),
        "rule": "MC1R E/e epistatic to ASIP A/a (Rieder 2001 / UC Davis VGL): e/e->chestnut, E_A_->bay, E_aa->black",
        "tier": "DEPLOYED-RULE INTEGRATION on a non-human OBSERVED label (clears the label-first gate IF "
                "colour is independently observed)",
        "anti_circular_gate": "colour MUST be independently observed, not genotype-derived",
    }
    if data_path is None or not data_path.exists():
        return {**base, "status": "VALIDATION_DATA_WALL",
                "note": ("no headlessly-fetchable INDEPENDENT-colour per-individual dataset: Dryad "
                         "10.5061/dryad.3q111 is genotype-DERIVED (circular) AND auth-gated (401); published "
                         "observed contingencies (Rieder 2001 / 709-horse Synergy / Noma CIE-Lab) are "
                         "PDF/paywalled. Provide --data <TSV: mc1r,asip,observed_colour> to score."),
                "data_path": str(data_path) if data_path else None}
    rows = list(csv.DictReader(data_path.open(encoding="utf-8"), delimiter="\t"))
    res = score_rows(rows)
    return {**base, "status": "SCORED" if res["n_scored"] else "NO_ROWS_SCORED", **res,
            "caveats": [
                "colour must be INDEPENDENTLY OBSERVED (caller asserts this) -- else the score is circular",
                "base-3 colours only (dark-bay/brown/grey shade-modifiers excluded)",
                "this RE-CONFIRMS a deployed rule (VGL/Rieder); near-100% concordance is expected -- the "
                "informative signal is any DISCORDANT cell (a rule-breaker)",
            ]}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=None, help="TSV: mc1r, asip, observed_colour (observed!)")
    a = ap.parse_args(argv)
    res = run(a.data)
    if res.get("status") == "SCORED":
        out = REPO / "wiki" / f"horse_coat_validation_{_date.today().isoformat()}.json"
        out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {out}]")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
