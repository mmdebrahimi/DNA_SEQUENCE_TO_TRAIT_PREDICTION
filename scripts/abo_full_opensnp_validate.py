"""Validate the full ABO O/A/B/AB decoder against the staged OpenSNP self-report labels.

Reuses the J3-ABO substrate (`data/j3_abo/j3_abo_substrate.json` — 395 samples with abo_group label +
rs8176719/rs8176746/rs8176747), so NO zip re-scan. Emits 4-way concordance + confusion + the O-vs-non-O
sub-metric. Label tier = self-reported ABO (near-independent, non-circular, ~15% noisy).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.abo_full import call_abo_full  # noqa: E402


def run(substrate_json: Path) -> dict:
    if not substrate_json.exists():
        return {"status": "NO_SUBSTRATE", "hint": "run scripts/abo_opensnp_ingest.py first"}
    samples = json.loads(substrate_json.read_text(encoding="utf-8"))["samples"]
    n_scored = n_correct = 0
    confusion: Counter = Counter()          # (label, pred)
    n_indet = n_no_label = 0
    o_tp = o_fp = o_tn = o_fn = 0
    for s in samples:
        label = s.get("abo_group")          # O/A/B/AB from self-report
        if label not in ("O", "A", "B", "AB"):
            n_no_label += 1
            continue
        pred = call_abo_full(s.get("rs8176719"), s.get("rs8176746"), s.get("rs8176747"))
        if pred == "Indeterminate":
            n_indet += 1
            continue
        n_scored += 1
        n_correct += int(pred == label)
        confusion[(label, pred)] += 1
        # O-vs-non-O sub-metric (O positive)
        li, pi = (label == "O"), (pred == "O")
        o_tp += int(li and pi); o_fp += int(not li and pi)
        o_tn += int(not li and not pi); o_fn += int(li and not pi)
    return {
        "status": "SCORED" if n_scored else "NO_SCORED",
        "schema": "abo-full-opensnp-validation-v1", "date": _date.today().isoformat(),
        "source": "OpenSNP archive dump 2017-12-08; staged substrate data/j3_abo/j3_abo_substrate.json",
        "rule": "3-variant deterministic ABO (rs8176719 + rs8176746 + rs8176747); UK-Biobank-style method",
        "label_tier": "self-reported ABO (near-independent, non-circular, ~15% noisy)",
        "n_scored": n_scored, "n_correct": n_correct,
        "accuracy_4way": round(n_correct / n_scored, 3) if n_scored else None,
        "n_indeterminate": n_indet, "n_no_label": n_no_label,
        "o_vs_nonO": {"TP": o_tp, "FP": o_fp, "TN": o_tn, "FN": o_fn,
                      "accuracy": round((o_tp + o_tn) / n_scored, 3) if n_scored else None},
        "confusion_label_pred": {f"{k[0]}->{k[1]}": v for k, v in sorted(confusion.items())},
        "caveats": ["self-report ~15% erroneous", "non-deletional O / A2 / cis-AB not captured",
                    "the A-vs-B extension of the O-vs-non-O cell (rs8176746/47 sourced, not fabricated)"],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--substrate", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_substrate.json")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / f"abo_full_opensnp_validation_{_date.today().isoformat()}.json")
    a = ap.parse_args(argv)
    res = run(a.substrate)
    if res.get("status") == "SCORED":
        a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"[wrote {a.out}]")
    print(json.dumps(res, indent=2))
    return 0 if res.get("status") == "SCORED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
