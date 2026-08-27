"""Score the v3 staleness run against the PRE-REGISTERED predictions and the hand adjudications.

Two scores, deliberately separate, because they answer different questions:

  1. AGGREGATE -- did the flag count / true positives / false positives clear the pre-committed bar?
     This is the number that says whether the tool got better.

  2. TARGETED -- did the three predicted flips actually happen, and did the must-hold true positive
     survive? This is the number that says whether it got better FOR THE REASON I DIAGNOSED.

The second is the one that can embarrass the diagnosis, which is why it exists. A run that clears the
aggregate bar while missing the targeted flips has improved by some other mechanism -- luck, a prompt
length change, a generation-order effect -- and the three-instance root-cause story would be wrong.

Ground truth for "true positive" is `staleness_adjudication.ADJUDICATIONS`, checked by hand against the
artifacts BEFORE v3 existed. It is not re-derived from this run.

Run: uv run python scripts/score_v3.py <results.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))


def _truth() -> dict[str, str]:
    """artifact -> the hand verdict, from the committed adjudication record."""
    from staleness_adjudication import ADJUDICATIONS
    # negative_results_map appears TWICE with opposite verdicts (a real stale bullet AND an accurate
    # historical one). Keying by artifact alone would silently drop one -- the shared-key overwrite trap.
    # Keep the TP so a lost real catch can never be hidden by the FP entry.
    out: dict[str, str] = {}
    for a in ADJUDICATIONS:
        if out.get(a.artifact) == "true_positive":
            continue
        out[a.artifact] = a.verdict
    return out


def score(results_path: Path) -> dict:
    from kaggle_staleness_auditor import parse_verdict
    from staleness_v3_preregistration import BAR, PREDICTIONS

    raw = json.load(open(results_path, encoding="utf-8"))
    verdicts = {r["artifact"]: parse_verdict(r["raw"]) for r in raw}
    flagged = {a for a, p in verdicts.items() if p["verdict"] == "stale"}
    truth = _truth()

    tp = sorted(a for a in flagged if truth.get(a) == "true_positive")
    fp = sorted(a for a in flagged if truth.get(a) != "true_positive")

    targeted = []
    for pr in PREDICTIONS:
        got = verdicts.get(pr.artifact, {}).get("verdict", "MISSING")
        targeted.append({"artifact": pr.artifact, "kind": pr.kind, "predicted": pr.predict,
                         "got": got, "hit": got == pr.predict})

    aggregate_pass = (len(flagged) <= BAR["max_flags"]
                      and len(tp) >= BAR["min_true_positives"]
                      and len(fp) <= BAR["max_false_positives"])
    core = [t for t in targeted if t["kind"] in {"finding", "correction", "historical", "must-hold-TP"}]
    return {"n_scored": len(raw), "flags": sorted(flagged), "true_positives": tp, "false_positives": fp,
            "aggregate_pass": aggregate_pass, "targeted": targeted,
            "targeted_hits": sum(t["hit"] for t in core), "targeted_total": len(core),
            "unparseable": sum(1 for p in verdicts.values() if not p["parse_ok"])}


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    r = score(Path(sys.argv[1]))
    print(f"scored {r['n_scored']} items, {r['unparseable']} unparseable\n")
    print(f"flags: {len(r['flags'])}   true positives: {len(r['true_positives'])}   "
          f"false positives: {len(r['false_positives'])}")
    print(f"AGGREGATE BAR: {'PASS' if r['aggregate_pass'] else 'FAIL'}")
    print(f"TARGETED (the diagnosis): {r['targeted_hits']}/{r['targeted_total']}\n")
    for t in r["targeted"]:
        mark = "HIT " if t["hit"] else "MISS"
        print(f"  [{mark}] {t['kind']:14} {t['artifact']}")
        print(f"         predicted {t['predicted']}, got {t['got']}")
    if r["flags"]:
        print("\nremaining flags:")
        for f in r["flags"]:
            print(f"  {'TP' if f in r['true_positives'] else 'FP'}  {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
