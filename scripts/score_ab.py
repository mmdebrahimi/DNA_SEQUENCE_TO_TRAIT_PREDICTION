"""Score the v2-vs-v3 A/B on the PRE-REPAIR snapshot, where the prompt is the only difference.

WHY THIS EXISTS. The first v2-vs-v3 comparison was uninterpretable because three things changed at once:
the prompt, the excerpt cap (my OOM fix), and the corpus text itself -- I had repaired the very claims the
v2 run found, which removed them from the test set. This run fixes all three:

  * corpus      -- reconstructed from git at the commit BEFORE the repairs (both known stale claims intact)
  * excerpt     -- 3000 chars for BOTH arms
  * prompt      -- the only variable

GROUND TRUTH on this snapshot. The pre-repair CLAUDE.md contains exactly two genuinely stale claims, both
hand-verified at the time and both since fixed in the live file:

  1. "8 reusable rejection GATES ... G1-G8"  -- the map has ten (G9/G10 added 2026-08-26)
  2. "8 unit tests at tests/test_models_cache.py" -- the file has 42

Claim 2 is cited in a bullet that names TWO files, so it yields two pairs (test_models_cache.py and
probe_nt_cache.py). Both carry the same stale count, so both are true positives -- counting only one would
understate recall. Everything else is treated as not-stale, which is an assumption inherited from the hand
adjudication of all 110, not a fresh judgement.

Run: uv run python scripts/score_ab.py <v2_results.json> <v3_results.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

SNAPSHOT = ROOT / "wiki" / "staleness_eval_snapshot_prerepair.json"

# The two stale claim texts, matched against the snapshot's claim field. Matching on TEXT rather than
# hardcoding item ids means the ground truth is derived from the corpus rather than asserted beside it --
# the hardcoded-list trap this project has hit three times.
STALE_MARKERS = ("8 reusable rejection GATES", "8 unit tests at")


def truth_items() -> set[str]:
    """item_ids of the genuinely-stale pairs, DERIVED from the snapshot text."""
    items = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    return {it["item_id"] for it in items if any(m in it["claim"] for m in STALE_MARKERS)}


def score_one(results_path: Path) -> dict:
    from kaggle_staleness_auditor import parse_verdict
    raw = json.load(open(results_path, encoding="utf-8"))
    truth = truth_items()
    per = {r["item_id"]: parse_verdict(r["raw"]) for r in raw}
    flagged = {i for i, p in per.items() if p["verdict"] == "stale"}
    tp = sorted(flagged & truth)
    fn = sorted(truth - flagged)
    fp = sorted(flagged - truth)
    n_neg = len(per) - len(truth)
    return {
        "n": len(raw), "unparseable": sum(1 for p in per.values() if not p["parse_ok"]),
        "flags": len(flagged), "tp": len(tp), "fn": len(fn), "fp": len(fp),
        "recall": len(tp) / len(truth) if truth else 0.0,
        "specificity": (n_neg - len(fp)) / n_neg if n_neg else 0.0,
        "precision": len(tp) / len(flagged) if flagged else 0.0,
        "missed": fn, "false_flags": fp,
    }


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    a, b = score_one(Path(sys.argv[1])), score_one(Path(sys.argv[2]))
    n_truth = len(truth_items())
    print(f"PRE-REPAIR snapshot: {a['n']} pairs, {n_truth} genuinely-stale claims\n")
    print(f"{'':14}{'v2':>10}{'v3':>10}")
    for k, label in (("flags", "flags"), ("tp", "true pos"), ("fn", "MISSED"), ("fp", "false pos"),
                     ("recall", "recall"), ("precision", "precision"), ("specificity", "specificity"),
                     ("unparseable", "unparseable")):
        fa, fb = a[k], b[k]
        fmt = (lambda v: f"{v:.3f}") if isinstance(fa, float) else (lambda v: f"{v}")
        print(f"  {label:12}{fmt(fa):>10}{fmt(fb):>10}")
    print()
    for name, r in (("v2", a), ("v3", b)):
        if r["missed"]:
            print(f"{name} MISSED (false negatives): {r['missed']}")
    print("\nOnly the prompt differed. Corpus, excerpt length, model and decoding were identical.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
