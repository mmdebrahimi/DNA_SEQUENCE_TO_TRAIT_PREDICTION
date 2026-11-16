#!/usr/bin/env python
"""Reconcile a pulled chained session against the checkpoint curves it left on disk.

WHAT THIS IS FOR. The bounded chain's whole claim is that progress accumulates across sessions, and
the only evidence for that is the per-epoch curve files. A run record is a SUMMARY of those files
written by the same process that wrote them, so reading the record alone cannot detect a session that
reported more progress than it made. This recomputes the totals from the curves INDEPENDENTLY and
refuses when they disagree.

IT ALSO ANSWERS THE ONE QUESTION THAT DECIDES WHETHER TO KEEP CHAINING: did this link ADD epochs? A
session that seeds correctly and adds nothing is a stall, not progress, and must stop the chain --
otherwise each further link costs a GPU session and a quota hour for no epochs. Pass `--against` a
previous link's reconciliation JSON to get that comparison as a verdict rather than an impression.

A BOUNDED STOP IS NOT A FAILURE and this script will not call it one. `run()` orders its phases
budget -> nulls -> diagnostic extension precisely so that a stop after the first two is gradeable,
and the scorer returns INDETERMINATE_BUDGET_EXHAUSTED until every seed reaches the frozen budget. So
the early links of a chain produce no verdict BY DESIGN; `gradeable` below reports whether the prefix
is complete, and `status` reports how the session ended. They are different questions.

  uv run python scripts/kaggle/reconcile_finetune_session.py --pulled ./out [--against prev.json]
      [--out wiki/finetune_arm_chain_link1.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scripts.finetune_arm_run import CKPT_NAME_RE  # noqa: E402  the SINGLE definition


def find_one(pulled: Path, name: str) -> Path | None:
    hits = sorted(p for p in pulled.rglob(name) if p.is_file())
    if len(hits) > 1:
        raise SystemExit("REFUSING: %d files named %r under %s (%s). Pull one session at a time, or "
                         "the reconciliation silently mixes two runs."
                         % (len(hits), name, pulled, ", ".join(str(h) for h in hits)))
    return hits[0] if hits else None


def read_curves(pulled: Path) -> dict[str, dict]:
    """Every accepted curve under the pulled tree, keyed by filename.

    Uses the same `CKPT_NAME_RE` the kernel and the uploader use, so a file none of them would accept
    cannot inflate the count here either. `smoke_*` is excluded, which matters: the rate measurement
    writes a curve shaped exactly like a real one.
    """
    out: dict[str, dict] = {}
    for p in sorted(pulled.rglob("*_curve.json")):
        if not (p.is_file() and CKPT_NAME_RE.match(p.name)):
            continue
        if p.name in out:
            raise SystemExit("REFUSING: two curves named %r under %s; totals would double-count."
                             % (p.name, pulled))
        try:
            out[p.name] = json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            raise SystemExit("REFUSING: %s is not readable JSON (%s). A truncated curve must not be "
                             "silently counted as zero epochs." % (p, e))
    return out


def reconcile(pulled: Path) -> dict:
    rec_path = find_one(pulled, "finetune_arm_run_kaggle.json")
    curves = read_curves(pulled)

    seeds = {k: v for k, v in curves.items() if k.startswith("seed_")}
    nulls = {k: v for k, v in curves.items() if k.startswith("null_")}

    def epochs(d: dict) -> int:
        c = d.get("val_spearman_by_epoch")
        return len(c) if isinstance(c, list) else int(d.get("epochs_done") or 0)

    from_curves = {
        "n_seed_curves": len(seeds),
        "n_null_curves": len(nulls),
        "epochs_by_seed": {k: epochs(v) for k, v in sorted(seeds.items())},
        "epochs_by_null": {k: epochs(v) for k, v in sorted(nulls.items())},
        "total_epochs": sum(epochs(v) for v in curves.values()),
    }

    out = {
        "schema": "finetune-arm-chain-link-v1",
        "pulled_from": str(pulled),
        "from_curves": from_curves,
        "run_artifact_present": rec_path is not None,
    }

    if rec_path is None:
        out["status"] = None
        out["verdict"] = "NO_RUN_ARTIFACT"
        out["why"] = (
            "No finetune_arm_run_kaggle.json under the pulled tree. The session did not reach its "
            "artifact write. Per-epoch curves may still exist and are counted above -- that is the "
            "point of writing them inside the epoch loop -- but this link produced no record, so "
            "nothing can be scored from it and the cause must be read from the kernel log.")
        return out

    rec = json.loads(rec_path.read_text(encoding="utf-8"))
    out["status"] = rec.get("status")
    out["why_status"] = rec.get("why_status")
    out["horizon_divergence_recorded"] = "horizon_divergence" in rec
    if out["horizon_divergence_recorded"]:
        out["horizon_epochs_supplied"] = rec["horizon_divergence"].get("horizon_epochs_supplied")

    # BOUNDEDNESS IS DERIVED, because the obvious field does not exist at this level. `run()`'s record
    # has no top-level `bounded_by_wall_clock`; that field is set by `train_one_seed` on a PER-SEED
    # record. The plan's verification criterion named the top-level field and was wrong about where it
    # lives. The property itself is recorded twice, so both are read and both are reported.
    phases = rec.get("phases") or []
    bounded_phases = [p.get("phase") for p in phases
                      if isinstance(p, dict) and p.get("bounded_here")]
    out["bounded_signals"] = {
        "status_is_bounded": rec.get("status") == "BOUNDED_BY_WALL_CLOCK",
        "phases_reporting_bounded_here": bounded_phases,
        "top_level_bounded_field_present": "bounded_by_wall_clock" in rec,
    }
    out["bounded"] = bool(out["bounded_signals"]["status_is_bounded"] or bounded_phases)

    # THE INDEPENDENT CHECK: the record's own per-seed curves against the files on disk.
    # The key is `per_seed`. An earlier version read `seeds`, which does not exist, so it compared
    # NOTHING and reported success -- the vacuity guard below exists because of that.
    rec_seeds = rec.get("per_seed")
    if not isinstance(rec_seeds, dict):
        rec_seeds = {}
    mismatches = []
    compared = 0
    for key, srec in sorted(rec_seeds.items()):
        if not isinstance(srec, dict):
            continue
        c = srec.get("val_spearman_by_epoch")
        if not isinstance(c, list):
            continue
        fname = "seed_%s_curve.json" % key
        disk = seeds.get(fname)
        if disk is None:
            mismatches.append("%s is in the record but has no curve file on disk" % fname)
            continue
        compared += 1
        if epochs(disk) != len(c):
            mismatches.append("%s: record says %d epoch(s), the file on disk has %d"
                              % (fname, len(c), epochs(disk)))

    # VACUITY IS CHECKED FIRST, and the order is load-bearing. A record whose per-seed section this
    # script cannot read ALSO leaves every disk curve "unmentioned", so the direction check below
    # would fire and report RECORD_DISAGREES_WITH_DISK -- the right refusal with the WRONG
    # diagnosis. The two need opposite responses: re-read the record's schema, versus distrust this
    # link's data. Applicability before measurement.
    out["n_record_seeds_compared"] = compared
    if compared == 0 and seeds:
        out["record_vs_disk_mismatches"] = mismatches
        out["reconciles"] = False
        out["verdict"] = "RECONCILIATION_VACUOUS"
        out["why"] = (
            "There are %d seed curve(s) on disk and the record yielded NOTHING to compare them "
            "against, so this reconciliation checked nothing. Reporting success here would "
            "manufacture confidence in the one claim the chain rests on. The cause is a record "
            "schema this script does not understand -- read its keys before trusting any link."
            % len(seeds))
        return out

    # The other direction: a curve on disk the record never mentions.
    rec_names = {"seed_%s_curve.json" % k for k in rec_seeds}
    for fname in sorted(set(seeds) - rec_names):
        mismatches.append("%s is on disk but the record does not mention that seed" % fname)

    out["record_vs_disk_mismatches"] = mismatches
    out["reconciles"] = (not mismatches) and compared > 0

    if mismatches:
        out["verdict"] = "RECORD_DISAGREES_WITH_DISK"
        out["why"] = ("The run record and the checkpoint files do not agree. The record is a summary "
                      "written by the process that wrote the curves, so a disagreement means one of "
                      "the two is not describing this session. Do not chain from this link.")
    else:
        # Set here, not in main(), so EVERY caller of reconcile() gets a complete verdict. An
        # earlier version left the happy path without one and only main() filled it in, which made
        # reconcile() raise KeyError for any other consumer.
        out["verdict"] = "LINK_OK"
        out["why"] = "The record reconciles with the curves on disk."
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pulled", required=True, help="directory a `pull` wrote")
    ap.add_argument("--against", default=None,
                    help="a previous link's reconciliation JSON, to test that epochs INCREASED")
    ap.add_argument("--out", default=None, help="write the reconciliation here")
    args = ap.parse_args(argv)

    pulled = Path(args.pulled)
    if not pulled.is_dir():
        ap.error("--pulled %s is not a directory" % pulled)

    out = reconcile(pulled)

    if args.against:
        prev = json.loads(Path(args.against).read_text(encoding="utf-8"))
        before = int(prev.get("from_curves", {}).get("total_epochs") or 0)
        now = int(out["from_curves"]["total_epochs"])
        out["previous_total_epochs"] = before
        out["epochs_added_by_this_link"] = now - before
        if now <= before:
            out["verdict"] = "STALLED_NO_EPOCHS_ADDED"
            out["why"] = (
                "This link added %d epoch(s) (was %d, now %d). A session that seeds correctly and "
                "adds nothing is a stall, not progress: each further link costs a GPU session and "
                "quota hours for no epochs, so the chain must STOP and the cause be diagnosed "
                "before another push." % (now - before, before, now))

    out.setdefault("verdict", "LINK_OK")
    out.setdefault("why", "The record reconciles with the curves on disk.")

    print(json.dumps(out, indent=2, sort_keys=True))
    if args.out:
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("WROTE %s" % p, file=sys.stderr)

    # Exit codes so a chaining loop can branch without parsing prose.
    return {"LINK_OK": 0, "NO_RUN_ARTIFACT": 2,
            "RECORD_DISAGREES_WITH_DISK": 3, "STALLED_NO_EPOCHS_ADDED": 4}.get(out["verdict"], 1)


if __name__ == "__main__":
    raise SystemExit(main())
