#!/usr/bin/env python
"""Stage a pulled checkpoint directory as a Kaggle dataset, so the next chained session resumes.

THE OTHER END OF A ROUND TRIP, and the reason this is a script rather than three CLI calls: it must
stage EXACTLY the files the kernel will accept. It imports `CKPT_NAME_RE` from
`scripts/finetune_arm_run.py` -- the module whose `_ckpt_path` and `_curve_path` produce the names --
so the two ends cannot disagree. Two copies of a selection rule across a round trip fail by PASSING:
the uploader stages a file the kernel ignores, the kernel seeds zero, and its refusal blames the
mount. That misdiagnosis cost two probe pushes on 2026-11-15 before the rule was single-sourced.

FILES ARE STAGED FLAT, and that is measured rather than tidy. Three things were established on
2026-11-15 by probing a real dataset:

  1. `kaggle datasets create` defaults to `--dir-mode skip`, so a staging folder containing a
     subdirectory uploads NOTHING and yields an empty dataset.
  2. `--dir-mode zip` uploads one archive per directory, and Kaggle then EXTRACTS it and FLATTENS
     the result -- a dataset built from `finetune_checkpoints/seed_99_curve.json` stores exactly
     `seed_99_curve.json` at its root.
  3. The mount path is not `/kaggle/input/<slug>`. A dataset attached via `dataset_sources`, with
     status `ready`, did not appear there; `/kaggle/input` held only `datasets/`.

Staging flat therefore sidesteps (1) and (2) entirely -- there are no directories to skip or zip --
and (3) is handled on the kernel side, which falls back to searching `/kaggle/input` recursively and
identifies checkpoints by NAME rather than by path.

  first link : uv run python scripts/kaggle/upload_finetune_checkpoints.py \\
                   --from ./pulled/finetune_checkpoints --slug ft-ckpt-s1 --new
  later links: ... --slug ft-ckpt-s2 --new          (a FRESH slug per link, see --new below)
  then       : build with --ckpt-dataset <slug> AND push with --dataset <slug>; both are required,
               the first so the kernel knows what to expect, the second so Kaggle attaches it.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from scripts.finetune_arm_run import CKPT_NAME_RE  # noqa: E402  the SINGLE definition


def _api():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api = KaggleApi()
    api.authenticate()
    return api


def classify(src: Path) -> tuple[dict[str, Path], dict[str, list[str]]]:
    """Partition every file under `src` into accepted checkpoints and reasons for exclusion.

    The exclusions are REPORTED rather than silent. A chained run that skips a file the operator
    believed was uploaded is the failure this whole module exists to prevent, and "I uploaded the
    directory" is not evidence about which files were taken.
    """
    taken: dict[str, Path] = {}
    skipped: dict[str, list[str]] = {"smoke (excluded by design)": [],
                                     "unfinished write (.tmp)": [],
                                     "not a checkpoint name": []}
    for p in sorted(src.rglob("*")):
        if not p.is_file():
            continue
        name = p.name
        if CKPT_NAME_RE.match(name):
            if name in taken:
                raise SystemExit(
                    "REFUSING: two files under %s are both named %r (%s and %s). The dataset is "
                    "flat and the kernel copies flat, so one of these two curves would silently "
                    "win. Resolve it before uploading." % (src, name, taken[name], p))
            taken[name] = p
            continue
        rel = p.relative_to(src).as_posix()
        if name.startswith("smoke_"):
            skipped["smoke (excluded by design)"].append(rel)
        elif name.endswith(".tmp"):
            skipped["unfinished write (.tmp)"].append(rel)
        else:
            skipped["not a checkpoint name"].append(rel)
    return taken, skipped


def summarize(taken: dict[str, Path]) -> dict:
    """What the dataset will let the next session skip, read from the curves themselves."""
    per_seed = {}
    for name, p in sorted(taken.items()):
        if not name.endswith("_curve.json"):
            continue
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            per_seed[name] = {"unreadable": True}
            continue
        per_seed[name] = {
            "epochs_done": rec.get("epochs_done"),
            "horizon_epochs": rec.get("horizon_epochs"),
            "complete": rec.get("complete"),
            "shuffled_labels": rec.get("shuffled_labels"),
        }
    weights = sorted(n for n in taken if n.endswith(".pt"))
    return {"curves": per_seed, "weight_files": weights,
            "total_epochs_on_disk": sum(v.get("epochs_done") or 0 for v in per_seed.values()
                                        if isinstance(v.get("epochs_done"), int))}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="src", required=True,
                    help="the pulled checkpoint directory (searched recursively)")
    ap.add_argument("--slug", required=True, help="dataset slug, lowercase with hyphens")
    ap.add_argument("--title", default=None, help="dataset title (defaults to the slug)")
    ap.add_argument("--new", action="store_true",
                    help="create a NEW dataset rather than versioning an existing one. Prefer a "
                         "fresh slug per link: a new VERSION of an existing dataset changes what a "
                         "previously-pushed kernel would mount, which makes an earlier session's "
                         "recorded input no longer reconstructible.")
    ap.add_argument("--notes", default="chained finetune-arm checkpoints",
                    help="version notes, when versioning rather than creating")
    ap.add_argument("--dry-run", action="store_true",
                    help="classify and report, upload nothing")
    args = ap.parse_args(argv)

    src = Path(args.src)
    if not src.is_dir():
        ap.error("--from %s is not a directory" % src)

    taken, skipped = classify(src)
    print("ACCEPTED %d file(s) matching %s" % (len(taken), CKPT_NAME_RE.pattern))
    total = 0
    for name, p in sorted(taken.items()):
        size = p.stat().st_size
        total += size
        print("  %-28s %9.2f MB   <- %s" % (name, size / 1e6, p))
    print("  total %.2f MB" % (total / 1e6))
    for reason, names in skipped.items():
        if names:
            print("SKIPPED (%s): %s%s"
                  % (reason, ", ".join(names[:10]), " ..." if len(names) > 10 else ""))

    summ = summarize(taken)
    print("RESUMABLE STATE: %d epoch(s) on disk across %d curve(s)"
          % (summ["total_epochs_on_disk"], len(summ["curves"])))
    for name, info in summ["curves"].items():
        print("  %-28s %s" % (name, info))

    if not taken:
        raise SystemExit(
            "REFUSING: nothing under %s matches %s, so the dataset would be empty and the next "
            "session would recompute every seed. A directory that LOOKS full of checkpoints but "
            "holds only smoke output or .tmp files produces exactly this; see the SKIPPED lines "
            "above for what was there instead." % (src, CKPT_NAME_RE.pattern))

    if args.dry_run:
        print("DRY RUN -- nothing uploaded")
        return 0

    api = _api()
    user = api.get_config_value("username") or "emanueleebrahimi"
    with tempfile.TemporaryDirectory() as td:
        stage = Path(td)
        for name, p in taken.items():
            shutil.copyfile(p, stage / name)   # FLAT, deliberately -- see the module docstring
        meta = {"title": args.title or args.slug,
                "id": "%s/%s" % (user, args.slug),
                "licenses": [{"name": "CC0-1.0"}]}
        (stage / "dataset-metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        print("staging %d flat file(s) in %s" % (len(taken), stage))
        if args.new:
            api.dataset_create_new(str(stage), dir_mode="skip")
        else:
            api.dataset_create_version(str(stage), args.notes, dir_mode="skip")

    print()
    print("dataset: %s/%s" % (user, args.slug))
    print("NEXT -- both flags are required, and for different reasons:")
    print("  build: --ckpt-dataset %s      (the kernel's expected mount, used as a HINT)" % args.slug)
    print("  push : --dataset %s           (Kaggle actually attaches it)" % args.slug)
    print("Wait for `kaggle datasets status %s/%s` to report ready before pushing: a kernel pushed "
          "against a not-yet-attachable dataset ran with NO mount at all on 2026-11-15."
          % (user, args.slug))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
