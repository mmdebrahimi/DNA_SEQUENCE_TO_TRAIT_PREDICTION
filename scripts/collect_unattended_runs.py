"""Roll up the three successor arms: launcher state, result artifact, metric — or an explicit refusal.

READ-ONLY. It reports; it records nothing in the ledger (that is `record_successor_arm_results.py`'s
job, and it must not happen on a refused arm).

THREE ARMS, THREE SCHEMAS, THREE DIFFERENT METRIC PATHS. This is the whole reason the collector
exists rather than a one-line read. The as-given arm keeps its metric at `rung.median_spearman`, the
heads arm at `arm_metric`, and the fine-tune arm at whatever its scorer emits. A reader that assumed
one shape would silently return `None` for two of three — and `None` flowing into a ledger as a
metric is exactly the failure that produces a confident wrong number with no symptom.

SO AN UNRECOGNISED SCHEMA REFUSES rather than falling through to a default extractor. A default is
the dangerous choice here: it would read `None` and report it as a measured absence of signal, which
is the same trap as the scorer's NaN-to-0.0 mapping one layer up.

NOT EVERY ARM HAS A DETACHED JOB, and requiring one would refuse a good result. The as-given arm was
run ATTENDED — there is no launcher, no heartbeat, no sentinel. Gating it on a `COMPLETE` launcher
state would refuse an arm whose result has been on disk and committed for days. So the job is
declared per-arm and is `None` where none exists; where a job IS declared, only `COMPLETE` licenses
reading the artifact as the arm's result.

EXIT CODE IS A READINESS CHECK. Non-zero while any arm is uncollected, zero only when all three are
in. That makes it usable as the one command to run on return: zero means record them.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "wiki"

#: The declared surface: arm name -> where its result lives, which schema it must carry, how to read
#: its metric, and which detached job (if any) gates reading it.
#:
#: `metric_path` is a dotted path, resolved generically. Declaring the path as DATA rather than
#: writing three bespoke readers is what makes a new arm a one-line addition and keeps the refusal
#: logic in one place.
ARMS = {
    "asgiven-position-resolved-promoteratlas-150bp": {
        "artifact": "wiki/asgiven_arm_result_2026-11-02.json",
        "schema": "asgiven-arm-result-v1",
        "metric_path": "rung.median_spearman",
        "interval_path": "rung.gain_bootstrap",
        "job": None,          # run ATTENDED -- no launcher exists, and that is not a defect
        "why_no_job": "this arm was scored attended, before the detached-launcher harness existed",
    },
    "shipped-expression-heads-zeroshot-150bp": {
        "artifact": "wiki/expression_heads_result_2026-11-02.json",
        "schema": "expression-heads-result-v1",
        "metric_path": "arm_metric",
        "interval_path": "best_primary.gain_bootstrap",
        "job": "expression_heads",
    },
    "finetuned-promoteratlas-encoder-150bp": {
        "artifact": "wiki/finetune_arm_result_2026-11-13.json",
        "schema": "finetune-arm-result-v1",
        "metric_path": "arm_metric",
        "interval_path": "scored.gain_bootstrap",
        "job": "finetune_arm",
    },
}

#: Launcher states that license reading an artifact as an arm's result. Only one.
TERMINAL_OK = ("COMPLETE",)


class CollectRefused(ValueError):
    """This arm cannot be read as a result. Carries the reason, never a metric."""


def dotted(obj, path: str):
    """Resolve a dotted path, or raise. Raising beats returning None: a missing metric must be a
    refusal with a reason, not a null that reads like a measurement."""
    cur = obj
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise CollectRefused("path %r does not resolve (stopped at %r)" % (path, part))
        cur = cur[part]
    return cur


def collect_arm(name: str, spec: dict, *, repo: Path | None = None,
                state_dir: Path | None = None) -> dict:
    """One arm's row: state, verdict, metric — or a refusal with its reason."""
    repo = repo or REPO
    row: dict = {"arm": name, "artifact": spec["artifact"], "job": spec["job"]}

    # 1. The launcher gate, where a launcher exists at all.
    if spec["job"]:
        from scripts.detached_launch import classify

        st = classify(spec["job"], state_dir=state_dir)
        row["job_state"] = st.state
        row["job_detail"] = st.reason
        if st.state not in TERMINAL_OK:
            row["collected"] = False
            row["refusal"] = ("job %r is %s, not COMPLETE -- an artifact is not readable as the "
                              "arm's result until its run terminated cleanly"
                              % (spec["job"], st.state))
            return row
    else:
        row["job_state"] = "n/a"
        row["job_detail"] = spec.get("why_no_job", "no detached job declared for this arm")

    # 2. The artifact.
    p = repo / spec["artifact"]
    if not p.is_file():
        row["collected"] = False
        row["refusal"] = "artifact absent: %s" % spec["artifact"]
        return row
    try:
        art = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        row["collected"] = False
        row["refusal"] = "artifact is not readable JSON (%s)" % e
        return row

    # 3. The schema. An UNRECOGNISED schema refuses -- never a default extractor.
    got = art.get("schema")
    if got != spec["schema"]:
        row["collected"] = False
        row["refusal"] = (
            "schema is %r, expected %r. Refusing rather than guessing an extractor: a default "
            "would read None and report it as a measured absence of signal." % (got, spec["schema"]))
        return row
    row["schema"] = got

    # 4. The metric and its interval, by declared path.
    try:
        row["metric"] = dotted(art, spec["metric_path"])
        row["metric_path"] = spec["metric_path"]
    except CollectRefused as e:
        row["collected"] = False
        row["refusal"] = "metric unreadable: %s" % e
        return row
    if row["metric"] is None:
        row["collected"] = False
        row["refusal"] = ("metric at %r is null. A null metric is not a result; it is a missing "
                          "one." % spec["metric_path"])
        return row

    try:
        row["interval"] = dotted(art, spec["interval_path"])
    except CollectRefused:
        row["interval"] = None          # reported absent, not fatal: not every arm carries one

    row["verdict"] = art.get("verdict")
    row["verdict_detail"] = art.get("verdict_detail")
    row["verdict_is_rederivable"] = art.get("verdict_is_rederivable_from_this_artifact")
    row["collected"] = True
    return row


def collect(*, repo: Path | None = None, state_dir: Path | None = None,
            arms: dict | None = None) -> dict:
    arms = arms or ARMS
    rows = {name: collect_arm(name, spec, repo=repo, state_dir=state_dir)
            for name, spec in arms.items()}
    n_ok = sum(1 for r in rows.values() if r.get("collected"))
    return {
        "schema": "unattended-runs-collection-v1",
        "n_arms": len(rows),
        "n_collected": n_ok,
        "n_refused": len(rows) - n_ok,
        "all_collected": n_ok == len(rows),
        "arms": rows,
        "what_this_is_not": (
            "a READ-ONLY roll-up. It records nothing in the arms ledger, and a refused arm must not "
            "be recorded -- recording a refusal as a metric is how a null becomes a number."),
        "terminal_states_that_license_reading": list(TERMINAL_OK),
    }


def render(rep: dict) -> str:
    out = ["%d of %d arms collected" % (rep["n_collected"], rep["n_arms"]), ""]
    for name, r in rep["arms"].items():
        if r.get("collected"):
            iv = r.get("interval") or {}
            ci = ""
            if isinstance(iv, dict) and "ci_lo" in iv:
                ci = "  CI [%+.4f, %+.4f]" % (iv["ci_lo"], iv["ci_hi"])
            out.append("  OK      %-48s %-18s metric %.4f%s"
                       % (name[:48], r.get("verdict") or "?", float(r["metric"]), ci))
        else:
            out.append("  REFUSED %-48s %s" % (name[:48], r["refusal"]))
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None,
                    help="write the collection JSON here (default: skip the write)")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)

    rep = collect()
    print(json.dumps(rep, indent=2, sort_keys=True) if a.json else render(rep))

    if a.out:
        p = Path(a.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(rep, indent=2, sort_keys=True), encoding="utf-8")
        print()
        print("-> %s" % p)

    # Readiness check: zero ONLY when every arm is in.
    return 0 if rep["all_collected"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
