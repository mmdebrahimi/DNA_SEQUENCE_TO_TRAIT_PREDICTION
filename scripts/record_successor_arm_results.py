"""Record the three successor arms' measured results in the arms ledger.

DIRECT DICT WRITES, NEVER `measure_arm`. This is not a style preference; it is a measured hazard.
`ArmsLedger.measure_arm` ends in `return Arm(**base)`, and `Arm` is a 7-field dataclass with no
field filter. All three successor arms carry `rationale` and `declared_from`, so that construction
raises `TypeError` -- AFTER `arms[name] = base` has already mutated the in-memory dict. The caller
then sees an exception over a ledger that has silently changed, which is the worst of both outcomes.
`kill_arm`, two methods below it, already filters with
`{k: v for k, v in base.items() if k in Arm.__dataclass_fields__}`; `measure_arm` does not. Writing
the dict directly sidesteps the whole question and keeps `rationale` byte-for-byte.

IT REFUSES MORE THAN IT WRITES, and the refusals are the point.
  - An arm the COLLECTOR refused is not recorded. The collector returns a row with no `metric` key
    at all when it cannot produce a trustworthy number, so recording one would mean inventing it.
  - An arm already `measured` with a DIFFERENT metric is refused, never overwritten. A published
    number that moves silently is the failure the ledger's own historical-metric guard exists for.
  - Re-running with the SAME metric is a no-op, so the script is safely idempotent.

THE VERDICT IS REPORTED, NOT ACTED ON. With no `planned` arm left, `ceiling_verdict` reads `met` --
which is a statement that the declared ARM SET is exhausted, NOT that the 0.941 replicate ceiling
has been reached. Those are different claims and conflating them would be the over-claim this whole
arm structure exists to prevent. Accepting `met` is a user authority call; this script prints it and
stops.

A LEADER CHANGE IS SURFACED LOUDLY. `best` is the argmax over point metrics, so an arm can take the
lead by a margin its own confidence interval says is not a gain. Where that happens the script says
so explicitly, because the ledger field cannot.

Usage:
    python scripts/record_successor_arm_results.py            # dry run, prints the diff
    python scripts/record_successor_arm_results.py --apply    # writes the ledger
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

from scripts.collect_unattended_runs import ARMS as COLLECTOR_ARMS
from scripts.collect_unattended_runs import collect
from scripts.declare_successor_arms import ceiling_criterion, verdict_now

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "soraya_runs/arms-ledger.json"

#: The arms this script records. Keys are ledger arm names AND collector keys -- deliberately the
#: same strings, so a rename cannot leave the two halves pointing at different arms. Asserted equal
#: to the collector's declared surface by test rather than trusted.
RECORD_ARMS = tuple(sorted(COLLECTOR_ARMS))

#: Fields this script is allowed to write. Anything else in an arm dict is preserved untouched.
WRITABLE = ("status", "metric", "evidence", "resolved_at")

STATUS_MEASURED = "measured"
STATUS_PLANNED = "planned"


class RecordRefused(ValueError):
    """This arm cannot be recorded. Carries the reason; never writes a number."""


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_ledger(path: Path | str = LEDGER) -> dict:
    p = Path(path)
    if not p.is_file():
        raise RecordRefused("arms ledger absent: %s" % p)
    return json.loads(p.read_text(encoding="utf-8"))


def save_ledger(data: dict, path: Path | str = LEDGER) -> Path:
    """Atomic durable write: tmp -> flush -> fsync -> replace.

    `write_text` returns once the data is in the page cache, so a power cut keeps the file's
    metadata and discards its contents -- which on this host produced a NUL-filled artifact that
    read as a malformed result rather than an absent one. The ledger is the last place that should
    happen.
    """
    final = Path(path)
    final.parent.mkdir(parents=True, exist_ok=True)
    tmp = final.with_suffix(final.suffix + ".tmp")
    blob = json.dumps(data, indent=1, sort_keys=True)
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(blob)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, final)
    return final


def evidence_for(name: str, row: dict) -> str:
    """An executed-evidence string naming the artifact a reader can open."""
    iv = row.get("interval") or {}
    ci = ""
    if isinstance(iv, dict) and "ci_lo" in iv:
        ci = "; CI [%+.4f, %+.4f]" % (iv["ci_lo"], iv["ci_hi"])
    return ("%s (%s; metric %.4f%s) collected by scripts/collect_unattended_runs.py"
            % (row["artifact"], row.get("verdict") or "?", float(row["metric"]), ci))


def plan_records(data: dict, report: dict) -> dict:
    """Decide, WITHOUT mutating, what each arm's recording would be.

    Returns {name: {"action": record|noop|refused, ...}}. Separating the decision from the write is
    what lets a dry run print exactly what an apply would do.
    """
    crit = ceiling_criterion(data)
    arms = data["criteria"][crit]
    out: dict = {}
    for name in RECORD_ARMS:
        row = (report.get("arms") or {}).get(name)
        if row is None:
            out[name] = {"action": "refused",
                         "reason": "the collector produced no row for this arm"}
            continue
        if not row.get("collected"):
            out[name] = {"action": "refused",
                         "reason": "the collector REFUSED this arm: %s" % row.get("refusal")}
            continue
        if "metric" not in row or row["metric"] is None:
            out[name] = {"action": "refused",
                         "reason": "the collected row carries no metric; a null is not a result"}
            continue

        metric = float(row["metric"])
        existing = arms.get(name)
        if existing is None:
            out[name] = {"action": "refused",
                         "reason": ("arm is not declared in the ledger. Recording an undeclared "
                                    "arm would create it as a side effect of measuring it; "
                                    "declare it first.")}
            continue
        if existing.get("status") == STATUS_MEASURED:
            old = existing.get("metric")
            if old is not None and abs(float(old) - metric) > 1e-12:
                out[name] = {"action": "refused",
                             "reason": ("already measured at %r and the collector now reads %r. A "
                                        "published metric does not move silently; resolve which is "
                                        "right rather than overwriting."
                                        % (old, metric))}
            else:
                out[name] = {"action": "noop", "metric": metric,
                             "reason": "already recorded with this metric"}
            continue

        out[name] = {"action": "record", "metric": metric,
                     "verdict": row.get("verdict"),
                     "evidence": evidence_for(name, row),
                     "from_status": existing.get("status")}
    return out


def apply_records(data: dict, decisions: dict, *, now: str | None = None) -> list:
    """Write the planned records by DIRECT dict assignment. Returns the names recorded."""
    crit = ceiling_criterion(data)
    arms = data["criteria"][crit]
    stamp = now or _now()
    recorded = []
    for name, d in sorted(decisions.items()):
        if d["action"] != "record":
            continue
        arm = arms[name]
        before = {k: v for k, v in arm.items() if k not in WRITABLE}
        arm["status"] = STATUS_MEASURED
        arm["metric"] = d["metric"]
        arm["evidence"] = d["evidence"]
        arm["resolved_at"] = stamp
        after = {k: v for k, v in arm.items() if k not in WRITABLE}
        if before != after:                      # defence, not decoration: direct writes are easy
            raise RecordRefused(               # to get wrong and `rationale` is unrecoverable
                "recording %r altered a non-writable field; refusing. before=%r after=%r"
                % (name, sorted(before), sorted(after)))
        recorded.append(name)
    return recorded


def leader_note(data: dict) -> dict:
    """Who leads on the POINT metric, and whether that lead is a demonstrated gain.

    `best` is a bare argmax, so a 0.0010 lead and a 0.1000 lead look identical in the field. Where
    the leader's own interval contains zero against the runner-up, that is said here, because
    nothing in the ledger can say it.
    """
    crit = ceiling_criterion(data)
    arms = data["criteria"][crit]
    measured = {n: a["metric"] for n, a in arms.items()
                if a.get("status") == STATUS_MEASURED and a.get("metric") is not None}
    if not measured:
        return {"leader": None, "detail": "no measured arm carries a metric"}
    ordered = sorted(measured.items(), key=lambda kv: -float(kv[1]))
    leader, lead_metric = ordered[0]
    runner, runner_metric = (ordered[1] if len(ordered) > 1 else (None, None))
    margin = (float(lead_metric) - float(runner_metric)) if runner else None
    return {
        "leader": leader, "metric": float(lead_metric),
        "runner_up": runner, "runner_up_metric": runner_metric,
        "margin_over_runner_up": None if margin is None else round(margin, 4),
        "caution": ("`best` is the argmax over POINT metrics. A margin smaller than the arms' own "
                    "confidence intervals is not a demonstrated gain, and this field cannot say "
                    "so -- read each arm's interval before calling the leader better."),
    }


def prospective(data: dict, decisions: dict) -> tuple:
    """What the verdict and the leader WOULD be, computed through the real skill code.

    A dry run that reported the CURRENT leader would hide exactly the event worth seeing: the
    as-given arm at 0.5086 overtakes the atlas arm at 0.5076 the moment it is recorded. So the
    preview applies the decisions to a deep copy, writes it into a throwaway cwd, and asks
    `verdict_now` about THAT -- no re-derivation of the verdict rule, and the live ledger is
    untouched.
    """
    import copy
    import tempfile

    trial = copy.deepcopy(data)
    apply_records(trial, decisions, now="(dry-run)")
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "soraya_runs" / LEDGER.name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(trial, indent=1, sort_keys=True), encoding="utf-8")
        return verdict_now(td), leader_note(trial)


def render(decisions: dict, before: dict, after: dict, leader: dict, *, applied: bool) -> str:
    L = ["ARMS RECORDING (%s)" % ("APPLIED" if applied else "DRY RUN -- nothing written"), ""]
    for name, d in sorted(decisions.items()):
        if d["action"] == "record":
            L.append("  RECORD  %-50s metric %.4f  (%s)"
                     % (name[:50], d["metric"], d.get("verdict") or "?"))
        elif d["action"] == "noop":
            L.append("  no-op   %-50s %s" % (name[:50], d["reason"]))
        else:
            L.append("  REFUSED %-50s %s" % (name[:50], d["reason"]))
    L.append("")
    label = "after: " if applied else "WOULD BE:"
    L.append("CEILING VERDICT  before:   %s" % before.get("verdict"))
    L.append("                 %s %s" % (label, after.get("verdict")))
    L.append("                 %s" % after.get("detail"))
    L.append("")
    L.append("WHAT `met` DOES AND DOES NOT MEAN")
    L.append("  `met` says the declared ARM SET is exhausted -- every declared arm has been run.")
    L.append("  It does NOT say the 0.941 replicate-derived ceiling has been reached. Accepting it")
    L.append("  as a resolution of the ceiling question is a USER AUTHORITY call, not this")
    L.append("  script's, and nothing here acts on it.")
    L.append("")
    if leader.get("leader"):
        L.append("LEADER  %s at %.4f" % (leader["leader"], leader["metric"]))
        if leader.get("runner_up"):
            L.append("        runner-up %s at %.4f (margin %+.4f)"
                     % (leader["runner_up"], float(leader["runner_up_metric"]),
                        leader["margin_over_runner_up"]))
        L.append("        %s" % leader["caution"])
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true",
                    help="write the ledger; the default is a dry run that prints the same decisions")
    ap.add_argument("--ledger", default=str(LEDGER))
    a = ap.parse_args(argv)

    try:
        data = load_ledger(a.ledger)
        report = collect()
        decisions = plan_records(data, report)
    except RecordRefused as e:
        print("REFUSED: %s" % e)
        return 2

    before = verdict_now()
    has_writes = any(d["action"] == "record" for d in decisions.values())
    if a.apply and has_writes:
        try:
            apply_records(data, decisions)
        except RecordRefused as e:
            print("REFUSED: %s" % e)
            return 2
        save_ledger(data, a.ledger)
        after, leader = verdict_now(), leader_note(data)
    elif has_writes:
        # A dry run must preview the PROSPECTIVE state, or it hides the leader change.
        try:
            after, leader = prospective(data, decisions)
        except RecordRefused as e:
            print("REFUSED: %s" % e)
            return 2
    else:
        after, leader = before, leader_note(data)

    print(render(decisions, before, after, leader, applied=bool(a.apply)))

    refused = [n for n, d in decisions.items() if d["action"] == "refused"]
    return 1 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
