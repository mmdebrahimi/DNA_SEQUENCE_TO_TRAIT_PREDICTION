"""Declare the three successor arms on the tile-chain ceiling criterion under R6. IDEMPOTENT.

The positional arm's execution ended with `ceiling_verdict` reading `met`, because recording it
removed the last `planned` arm. That is a statement about the ARM SET being exhausted, not about the
ceiling having been reached -- `best` sits at 0.5076 against a 0.941 replicate-derived ceiling, 54%
of it. The user's decision (2026-11-02) was the elaborate branch: declare and run all three named
successor arms rather than accept the terminal.

**Declaring them is what makes them unable to evaporate, which is the whole point of R6.** An arm
declared only once its cheap sibling succeeds is an arm that quietly never happens.

**Why a committed script rather than a hand edit.** The rationales below are the arms' provenance,
and a hand edit leaves no re-runnable record of what was declared or why. Re-running this script is
a no-op, so it is safe to use as the audit trail.

**Why a DIRECT dict write rather than `ArmsLedger.declare_arm`.** `declare_arm` returns
`Arm(**arms[name])` for an arm that already exists, and `Arm` is a 7-field dataclass with no field
filter -- unlike `kill_arm`, which filters to `Arm.__dataclass_fields__`. So the moment an arm
carries an annotation key (and all three of these carry `rationale` and `declared_from`), the API
path raises `TypeError` on the RETURN construction, after having already mutated the dict. A direct
write preserves every annotation key by construction and never half-applies.

**`evidence` stays EMPTY and that is enforced elsewhere.** `evidence` means EXECUTED evidence under
R6/R7. A rationale, a substrate fact or a wall classification is not that, so none of them may be
parked there; `test_every_PLANNED_arm_still_has_an_empty_evidence_string` is the guard.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "soraya_runs/arms-ledger.json"
CRITERION_PREFIX = "ceiling:"
PLAN = ("plans/Declare_And_Run_All_Three_Successor_Arms_And_Convert_The_Guards_From_State_Pins_"
        "To_Invariants_Plan/technical-plan.md")

#: Arm fields on the `Arm` dataclass. Everything else is an annotation key that a direct write
#: preserves and the API path would raise on.
ARM_FIELDS = ("name", "cost_tier", "status", "metric", "evidence", "declared_at", "resolved_at")

ARMS: dict[str, dict] = {
    "asgiven-position-resolved-promoteratlas-150bp": {
        "cost_tier": "cheap",
        "rationale": (
            "Declared 2026-11-02 under R6 after the positional arm measured. **Motivated by a "
            "MEASUREMENT rather than a preference, which is the only reason it is worth an arm:** "
            "the positional arm's orientation diagnostic resolved AGAINST orienting -- as-given "
            "0.5086 against oriented 0.4893 -- so reversing the bin axis on minus-strand tiles "
            "DESTROYS about 0.0193 of signal rather than aligning the cohorts. The positional arm's "
            "declared co-primary was the ORIENTED rung, so its eligible primary was handicapped by "
            "the transform that same run showed to be wrong. This arm changes the DECLARED PRIMARY "
            "to the as-given reduction, which is what makes it a distinct arm under R6 and NOT a "
            "re-read of the positional arm. **Pre-committed honesty: it is expected to TIE, not "
            "win.** 0.5086 sits above the incumbent 0.5076 but its paired CI [-0.0079, +0.0017] "
            "CONTAINS zero, so the bar predicts a tie and the arm's job is not to re-report 0.5086 "
            "as a gain. **And its bar is NOT a blind pre-registration** -- the number was already "
            "visible as the positional run's INELIGIBLE secondary, which is a real weakness stated "
            "rather than hidden. Cheap because the features are already extracted and committed on "
            "D:; it needs no second forward pass."
        ),
    },
    "shipped-expression-heads-zeroshot-150bp": {
        "cost_tier": "medium",
        "rationale": (
            "Declared 2026-11-02 under R6. The PromoterAtlas clone ships five expression-prediction "
            "heads and a segmentation head beside the frozen base, and none has ever been scored on "
            "this rung. **Measured, not assumed: all four transcription heads carry base weights "
            "BYTE-IDENTICAL to promoteratlas-base.pt -- 124/124 tensors, max |dev| 0.000e+00 -- so "
            "this arm measures FROZEN encoder + trained head, not a fine-tuned encoder**, and all "
            "four can share one forward pass. **It is NOT a drop-in, and the obstruction is "
            "arithmetic:** every head's fc1.weight is (64, 11008) where 11008 = 128 x 86, so the "
            "heads are geometry-locked to 86 bp while this substrate is 150 bp (which flattens to "
            "128 x 150 = 19200). The arm therefore needs a windowing decision frozen in its bar "
            "BEFORE any head output is scored, and measured against a true per-window recompute "
            "because conv padding='same' makes the cheap full-latent reduction differ at the edges. "
            "**Zero-shot TRANSFER framing, pre-committed:** these heads were trained on OTHER "
            "libraries, so a null bounds transfer to the tile rung rather than bounding the heads. "
            "**Circularity scope, narrower than it looks:** the urtecho2018 head IS circular on the "
            "sigma70 rung (GSE108535) and must never be a comparator there, but this substrate is "
            "GSE144621 peak tiles and the repo's own withdrawn lead established that no "
            "genome-wide/tiling Urtecho head exists -- so it is not circular HERE. Recorded so "
            "nobody reuses it on the wrong rung."
        ),
    },
    "finetuned-promoteratlas-encoder-150bp": {
        "cost_tier": "large",
        "rationale": (
            "Declared 2026-11-02 under R6, and it is the arm the positional arm's own limits point "
            "at: that arm froze the encoder, so its NO_GAIN bounds nothing about a fine-tuned one. "
            "**Tractable here precisely because the encoder is small -- 1,847,428 parameters, not a "
            "1.2B LM -- so this is NOT comparable to the project's blocked large-LM fine-tunes.** "
            "Local compute is CPU-only (cuda False, 4 threads) at a measured 57.5 seq/s forward; "
            "the established GPU route is a Kaggle T4 kernel and scripts/kaggle/ already holds "
            "seven builders, so the route is chosen from a MEASURED probe rather than estimated. "
            "**The outcome that must not be over-read, pre-committed:** beating 0.5076 shows the "
            "frozen representation was the limit, while failing to shows either the SUBSTRATE or "
            "the TRAINING BUDGET was -- and a budget-limited null is not a representation-limited "
            "one. The probe's measured epoch-to-early-stop curve is what separates them, and "
            "without that curve this arm should not be run at all."
        ),
    },
}


class DeclareError(RuntimeError):
    """The ledger is not in a state where these arms can be declared."""


def ceiling_criterion(data: dict) -> str:
    crits = [k for k in data.get("criteria", {}) if k.startswith(CRITERION_PREFIX)]
    if len(crits) != 1:
        raise DeclareError(
            "expected exactly one %r criterion, found %d: %r" % (CRITERION_PREFIX, len(crits), crits)
        )
    return crits[0]


def verdict_now(cwd: str | None = None) -> dict:
    """Recompute `ceiling_verdict` through the real skill code, from DISK."""
    sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")
    try:
        from arms_ledger import ArmsLedger  # type: ignore
    finally:
        sys.path.pop(0)
    led = ArmsLedger(cwd or os.getcwd())
    return led.ceiling_verdict(ceiling_criterion(led._data)).to_dict()


def declare(data: dict, *, now: str | None = None) -> tuple[list[str], list[str]]:
    """Add any missing arm. Returns (added, already_present). Never touches an existing arm."""
    crit = ceiling_criterion(data)
    arms = data["criteria"][crit]
    stamp = now or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    added, present = [], []
    for name, spec in ARMS.items():
        if name in arms:
            present.append(name)
            continue
        arms[name] = {
            "name": name,
            "cost_tier": spec["cost_tier"],
            "status": "planned",
            "metric": None,
            "evidence": "",
            "declared_at": stamp,
            "resolved_at": "",
            "rationale": spec["rationale"],
            "declared_from": PLAN,
        }
        added.append(name)
    return sorted(added), sorted(present)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    before = verdict_now()
    print("BEFORE  verdict %s | measured %d / planned %d / killed %d | best %s"
          % (before["verdict"], len(before["measured"]), len(before["planned"]),
             len(before["killed"]), before["best"]))

    data = json.loads(LEDGER.read_text(encoding="utf-8"))
    crit = ceiling_criterion(data)
    pre = {n: (a.get("status"), a.get("metric")) for n, a in data["criteria"][crit].items()}
    n_before = len(data["criteria"][crit])

    added, present = declare(data)
    if present:
        print("already declared (no-op): %s" % ", ".join(present))
    if not added:
        print("nothing to add; the ledger is unchanged")
        return 0
    if args.dry_run:
        print("DRY RUN would add: %s" % ", ".join(added))
        return 0

    LEDGER.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")

    # --- semantic diff: no pre-existing arm may have moved -------------------------------------
    back = json.loads(LEDGER.read_text(encoding="utf-8"))
    post_arms = back["criteria"][crit]
    moved = {n: (pre[n], (post_arms[n].get("status"), post_arms[n].get("metric")))
             for n in pre if (post_arms[n].get("status"), post_arms[n].get("metric")) != pre[n]}
    if moved:
        raise DeclareError("pre-existing arm(s) changed, which this script must never do: %r" % moved)
    lost = sorted(set(pre) - set(post_arms))
    if lost:
        raise DeclareError("pre-existing arm(s) vanished: %r" % lost)

    print("added: %s" % ", ".join(added))
    print("arms  %d -> %d  | no pre-existing arm's status or metric moved (semantic diff clean)"
          % (n_before, len(post_arms)))

    after = verdict_now()
    print("AFTER   verdict %s | measured %d / planned %d / killed %d | best %s"
          % (after["verdict"], len(after["measured"]), len(after["planned"]),
             len(after["killed"]), after["best"]))
    print("        planned: %s" % ", ".join(after["planned"]))
    if after["verdict"] != "blocked:untried-arm":
        raise DeclareError(
            "declaring three arms should have returned the verdict to blocked:untried-arm, but it "
            "reads %r. The whole point of the declaration is that an untried arm blocks the "
            "criterion." % after["verdict"]
        )
    if after["best"] != before["best"]:
        raise DeclareError(
            "`best` moved from %r to %r. Declaring an arm must never move `best` -- a planned arm "
            "has no metric." % (before["best"], after["best"])
        )
    print()
    print("the verdict returned to blocked:untried-arm, which is what converts the positional "
          "plan's guard ESCALATION into an ordinary re-pin")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
