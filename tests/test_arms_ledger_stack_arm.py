"""Guards for the declared orthogonal-modality stack arm.

RE-PINNED A FIFTH TIME ON 2026-11-01, and this time the PLANNED arm is a different arm rather than
a shorter list. `substrate-change-constructed-eukaryotic-or-human-panel` RAN locally over 10
published folds of HepG2 (n=61,463, join rate 1.0000) against a bar frozen beforehand at
`wiki/eukaryotic_arm_acceptance_bar.json`, and returned SUBSTRATE_LIMITED: a retrained conv encoder
reaches 0.4853 against a GC-only ridge at 0.2362, a gain of +0.2471 with a paired-bootstrap CI
[+0.2399, +0.2536] over folds, against the BACTERIAL anchor of +0.0122. Legitimate under R6 -- the
arm resolved on its OWN executed evidence. In the same breath, G-A candidate 8's
`positional-non-mean-pooled-promoteratlas-150bp` was DECLARED under R6, because an arm declared only
once its cheap sibling succeeds is an arm that evaporates. So the set is now 11 measured / 1 planned
/ 1 KILLED over 13.

What is UNCHANGED, and it is what keeps this a re-pin rather than an escalation: the ceiling verdict
is STILL `blocked:untried-arm` and STILL `not met`, and `best` is STILL the atlas arm at 0.5076. By
this file's own standing rule, had re-pinning flipped the verdict the correct action would have been
to stop and escalate instead.

What this re-pin makes STRONGER, and it is a new KIND of invariant here:
`test_the_human_arms_metric_is_not_read_as_best` pins that the eukaryotic arm's metric -- a HUMAN
lentiMPRA Spearman sitting in a criterion whose other ten are BACTERIAL -- cannot silently become
`best`. `ceiling_verdict` picks `best` as a max over metrics, so the moment a human arm outscores
0.5076 that max would be comparing across substrates with different noise floors, which is exactly
what the arm's own bar forbids. It is harmless today (0.4853 < 0.5076) and the guard exists so it
cannot become harmful quietly; the fix if it fires is to move the arm to its own criterion, which is
a scope change and therefore the user's call, NOT to relax the assertion.

Two test names were also made COUNT-FREE (`test_the_measured_arms_are_unchanged`,
`test_the_arm_set_is_exactly_the_pinned_partition`). They read "ten" while asserting over eleven
arms, having carried a count through four re-pins -- a guard whose name misstates what it guards is
the same drift this project keeps finding in orientation strings.

RE-PINNED A FOURTH TIME ON 2026-11-01, and this time the BEST arm moves, which is a substantive
change rather than a count. `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` RAN on a
Kaggle T4 over 44,106 peak tiles x 10 seeds and closed `ATLAS_BEATS_THERMODYNAMIC` at median
Spearman 0.5076, against the previous best overall -- this very stack arm -- at 0.4404. Every gate
passed, including a shuffled-label null pooled to 270 seeds at median -0.0014, ci95
[-0.0114, +0.0098]; the 10-seed null had failed at 0.0203 and the run was NOT graded on it. So the
set is now 10 measured / 1 planned / 1 KILLED over 12. Legitimate under R6: the arm resolved on its
OWN executed evidence against a bar frozen before the run.

What is UNCHANGED, and it is what keeps this a re-pin: the ceiling verdict is STILL
`blocked:untried-arm` and STILL `not met`, because one declared arm has still never run. A new best
arm does not clear a criterion that an untried arm blocks.

What CHANGED and is asserted explicitly rather than quietly updated: `verdict.best` is no longer
NEW_ARM. The previous re-pin added `assert verdict.best == NEW_ARM` precisely so that a shift in the
leader could not pass silently, and it has now fired as designed.

RE-PINNED A THIRD TIME ON 2026-10-31, for the same reason and by the precedent this docstring
already records. `pyramid-4scale-2.22decades-499bp-fragments` RAN on 2026-10-28; its verdict was
void under INVALID_NULL_SCORED until the 270-seed pooled null PASSED the validity gate on
2026-10-31 under the pre-registered Branch A of `wiki/fragment_null_270seed_prereg_2026-10-30.md`,
which made the arm's verdict computable as NO_GAIN (paired median -0.0127 [-0.0205,-0.0040] against
a 4-mer ridge, 2 of 10 seeds winning against a pre-registered >=8). So the set is now 9 measured /
2 planned / 1 KILLED over 12 arms. Legitimate under R6 -- the arm resolved on its OWN executed
evidence, not a sibling's -- and the gate it passed was frozen before any batch ran.

What is DELIBERATELY UNCHANGED, because that is what makes this a re-pin and not a rewrite: the
ceiling verdict is STILL `blocked:untried-arm` and STILL `not met`. Retiring this arm did not clear
the criterion; it only shortened the list of reasons it is blocked. If re-pinning had flipped the
verdict, the correct action would have been to stop and escalate instead.
`test_the_successor_arm_records_WHY_it_is_worth_running` is STRONGER than before: it still asserts
the rationale is intact and that it has NOT leaked into the executed-evidence field, and now also
asserts that field holds real executed evidence.

RE-PINNED TWICE ON 2026-10-27. Second time: the 4-decade arm was RETIRED on its own executed
evidence (`scripts/pyramid_span_feasibility_probe.py` measured the maximum HONEST span at 2.2210
decades against the 4.00 it declares -- 60.1x short at the encoder's 3 bp floor) and the reachable
successor `pyramid-4scale-2.22decades-499bp-fragments` was declared in its place. So the set is now
8 measured / 3 planned / 1 KILLED over 12 arms. That is a legitimate move under R6 -- an arm dying by
its OWN result, not a sibling's -- and the rails follow the ledger rather than freezing a stale set.

First time: the declared stack arm RAN. It was `planned` when these guards were
written and is now `measured` at 0.4404 on its own evidence, so four of them went red
because their SUBJECT legitimately moved -- not because the ledger was tampered with.
Re-pinning to the state that exists is the answer here; the invariants are UNCHANGED and
two are now stronger, asserted over EVERY planned arm rather than over this one.

The EIGHT measured arms are pinned by metric so a future edit cannot quietly restate a
result, and the refusal `kill_arm` makes on a redundancy reason is pinned too -- because
withdrawing the kill option is only safe while that refusal actually fires.

Every mutating test runs against a COPY of the ledger in `tmp_path`. A test that
killed an arm in the committed ledger would be a destructive test.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, r"C:/my_skills/skills/soraya/scripts")

import arms_ledger as al  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
LEDGER = REPO / "soraya_runs" / "arms-ledger.json"
CRITERION_PREFIX = "ceiling: best achievable prediction of bacterial promoter"
NEW_ARM = "orthogonal-modality-stack-promcalc-plus-N-channels"

#: The reachable re-scope that replaced multiresolution-pyramid-4-decades-UNRUN. PLANNED until
#: 2026-10-31, now MEASURED -- so it is pinned in MEASURED and named here for the kill-reason
#: cross-reference, which must keep pointing at it whatever its status.
SUCCESSOR_MEASURED = "pyramid-4scale-2.22decades-499bp-fragments"

#: Ran 2026-11-01 and closed ATLAS_BEATS_THERMODYNAMIC. It was `best` on this criterion until
#: 2026-11-04, which is why BEST_ARM below is a separate constant from NEW_ARM.
ATLAS_ARM = "promoteratlas-base-embeddings-frozen-150bp-padded-to-200"

#: RE-PINNED 2026-11-04, third firing of the leader guard, and the reason matters more than the move.
#:
#: The as-given arm was recorded MEASURED at 0.5086 against ATLAS_ARM's 0.5076, so the argmax moved
#: by **+0.0010**. Per this guard's own instruction the constant is re-pinned naming the new leader
#: rather than the assertion relaxed.
#:
#: **BUT THE LEAD IS NOT A DEMONSTRATED GAIN, AND NOTHING IN THE LEDGER CAN SAY SO.** The arm's own
#: verdict is `ASGIVEN_TIES` with a paired CI of [-0.0068, +0.0004] -- an interval CONTAINING ZERO,
#: and roughly seven times the margin by which it leads. The arm was pre-committed as expected to
#: TIE, in those words, before it ran. So `best` here has become a bookkeeping argmax over point
#: estimates that crosses two arms the evidence cannot separate: reading "the as-given arm beat the
#: atlas arm" off this field would be exactly the over-claim the arm's own bar was written to
#: prevent. Both are, on the measured evidence, the joint leaders at ~0.508.
ASGIVEN_ARM = "asgiven-position-resolved-promoteratlas-150bp"
BEST_ARM = ASGIVEN_ARM

#: Recorded MEASURED 2026-11-04 alongside ASGIVEN_ARM. A clean interpretable negative, not a null
#: result -- which is a distinction the metric alone does not carry; see its MEASURED note.
HEADS_ARM = "shipped-expression-heads-zeroshot-150bp"

#: RAN 2026-11-01 locally over 10 published folds of HepG2. Its metric is a HUMAN lentiMPRA
#: Spearman while every sibling metric is a BACTERIAL tile Spearman, so it is pinned here for
#: completeness and `test_the_human_arms_metric_is_not_read_as_best` pins that the incommensurability
#: cannot turn into a silent cross-substrate `best`.
EUK_ARM = "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"

#: Declared 2026-11-01 under R6 from G-A candidate 8, AT THE MOMENT its sibling ran -- which is the
#: point of R6: an arm declared only once its cheap sibling succeeds is an arm that evaporates.
POSITIONAL_ARM = "positional-non-mean-pooled-promoteratlas-150bp"

# The measured arms, pinned by metric. Restating a measured number elsewhere
# would now have to change this line too.
MEASURED = {
    "gc-single-feature": 0.3112,
    "conv-3width-13k-params-150bp": 0.3234,
    "pyramid-4scale-1.70decades-150bp": 0.3167,
    "promoter-calculator-biophysical-346params-fwd_max_tx": 0.3968,
    "promcalc-full-decomposition-ridge-150bp": 0.3925,
    "fragment-rung-321k-rows-size-lever-LB-M9": 0.185,
    "promcalc-plus-conv-rank-average-hybrid-150bp": 0.4254,
    # Ran 2026-10-26: WEAK_STACK, +0.0145 over the 2-channel hybrid at 10/10 seeds,
    # CI95 [0.0108, 0.0172] excluding zero, short only of the INHERITED 0.02 margin.
    NEW_ARM: 0.4404,
    # Ran 2026-10-28, GRADEABLE 2026-10-31 once the 270-seed null cleared the validity gate:
    # NO_GAIN. It LOSES to a 4-mer-count ridge (0.1001 vs 0.1120) on the fragment rung, which is
    # why retiring it shortens the blocked list without improving `best`.
    SUCCESSOR_MEASURED: 0.1001,
    # Ran 2026-11-01 on a Kaggle T4: ATLAS_BEATS_THERMODYNAMIC. Frozen cross-species encoder,
    # masked mean-pooled to 128 dims at NATIVE 150 bp (the declared padded-200 geometry scored
    # 0.4983, so padding costs 0.0093). Paired +0.1956 [0.1704, 0.2182] over conv at 10/10 seeds on
    # byte-identical test sets. Read as an UPPER BOUND: corpus-selection exposure is not excluded.
    ATLAS_ARM: 0.5076,
    # Ran 2026-11-01 LOCALLY (D: reattached), 10 published folds of HepG2, n=61,463, join rate
    # 1.0000: SUBSTRATE_LIMITED. conv 0.4853 vs a GC-only ridge 0.2362 -> gain +0.2471 with a
    # paired-bootstrap CI [+0.2399, +0.2536] against the BACTERIAL anchor +0.0122 (conv 0.3234 - gc
    # 0.3112). NOT COMMENSURABLE with the arms above: this is a human Spearman on a different assay
    # with a different noise floor, which is why the arm's contribution is the CONTRAST and not the
    # level, and why it must never be read as `best`.
    EUK_ARM: 0.4853,
}

#: THE SUCCESSOR ARMS ARE DELIBERATELY NOT IN `MEASURED`, and this note exists because adding them
#: is the obvious-looking move that breaks the non-vacuity harness.
#:
#: `MEASURED` is a HISTORICAL pin -- its contract is "a published number cannot change", not "this
#: is every measured arm". `test_no_historical_measured_metric_has_moved` was converted to a SUBSET
#: form precisely so the three successor arms could be measured without three more re-pins, and
#: `tests/test_arms_invariant_guards_non_vacuous.py` proves the guards still fire by building
#: synthetic states in which those same arms are `planned` with `metric=None`. An arm cannot be
#: both a historical pin and a lever the harness flips to planned, so the two leader comparisons
#: below read their metrics from the LIVE LEDGER instead -- which is also the stronger reading,
#: since a hand-maintained constant can drift from the file it describes.
#:
#: Recorded 2026-11-04, for the record rather than as pins:
#:   ASGIVEN_ARM 0.5086 -- ASGIVEN_TIES, paired CI [-0.0068, +0.0004] CONTAINING ZERO, ~7x the
#:     +0.0010 by which it now leads. Pre-committed as expected to tie. Quote it as a TIE.
#:   HEADS_ARM   0.1178 -- HEADS_NO_GAIN, gain -0.3922, CI [-0.4131, -0.3861], 0/10 seed wins,
#:     3 of 5 heads ANTI-correlated. A clean INTERPRETABLE negative, not a null result: the
#:     predictor-side permutation null passed at max|null| 0.0205 = 0.174x the real signal.

#: HISTORICAL, and deliberately no longer asserted anywhere (converted 2026-11-02). This was the
#: exact planned-set pin; three guards compared against it and all three went red every time an arm
#: legitimately resolved or a new one was declared. They now assert the INVARIANT the pin was
#: approximating -- `planned` non-empty iff the verdict is `blocked:untried-arm` -- so this set is
#: kept only as the record of what one day's arm set was. **Do not re-introduce it into an
#: assertion**: that is the drift the conversion exists to end.
HISTORICAL_PLANNED_2026_11_01 = {
    POSITIONAL_ARM,
}

#: Retired 2026-10-27 on its own executed evidence, NOT on a sibling's result.
KILLED_ARM = "multiresolution-pyramid-4-decades-UNRUN"
SUCCESSOR = "pyramid-4scale-2.22decades-499bp-fragments"
SPAN_PROBE = "wiki/pyramid_span_feasibility_2026-10-27.json"


def _live_ledger() -> al.ArmsLedger:
    return al.ArmsLedger(cwd=str(REPO))


def _criterion(ledger: al.ArmsLedger) -> str:
    return next(k for k in ledger._data["criteria"] if k.startswith(CRITERION_PREFIX))


@pytest.fixture
def sandbox(tmp_path: Path) -> Path:
    """A throwaway copy of the real ledger -- mutating tests must not touch the real one."""
    dest = tmp_path / "soraya_runs"
    dest.mkdir(parents=True)
    shutil.copy2(LEDGER, dest / "arms-ledger.json")
    return tmp_path


def test_the_stack_arm_RAN_and_carries_its_own_evidence():
    """Was `planned`; ran 2026-10-26. R6 is satisfied by its OWN executed evidence.

    The arm could not be killed on a sibling's result (see the withdrawal guards below) and
    it was not: it was RUN, and that is the only thing that legitimately moves an arm off
    `planned`. The metric is pinned in MEASURED, so a restatement has to edit that line too.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    assert NEW_ARM in arms
    arm = arms[NEW_ARM]
    assert arm["status"] == al.STATUS_MEASURED
    assert arm["metric"] == 0.4404
    assert arm["cost_tier"] == "cheap", "running it did not change what it cost"
    assert arm["evidence"].strip(), "a measured arm must carry its evidence"


def test_every_PLANNED_arm_still_has_an_empty_evidence_string():
    """Generalised from the stack arm, which used to be the one carrying this invariant.

    `evidence` means EXECUTED evidence (R6/R7). A substrate fact, a wall classification or a
    rationale is NOT that, so none of them may be parked in this field -- which is why the
    three untried arms' walls live in `wiki/ceiling_arm_wall_classification_2026-10-26.md`
    and not in the ledger. Asserting it over every planned arm is strictly stronger than
    asserting it over one.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    planned = {n: a for n, a in arms.items() if a["status"] == al.STATUS_PLANNED}
    assert planned, "vacuous: no planned arms left to check"
    for name, arm in planned.items():
        assert arm["evidence"] == "", "%s: a planned arm has no executed evidence" % name
        assert arm["metric"] is None, "%s: a planned arm has no metric" % name


def test_the_stack_arm_records_why_it_is_worth_running():
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][NEW_ARM]
    rationale = arm["rationale"]
    # The declaration rests on the shuffled-conv control, so it must be named.
    assert "0.4254" in rationale
    assert "SHUFFLED-CONV" in rationale.upper()
    assert "-0.1232" in rationale


def test_the_withdrawn_kill_option_is_recorded_not_merely_dropped():
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][NEW_ARM]
    note = arm["kill_option_withdrawn"]
    assert "redundancy" in note.lower()
    assert "authority" in note.lower()
    assert "scope reduction" in arm["verdict_consequence"].lower()


def test_no_historical_measured_metric_has_moved():
    """SUBSET pin, converted from an exact-equality pin on 2026-11-02.

    **What the exact pin actually protected and what it cost.** `measured == MEASURED` kept a
    published number from drifting, which is the real guard and is preserved below. But it ALSO
    asserted that no arm may ever be added to the measured set, so it went red on every legitimate
    measurement -- and its failure message said *"a measured arm's metric moved"* when in fact the
    measured SET had gained a member and no historical metric had moved at all. A guard whose
    message misdescribes its own failure sends the next reader looking for the wrong defect.

    The subset form keeps the protection that matters (a historical number cannot change) and
    survives the three successor arms declared 2026-11-02 being measured, instead of needing three
    more re-pins.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    measured = {n: a["metric"] for n, a in arms.items() if a["status"] == al.STATUS_MEASURED}

    for name, metric in MEASURED.items():
        assert name in measured, f"{name} is no longer measured"
        assert measured[name] == metric, (
            f"{name}'s published metric moved from {metric} to {measured[name]}")
        assert arms[name]["evidence"].strip(), f"{name} lost its evidence string"

    # new measured arms are permitted, and every one must carry a metric and its evidence
    for name in sorted(set(measured) - set(MEASURED)):
        assert measured[name] is not None, f"{name} is measured with no metric"
        assert arms[name]["evidence"].strip(), f"{name} is measured with no evidence"


def test_the_arm_partition_reconciles_with_the_whole():
    """INVARIANT pin, converted from an exact-set pin on 2026-11-02.

    The exact form asserted `planned == EXPECTED_PLANNED` and `len(arms) == 13`, both of which are
    properties of one day rather than of the ledger. What is actually invariant -- and what the
    original comment said it was reaching for -- is that the three buckets the verdict partitions on
    reconcile with the declared total, with nothing falling outside them. Reading the three sets
    cannot catch an arm whose status is a typo; the arithmetic can.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]

    buckets = {al.STATUS_MEASURED: set(), al.STATUS_PLANNED: set(), al.STATUS_KILLED: set()}
    orphans = []
    for name, a in arms.items():
        if a.get("status") in buckets:
            buckets[a["status"]].add(name)
        else:
            orphans.append((name, a.get("status")))

    assert not orphans, f"arm(s) in no bucket the verdict partitions on: {orphans!r}"
    assert sum(len(v) for v in buckets.values()) == len(arms)
    assert buckets[al.STATUS_KILLED] == {KILLED_ARM}, "exactly one arm has ever been retired here"
    assert MEASURED.keys() <= buckets[al.STATUS_MEASURED], "a historical measured arm left the set"
    # a FLOOR, not an equality: arms are appended over time and the count is not the guard.
    assert len(arms) >= len(MEASURED) + 1, "the declared total fell below its historical floor"


def test_the_verdict_is_blocked_IF_AND_ONLY_IF_an_arm_is_untried():
    """BICONDITIONAL pin, converted from an exact-state pin on 2026-11-02.

    **The property R5/R6 actually rests on is an IMPLICATION, not a state:** while any arm is still
    `planned` the ceiling question is OPEN, so a cheap arm clearing a threshold cannot settle it.
    The old form asserted `verdict == "blocked:untried-arm"` together with
    `set(planned) == EXPECTED_PLANNED` and `len(planned) == 1`, which is that implication PLUS a
    snapshot of one day's arm set. The snapshot is what went red when the positional arm was
    recorded on 2026-11-01 and again when three successors were declared on 2026-11-02; the
    implication never moved and still has not.

    Asserted in BOTH directions, so neither a `met` verdict over an untried arm nor a `blocked`
    verdict over a fully-measured set can pass.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    verdict = ledger.ceiling_verdict(_criterion(ledger))
    planned = {n for n, a in arms.items() if a["status"] == al.STATUS_PLANNED}

    assert set(verdict.planned) == planned, "the verdict and the file disagree about what is planned"
    if planned:
        assert verdict.verdict == "blocked:untried-arm", (
            f"{len(planned)} arm(s) are untried ({sorted(planned)}) but the verdict reads "
            f"{verdict.verdict!r}; an untried arm must block the criterion")
        assert not verdict.met
    else:
        assert verdict.verdict in ("met", "blocked:single-arm"), (
            f"no arm is planned but the verdict reads {verdict.verdict!r}")
        assert verdict.met == (verdict.verdict == "met")


def test_best_is_the_argmax_over_measured_and_a_leader_change_must_be_EXPLICIT():
    """`best` stays an EXACT pin on purpose -- this one is not converted.

    A leader change is a substantive event, not bookkeeping, and the previous re-pin added
    `best == NEW_ARM` precisely so a shift could not pass silently. It has fired twice as designed
    (the stack arm at 0.4404 was displaced by the atlas arm at 0.5076). **If a successor arm
    out-scores 0.5076 this will fire again, and the correct response is to re-pin BEST_ARM naming
    the new leader -- never to relax the assertion.**

    The argmax clause beside it is the invariant: it catches a `best` that is not actually the
    maximum, which an exact pin alone cannot.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    verdict = ledger.ceiling_verdict(_criterion(ledger))

    measured = {n: a["metric"] for n, a in arms.items()
                if a["status"] == al.STATUS_MEASURED and a["metric"] is not None}
    assert measured
    assert verdict.best in measured
    assert measured[verdict.best] == max(measured.values()), "`best` is not the argmax"

    assert verdict.best == BEST_ARM, (
        f"the leader moved from {BEST_ARM} to {verdict.best} "
        f"({measured.get(BEST_ARM)} -> {measured[verdict.best]}). This guard exists to make that "
        f"EXPLICIT: re-pin BEST_ARM with the reason, do not relax this assertion.")
    assert verdict.best != NEW_ARM, "the leader moved; this guard exists to make that explicit"
    # Read from the LIVE ledger, not from MEASURED: BEST_ARM is now a successor arm, which is
    # deliberately absent from that historical pin (see the note under MEASURED).
    assert measured[BEST_ARM] > measured[NEW_ARM]
    assert SUCCESSOR_MEASURED in verdict.measured
    assert ATLAS_ARM in verdict.measured


def test_the_retired_arm_cannot_be_READ_as_tested_and_failed():
    """`killed` is the only status the enum offers for "cannot be instantiated", so the REASON
    carries the whole distinction. Without these clauses a later reader concludes the 4-decade
    hypothesis was tested and lost, when the arm was never run at all.
    """
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][KILLED_ARM]
    assert arm["status"] == al.STATUS_KILLED
    assert arm["evidence"] == SPAN_PROBE, "the kill must cite the probe that grounds it"
    reason = arm["reason"]
    low = reason.lower()
    # R7: scope narrow enough that the evidence supports it.
    assert "scope" in low and "this host" in low, "the kill must bound itself to this substrate"
    assert "not a refutation" in low, "it must refuse the tested-and-failed reading outright"
    assert "predictive merit" in low, "it must disclaim saying anything about the arm's merit"
    # R7: the observation that re-opens it.
    assert "revival" in low
    assert "30000" in reason or "30,000" in reason
    # It must name the successor, or the re-scope is lost with the arm.
    assert SUCCESSOR in reason


def test_the_retired_arm_was_NOT_killed_for_redundancy():
    """The exact defect R6 exists to prevent: dying because a sibling cleared the floor."""
    ledger = _live_ledger()
    reason = ledger._data["criteria"][_criterion(ledger)][KILLED_ARM]["reason"]
    assert not al._REDUNDANCY_REASON.search(reason), \
        "the stored reason matches the redundancy pattern kill_arm refuses"
    assert "substrate-gated" in reason.lower()


def test_the_successor_arm_records_WHY_it_is_worth_running():
    """Planned arms carried no reason at all, which is why `blocked:untried-arm` could only ever
    name them. The successor sets `rationale` directly, since `declare_arm` takes no such argument.
    """
    ledger = _live_ledger()
    arm = ledger._data["criteria"][_criterion(ledger)][SUCCESSOR]
    assert arm["status"] == al.STATUS_MEASURED, "measured 2026-10-31; see this module's docstring"
    r = arm["rationale"]
    assert "2.22" in r and "1.70" in r, "it must state the span it buys over the measured sibling"
    assert "0.7928" in r, "it must name the FRAGMENT rung's own ceiling"
    assert "0.1313" in r, "it must name the fragment-rung incumbent it has to beat"
    assert "never be read against" in r.lower(), \
        "it must refuse the cross-rung comparison the census ceiling_rule forbids"
    assert "pad/crop" in r, "the variable-length design point must be declared pre-run"
    # And the rationale must NOT have leaked into the field reserved for executed evidence.
    # STRENGTHENED 2026-10-31: the arm is now measured, so `evidence` must hold EXECUTED evidence
    # -- asserting emptiness would no longer protect anything. The no-leak invariant is kept as a
    # direct check, which is a tighter statement than the old proxy of an empty string.
    ev = arm["evidence"]
    assert ev.strip(), "a measured arm must carry executed evidence"
    assert "glm_fragment_pyramid_arm_kaggle_full_2026-10-28.json" in ev,         "the evidence must cite the run that produced the number"
    assert "fragment_null_270seed_result_2026-10-31" in ev,         "and the null resolution that made the verdict gradeable at all"
    assert r not in ev, "the rationale must not have been copied into the evidence field"
    assert "2.22" not in ev or "NO_GAIN" in ev, "evidence must read as a result, not a plan"


def test_the_eukaryotic_arm_is_still_present():
    """Step 2's north-star step-4 cross-reference must not dangle."""
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    assert "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN" in arms


def test_kill_arm_still_refuses_a_redundancy_reason(sandbox):
    """The regression guard the withdrawal depends on.

    If this refusal ever stops firing, the kill option silently becomes available
    again and a sibling clearing the floor could cancel an unrun arm.
    """
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    with pytest.raises(al.ArmsError, match="REFUSED"):
        ledger.kill_arm(
            crit, "multiresolution-pyramid-4-decades-UNRUN",
            reason="redundant now that the hybrid arm cleared the floor",
            evidence="the hybrid arm measured 0.4254",
        )


def test_kill_arm_still_refuses_without_evidence(sandbox):
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    with pytest.raises(al.ArmsError, match="evidence"):
        ledger.kill_arm(crit, "promoteratlas-base-embeddings-frozen-150bp-padded-to-200",
                        reason="not worth it", evidence="")


def test_the_human_arms_metric_is_not_read_as_best():
    """The eukaryotic arm's metric is a HUMAN Spearman in a criterion whose other metrics are
    BACTERIAL, and `ceiling_verdict` picks `best` as a max over metrics. So `best` is only
    meaningful here while no human arm outscores the bacterial leader.

    Harmless TODAY -- 0.4853 < 0.5076 -- and this guard exists so it cannot become harmful
    silently. If a future human arm outscores the bacterial best, this fires, and the answer is NOT
    to relax it: it is to move the arm to its own criterion, which is a scope change and therefore
    the user's call. The structural issue is annotated on the arm itself rather than papered over.
    """
    ledger = _live_ledger()
    arms = ledger._data["criteria"][_criterion(ledger)]
    verdict = ledger.ceiling_verdict(_criterion(ledger))

    assert verdict.best == BEST_ARM, "`best` must stay the bacterial leader"
    # Both metrics read from the LIVE ledger: BEST_ARM is now a successor arm and so is
    # deliberately absent from the historical MEASURED pin (see the note under it).
    assert arms[EUK_ARM]["metric"] < arms[verdict.best]["metric"], (
        "a human-substrate metric now exceeds the bacterial best, so `best` is comparing across "
        "substrates; move the arm to its own criterion rather than relaxing this")
    # The incommensurability is RECORDED on the arm, not merely known here.
    note = arms[EUK_ARM].get("metric_not_commensurable_with_siblings", "")
    assert note.strip(), "the arm lost the annotation that its metric is not commensurable"
    assert "bacterial" in note.lower() and "best" in note.lower()


def test_killing_every_unrun_arm_would_change_the_verdict(sandbox):
    """Shows the withdrawal is load-bearing, not decorative.

    With evidence and a non-redundancy reason the kills are permitted -- which is
    exactly why a no-evidence kill has to route to user authority rather than being
    something this plan can do on its own.
    """
    ledger = al.ArmsLedger(cwd=str(sandbox))
    crit = _criterion(ledger)
    # DERIVED from the sandbox rather than from a pinned set, converted 2026-11-02: the old form
    # killed only EXPECTED_PLANNED, so once three successor arms were declared it killed one of
    # four and the verdict correctly stayed blocked -- the test failed while the property it
    # describes still held.
    unrun = sorted(n for n, a in ledger._data["criteria"][crit].items()
                   if a.get("status") == al.STATUS_PLANNED)
    if not unrun:
        pytest.skip("no arm is planned, so the withdrawal has no subject to be load-bearing for")
    for name in unrun:
        ledger.kill_arm(crit, name, reason="measured failure on its own run",
                        evidence="synthetic executed evidence (test only)")
    assert ledger.ceiling_verdict(crit).verdict != "blocked:untried-arm"

    # And the real ledger is untouched -- asserted as the biconditional, not a fixed string.
    live = _live_ledger()
    live_arms = live._data["criteria"][_criterion(live)]
    live_planned = {n for n, a in live_arms.items() if a["status"] == al.STATUS_PLANNED}
    live_verdict = live.ceiling_verdict(_criterion(live)).verdict
    assert (live_verdict == "blocked:untried-arm") == bool(live_planned)
    assert live_planned, "the sandbox kills must not have reached the live ledger"
