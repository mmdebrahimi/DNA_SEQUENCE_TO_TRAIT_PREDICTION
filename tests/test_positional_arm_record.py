"""Tests for recording the positional arm and rendering its memo.

Three things are pinned here. The ledger write must PRESERVE the arm's off-schema annotation keys,
which is the whole reason it is a direct edit rather than a `measure_arm` call. The memo generator
must REFUSE an `INDETERMINATE_*` verdict, because a memo on a withheld verdict reads as a finding.
And the arms partition must reconcile with the whole -- measured plus planned plus killed must equal
the declared arm count, checked by arithmetic rather than by reading the three lists.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "soraya_runs/arms-ledger.json"
RESULT = REPO / "wiki/positional_arm_result_2026-11-01.json"
BAR = REPO / "wiki/positional_arm_acceptance_bar.json"
MEMO = REPO / "wiki/positional_arm_result_2026-11-01.md"

ARM = "positional-non-mean-pooled-promoteratlas-150bp"
SORAYA_SCRIPTS = Path(r"C:/my_skills/skills/soraya/scripts")
ARM_FIELDS = {"name", "cost_tier", "status", "metric", "evidence", "declared_at", "resolved_at"}

#: The keys that are NOT on the Arm dataclass and must survive every write.
OFF_SCHEMA = ("declared_from", "rationale")


def _ledger() -> dict:
    return json.loads(LEDGER.read_text(encoding="utf-8"))


def _crit(d: dict) -> str:
    crits = [k for k in d["criteria"] if k.startswith("ceiling:")]
    assert len(crits) == 1, "expected exactly one ceiling criterion, found %r" % crits
    return crits[0]


def _arms(d: dict) -> dict:
    return d["criteria"][_crit(d)]


# --------------------------------------------------------------------------------------------------
# the arm's annotation keys must survive
# --------------------------------------------------------------------------------------------------
def test_the_arm_is_declared_on_the_single_ceiling_criterion():
    arms = _arms(_ledger())
    assert ARM in arms
    assert arms[ARM]["status"] in ("planned", "measured")


def test_the_arm_carries_off_schema_annotation_keys():
    """This is WHY the ledger write is a direct edit; if it ever stops being true, say so."""
    a = _arms(_ledger())[ARM]
    extra = sorted(set(a) - ARM_FIELDS)
    assert extra, "the arm has no off-schema keys, so the direct-edit rationale no longer holds"
    for k in OFF_SCHEMA:
        assert k in a, "%r was expected among the arm's annotation keys; found %r" % (k, extra)


@pytest.mark.skipif(not (SORAYA_SCRIPTS / "arms_ledger.py").is_file(),
                    reason="soraya arms_ledger.py not on this host")
def test_measure_arm_would_raise_on_this_arm_which_is_why_the_write_is_direct():
    """Non-vacuity for the design choice: the API path genuinely fails on these keys.

    `measure_arm` ends in `Arm(**base)` with no field filter, unlike `kill_arm`, which filters to
    `Arm.__dataclass_fields__`. So an arm carrying annotation keys raises on the RETURN construction
    -- after the dict has already been mutated, which is the other reason not to route through it.
    """
    sys.path.insert(0, str(SORAYA_SCRIPTS))
    try:
        import arms_ledger as mod  # type: ignore
    finally:
        sys.path.pop(0)

    a = dict(_arms(_ledger())[ARM])
    a.update(status="measured", metric=0.5, evidence="x", resolved_at="now")
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        mod.Arm(**a)
    # and the filtered construction -- what kill_arm does and what a direct edit preserves -- works
    assert mod.Arm(**{k: v for k, v in a.items() if k in mod.Arm.__dataclass_fields__})


def test_a_direct_edit_round_trips_the_annotation_keys(tmp_path):
    """Simulate the Step 6 write on a COPY and assert nothing is lost."""
    d = _ledger()
    arms = _arms(d)
    before = dict(arms[ARM])

    arms[ARM].update(status="measured", metric=0.4242, evidence="wiki/some_artifact.json",
                     resolved_at="2026-11-02T00:00:00Z")
    p = tmp_path / "ledger.json"
    p.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    back = json.loads(p.read_text(encoding="utf-8"))
    after = back["criteria"][_crit(back)][ARM]

    for k in OFF_SCHEMA:
        assert after[k] == before[k], "annotation key %r was altered by the write" % k
    assert set(before) <= set(after), "the write must never DROP a key"
    assert after["status"] == "measured" and after["metric"] == 0.4242
    assert after["evidence"]


# --------------------------------------------------------------------------------------------------
# the partition reconciles with the whole
# --------------------------------------------------------------------------------------------------
def test_the_arms_partition_reconciles_with_the_declared_arm_count():
    """Checked by ARITHMETIC. Reading the three lists would not catch an arm in none of them."""
    arms = _arms(_ledger())
    buckets = {"measured": 0, "planned": 0, "killed": 0}
    unknown = []
    for name, a in arms.items():
        s = a.get("status")
        if s in buckets:
            buckets[s] += 1
        else:
            unknown.append((name, s))
    assert not unknown, "arm(s) in no bucket the verdict partitions on: %r" % unknown
    assert sum(buckets.values()) == len(arms), (
        "partition %r sums to %d but %d arms are declared"
        % (buckets, sum(buckets.values()), len(arms)))
    assert buckets["measured"] >= 2, "a ceiling needs >=2 measured arms for `best` to mean anything"


def test_every_measured_arm_carries_a_metric_and_evidence():
    for name, a in _arms(_ledger()).items():
        if a.get("status") != "measured":
            continue
        assert a.get("metric") is not None, "%s is measured with no metric" % name
        assert str(a.get("evidence", "")).strip(), "%s is measured with no evidence" % name


def test_best_is_the_max_over_measured_metrics():
    arms = _arms(_ledger())
    measured = {n: a["metric"] for n, a in arms.items() if a.get("status") == "measured"}
    assert measured
    best = max(measured, key=lambda n: measured[n])
    assert measured[best] == max(measured.values())


@pytest.mark.skipif(not (SORAYA_SCRIPTS / "arms_ledger.py").is_file(),
                    reason="soraya arms_ledger.py not on this host")
def test_the_live_ceiling_verdict_agrees_with_the_partition_it_decomposes():
    """A verdict that disagreed with its own arm counts would be reading a different ledger.

    This is the check that caught a `blocked:no-arms` over twelve declared arms once: the loader had
    silently resolved a nonexistent path, and a nonexistent path and an empty ledger return the
    IDENTICAL shape.
    """
    sys.path.insert(0, str(SORAYA_SCRIPTS))
    try:
        from arms_ledger import ArmsLedger  # type: ignore
    finally:
        sys.path.pop(0)

    led = ArmsLedger(os.getcwd())
    crit = _crit(led._data)
    v = led.ceiling_verdict(crit).to_dict()
    arms = _arms(led._data)

    assert len(arms) > 0, "the loader resolved an EMPTY arm set, which is the git-bash path trap"
    assert len(v["measured"]) + len(v["planned"]) + len(v["killed"]) == len(arms)
    assert set(v["measured"]).isdisjoint(v["planned"])
    assert set(v["measured"]).isdisjoint(v["killed"])
    if v["planned"]:
        assert v["verdict"] == "blocked:untried-arm"
    else:
        assert v["verdict"] in ("met", "blocked:single-arm")
    if v["measured"]:
        assert v["best"] in v["measured"]


# --------------------------------------------------------------------------------------------------
# the memo generator's refusal
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("verdict", [
    "INDETERMINATE_PARITY_FAIL", "INDETERMINATE_SPLIT_MISMATCH", "INDETERMINATE_NULL_FAIL",
])
def test_the_memo_refuses_to_render_a_withheld_verdict(verdict, tmp_path, capsys):
    from scripts.positional_arm_memo import main

    stub = tmp_path / "res.json"
    stub.write_text(json.dumps({"verdict": verdict, "date": "2026-11-01"}), encoding="utf-8")
    out = tmp_path / "memo.md"
    assert main(["--result", str(stub), "--out", str(out)]) == 2
    assert "REFUSING" in capsys.readouterr().out
    assert not out.exists(), "a refused memo must not leave a file behind"


def test_the_memo_refuses_an_absent_result_artifact(tmp_path, capsys):
    from scripts.positional_arm_memo import main

    assert main(["--result", str(tmp_path / "nope.json"), "--out", str(tmp_path / "m.md")]) == 2
    assert "run the scorer first" in capsys.readouterr().out


@pytest.mark.skipif(not RESULT.is_file(), reason="result artifact absent; run the scorer")
def test_the_memo_renders_every_figure_from_the_artifact():
    """Spot-check that the headline numbers in the memo are the artifact's, not restated."""
    from scripts.positional_arm_memo import render

    res = json.loads(RESULT.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    md = render(res, bar, None)

    assert res["verdict"] in md
    assert "%.4f" % bar["inherited_not_chosen"]["comparator"]["atlas_median"] in md
    for name, r in res["rungs"].items():
        assert name in md
        if not r.get("skipped"):
            assert "%.4f" % r["median_spearman"] in md


@pytest.mark.skipif(not RESULT.is_file(), reason="result artifact absent; run the scorer")
def test_the_memo_is_ordered_by_required_reader_action():
    from scripts.positional_arm_memo import render

    md = render(json.loads(RESULT.read_text(encoding="utf-8")),
                json.loads(BAR.read_text(encoding="utf-8")), None)
    i_dec = md.index("## DECISIONS NEEDED")
    i_cor = md.index("## CORRECTIONS")
    i_res = md.index("## The result")
    assert i_dec < i_cor < i_res, "decisions, then corrections, then findings"


@pytest.mark.skipif(not RESULT.is_file(), reason="result artifact absent; run the scorer")
def test_the_memo_carries_both_scoping_corrections_and_the_honest_limits():
    from scripts.positional_arm_memo import render

    res = json.loads(RESULT.read_text(encoding="utf-8"))
    md = render(res, json.loads(BAR.read_text(encoding="utf-8")), None)
    assert len(res["scoping_corrections"]) == 2
    for c in res["scoping_corrections"]:
        assert c[:60] in md
    for lim in res["honest_limits"]:
        assert lim[:60] in md


@pytest.mark.skipif(not MEMO.is_file(), reason="memo not rendered yet")
def test_the_rendered_memo_on_disk_matches_a_fresh_render():
    """No figure may be hand-edited into the memo after generation."""
    from scripts.positional_arm_memo import render

    res = json.loads(RESULT.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    on_disk = MEMO.read_text(encoding="utf-8")
    # the DECISIONS section depends on the live ceiling verdict, so compare the parts that do not
    fresh_tail = render(res, bar, None).split("## The result", 1)[1]
    assert on_disk.split("## The result", 1)[1].strip() == fresh_tail.strip()
