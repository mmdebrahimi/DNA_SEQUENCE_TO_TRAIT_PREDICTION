"""Guards for the chained-session reconciler.

WHY IT EXISTS. A run record is a summary written by the same process that wrote the per-epoch curve
files, so reading the record alone cannot catch a session that reported more progress than it made.
The reconciler recomputes the totals from the curves independently. These tests check that it
actually refuses the disagreements it exists to find -- a reconciler that always says LINK_OK is
indistinguishable from no reconciler at all.

THE STALL CASE IS THE ONE WITH TEETH. A session that seeds correctly and adds no epochs costs a GPU
session and quota hours for nothing, and it looks exactly like a healthy bounded stop unless the
totals are compared across links. The weekly GPU quota is the binding constraint on this whole arm
(30.00 h against a 32.26 h prefix, measured 2026-11-15), so a stall that is not caught is paid for
in the resource there is least of.

Pure filesystem checks against tmp_path: no network, no Kaggle, no torch.
"""
from __future__ import annotations

import ast
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RECON = ROOT / "scripts" / "kaggle" / "reconcile_finetune_session.py"
RUNNER = ROOT / "scripts" / "finetune_arm_run.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_recon", RECON)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _link(base: Path, *, seeds: dict[int, int], nulls: dict[int, int] | None = None,
          record_seeds: dict[int, int] | None = None, status: str = "BOUNDED_BY_WALL_CLOCK",
          with_record: bool = True, smoke_epochs: int = 0) -> Path:
    """A synthetic pulled tree. `seeds` maps seed -> epochs on disk."""
    d = base / "finetune_checkpoints"
    d.mkdir(parents=True, exist_ok=True)
    for s, n in seeds.items():
        (d / ("seed_%d_curve.json" % s)).write_text(json.dumps(
            {"val_spearman_by_epoch": [0.1] * n, "epochs_done": n, "horizon_epochs": 38}),
            encoding="utf-8")
    for s, n in (nulls or {}).items():
        (d / ("null_%d_curve.json" % s)).write_text(json.dumps(
            {"val_spearman_by_epoch": [0.0] * n, "epochs_done": n, "shuffled_labels": True}),
            encoding="utf-8")
    if smoke_epochs:
        (d / "smoke_0_curve.json").write_text(json.dumps(
            {"val_spearman_by_epoch": [0.0] * smoke_epochs, "epochs_done": smoke_epochs}),
            encoding="utf-8")
    if with_record:
        rs = record_seeds if record_seeds is not None else seeds
        # `per_seed`, NOT `seeds`. An earlier version of this fixture invented `seeds`, which made
        # every test here pass against a schema no producer emits while the reconciler's central
        # comparison was dead on real input. The key name is the contract.
        (base / "finetune_arm_run_kaggle.json").write_text(json.dumps(
            {"status": status, "why_status": "synthetic",
             "phases": [{"phase": "1-seeds-to-frozen-budget", "bounded_here": True}],
             "horizon_divergence": {"horizon_epochs_supplied": 38},
             "per_seed": {str(s): {"val_spearman_by_epoch": [0.1] * n} for s, n in rs.items()}}),
            encoding="utf-8")
    return base


def test_the_reconciler_uses_the_SINGLE_pattern_definition():
    src = RECON.read_text(encoding="utf-8")
    assert "from scripts.finetune_arm_run import CKPT_NAME_RE" in src
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Assign):
            assert "CKPT_NAME_RE" not in {t.id for t in node.targets if isinstance(t, ast.Name)}
    assert RUNNER.read_text(encoding="utf-8").count("CKPT_NAME_RE = re.compile(") == 1


def test_a_healthy_link_reconciles(mod, tmp_path):
    out = mod.reconcile(_link(tmp_path, seeds={0: 3, 1: 2}))
    assert out["reconciles"] is True
    assert out["from_curves"]["total_epochs"] == 5
    assert out["status"] == "BOUNDED_BY_WALL_CLOCK"
    assert out["horizon_divergence_recorded"] is True
    assert out["horizon_epochs_supplied"] == 38
    # The comparison must have HAPPENED. `reconciles` was previously computable from an empty
    # comparison, which is how a dead reader reported success on the first real record.
    assert out["n_record_seeds_compared"] == 2, (
        "reconciles=True with %r seeds compared is a vacuous pass"
        % out.get("n_record_seeds_compared"))


def test_the_smoke_curve_is_NOT_counted_as_progress(mod, tmp_path):
    """The rate measurement writes a curve shaped exactly like a real one.

    Counting it would credit the chain with an epoch that was deliberately thrown away, and would
    make a stalled link look like it advanced.
    """
    out = mod.reconcile(_link(tmp_path, seeds={0: 3}, smoke_epochs=1))
    assert out["from_curves"]["total_epochs"] == 3, "the smoke curve was counted"
    assert "smoke_0_curve.json" not in out["from_curves"]["epochs_by_seed"]


def test_nulls_are_counted_but_kept_SEPARATE(mod, tmp_path):
    """Null epochs cost the same quota as real ones, so they count toward the prefix -- but they are
    not seed progress, and a reader has to be able to tell which is which."""
    out = mod.reconcile(_link(tmp_path, seeds={0: 7}, nulls={0: 7, 1: 7}))
    fc = out["from_curves"]
    assert fc["total_epochs"] == 21
    assert fc["n_seed_curves"] == 1 and fc["n_null_curves"] == 2
    assert set(fc["epochs_by_null"]) == {"null_0_curve.json", "null_1_curve.json"}


def test_a_record_claiming_more_than_disk_is_REFUSED(mod, tmp_path):
    out = mod.reconcile(_link(tmp_path, seeds={0: 1}, record_seeds={0: 7}))
    assert out["verdict"] == "RECORD_DISAGREES_WITH_DISK"
    assert out["reconciles"] is False
    assert "record says 7" in out["record_vs_disk_mismatches"][0]


def test_a_record_naming_a_seed_with_no_curve_file_is_REFUSED(mod, tmp_path):
    out = mod.reconcile(_link(tmp_path, seeds={0: 3}, record_seeds={0: 3, 5: 4}))
    assert out["verdict"] == "RECORD_DISAGREES_WITH_DISK"
    assert any("seed_5_curve.json" in m and "no curve file" in m
               for m in out["record_vs_disk_mismatches"]), out["record_vs_disk_mismatches"]


def test_a_missing_run_artifact_still_counts_the_epochs_that_survived(mod, tmp_path):
    """The point of writing curves inside the epoch loop: a session that dies late loses its record,
    not its work. The verdict must say the record is missing without discarding the epochs."""
    out = mod.reconcile(_link(tmp_path, seeds={0: 4}, with_record=False))
    assert out["verdict"] == "NO_RUN_ARTIFACT"
    assert out["from_curves"]["total_epochs"] == 4
    assert out["run_artifact_present"] is False
    assert "no record" in out["why"]


def test_two_curves_of_the_same_name_REFUSE_rather_than_double_count(mod, tmp_path):
    a = tmp_path / "a" / "finetune_checkpoints"
    b = tmp_path / "b" / "finetune_checkpoints"
    a.mkdir(parents=True)
    b.mkdir(parents=True)
    for d in (a, b):
        (d / "seed_0_curve.json").write_text('{"val_spearman_by_epoch":[0.1]}', encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        mod.read_curves(tmp_path)
    assert "double-count" in str(e.value)


def test_two_run_artifacts_REFUSE_rather_than_mixing_two_sessions(mod, tmp_path):
    (tmp_path / "x").mkdir()
    (tmp_path / "y").mkdir()
    for sub in ("x", "y"):
        (tmp_path / sub / "finetune_arm_run_kaggle.json").write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        mod.find_one(tmp_path, "finetune_arm_run_kaggle.json")
    assert "one session at a time" in str(e.value)


def test_a_truncated_curve_REFUSES_rather_than_counting_zero(mod, tmp_path):
    d = tmp_path / "finetune_checkpoints"
    d.mkdir(parents=True)
    (d / "seed_0_curve.json").write_text("{not json", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        mod.read_curves(tmp_path)
    assert "must not be" in str(e.value) and "zero epochs" in str(e.value)


# ---------------------------------------------------------------------------------------------------
# the stall check, and its exit codes
# ---------------------------------------------------------------------------------------------------
def test_a_stalled_link_is_caught_and_stops_the_chain(mod, tmp_path, capsys):
    prev = _link(tmp_path / "l1", seeds={0: 3, 1: 2})
    prev_json = tmp_path / "l1.json"
    rc = mod.main(["--pulled", str(prev), "--out", str(prev_json)])
    assert rc == 0
    capsys.readouterr()

    same = _link(tmp_path / "l2", seeds={0: 3, 1: 2})
    rc = mod.main(["--pulled", str(same), "--against", str(prev_json)])
    assert rc == 4, "a stalled link must not exit 0"
    out = json.loads(capsys.readouterr().out)
    assert out["verdict"] == "STALLED_NO_EPOCHS_ADDED"
    assert out["epochs_added_by_this_link"] == 0
    assert "must STOP" in out["why"]


def test_a_progressing_link_passes_the_stall_check(mod, tmp_path, capsys):
    prev = _link(tmp_path / "l1", seeds={0: 3})
    prev_json = tmp_path / "l1.json"
    mod.main(["--pulled", str(prev), "--out", str(prev_json)])
    capsys.readouterr()
    more = _link(tmp_path / "l2", seeds={0: 7, 1: 4})
    rc = mod.main(["--pulled", str(more), "--against", str(prev_json)])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["verdict"] == "LINK_OK"
    assert out["epochs_added_by_this_link"] == 8


def test_the_exit_codes_are_distinct_so_a_loop_can_branch_without_parsing_prose(mod, tmp_path,
                                                                               capsys):
    codes = {}
    codes["ok"] = mod.main(["--pulled", str(_link(tmp_path / "a", seeds={0: 1}))])
    capsys.readouterr()
    codes["no_artifact"] = mod.main(
        ["--pulled", str(_link(tmp_path / "b", seeds={0: 1}, with_record=False))])
    capsys.readouterr()
    codes["mismatch"] = mod.main(
        ["--pulled", str(_link(tmp_path / "c", seeds={0: 1}, record_seeds={0: 9}))])
    capsys.readouterr()
    assert codes == {"ok": 0, "no_artifact": 2, "mismatch": 3}, codes
    assert len(set(codes.values())) == 3, "the exit codes collapsed; a loop could not branch on them"


# ===================================================================================================
# the vacuity guard, and the real record -- added after the first real session exposed a dead reader
# ===================================================================================================
REAL_RECORD = ROOT / "wiki" / "finetune_arm_session1_record_2026-11-15.json"
REAL_CURVE = ROOT / "wiki" / "finetune_arm_session1_seed_0_curve_2026-11-15.json"


def test_a_record_whose_per_seed_section_cannot_be_READ_is_VACUOUS_not_OK(mod, tmp_path):
    """The defect this guard exists for, reproduced deliberately.

    A record carrying the per-seed data under a key this script does not know yields nothing to
    compare. Reporting LINK_OK there manufactures confidence in the one claim the chain rests on --
    that progress accumulates -- so it must refuse instead.
    """
    base = _link(tmp_path, seeds={0: 3}, with_record=False)
    (base / "finetune_arm_run_kaggle.json").write_text(json.dumps(
        {"status": "BOUNDED_BY_WALL_CLOCK",
         "seeds": {"0": {"val_spearman_by_epoch": [0.1, 0.2, 0.3]}}}),  # the WRONG key, on purpose
        encoding="utf-8")
    out = mod.reconcile(base)
    assert out["verdict"] == "RECONCILIATION_VACUOUS", out
    assert out["reconciles"] is False, "a vacuous reconciliation must not report reconciles=True"
    assert out["n_record_seeds_compared"] == 0
    assert "checked nothing" in out["why"]


def test_a_curve_on_disk_the_record_never_mentions_is_a_MISMATCH(mod, tmp_path):
    """The other direction. Only checking record -> disk misses a seed the record forgot."""
    out = mod.reconcile(_link(tmp_path, seeds={0: 3, 1: 2}, record_seeds={0: 3}))
    assert out["verdict"] == "RECORD_DISAGREES_WITH_DISK"
    assert any("seed_1_curve.json" in m and "does not mention" in m
               for m in out["record_vs_disk_mismatches"]), out["record_vs_disk_mismatches"]


def test_boundedness_is_derived_because_the_obvious_field_does_not_EXIST(mod, tmp_path):
    """`run()` has no top-level `bounded_by_wall_clock`; the plan's criterion named the wrong level.

    `train_one_seed` sets that field on a PER-SEED record. At the run level the property lives in
    `status` and in `phases[i].bounded_here`, so both are read and the absence of the top-level field
    is REPORTED rather than silently rendered as null.
    """
    out = mod.reconcile(_link(tmp_path, seeds={0: 3}))
    assert out["bounded"] is True
    sig = out["bounded_signals"]
    assert sig["status_is_bounded"] is True
    assert sig["phases_reporting_bounded_here"] == ["1-seeds-to-frozen-budget"]
    assert sig["top_level_bounded_field_present"] is False, (
        "if run() ever gains a top-level bounded_by_wall_clock, read it directly and simplify this")


# ---------------------------------------------------------------------------------------------------
# the only test here whose input I did not invent
# ---------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not REAL_RECORD.is_file(), reason="the committed session-1 record is absent")
def test_the_reader_works_on_the_REAL_committed_session_record(mod, tmp_path):
    """Reads the verbatim artifact pulled from the first Kaggle run of `run()`.

    This is the test that would have caught the original defect. Every other case here is built from
    a fixture I wrote, and a fixture that speaks my invented schema cannot tell me the reader is
    broken -- it tells me the reader matches my invention.
    """
    import shutil
    d = tmp_path / "finetune_checkpoints"
    d.mkdir(parents=True)
    shutil.copyfile(REAL_CURVE, d / "seed_0_curve.json")
    shutil.copyfile(REAL_RECORD, tmp_path / "finetune_arm_run_kaggle.json")

    out = mod.reconcile(tmp_path)
    assert out["verdict"] == "LINK_OK", out
    assert out["n_record_seeds_compared"] == 1, (
        "the reader compared %r seeds against the REAL record; a zero here is the original defect"
        % out.get("n_record_seeds_compared"))
    assert out["reconciles"] is True
    assert out["record_vs_disk_mismatches"] == []
    assert out["from_curves"]["total_epochs"] == 4
    assert out["status"] == "BOUNDED_BY_WALL_CLOCK"
    assert out["bounded"] is True
    assert out["horizon_divergence_recorded"] is True
    assert out["horizon_epochs_supplied"] == 38


@pytest.mark.skipif(not REAL_RECORD.is_file(), reason="the committed session-1 record is absent")
def test_the_real_record_uses_per_seed_and_has_no_top_level_bounded_field():
    """Pins the two schema facts the reader depends on, read off the producer's own output."""
    rec = json.loads(REAL_RECORD.read_text(encoding="utf-8"))
    assert "per_seed" in rec, "the record key the reader depends on is gone; re-read before fixing"
    assert "seeds" not in rec, (
        "a top-level `seeds` key appeared; if the schema now carries both, the reader must say which "
        "one it trusts rather than silently preferring one")
    assert "bounded_by_wall_clock" not in rec, (
        "run() now has a top-level bounded field; read it directly and simplify the derivation")
    assert rec["horizon_divergence"]["horizon_epochs_supplied"] == 38
    # And the per-seed record carries the curve the reconciliation compares.
    assert isinstance(rec["per_seed"]["0"]["val_spearman_by_epoch"], list)
