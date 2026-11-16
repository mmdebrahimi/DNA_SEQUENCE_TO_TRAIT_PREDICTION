"""Guards for the checkpoint-dataset uploader -- the other end of the chained-session round trip.

THE FAILURE THIS PREVENTS IS SILENT. The uploader has to stage exactly the files the Kaggle kernel
will accept. If the two disagree, nothing errors: the uploader reports success, the kernel seeds
zero, and the kernel's refusal blames the mount. That misdiagnosis is not hypothetical -- it cost two
probe pushes on 2026-11-15 before the selection rule was single-sourced into
`scripts/finetune_arm_run.py`, the module whose `_ckpt_path` and `_curve_path` produce the names.

So the first test here is the one that matters: the uploader must not own a copy of the pattern.

Pure filesystem + structural checks against tmp_path: no network, no Kaggle, no torch.
"""
from __future__ import annotations

import ast
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
UPLOADER = ROOT / "scripts" / "kaggle" / "upload_finetune_checkpoints.py"
RUNNER = ROOT / "scripts" / "finetune_arm_run.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_upl", UPLOADER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture
def ckpt_dir(tmp_path: Path) -> Path:
    """A realistic pulled directory: two real seeds, a null, and three things to exclude."""
    d = tmp_path / "pulled" / "finetune_checkpoints"
    d.mkdir(parents=True)
    (d / "seed_0_curve.json").write_text(json.dumps(
        {"epochs_done": 3, "horizon_epochs": 38, "complete": False,
         "shuffled_labels": False, "val_spearman_by_epoch": [0.1, 0.2, 0.3]}), encoding="utf-8")
    (d / "seed_0.pt").write_bytes(b"weights")
    (d / "null_0_curve.json").write_text(json.dumps(
        {"epochs_done": 7, "horizon_epochs": 7, "complete": True,
         "shuffled_labels": True}), encoding="utf-8")
    (d / "smoke_0.pt").write_bytes(b"throwaway")
    (d / "seed_9.pt.tmp").write_bytes(b"half-written")
    (d / "README.txt").write_text("notes", encoding="utf-8")
    return tmp_path / "pulled"


# ---------------------------------------------------------------------------------------------------
# the single-source guard -- the reason this file leads with it
# ---------------------------------------------------------------------------------------------------
def test_the_uploader_imports_the_pattern_and_never_defines_its_own():
    src = UPLOADER.read_text(encoding="utf-8")
    assert "from scripts.finetune_arm_run import CKPT_NAME_RE" in src, (
        "the uploader must import the single definition from the module that owns the filenames")
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            names = {t.id for t in node.targets if isinstance(t, ast.Name)}
            assert "CKPT_NAME_RE" not in names, (
                "the uploader defines its own copy of the checkpoint-name pattern; two copies "
                "across a round trip fail by PASSING")
    # And it must be the same object the runner defines, not a same-looking literal.
    assert RUNNER.read_text(encoding="utf-8").count("CKPT_NAME_RE = re.compile(") == 1


def test_the_uploaders_pattern_IS_the_runners(mod):
    runner_src = RUNNER.read_text(encoding="utf-8")
    pat = None
    for node in ast.walk(ast.parse(runner_src)):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "CKPT_NAME_RE" for t in node.targets)):
            pat = ast.literal_eval(node.value.args[0])
    assert pat, "could not read the runner's pattern"
    assert mod.CKPT_NAME_RE.pattern == pat


# ---------------------------------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------------------------------
def test_classify_accepts_real_checkpoints_and_reports_every_exclusion(mod, ckpt_dir):
    taken, skipped = mod.classify(ckpt_dir)
    assert set(taken) == {"seed_0_curve.json", "seed_0.pt", "null_0_curve.json"}
    assert skipped["smoke (excluded by design)"] == ["finetune_checkpoints/smoke_0.pt"]
    assert skipped["unfinished write (.tmp)"] == ["finetune_checkpoints/seed_9.pt.tmp"]
    assert skipped["not a checkpoint name"] == ["finetune_checkpoints/README.txt"]


def test_the_smoke_exclusion_is_the_load_bearing_one(mod, tmp_path):
    """A directory holding ONLY smoke output must yield nothing, and say why.

    The rate measurement trains a throwaway epoch at seed 0 with the real hyperparameters, split and
    builder, so its curve is indistinguishable in shape from a real one -- and `advance()` reads any
    curve that already reaches the requested depth. Uploading smoke output would make a discarded
    measurement epoch the arm's trajectory prefix.
    """
    d = tmp_path / "smoke_only"
    d.mkdir()
    (d / "smoke_0_curve.json").write_text("{}", encoding="utf-8")
    (d / "smoke_0.pt").write_bytes(b"x")
    taken, skipped = mod.classify(d)
    assert taken == {}
    assert len(skipped["smoke (excluded by design)"]) == 2


def test_classify_searches_recursively(mod, ckpt_dir):
    """The pulled tree nests the checkpoints one level down, so a non-recursive scan finds nothing."""
    taken, _ = mod.classify(ckpt_dir)
    assert taken, "nothing found; the search is not recursing into the pulled subdirectory"
    assert all(p.parent.name == "finetune_checkpoints" for p in taken.values())


def test_duplicate_basenames_across_subdirectories_REFUSE(mod, tmp_path):
    """The dataset is flat and the kernel copies flat, so one of two same-named curves would win."""
    a = tmp_path / "a" / "finetune_checkpoints"
    b = tmp_path / "b" / "finetune_checkpoints"
    a.mkdir(parents=True)
    b.mkdir(parents=True)
    (a / "seed_0_curve.json").write_text('{"epochs_done": 3}', encoding="utf-8")
    (b / "seed_0_curve.json").write_text('{"epochs_done": 7}', encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        mod.classify(tmp_path)
    msg = str(e.value)
    assert "both named" in msg and "seed_0_curve.json" in msg
    assert "silently win" in msg


# ---------------------------------------------------------------------------------------------------
# the resumable-state summary
# ---------------------------------------------------------------------------------------------------
def test_summarize_reports_what_the_next_session_can_skip(mod, ckpt_dir):
    taken, _ = mod.classify(ckpt_dir)
    s = mod.summarize(taken)
    assert s["total_epochs_on_disk"] == 10, s
    assert s["weight_files"] == ["seed_0.pt"]
    assert s["curves"]["seed_0_curve.json"]["epochs_done"] == 3
    assert s["curves"]["seed_0_curve.json"]["complete"] is False
    assert s["curves"]["null_0_curve.json"]["complete"] is True
    assert s["curves"]["null_0_curve.json"]["shuffled_labels"] is True, (
        "the null's shuffled-label flag must survive; a null mistaken for a real seed would be a "
        "silent corruption of the arm's result")


def test_summarize_survives_an_unreadable_curve(mod, tmp_path):
    """A truncated JSON must be reported, never crash the upload or pass as zero epochs."""
    d = tmp_path / "finetune_checkpoints"
    d.mkdir(parents=True)
    (d / "seed_0_curve.json").write_text("{not json", encoding="utf-8")
    taken, _ = mod.classify(tmp_path)
    s = mod.summarize(taken)
    assert s["curves"]["seed_0_curve.json"] == {"unreadable": True}
    assert s["total_epochs_on_disk"] == 0


# ---------------------------------------------------------------------------------------------------
# the refusal
# ---------------------------------------------------------------------------------------------------
def test_an_empty_selection_REFUSES_rather_than_uploading_nothing(mod, tmp_path, capsys):
    """An empty dataset makes the next session recompute every seed while looking normal."""
    d = tmp_path / "finetune_checkpoints"
    d.mkdir(parents=True)
    (d / "smoke_0.pt").write_bytes(b"x")
    with pytest.raises(SystemExit) as e:
        mod.main(["--from", str(tmp_path), "--slug", "ft-ckpt-x", "--dry-run"])
    msg = str(e.value)
    assert "REFUSING" in msg
    assert "recompute every seed" in msg
    # The refusal must point at the SKIPPED report, which is the only thing that says what WAS there.
    assert "SKIPPED" in msg
    out = capsys.readouterr().out
    assert "smoke (excluded by design)" in out


def test_a_dry_run_uploads_nothing_and_still_reports(mod, ckpt_dir, capsys):
    rc = mod.main(["--from", str(ckpt_dir), "--slug", "ft-ckpt-x", "--dry-run"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "DRY RUN -- nothing uploaded" in out
    assert "ACCEPTED 3 file(s)" in out
    assert "RESUMABLE STATE: 10 epoch(s)" in out


def test_the_next_step_instructions_name_BOTH_required_flags(mod, ckpt_dir, capsys):
    """One flag tells the kernel what to expect; the other makes Kaggle attach it. Both or neither.

    A build with --ckpt-dataset and a push without --dataset produces the exact failure that was
    measured: the kernel expects a mount, nothing is attached, and it refuses.
    """
    src = UPLOADER.read_text(encoding="utf-8")
    assert "--ckpt-dataset" in src and "--dataset" in src
    assert "both flags are required" in src.lower() or "both are required" in src.lower()
    # And it must warn about pushing before the dataset is attachable, which was measured.
    assert "ready before pushing" in src
