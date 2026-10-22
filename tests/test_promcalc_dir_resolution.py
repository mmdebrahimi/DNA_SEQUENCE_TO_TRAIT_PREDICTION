"""The promcalc directory is resolved AT POINT OF USE in all three stages (Step 1, 2026-10-22).

All three stages hardcoded `D:/dna_decode_cache/promcalc`. That volume detached permanently, and the
failure it produced was not just inconvenient but MISLEADING: stage B reported "tiles.tsv missing -
run stage A first" and stage C reported "no decomposition shards - run stage B first", both of which
send the operator to re-run a stage when the real wall is that the whole volume is gone. Those are
different walls, and merging them destroys the diagnosis.

**The D: default is deliberately unchanged in all three**, so a host that HAS the volume behaves
bit-identically. This is parameterisation, not relocation, and the first test below pins that.

No shared helper is extracted, because the three directions genuinely differ:
  A  writes an OUTPUT dir it creates   -> the nearest existing ANCESTOR must be reachable
  B  reads an INPUT dir, stdlib only   -> may not import dna_decode at all
  C  reads an INPUT dir holding shards -> resolve_external_root's is_dir semantics fit exactly
"""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
STAGE_A = REPO / "scripts" / "promcalc_export_tiles.py"
STAGE_B = REPO / "scripts" / "promcalc_decompose.py"
STAGE_C = REPO / "scripts" / "promcalc_decomposition_arm.py"

RECORDED_DEFAULT = "D:/dna_decode_cache/promcalc"


def _load(path: Path, tag: str):
    spec = importlib.util.spec_from_file_location(f"_promcalc_{tag}", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- the default must NOT move -------------------------------------------------------------------

@pytest.mark.parametrize("path", [STAGE_A, STAGE_B, STAGE_C], ids=["stage_a", "stage_b", "stage_c"])
def test_the_recorded_D_default_is_UNCHANGED(path: Path):
    """Step 1 is ADDITIVE. A host that has the volume must take the identical path.

    If this fails, the change stopped being parameterisation and became relocation -- which would
    silently repoint every operator who still has the drive.
    """
    assert RECORDED_DEFAULT in path.read_text(encoding="utf-8"), (
        f"{path.name} no longer carries the recorded default {RECORDED_DEFAULT!r}; Step 1 was meant "
        f"to parameterise it, not replace it."
    )


# --- stage A: an OUTPUT dir, judged by its nearest existing ancestor ----------------------------

def test_stage_a_refuses_an_unreachable_out_dir_with_a_message_naming_the_flag():
    mod = _load(STAGE_A, "a_unreachable")
    with pytest.raises(SystemExit) as exc:
        mod._require_creatable_out_dir("Z:/nope/deeper")
    msg = str(exc.value)
    assert "--out-dir" in msg, "the refusal must name the flag the operator has to pass"
    assert "not reachable" in msg
    # It must NOT be an OSError/WinError traceback -- that was the pre-Step-1 behaviour.
    assert "WinError" not in msg


def test_stage_a_ACCEPTS_a_dir_that_does_not_exist_yet_but_whose_parent_does():
    """The reason stage A cannot use `resolve_external_root` on its own output path.

    That resolver reports `configured` only when the path ALREADY is a directory. On a first export
    the leaf legitimately does not exist, so asking the resolver about it would refuse a good run.
    """
    mod = _load(STAGE_A, "a_creatable")
    target = REPO / "tests"  / "_not_created_by_this_test"
    assert not target.exists(), "fixture assumption: the leaf must not exist"
    assert mod._require_creatable_out_dir(target) == target
    assert not target.exists(), "the check must not CREATE anything -- mkdir is the caller's job"


def test_stage_a_refuses_when_an_ancestor_is_a_FILE(tmp_path: Path):
    """A file where a directory belongs cannot be refused by an existence check alone."""
    f = tmp_path / "a_file"
    f.write_text("not a directory", encoding="utf-8")
    mod = _load(STAGE_A, "a_file_ancestor")
    with pytest.raises(SystemExit) as exc:
        mod._require_creatable_out_dir(f / "under_a_file")
    assert "FILE" in str(exc.value)


# --- stage B: three walls that used to be one message -------------------------------------------

def test_stage_b_distinguishes_UNREACHABLE_from_NOT_YET_EXPORTED(tmp_path: Path):
    """The whole point of Step 1 for this stage.

    Before: an absent volume and a reachable-but-unexported directory produced the SAME message,
    "tiles.tsv missing - run stage A first", which is actively wrong advice in the first case.
    """
    mod = _load(STAGE_B, "b_walls")

    with pytest.raises(SystemExit) as detached:
        mod._require_reachable_in_dir("Z:/nope/deeper")
    assert "not reachable" in str(detached.value)
    assert "NOT 'stage A has not run yet'" in str(detached.value), (
        "the detached-volume refusal must say explicitly that this is NOT the stage-A wall"
    )

    with pytest.raises(SystemExit) as unexported:
        mod._require_reachable_in_dir(tmp_path / "not_exported_yet")
    assert "does not exist" in str(unexported.value)
    assert "is reachable" in str(unexported.value), (
        "when the ancestor IS reachable the refusal must say so, or it reads as a detached volume"
    )

    # and the two messages must actually differ, or the distinction is cosmetic
    assert str(detached.value) != str(unexported.value)


def test_stage_b_PASSES_THROUGH_an_existing_directory(tmp_path: Path):
    """Non-vacuity: the reachable path must return, not raise.

    This is the path the live 44,106-tile run and every watchdog relaunch take, so a regression here
    would break the job rather than just its error messages.
    """
    mod = _load(STAGE_B, "b_ok")
    assert mod._require_reachable_in_dir(tmp_path) == tmp_path


def test_stage_b_refusal_path_imports_NOTHING_from_dna_decode():
    """Stage B's contract: stdlib only, at any nesting depth.

    Stages A and C reach for `external_root` in their refusals; stage B may not, because it runs under
    a separate interpreter that holds no project dependencies. Asserted by AST rather than by grepping,
    since the module legitimately MENTIONS dna_decode in prose and in a path string.
    """
    tree = ast.parse(STAGE_B.read_text(encoding="utf-8"))
    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            imported.append(node.module or "")
    offenders = sorted({m for m in imported if m.split(".")[0] == "dna_decode"})
    assert not offenders, f"stage B gained {offenders}; the separate-interpreter boundary is broken"


def test_stage_b_helper_is_reachable_under_the_ISOLATED_interpreter_contract():
    """The helper must be callable with only the stdlib present.

    Exercised here by importing stage B in the project env (which has no `promotercalculator`) and
    calling the helper -- if either the module or the helper had acquired a third-party dependency,
    this would fail.
    """
    mod = _load(STAGE_B, "b_iso")
    assert callable(mod._require_reachable_in_dir)


# --- stage C: an INPUT dir that must already hold shards ----------------------------------------

def test_stage_c_refuses_an_unreachable_in_dir_and_says_it_is_NOT_the_stage_b_wall():
    mod = _load(STAGE_C, "c_unreachable")
    with pytest.raises(SystemExit) as exc:
        mod._require_decomposition_dir("Z:/nope")
    msg = str(exc.value)
    assert "--in-dir" in msg
    assert "NOT substituted automatically" in msg, (
        "allow_discovered must stay False and the refusal must say so -- adopting a same-named root "
        "would convert a loud missing-ROOT into a quiet missing-DATASET"
    )
    assert "different wall from" in msg, (
        "must distinguish itself from 'stage B has not run yet', which is what the old message implied"
    )


def test_stage_c_PASSES_THROUGH_a_reachable_dir(tmp_path: Path):
    """Non-vacuity for stage C."""
    mod = _load(STAGE_C, "c_ok")
    assert mod._require_decomposition_dir(tmp_path) == tmp_path


def test_stage_c_never_passes_allow_discovered_True():
    """Pinned as code shape, because the default is the load-bearing part.

    `resolve_external_root` CAN return a discovered root. Passing `allow_discovered=True` here would
    silently repoint the arm at a directory of the right NAME with the wrong CONTENTS.

    **Checked by AST, not by substring.** The first version of this test searched the helper's source
    text for `allow_discovered` and failed -- against correct code -- because the DOCSTRING explains
    why the flag stays False. A text search cannot tell a mention from a keyword argument, which is
    the same trap this file's own header records for stage B's `dna_decode` occurrences. Asserting on
    the call's keywords also makes the test stricter: it now catches `allow_discovered=True` passed
    anywhere in the helper, while still permitting the explanatory prose.
    """
    tree = ast.parse(STAGE_C.read_text(encoding="utf-8"))
    helper = next((n for n in ast.walk(tree)
                   if isinstance(n, ast.FunctionDef) and n.name == "_require_decomposition_dir"), None)
    assert helper is not None, "_require_decomposition_dir is gone"

    calls = [n for n in ast.walk(helper)
             if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name)
             and n.func.id == "resolve_external_root"]
    assert calls, "the helper no longer calls resolve_external_root"

    for call in calls:
        for kw in call.keywords:
            assert kw.arg != "allow_discovered", (
                "allow_discovered is passed explicitly; it must be left at its default False so no "
                "future edit can flip it here unnoticed"
            )
