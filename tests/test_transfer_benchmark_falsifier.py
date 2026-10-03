"""Guards for the three-arm F1 falsifier (Step 6).

Arms B and C are fully offline; A1 needs `D:` and is skipped with an explicit reason when absent --
and its skip must FORCE indeterminacy rather than letting two arms look like a pass.
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import scripts.transfer_benchmark_falsifier as F  # noqa: E402

SRC = ROOT / "scripts" / "transfer_benchmark_falsifier.py"
FROZEN = (
    "dna_decode/eval/amr_rules.py",
    "dna_decode/data/calibrated_amr_rules.json",
    "dna_decode/data/mic_tiers.py",
    "dna_decode/data/shipped_decoder_surface.py",
    "dna_decode/eval/cohort_manifest.py",
)


# --- the bar is frozen above the data ---------------------------------------------------------------

def test_the_preregistered_bar_appears_BEFORE_the_first_substrate_path_in_source_order():
    """Pre-committing the verdict-conditional response before results land is the discipline; a bar
    defined after the data is read is not pre-registered in any meaningful sense."""
    tree = ast.parse(SRC.read_text(encoding="utf-8"))
    line_of = {}
    for node in tree.body:
        for tgt in getattr(node, "targets", []):
            if isinstance(tgt, ast.Name):
                line_of[tgt.id] = node.lineno
    assert "PREREGISTERED" in line_of, "the bar is not a module-level constant"
    for path_const in ("BLOOM_GENO", "BLOOM_PHENO", "G2_FIXTURE"):
        assert line_of["PREREGISTERED"] < line_of[path_const], \
            f"{path_const} is defined before the bar"


def test_the_bar_names_only_A1_and_B_as_pass_requirements():
    assert F.PREREGISTERED["pass_requires"] == ["A1_positive", "B_negative"]
    assert F.PREREGISTERED["C_counts_toward_pass"] is False


# --- the dropped cross-panel control cannot be reinstated silently --------------------------------

def test_there_is_NO_cross_panel_flag():
    """The v2 draft offered `--incompatible-feature-space-control`. It was dropped because its result
    is knowable without running it; the CLI must reject it so it cannot creep back."""
    with pytest.raises(SystemExit):
        F.main(["--incompatible-feature-space-control"])
    src = SRC.read_text(encoding="utf-8")
    assert "add_argument(\"--incompatible-feature-space-control\"" not in src


# --- the offline arms ------------------------------------------------------------------------------

def test_arm_B_is_negative_with_the_machinery_actually_exercised():
    b = F.arm_b_synthetic_negative()
    assert b["negative"] is True
    assert b["machinery_exercised"] is True
    assert b["verdict"] != "BEATS_ALL_CONTROLS"
    assert b["n_folds"] == 4 and b["n_leakage_audits_clean"] == 4
    assert b["report"]["n_controls_run"] == 5


def test_arm_C_is_non_positive_on_every_seed_DESPITE_winning_global_r2():
    """THE HARD CASE. The candidate beats the baseline on one metric; the verdict must still refuse,
    driven by the de-confounded cell rather than by whichever metric the candidate happens to win."""
    c = F.arm_c_verdict_replay()
    assert c["all_seeds_non_positive"] is True
    assert c["won_global_r2_on_every_seed"] is True, \
        "if the embedding no longer wins global r2, Arm C is not testing what it claims"
    assert c["counts_toward_pass"] is False
    assert c["machinery_exercised"] is False
    assert len(c["per_seed"]) == 3
    for s in c["per_seed"]:
        assert s["within_group_r2"] < 0


def test_arm_C_alone_can_never_produce_a_PASS():
    """The loophole: replaying a historical conclusion while the live path is only synthetic-tested."""
    res = F.run(offline_only=True)
    assert res["arms"]["C"]["all_seeds_non_positive"] is True
    assert res["status"] == F.INDETERMINATE and res["exit_code"] == 3


def test_skipping_A1_FORCES_indeterminate_rather_than_a_two_arm_pass():
    res = F.run(offline_only=True)
    assert "A1" not in res["arms"]
    assert res["status"] != F.DISCRIMINATED
    assert any("A1" in r for r in res["refusals"])


def test_a_missing_arm_C_fixture_does_NOT_block_a_pass(monkeypatch):
    """`cannot run` and `ran and broke` are different facts. A missing regression fixture is recorded
    but must not be load-bearing for the pass; only a C that RAN and came out positive blocks."""
    monkeypatch.setattr(F, "G2_FIXTURE", ROOT / "wiki" / "does_not_exist_g2.json")
    monkeypatch.setattr(F, "arm_a1_scoring_core",
                        lambda *a, **k: {"arm": "A1_scoring_core_control", "positive": True,
                                         "is_transfer_evidence": False})
    res = F.run()
    assert "C" not in res["arms"]
    assert any(r.startswith("C:") for r in res["refusals"])
    assert res["status"] == F.DISCRIMINATED and res["exit_code"] == 0


def test_a_REGRESSED_arm_C_blocks_the_pass(monkeypatch):
    """Non-vacuity for the asymmetry above: the other branch must actually bite."""
    monkeypatch.setattr(F, "arm_a1_scoring_core",
                        lambda *a, **k: {"arm": "A1_scoring_core_control", "positive": True,
                                         "is_transfer_evidence": False})
    monkeypatch.setattr(F, "arm_c_verdict_replay",
                        lambda: {"arm": "C_verdict_replay_regression",
                                 "all_seeds_non_positive": False, "counts_toward_pass": False})
    res = F.run()
    assert res["status"] == F.BAR_NOT_MET and res["exit_code"] == 1


# --- the unmeasured verdict is a VALUE ---------------------------------------------------------------

def test_transfer_positive_unmeasured_is_emitted_on_every_run_as_a_value():
    for kwargs in ({"offline_only": True},):
        res = F.run(**kwargs)
        assert res["transfer_positive"] == F.TRANSFER_POSITIVE_UNMEASURED
        assert "single constructed population" in res["transfer_positive_reason"]


# --- the genotype integrity gate ---------------------------------------------------------------------

def test_a_truncated_genotype_matrix_is_REFUSED_not_scored(tmp_path):
    """The plan named BYxRM_GenoData.txt for Arm A1; that file is a truncated download (exactly 1 MiB,
    511 of 11,623 rows, final row 626 fields of 1009). Building against it would have scored a ragged
    fifth of the panel. The check is RAGGEDNESS, not size -- a size threshold would be a guess."""
    good = tmp_path / "ok.txt"
    good.write_text("m\ta\tb\n1\tB\tR\n2\tR\tB\n", encoding="utf-8")
    assert F.assert_genotype_file_intact(good)["intact"] is True

    cut = tmp_path / "cut.txt"
    cut.write_text("m\ta\tb\n1\tB\tR\n2\tR\n", encoding="utf-8")       # last row short
    with pytest.raises(F.ArmRefused) as exc:
        F.assert_genotype_file_intact(cut)
    assert "TRUNCATED" in str(exc.value)

    with pytest.raises(F.ArmRefused):
        F.assert_genotype_file_intact(tmp_path / "absent.txt")


def test_the_real_truncated_bloom_file_is_refused_if_present():
    """Scoped so CI without D: still passes, but on this host it pins the actual corrupt artifact."""
    bad = Path("D:/dna_decode_cache/bloom/BYxRM_GenoData.txt")
    if not bad.exists():
        pytest.skip("D: not mounted")
    with pytest.raises(F.ArmRefused):
        F.assert_genotype_file_intact(bad)


def test_the_falsifier_points_at_the_INTACT_genotype_file():
    assert F.BLOOM_GENO.name == "geno_v2.txt", \
        "BYxRM_GenoData.txt is the truncated one; geno_v2.txt is the full 11,623-marker panel"


# --- artifact + frozen surface -----------------------------------------------------------------------

def test_the_emitted_json_is_reparseable_and_the_md_is_written(tmp_path):
    """A hand-edited or malformed artifact has gone unnoticed in this repo before."""
    rc = F.main(["--offline-only", "--out-dir", str(tmp_path)])
    assert rc == 3
    jsons = list(tmp_path.glob("transfer_benchmark_falsifier_*.json"))
    mds = list(tmp_path.glob("transfer_benchmark_falsifier_*.md"))
    assert len(jsons) == 1 and len(mds) == 1
    d = json.loads(jsons[0].read_text(encoding="utf-8"))
    assert d["schema"] == "transfer-benchmark-falsifier-v1"
    assert d["transfer_positive"] == F.TRANSFER_POSITIVE_UNMEASURED
    assert F.TRANSFER_POSITIVE_UNMEASURED in mds[0].read_text(encoding="utf-8")


def test_the_five_frozen_files_EXIST_before_anything_claims_they_are_unchanged():
    """`git diff --quiet` on a missing path exits 0 and tests nothing -- a guard on a nonexistent file
    passes vacuously, which this repo has already paid for."""
    missing = [p for p in FROZEN if not (ROOT / p).exists()]
    assert not missing, f"frozen paths do not resolve, so any diff check on them is vacuous: {missing}"


def test_the_falsifier_does_not_import_any_frozen_module():
    imported: set[str] = set()
    for node in ast.walk(ast.parse(SRC.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    for bad in ("amr_rules", "mic_tiers", "shipped_decoder_surface", "cohort_manifest"):
        assert not any(bad in m for m in imported), f"{bad} imported: {imported}"
    assert any("transfer_gauntlet" in m for m in imported), \
        f"the scanner is not reading imports: {imported}"
