"""The probe must not mistake a corrupt input, or a useless flip, for a finding.

Two failure modes are pinned because both actually occurred while writing it: a head-only gzip check
that passes truncated files, and a verdict that counted flips without asking whether any flip helped.
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "pneumo_selection_rule_probe.py"
ARTIFACT = ROOT / "wiki" / "pneumo_selection_rule_probe_2026-09-04.json"


def test_a_truncated_gzip_is_rejected_not_silently_used(tmp_path):
    """A truncated stream decompresses its FIRST blocks fine and raises only at the end, so a
    head-only check passes exactly the files that later fail. Full decompression is the guard."""
    good = gzip.compress(b">contig1\n" + b"ACGT" * 20000)
    truncated = good[: len(good) // 2]
    p = tmp_path / "t.fa.gz"
    p.write_bytes(truncated)
    with pytest.raises(Exception):
        gzip.decompress(p.read_bytes())
    # and a head-only read does NOT raise -- which is precisely why it is the wrong guard
    with gzip.open(p, "rb") as fh:
        head = fh.read(64)
    assert head.startswith(b">"), "head-only read succeeds on a truncated file (the trap)"


def test_probe_refuses_when_no_assembly_produced_a_call(tmp_path):
    """Zero calls means the rule was never exercised; a zero-flip verdict there would be hollow."""
    asm = tmp_path / "asm"
    asm.mkdir()
    (asm / "ERS999999.fa.gz").write_bytes(b"<!DOCTYPE HTML>\n403 Forbidden\n")  # the real cache shape
    out = tmp_path / "o.json"
    r = subprocess.run([sys.executable, str(SCRIPT), "--asm-dir", str(asm),
                        "--cohort", str(tmp_path / "missing.tsv"), "--out", str(out)],
                       capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 3, r.stdout + r.stderr
    assert "REFUSING" in r.stderr
    assert not out.exists(), "a refused run must not leave an artifact"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="probe artifact absent")
def test_the_committed_verdict_scores_flips_against_the_label_not_just_counts_them():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    fo = d["flip_outcomes_vs_measured_serogroup"]
    assert set(fo) == {"improved", "worsened", "neither_right"}
    # The finding: the rule flips a call but never turns a miss into a hit.
    assert d["n_flipped"] >= 1, "a zero-flip artifact would make the verdict a different claim"
    assert fo["improved"] == 0
    assert d["verdict"] == "RULE_FLIPS_BUT_NEVER_IMPROVES_ON_SAMPLE"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="probe artifact absent")
def test_the_artifact_records_the_cache_corruption():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ci = d["cache_integrity"]
    assert ci["n_corrupt"] > 0, "the corrupt cached cohort must stay recorded, not be quietly dropped"
    assert any("403" in str(c) or "EOFError" in str(c) or "too small" in str(c)
               for c in ci["corrupt"]), "the corruption REASONS must be auditable"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="probe artifact absent")
def test_small_sample_limit_is_stated():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_cached_scored"] < 260
    assert any("SMALL SAMPLE" in s for s in d["honest_limits"])
