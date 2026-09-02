"""Pins for the PEAR forward-replication run: coordinates, tie handling, and the refuse-guard.

The correlation itself is a data result and lives in the artifact. What is pinned here is the machinery
that makes the result trustworthy -- each of these silently produced a WRONG-BUT-PLAUSIBLE number first.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "pear_forward_replication.py"
ARTIFACT = ROOT / "wiki" / "pear_forward_replication_esm2_2026-09-01.json"
REF = Path("D:/dna_decode_cache/pear/CTXM-14/Genotype_barcode_calling/Ref_CTXM.fasta")
TABLE = Path("D:/dna_decode_cache/pear/extracted/Figure3.A__data.tsv")


def _mod():
    spec = importlib.util.spec_from_file_location("pear_forward_replication", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --- tie-safe rank correlation (the documented tie trap) ------------------------------------------

def test_spearman_is_exact_on_a_perfect_monotone_relation():
    m = _mod()
    assert m.spearman([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]) == pytest.approx(1.0)
    assert m.spearman([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]) == pytest.approx(-1.0)


def test_spearman_uses_mid_ranks_for_ties():
    """sorted()-order tie-breaking silently shifts correlations; mid-ranks are the fix."""
    m = _mod()
    # y is constant within the tied x block -> mid-ranks make this exactly 0, order-independent
    assert m.spearman([1.0, 1.0, 2.0, 2.0], [1.0, 2.0, 1.0, 2.0]) == pytest.approx(0.0)


def test_spearman_returns_none_on_a_degenerate_input():
    m = _mod()
    assert m.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None   # zero variance in x
    assert m.spearman([1.0], [2.0]) is None                        # too few points


# --- the artifact's own claims stay self-consistent ------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.is_file(), reason="run artifact not present")
def test_the_committed_run_passed_its_own_coordinate_sanity_check():
    """Nonsense must be less fit than silent on BOTH drugs. If this is not True, the mapping is wrong
    and the correlation in the same artifact means nothing."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    chk = d["sanity_check"]
    assert chk["passed"] is True
    assert chk["median_nonsense_ctx"] < chk["median_silent_ctx"]
    assert chk["median_nonsense_caz"] < chk["median_silent_caz"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="run artifact not present")
def test_the_committed_run_actually_scored_missense_variants():
    """The failure this guards: an uncoerced ESM table makes every missense variant a SKIP, leaving a
    run that reports 'missense n=0' and a silent-only table instead of failing."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["by_consequence"].get("missense", 0) > 1000
    assert d["n_skipped"] == 0
    assert d["correlation"]["ctx"]["spearman_score_vs_fitness"] is not None


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="run artifact not present")
def test_the_artifact_states_its_sign_convention_and_benchmark():
    """A bare correlation is unreadable without the sign convention, and the TEM-1 number must be
    labelled as a different protein and drug rather than a target."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert "MORE NEGATIVE = more damaging" in d["sign_convention"]
    assert d["benchmark"]["tem1_ampicillin_genome_edit_spearman"] == 0.7611
    assert "DIFFERENT" in d["benchmark"]["note"]


# --- the coordinate convention, re-derived from the real files when they are present ---------------

@pytest.mark.skipif(not (REF.is_file() and TABLE.is_file()), reason="PEAR extract not on this host")
def test_pear_coordinates_are_reference_based_not_offset():
    """2114/2114 on ref coords vs 452/2114 on the +81 Figure-2 axis convention. Re-derived, not recalled."""
    import csv
    import re
    m = _mod()
    seq = m.load_ref(REF)
    pat = re.compile(r"^([ACGT])(\d+)([ACGT])$")
    ok_ref = ok_off = total = 0
    for r in csv.DictReader(TABLE.open(encoding="utf-8"), delimiter="\t"):
        g = pat.match(r["genotype"])
        if not g:
            continue
        wt, p = g.group(1), int(g.group(2))
        total += 1
        ok_ref += seq[p - 1] == wt if p <= len(seq) else 0
        ok_off += seq[p - 82] == wt if 0 < p - 81 <= len(seq) else 0
    assert total == 2114
    assert ok_ref == total, "reference-coordinate reading must match every variant"
    assert ok_off < total // 2, "the offset reading must NOT match (it is the wrong convention)"
