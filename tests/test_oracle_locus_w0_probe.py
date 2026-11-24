"""Tests for the W0 probe's PURE logic. No files, no network, no Excel.

The probe's whole job is to refuse rather than return a number that reads as a clean result, so
these tests are mostly about the refusals firing and about the start-window arithmetic being
right -- that arithmetic is the one place a silent error would change every null in the artifact.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    p = ROOT / "scripts" / "oracle_locus_w0_probe.py"
    spec = importlib.util.spec_from_file_location("oracle_locus_w0_probe", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["oracle_locus_w0_probe"] = mod
    spec.loader.exec_module(mod)
    return mod


W0 = _load()


# ---- chromosome naming: the three cells disagree, so one normaliser both ways ----------------

@pytest.mark.parametrize("raw,want", [("1", "chr1"), ("chr1", "chr1"), ("X", "chrX"),
                                      ("chrX", "chrX"), (" 10 ", "chr10")])
def test_to_chr_is_idempotent_and_handles_bare(raw, want):
    assert W0.to_chr(raw) == want
    assert W0.to_chr(W0.to_chr(raw)) == want


@pytest.mark.parametrize("raw,want", [("chr1", "1"), ("1", "1"), ("chrX", "X")])
def test_strip_chr_is_the_inverse(raw, want):
    assert W0.strip_chr(raw) == want


def test_wtc11_bare_and_hepg2_prefixed_join_after_normalising():
    """WTC11 writes `1`; HepG2 writes `chr1`. Unnormalised they never join."""
    assert "1" != "chr1"
    assert W0.to_chr("1") == W0.to_chr("chr1")


# ---- the start-window union: the arithmetic every null depends on ----------------------------

def test_start_window_union_on_a_hand_computed_overlapping_pair():
    """Two runs 100 bp apart -- CLOSER than the 200 bp element, so a sum double-counts.

    Each run's start-window at k=1 is [a-200+1, b-1+1): (1000,1200) -> [801,1200) = 399 starts,
    (1300,1500) -> [1101,1500) = 399. They overlap from 1101, so the sum is 798 and the true
    union is [801,1500) = 699.
    """
    by = {"chr1": [(1000, 1200), (1300, 1500)]}
    assert W0.naive_window_sum(by, 1) == 798
    assert W0.start_window_union(by, 1) == 699


def test_far_apart_runs_are_not_double_counted_away():
    """The correction must only remove REAL overlap: disjoint windows give union == sum."""
    far = {"chr1": [(0, 200), (10_000, 10_200)]}
    assert W0.start_window_union(far, 1) == W0.naive_window_sum(far, 1)


def test_identical_interval_tier_admits_exactly_one_start_per_run():
    """At k=200 only a start exactly at `a` fully covers a 200 bp run."""
    by = {"chr1": [(1000, 1200), (1300, 1500)]}
    assert W0.start_window_union(by, 200) == 2


def test_union_never_exceeds_the_naive_sum_at_any_tier():
    by = {"chr1": [(0, 200), (150, 400), (1000, 1200)], "chr2": [(5, 205)]}
    for _name, k in W0.TIERS:
        assert W0.start_window_union(by, k) <= W0.naive_window_sum(by, k)


def test_chromosomes_are_unioned_independently():
    """A run on chr1 must not cancel opportunity on chr2."""
    one = {"chr1": [(1000, 1200)]}
    two = {"chr1": [(1000, 1200)], "chr2": [(1000, 1200)]}
    assert W0.start_window_union(two, 1) == 2 * W0.start_window_union(one, 1)


# ---- merge integrity: an unmerged footprint invalidates every null ---------------------------

def test_merge_integrity_detects_each_defect_separately():
    assert W0.merge_integrity({"c": [(0, 100), (100, 200)]})["abutting"] == 1
    assert W0.merge_integrity({"c": [(0, 100), (50, 200)]})["overlapping"] == 1
    assert W0.merge_integrity({"c": [(100, 200), (0, 50)]})["unsorted"] == 1


def test_merge_integrity_passes_a_genuine_merge():
    assert not any(W0.merge_integrity({"c": [(0, 100), (200, 300)]}).values())


# ---- the refusals -- a probe that returns a number on broken input is the whole hazard -------

def test_find_header_locates_a_moving_header_row():
    """Measured: the header sits at row 4 in HepG2, 2 in K562, 1 in WTC11."""
    rows = [("Elements were designed",), (None,), (None,), ("name", "category", "chr.hg38")]
    assert W0.find_header(rows) == 3
    assert W0.find_header([("name", "category")]) == 0


def test_find_header_refuses_rather_than_guessing():
    with pytest.raises(W0.ProbeRefusal):
        W0.find_header([("Elements were designed",), (None,), ("something else",)])


def test_classify_columns_refuses_without_a_complete_hg38_trio():
    """An absent hg38 coordinate is reported missing, NEVER filled from hg19."""
    with pytest.raises(W0.ProbeRefusal):
        W0.classify_columns(("name", "category", "chr.hg19", "start.hg19", "stop.hg19"))


def test_classify_columns_finds_both_assemblies_when_present():
    hdr = ("name", "category", "chr.hg19", "start.hg19", "stop.hg19", "str.hg19",
           "chr.hg38", "start.hg38", "stop.hg38")
    got = W0.classify_columns(hdr)
    assert sorted(got["assemblies"]) == ["hg19", "hg38"]
    assert got["assemblies"]["hg38"] == {"chr": 6, "start": 7, "stop": 8}


def test_the_frozen_lists_are_non_empty():
    """A guard looping over an empty list passes vacuously forever."""
    assert W0.TIERS and W0.OUTCOMES
    assert "REFUSED_POSITIVE_CONTROL_FAILED" in W0.OUTCOMES


def test_tiers_are_ordered_loosest_to_strictest_and_cover_the_element():
    ks = [k for _n, k in W0.TIERS]
    assert ks == sorted(ks)
    assert ks[-1] == W0.ELEMENT_BP


def test_split_chromosome_sets_are_disjoint():
    """The plumbing gate is only meaningful if the published split rule partitions."""
    val, test = W0.SPLIT_CHROMS["val"], W0.SPLIT_CHROMS["test"]
    assert not (val & test)


def test_self_check_passes():
    assert W0._self_check() == 0
