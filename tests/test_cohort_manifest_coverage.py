"""The leakage registry's DISCOVERY blind spot — measured, pinned, and deliberately NOT fixed.

`dna_decode/eval/cohort_manifest.py` globs **exactly** `data/raw/*/selected.tsv`. Every AR-Bank and
external-cohort arm writes `selected_strict.tsv` / `selected_relaxed.tsv`, so none of them is visible to
a guard whose whole job is to answer "was this accession used before?". The guard fails OPEN in the one
direction it exists to prevent.

It is NOT fixed in-place because `cohort_manifest.py` is sha256-pinned by the ACTIVE prospective lock.
Broadening the glob would retire every prospective number scored against that lock. See
`scripts/cohort_manifest_coverage_audit.py` for the full reasoning.

These tests pin STRUCTURE, not counts. An exact accession count against a live `data/raw/` tree fires on
ordinary cohort accrual rather than on the regression — the trap this repo has on record.
"""
from __future__ import annotations

import glob
from pathlib import Path

import pytest

from scripts.cohort_manifest_coverage_audit import (
    BROAD_PATTERN,
    FROZEN_PATTERN,
    _accessions,
    coverage,
)

REPO = Path(__file__).resolve().parent.parent


def test_the_audit_is_non_vacuous_there_really_is_a_blind_spot():
    """If this ever reports zero missed cohorts, either the frozen glob was broadened (see the
    self-retiring test below) or the scan is broken. A clean bill here is a claim, not a default."""
    rep = coverage()
    assert rep["n_cohorts_on_disk"] > 0, "scan found no cohorts at all — broken, not clean"
    assert rep["n_cohorts_missed"] > 0, (
        "no cohort is missed: if cohort_manifest.py was broadened, the prospective lock needs "
        "re-minting and test_the_frozen_glob_is_still_narrow should have caught it")
    assert rep["n_accessions_missed"] > 0


def test_every_suffixed_cohort_dir_is_reported_missed():
    """The structural invariant: a directory carrying `selected_*.tsv` but no bare `selected.tsv`
    CANNOT be seen by the frozen glob, so it must appear in `missed_cohorts`. Holds at any corpus
    size, which an accession count does not."""
    raw = REPO / "data" / "raw"
    with_broad = {Path(p).parent.name for p in glob.glob(f"{raw}/*/{BROAD_PATTERN}")}
    with_exact = {Path(p).parent.name for p in glob.glob(f"{raw}/*/{FROZEN_PATTERN}")}
    expected_missed = with_broad - with_exact
    reported = set(coverage()["missed_cohorts"])
    assert reported == expected_missed, (
        f"audit disagrees with the filesystem: only-in-report={reported - expected_missed}, "
        f"only-on-disk={expected_missed - reported}")


def test_the_audit_agrees_with_the_LIVE_frozen_manifest():
    """Cross-check against the real `build_manifest()` rather than trusting our own glob. If these
    diverge, the audit is measuring something other than what the guard actually does."""
    from dna_decode.eval.cohort_manifest import build_manifest
    m = build_manifest()
    live_tsv_cohorts = {c.name for c in m.cohorts if c.source == "selected_tsv"}
    rep = coverage()
    assert rep["n_cohorts_seen"] == len(live_tsv_cohorts), (
        f"audit says the frozen glob sees {rep['n_cohorts_seen']} selected-tsv cohorts; "
        f"build_manifest() actually reports {len(live_tsv_cohorts)}")
    # and no cohort we call MISSED may actually be present in the live manifest
    assert not (set(rep["missed_cohorts"]) & live_tsv_cohorts), (
        "a cohort reported as invisible is in fact registered — the audit is over-claiming")


def test_the_frozen_glob_is_still_narrow_SELF_RETIRING():
    """Self-retiring premise check. If someone broadens the glob, this fails ON PURPOSE and says what
    else must happen: the prospective lock pins this file's sha256, so the fix REQUIRES re-minting the
    lock and restarting the accrual clock. Failing loudly here is cheaper than discovering it when
    `verify_lock` hard-fails mid-scoring."""
    src = (REPO / "dna_decode" / "eval" / "cohort_manifest.py").read_text(encoding="utf-8")
    assert src.strip(), "the frozen manifest module is missing or empty — assert existence first"
    assert '/selected.tsv"' in src or "/selected.tsv'" in src, (
        "cohort_manifest.py no longer globs the narrow `selected.tsv`. If that was intentional, the "
        "ACTIVE prospective lock pins this file and must be re-minted (a NEW lock date), which retires "
        "every prospective number scored against the old one. Update this test only alongside that.")
    assert "selected*.tsv" not in src, (
        "the broad pattern appears in the frozen module — see the message above")


def test_accession_reader_is_not_vacuous_on_a_real_cohort_file():
    """The whole measurement rests on reading column 0. Prove it returns real accessions on a real
    file rather than silently returning an empty set (the shape that made an earlier overlap check
    report a flattering zero)."""
    found = sorted(glob.glob(f"{REPO / 'data' / 'raw'}/*/{BROAD_PATTERN}"))
    if not found:
        pytest.skip("no cohort TSVs on disk")
    non_empty = [p for p in found if _accessions(Path(p))]
    assert non_empty, "every cohort TSV parsed to an EMPTY accession set — the reader is broken"
    sample = _accessions(Path(non_empty[0]))
    assert any(a[:4] in ("SAMN", "SAME", "SAMD", "GCA_", "GCF_") for a in sample), (
        f"parsed tokens do not look like accessions: {sorted(sample)[:4]}")
