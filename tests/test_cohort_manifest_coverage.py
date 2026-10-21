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
    GAP_MEASURED,
    NO_GAP,
    UNTESTABLE_NO_FILES,
    UNTESTABLE_VACUOUS,
    _accessions,
    classify,
    coverage,
)

REPO = Path(__file__).resolve().parent.parent


def test_the_audit_is_non_vacuous_there_really_is_a_blind_spot():
    """A clean bill here is a claim, not a default — but there are THREE ways to get a zero, and the
    original version of this test named only two of them.

    It asserted `n_cohorts_missed > 0` and, on failure, said the frozen glob had probably been
    broadened and the prospective lock therefore needed re-minting. On this laptop (2026-10-21) the
    real cause is the third one: every cohort file is named exactly `selected.tsv` because the
    external/AR-Bank arms live on the detached data volume, so the broad glob adds nothing and the
    zero is an IDENTITY. Acting on the old message would have meant re-minting a lock — retiring
    every prospective number — to fix a filesystem that was merely unplugged.

    So: skip on a PROVABLE absence (no variant-named file exists), and keep the hard failure for the
    case where variant files are present and still unmissed, which is the only state that really does
    implicate the glob or the scan.
    """
    rep = coverage()
    assert rep["n_cohorts_on_disk"] > 0, "scan found no cohorts at all — broken, not clean"
    if rep["n_variant_files_on_disk"] == 0:
        pytest.skip(
            f"UNTESTABLE here, not clean: {rep['verdict']}. No `selected_*.tsv` variant file is on "
            "disk, so the frozen and broad globs are the same query and MISSED 0 is true by "
            "construction. Attach the external-cohort volume to measure the blind spot. This is NOT "
            "evidence that cohort_manifest.py's glob was broadened — do not re-mint the lock for it.")
    assert rep["n_cohorts_missed"] > 0, (
        f"{rep['n_variant_files_on_disk']} variant-named cohort file(s) are on disk "
        f"({rep['variant_basenames_seen']}) yet NONE is reported missed. The corpus is present, so "
        "this implicates the audit or the frozen glob: if cohort_manifest.py was broadened, the "
        "prospective lock pins it and needs re-minting, and test_the_frozen_glob_is_still_narrow "
        "should have fired too.")
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


# --------------------------------------------------------------------------------------------- #
# A zero has three causes, and the artifact must say WHICH (added 2026-10-21)
# --------------------------------------------------------------------------------------------- #
# The audit shipped a report that read as both "coverage is complete" (missed 0 of 490) and "there is
# a measured gap" (hardcoded consequence/limit prose) on a host where the comparison could not run at
# all. `BROAD_PATTERN` is a strict superset of `FROZEN_PATTERN`, so with no variant-named file present
# the two globs are the same query and a zero is an identity. These pin the distinction.
def test_classify_names_the_vacuous_case_rather_than_calling_it_clean():
    """No variant-named file => the comparison could not run, whatever the missed count says."""
    assert classify(n_all_files=10, n_variant_files=0, n_missed_accs=0) == UNTESTABLE_VACUOUS


def test_classify_separates_an_EMPTY_scan_from_a_vacuous_one():
    assert classify(n_all_files=0, n_variant_files=0, n_missed_accs=0) == UNTESTABLE_NO_FILES


def test_classify_reports_a_REAL_gap_when_variants_exist_and_are_missed():
    assert classify(n_all_files=12, n_variant_files=2, n_missed_accs=461) == GAP_MEASURED


def test_classify_reports_a_GENUINE_clean_bill_only_when_variants_were_present():
    """The one state in which 'no gap' is a finding: the files that could be missed exist, and aren't."""
    assert classify(n_all_files=12, n_variant_files=2, n_missed_accs=0) == NO_GAP


def _cohort(root, dirname, basename, accs):
    d = root / dirname
    d.mkdir(parents=True, exist_ok=True)
    (d / basename).write_text("".join(f"{a}\tR\n" for a in accs), encoding="utf-8")


def test_a_VACUOUS_run_REFUSES_to_ship_the_measured_gap_prose(tmp_path):
    """The defect itself: gap prose attached to a run that measured no gap and could not."""
    _cohort(tmp_path, "cohort_a", "selected.tsv", ["SAMN001", "SAMN002"])
    rep = coverage(tmp_path)
    assert rep["verdict"] == UNTESTABLE_VACUOUS
    assert rep["comparison_was_testable"] is False
    assert rep["n_accessions_missed"] == 0
    assert "NOT MEASURED BY THIS RUN" in rep["consequence_if_unfixed"]
    assert "IDENTITY, NOT A RESULT" in rep["honest_limit"]
    # and it must not speak of accessions it never found
    assert "the missed accessions" not in rep["honest_limit"]


def test_a_TESTABLE_run_DOES_ship_the_measured_gap_prose(tmp_path):
    """Non-vacuity of the branch above: with a variant file present the real finding is reported."""
    _cohort(tmp_path, "cohort_a", "selected.tsv", ["SAMN001"])
    _cohort(tmp_path, "arbank_kleb", "selected_strict.tsv", ["SAMN900", "SAMN901"])
    rep = coverage(tmp_path)
    assert rep["verdict"] == GAP_MEASURED
    assert rep["comparison_was_testable"] is True
    assert rep["missed_cohorts"] == ["arbank_kleb"]
    assert rep["n_accessions_missed"] == 2
    assert rep["variant_basenames_seen"] == ["selected_strict.tsv"]
    assert "fails OPEN" in rep["consequence_if_unfixed"]
    assert "DISCOVERY coverage only" in rep["honest_limit"]


def test_a_variant_file_BESIDE_a_bare_one_is_still_counted_as_testable(tmp_path):
    """The subtle case: the directory IS visible to the frozen glob, so it is not a missed COHORT --
    but its extra accessions are invisible, and the run is testable either way."""
    d = tmp_path / "both_shapes"
    d.mkdir()
    (d / "selected.tsv").write_text("SAMN001\tR\n", encoding="utf-8")
    (d / "selected_relaxed.tsv").write_text("SAMN001\tR\nSAMN777\tS\n", encoding="utf-8")
    rep = coverage(tmp_path)
    assert rep["comparison_was_testable"] is True
    assert rep["missed_cohorts"] == [], "the dir is reachable by the frozen glob"
    assert rep["n_accessions_missed"] == 1, "but SAMN777 is not, and that is the fail-open direction"
    assert rep["verdict"] == GAP_MEASURED


def test_the_LIVE_artifact_carries_a_verdict_consistent_with_the_disk():
    """Cross-check the shipped fields against the filesystem rather than against each other."""
    raw = REPO / "data" / "raw"
    variants = [p for p in glob.glob(f"{raw}/*/{BROAD_PATTERN}") if Path(p).name != FROZEN_PATTERN]
    rep = coverage()
    assert rep["n_variant_files_on_disk"] == len(variants)
    assert rep["comparison_was_testable"] is bool(variants)
    if not variants:
        assert rep["verdict"] == UNTESTABLE_VACUOUS
        assert rep["n_accessions_missed"] == 0, "must be 0 by construction when the globs coincide"
