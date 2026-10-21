"""Measure what the leakage manifest CANNOT SEE — read-only, NON-FROZEN, exit 0 always.

WHY THIS IS A SEPARATE SCRIPT AND NOT A ONE-LINE FIX
----------------------------------------------------
`dna_decode/eval/cohort_manifest.py` discovers cohorts with `glob("data/raw/*/selected.tsv")` —
**exactly** `selected.tsv`. Every AR-Bank and external-cohort arm writes `selected_strict.tsv` /
`selected_relaxed.tsv` instead, so all of them are INVISIBLE to the leakage registry. Measured
2026-10-03, **on a host where those directories were present**: 23 cohort directories and 461
accessions absent from a fail-closed guard whose entire job is to answer "was this isolate used
before?".

THAT FIGURE IS NOT REPRODUCIBLE ON EVERY HOST, AND A ZERO HERE IS NOT A CLEAN BILL
----------------------------------------------------------------------------------
`BROAD_PATTERN` is a strict superset of `FROZEN_PATTERN`, so the only files it can add are those whose
basename is NOT `selected.tsv`. When no such file is on disk the two globs are **the same query** and
`n_accessions_missed == 0` is true BY CONSTRUCTION — the audit compared a set to itself and measured
nothing. Measured 2026-10-21 on this laptop: all 10 cohort files are named exactly `selected.tsv` (the
external/AR-Bank arms live on the detached data volume), so the run is vacuous in precisely that way.

So the report carries a `verdict`, and the gap prose is CONDITIONAL on it. Without that, one artifact
reads as both "coverage is complete" (the counts) and "there is a measured gap" (the prose), and
neither is a claim such a run can support. This is the repo's own recorded trap — a restriction that
removes nothing is not a control — applied to this audit itself.

The obvious fix (broaden the glob to `selected*.tsv`, union per directory) is **three lines** and would
make the manifest strictly MORE fail-closed. It is deliberately NOT applied here, because
`cohort_manifest.py` is one of the FIVE sha256-pinned files in the active prospective lock
(`wiki/prospective_lock_manifest_2026-08-31.json`). Editing it would:

  * invalidate the v2 lock,
  * RETIRE every prospective number scored against it — including the 2026-09-28 Campylobacter cipro
    cell (acc 0.998) and the 2026-10-03 Klebsiella meropenem cell,
  * and force a third restart of the accrual clock.

That is an AUTHORITY decision about the frozen surface, not an executor one. So this script delivers
the VISIBILITY (which is the actual value) with zero lock cost, and the frozen fix stays a user call.

    uv run python scripts/cohort_manifest_coverage_audit.py
"""
from __future__ import annotations

import glob
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# The pattern the FROZEN manifest actually uses, and the one that would cover the real corpus.
FROZEN_PATTERN = "selected.tsv"
BROAD_PATTERN = "selected*.tsv"

# Verdicts. The first two are NOT findings about the glob — they are findings about the corpus
# this run could see, and they must never be read as a clean bill.
UNTESTABLE_NO_FILES = "UNTESTABLE_NO_COHORT_FILES_ON_DISK"
UNTESTABLE_VACUOUS = "UNTESTABLE_COMPARISON_VACUOUS_NO_VARIANT_FILES_PRESENT"
GAP_MEASURED = "GAP_MEASURED"
NO_GAP = "NO_GAP_VARIANT_FILES_PRESENT_AND_COVERED"


def classify(n_all_files: int, n_variant_files: int, n_missed_accs: int) -> str:
    """Why a zero is a zero. Pure, so the three states are testable without a filesystem.

    `n_variant_files` is the count of broad-glob hits whose basename is NOT `FROZEN_PATTERN` — i.e.
    the only files that could ever be missed. At zero, the frozen and broad globs resolve to the
    same set and a zero-missed result is an identity, not a measurement.
    """
    if n_all_files == 0:
        return UNTESTABLE_NO_FILES
    if n_variant_files == 0:
        return UNTESTABLE_VACUOUS
    return GAP_MEASURED if n_missed_accs > 0 else NO_GAP


def _accessions(path: Path) -> set[str]:
    """Column 0 of a headerless 2-col selected TSV. Same shape the frozen loader consumes."""
    out = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        tok = line.split("\t")[0].strip()
        if tok:
            out.add(tok)
    return out


def coverage(data_raw: Path | None = None) -> dict:
    """What the frozen glob sees vs what is on disk. Pure + offline; no network, no Docker."""
    raw = Path(data_raw) if data_raw else REPO / "data" / "raw"
    seen_files = sorted(Path(p) for p in glob.glob(f"{raw}/*/{FROZEN_PATTERN}"))
    all_files = sorted(Path(p) for p in glob.glob(f"{raw}/*/{BROAD_PATTERN}"))

    seen_dirs = {p.parent.name for p in seen_files}
    all_dirs = {p.parent.name for p in all_files}

    seen_accs: set[str] = set()
    for p in seen_files:
        seen_accs |= _accessions(p)
    all_accs: set[str] = set()
    for p in all_files:
        all_accs |= _accessions(p)

    missed_dirs = sorted(all_dirs - seen_dirs)
    missed_accs = all_accs - seen_accs

    # The only files the broad glob can add. At zero, this run cannot test the blind spot at all.
    variant_files = [p for p in all_files if p.name != FROZEN_PATTERN]
    verdict = classify(len(all_files), len(variant_files), len(missed_accs))
    testable = verdict in (GAP_MEASURED, NO_GAP)

    report = {
        "artifact": "cohort_manifest_coverage_audit",
        "schema": "cohort-manifest-coverage-v1",
        "generated": date.today().isoformat(),
        "frozen_pattern": FROZEN_PATTERN,
        "broad_pattern": BROAD_PATTERN,
        "n_files_seen": len(seen_files),
        "n_files_on_disk": len(all_files),
        "n_cohorts_seen": len(seen_dirs),
        "n_cohorts_on_disk": len(all_dirs),
        "n_cohorts_missed": len(missed_dirs),
        "n_accessions_seen": len(seen_accs),
        "n_accessions_on_disk": len(all_accs),
        "n_accessions_missed": len(missed_accs),
        "missed_fraction": (round(len(missed_accs) / len(all_accs), 4) if all_accs else None),
        "missed_cohorts": missed_dirs,
        "verdict": verdict,
        "comparison_was_testable": testable,
        "n_variant_files_on_disk": len(variant_files),
        "variant_basenames_seen": sorted({p.name for p in variant_files}),
        "why_not_fixed_here": (
            "dna_decode/eval/cohort_manifest.py is one of the five sha256-pinned files in the ACTIVE "
            "prospective lock (wiki/prospective_lock_manifest_2026-08-31.json). Broadening its glob "
            "would invalidate that lock and retire every prospective number scored against it. That is "
            "a user AUTHORITY decision about the frozen surface, not an executor one."
        ),
    }

    if testable:
        report["consequence_if_unfixed"] = (
            "prior_accessions() under-reports the exclusion pool, so a NEW cohort can silently reuse an "
            "accession already consumed by an AR-Bank or external arm. The guard fails OPEN in exactly "
            "the direction it exists to prevent."
        )
        report["honest_limit"] = (
            "This measures DISCOVERY coverage only. It does not claim the missed accessions were "
            "actually reused anywhere — only that the registry cannot answer the question for them."
        )
    else:
        # The gap may well be real; this run simply did not and could not measure it. Shipping the
        # measured-gap prose here would assert a finding off an identity.
        report["consequence_if_unfixed"] = (
            "NOT MEASURED BY THIS RUN. The narrow glob remains a LATENT fail-open risk — if an "
            "AR-Bank or external arm writes `selected_strict.tsv` / `selected_relaxed.tsv`, "
            "prior_accessions() cannot see it — but no such file was on disk, so this run offers no "
            "evidence either way. Re-run with the external-cohort volume attached to measure it."
        )
        report["honest_limit"] = (
            "THE ZERO ABOVE IS AN IDENTITY, NOT A RESULT. `BROAD_PATTERN` is a superset of "
            "`FROZEN_PATTERN`, so with no variant-named file present the two globs are the same "
            "query and `n_accessions_missed` is 0 by construction. This run says nothing about "
            "whether the frozen glob is adequate; it says the corpus that would test it is absent."
        )
    return report


def main() -> int:
    rep = coverage()
    out = REPO / "wiki" / f"cohort_manifest_coverage_{rep['generated']}.json"
    out.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    print(f"cohorts  seen {rep['n_cohorts_seen']:3d} / on-disk {rep['n_cohorts_on_disk']:3d}"
          f"   MISSED {rep['n_cohorts_missed']}")
    print(f"accessions seen {rep['n_accessions_seen']:4d} / on-disk {rep['n_accessions_on_disk']:4d}"
          f"   MISSED {rep['n_accessions_missed']} ({rep['missed_fraction']})")
    print(f"\nverdict: {rep['verdict']}")
    if not rep["comparison_was_testable"]:
        # Say this where a human will read it. A quiet zero is how a vacuous run gets quoted as
        # a clean bill of health.
        print("  The frozen and broad globs resolved to the SAME file set "
              f"({rep['n_variant_files_on_disk']} variant-named files on disk), so MISSED 0 is true "
              "by construction.\n  This run did NOT test the blind spot. Attach the "
              "external-cohort volume and re-run to measure it.")
    if rep["n_cohorts_missed"]:
        print("\nINVISIBLE to the leakage registry:")
        for d in rep["missed_cohorts"]:
            print(f"    {d}")
        print("\nNOT fixed here: cohort_manifest.py is sha256-pinned by the active prospective lock.")
        print("Broadening its glob is a FROZEN-SURFACE edit -> an authority call (see --help).")
    print(f"\nwrote {out.relative_to(REPO).as_posix()}")
    return 0  # a REPORT, never a gate


if __name__ == "__main__":
    raise SystemExit(main())
