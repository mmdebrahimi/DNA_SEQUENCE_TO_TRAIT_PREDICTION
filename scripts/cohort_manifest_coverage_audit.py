"""Measure what the leakage manifest CANNOT SEE — read-only, NON-FROZEN, exit 0 always.

WHY THIS IS A SEPARATE SCRIPT AND NOT A ONE-LINE FIX
----------------------------------------------------
`dna_decode/eval/cohort_manifest.py` discovers cohorts with `glob("data/raw/*/selected.tsv")` —
**exactly** `selected.tsv`. Every AR-Bank and external-cohort arm writes `selected_strict.tsv` /
`selected_relaxed.tsv` instead, so all of them are INVISIBLE to the leakage registry. Measured
2026-10-03: **23 cohort directories and 461 accessions** absent from a fail-closed guard whose entire
job is to answer "was this isolate used before?".

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
    return {
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
        "why_not_fixed_here": (
            "dna_decode/eval/cohort_manifest.py is one of the five sha256-pinned files in the ACTIVE "
            "prospective lock (wiki/prospective_lock_manifest_2026-08-31.json). Broadening its glob "
            "would invalidate that lock and retire every prospective number scored against it. That is "
            "a user AUTHORITY decision about the frozen surface, not an executor one."
        ),
        "consequence_if_unfixed": (
            "prior_accessions() under-reports the exclusion pool, so a NEW cohort can silently reuse an "
            "accession already consumed by an AR-Bank or external arm. The guard fails OPEN in exactly "
            "the direction it exists to prevent."
        ),
        "honest_limit": (
            "This measures DISCOVERY coverage only. It does not claim the missed accessions were "
            "actually reused anywhere — only that the registry cannot answer the question for them."
        ),
    }


def main() -> int:
    rep = coverage()
    out = REPO / "wiki" / f"cohort_manifest_coverage_{rep['generated']}.json"
    out.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    print(f"cohorts  seen {rep['n_cohorts_seen']:3d} / on-disk {rep['n_cohorts_on_disk']:3d}"
          f"   MISSED {rep['n_cohorts_missed']}")
    print(f"accessions seen {rep['n_accessions_seen']:4d} / on-disk {rep['n_accessions_on_disk']:4d}"
          f"   MISSED {rep['n_accessions_missed']} ({rep['missed_fraction']})")
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
