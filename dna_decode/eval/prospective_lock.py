"""Prospective-lock validation layer for the frozen deterministic AMR decoder (LA-3, 2026-06-22).

The reproducibility freeze (2026-06-13, `wiki/reproducibility_freeze_2026-06-13.md`) named TWO non-foreclosed
forward paths; this is path #2: **prospective-lock**. The idea defeats every circularity / leakage gate
(G1–G8 in the negative-results map) NOT by acquiring a new label, but by TIME:

  Pin the exact frozen decoder NOW (sha256 of its load-bearing files) + stamp a lock DATE. Any isolate whose
  genome first became public STRICTLY AFTER the lock date is a leakage-free test case BY CONSTRUCTION — the
  frozen decoder provably could not have been tuned to data that did not yet exist. No new label acquisition
  required; the validation NUMBER accrues as post-lock measured-AST genomes appear (that is what
  "prospective" means — commit now, score later).

This module is PURE + offline (no network/Docker). It owns three guarantees:
  1. `verify_lock(manifest)` — re-hash the frozen surface; the thing being scored is PROVABLY the locked
     decoder (tamper-evident). Mirrors the TB leak-guard sha256 pin.
  2. `is_prospective_eligible(first_public_date, lock_date)` — an isolate is eligible iff its EARLIEST
     possible public date is strictly after the lock (fail-closed: a missing / year-only / pre-lock date is
     INELIGIBLE, never silently admitted — an un-provable date could be leakage).
  3. `compute_lock_manifest(...)` — the committed, timestamped commitment artifact.

The scorer (`scripts/prospective_lock_validate.py`) reuses the frozen `call_resistance` + the external-cohort
arm's `_conf`; the only NEW gate vs the provenance-disjoint / external arms is the TEMPORAL eligibility
filter here (instead of accession / BioSample overlap).
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from datetime import date
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent.parent

# The load-bearing units that DETERMINE a call (the frozen decoder surface). A prospective result is only
# meaningful if THESE are byte-identical to the frozen decoder — verify_lock proves it. Outputs (report
# card, provdisjoint JSONs) are NOT pinned here: they are products of the decoder, not the decoder.
FROZEN_SURFACE_FILES: tuple[str, ...] = (
    "dna_decode/eval/amr_rules.py",                 # the call_resistance engine
    "dna_decode/data/calibrated_amr_rules.json",    # the deployed per-(organism,drug) config
    "dna_decode/data/mic_tiers.py",                 # CLSI/EUCAST breakpoints + mechanism→loci catalogs
    "dna_decode/data/shipped_decoder_surface.py",   # the authoritative deployed-claim grid
    "dna_decode/eval/cohort_manifest.py",           # the leakage registry (fail-closed discipline)
)

LOCK_SCHEMA = "prospective-lock-manifest-v1"

# The 10 SCORED cells the freeze validated (the grid a prospective test extends). organism keys are the
# call_resistance registry organism; drug is the CLI drug name.
SCORED_CELLS: tuple[tuple[str, str], ...] = (
    ("Campylobacter", "ciprofloxacin"),
    ("Escherichia_coli_Shigella", "ciprofloxacin"),
    ("Escherichia_coli_Shigella", "ceftriaxone"),
    ("Escherichia_coli_Shigella", "gentamicin"),
    ("Escherichia_coli_Shigella", "tetracycline"),
    ("Klebsiella", "ciprofloxacin"),
    ("Klebsiella", "ceftriaxone"),
    ("Klebsiella", "gentamicin"),
    ("Klebsiella", "meropenem"),
    ("Klebsiella", "tetracycline"),
)


def _sha256_file(rel: str, repo: Path = _REPO) -> str:
    return hashlib.sha256((repo / rel).read_bytes()).hexdigest()


def surface_hashes(repo: Path = _REPO) -> dict[str, str]:
    """Current sha256 of every frozen-surface file (the live decoder state)."""
    return {rel: _sha256_file(rel, repo) for rel in FROZEN_SURFACE_FILES}


def compute_lock_manifest(lock_date: str, frozen_commit: str, lock_commit: str,
                          repo: Path = _REPO) -> dict:
    """The timestamped, hash-pinned prospective-lock commitment.

    `lock_date` (YYYY-MM-DD) is the cutoff: an isolate is an eligible prospective test case iff it became
    public strictly after this date. `frozen_commit` = the reproducibility-freeze commit (b3761c8);
    `lock_commit` = the commit at which this manifest was created (the decoder files are byte-identical to
    the freeze, which the hashes prove regardless of commit)."""
    return {
        "schema": LOCK_SCHEMA,
        "lock_date": lock_date,
        "frozen_commit": frozen_commit,
        "lock_commit": lock_commit,
        "eligibility_rule": ("an isolate is a leakage-free prospective test case IFF its earliest possible "
                             "public date (INSDC first_public / release_date) is STRICTLY AFTER lock_date; "
                             "missing / year-only / pre-lock dates are INELIGIBLE (fail-closed)"),
        "protocol": ("verify_lock (decoder == frozen) -> filter cohort to prospective-eligible -> frozen "
                     "call_resistance(organism, drug) -> independent_cohort_validate._conf; the result is "
                     "stamped with this manifest's surface hashes + lock_date so it is provably prospective"),
        "surface_sha256": surface_hashes(repo),
        "scored_cells": [{"organism": o, "drug": d} for o, d in SCORED_CELLS],
        "what_this_is_not": ("NOT a validation number today — the number accrues as post-lock measured-AST "
                             "genomes appear. This artifact is the tamper-evident commitment + eligibility "
                             "rule that makes any later score PROVABLY prospective (leakage-free by time)."),
    }


@dataclass(frozen=True)
class LockVerification:
    ok: bool
    drifted: list[str]          # files whose live sha256 != the manifest
    missing: list[str]          # files in the manifest absent on disk


def verify_lock(manifest: dict, repo: Path = _REPO) -> LockVerification:
    """Re-hash the frozen surface; confirm the LIVE decoder is byte-identical to the locked one.

    A prospective score is only honest if the decoder scoring the post-lock data is the SAME decoder that
    was committed at lock time. Any drift / missing file => the lock is broken; the scorer MUST refuse."""
    pinned: dict[str, str] = manifest.get("surface_sha256", {})
    drifted, missing = [], []
    for rel, want in pinned.items():
        p = repo / rel
        if not p.exists():
            missing.append(rel)
            continue
        if _sha256_file(rel, repo) != want:
            drifted.append(rel)
    return LockVerification(ok=not drifted and not missing, drifted=sorted(drifted), missing=sorted(missing))


def _earliest_possible_date(date_str: str | None) -> date | None:
    """Parse an INSDC-style date to the EARLIEST day it could denote, or None if unusable.

    Full `YYYY-MM-DD` -> that day. `YYYY-MM` -> first of month. Bare `YYYY` -> None (too coarse to PROVE
    post-lock — fail-closed). Anything unparseable -> None. Conservative on purpose: we admit an isolate
    only when even its earliest possible interpretation is after the lock."""
    if not date_str:
        return None
    s = str(date_str).strip()[:10]
    parts = s.replace("/", "-").split("-")
    try:
        if len(parts) >= 3 and all(parts[:3]):
            return date(int(parts[0]), int(parts[1]), int(parts[2]))
        if len(parts) == 2 and all(parts):
            return date(int(parts[0]), int(parts[1]), 1)
    except (ValueError, TypeError):
        return None
    return None  # bare year (or junk) -> cannot prove post-lock


@dataclass(frozen=True)
class EligibilityVerdict:
    eligible: bool
    reason: str                 # post_lock / pre_or_equal_lock / undatable_fail_closed


def is_prospective_eligible(first_public_date: str | None, lock_date: str) -> EligibilityVerdict:
    """True iff the isolate provably became public AFTER lock_date (leakage-free by construction).

    Fail-closed: an undatable / year-only / on-or-before-lock isolate is INELIGIBLE — it cannot be PROVEN to
    postdate the frozen decoder, so admitting it would risk leakage (the exact trap prospective-lock avoids)."""
    lock = _earliest_possible_date(lock_date)
    if lock is None:
        raise ValueError(f"lock_date must be a full YYYY-MM-DD; got {lock_date!r}")
    earliest = _earliest_possible_date(first_public_date)
    if earliest is None:
        return EligibilityVerdict(False, "undatable_fail_closed")
    if earliest > lock:
        return EligibilityVerdict(True, "post_lock")
    return EligibilityVerdict(False, "pre_or_equal_lock")
