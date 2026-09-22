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

import json
import hashlib
from dataclasses import dataclass, field
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
                          repo: Path = _REPO, manifest_created: str = "",
                          cutoff_justification: str = "") -> dict:
    """The timestamped, hash-pinned prospective-lock commitment.

    `lock_date` (YYYY-MM-DD) is the eligibility CUTOFF: an isolate is an eligible prospective test case iff it
    became public strictly after this date. It is the decoder-FREEZE date (when the decoder state became
    immutable), NOT necessarily the manifest-creation date — the manifest is a tamper-evident RECORD of a
    state that already existed, so an isolate postdating the freeze is provably unseen regardless of when the
    record was written. `frozen_commit` = the reproducibility-freeze commit (b3761c8); `lock_commit` = the
    commit at which this manifest was created; `manifest_created` = that creation date; `cutoff_justification`
    = why `lock_date` is honest (e.g. byte-identity of the frozen surface to `frozen_commit`)."""
    return {
        "schema": LOCK_SCHEMA,
        "lock_date": lock_date,
        "manifest_created": manifest_created or lock_date,
        "cutoff_justification": cutoff_justification,
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
    # Files the manifest FAILED TO PIN (or pinned that are not part of the frozen surface). Distinct
    # from `missing`, which means "pinned but absent from disk": this means the manifest never made
    # the commitment at all, so there is nothing to verify rather than something that failed.
    incomplete_pin: list[str] = field(default_factory=list)


def verify_lock(manifest: dict, repo: Path = _REPO) -> LockVerification:
    """Re-hash the frozen surface; confirm the LIVE decoder is byte-identical to the locked one.

    A prospective score is only honest if the decoder scoring the post-lock data is the SAME decoder that
    was committed at lock time. Any drift / missing file => the lock is broken; the scorer MUST refuse."""
    pinned: dict[str, str] = manifest.get("surface_sha256", {})
    # A manifest that pins NOTHING must not verify. Iterating an empty dict finds no drift and no
    # missing file, so the original returned ok=True for `surface_sha256: {}`, for a missing key, and
    # for a single-file subset -- a lock that checks nothing, reported as a lock that holds. Harmless
    # while every caller named one trusted manifest; NOT harmless once `resolve_active_lock` globs every
    # committed manifest and picks a verifying one, because a vacuous manifest would be selected
    # EXACTLY when the real lock stops verifying -- i.e. it would defeat fail-closed at the only moment
    # fail-closed matters. Same shape as the vacuous filters this project has caught three times before.
    if set(pinned) != set(FROZEN_SURFACE_FILES):
        return LockVerification(
            ok=False,
            drifted=[],
            missing=sorted(set(FROZEN_SURFACE_FILES) - set(pinned)),
            incomplete_pin=sorted(set(pinned) ^ set(FROZEN_SURFACE_FILES)),
        )
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


class NoActiveLock(RuntimeError):
    """No committed manifest describes the LIVE frozen surface."""


_HEX64 = 64


def is_well_formed_manifest(manifest: dict) -> list[str]:
    """Reasons this is not a valid prospective-lock COMMITMENT. Empty list = well formed.

    Hash agreement alone does not make a file a lock. A manifest carrying the correct hashes but no
    schema and `lock_date: "not-a-date"` was selectable as the active lock (verified 2026-09-22), and
    that bogus cutoff would then have flowed straight into the eligibility filter. Verifying answers
    "does this describe the live decoder"; this answers "is this a commitment at all".

    DELIBERATELY NOT CHECKED: agreement between the filename date and `lock_date`. It is tempting and
    it is FALSE here -- `prospective_lock_manifest_2026-06-22.json` legitimately carries lock_date
    2026-06-13, because the manifest RECORDS a freeze that already happened rather than creating one.
    Enforcing it would reject the repo's own real manifest.
    """
    problems: list[str] = []
    if manifest.get("schema") != LOCK_SCHEMA:
        problems.append(f"schema is {manifest.get('schema')!r}, expected {LOCK_SCHEMA!r}")
    raw = manifest.get("lock_date")
    try:
        date.fromisoformat(str(raw))
    except (TypeError, ValueError):
        problems.append(f"lock_date {raw!r} is not an ISO YYYY-MM-DD date")
    pinned = manifest.get("surface_sha256")
    if not isinstance(pinned, dict):
        problems.append("surface_sha256 is missing or not an object")
    else:
        bad = sorted(f for f, h in pinned.items()
                     if not (isinstance(h, str) and len(h) == _HEX64
                             and all(c in "0123456789abcdef" for c in h.lower())))
        if bad:
            problems.append(f"surface_sha256 entries are not 64-hex digests: {bad}")
    return problems


def resolve_active_lock(repo: Path = _REPO) -> tuple[Path, dict]:
    """The lock manifest that describes the surface as it is RIGHT NOW -- the single source of truth
    for the eligibility cutoff.

    WHY THIS EXISTS (found 2026-09-22). `scripts/fetch_prospective_cohort.py` carried
    `LOCK_DATE = "2026-06-13"` as a hardcoded module constant and used it as the `--lock-date` default.
    The v2 gentamicin lock moved the cutoff to 2026-08-31 on 2026-08-31, and the constant did not follow,
    so the sweep would have collected isolates dated after the RETIRED v1 cutoff and labelled them
    prospective. For the deployed rule those are PRE-lock -- a false prospective claim, which is exactly
    the leakage the lock exists to prevent. Same class as the salmserovar CLI that ran a threshold
    replaced five days earlier ("a shipped default is a promise").

    The resolution rule is definitional rather than conventional: a manifest pins sha256 hashes of the
    frozen surface, so the manifest that still VERIFIES is the one describing the live decoder, and one
    that no longer verifies describes a retired one. Nothing has to be kept in sync by hand -- retiring a
    lock is precisely the act of making its manifest stop verifying.

    Raises NoActiveLock when none verifies. That is FAIL-CLOSED on purpose: falling back to the newest
    file, or to a default, is how a stale cutoff silently becomes a prospective claim.
    """
    candidates = sorted((repo / "wiki").glob("prospective_lock_manifest_*.json"))
    verified: list[tuple[Path, dict]] = []
    malformed: list[str] = []
    for p in candidates:
        try:
            m = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        # A file must be a well-formed COMMITMENT before its hashes are allowed to mean anything.
        # Malformed manifests are collected and reported rather than silently skipped: one sitting
        # unnoticed in wiki/ is exactly how a bad cutoff would enter later.
        problems = is_well_formed_manifest(m)
        if problems:
            malformed.append(f"{p.name}: {'; '.join(problems)}")
            continue
        if verify_lock(m, repo=repo).ok:
            verified.append((p, m))
    if not verified:
        detail = f" Malformed manifests ignored: {malformed}." if malformed else ""
        raise NoActiveLock(detail + 
            "no committed prospective-lock manifest verifies against the live frozen surface "
            f"(checked {len(candidates)}). The surface has drifted from every recorded lock, so there "
            "is no honest eligibility cutoff -- re-lock before sweeping for prospective isolates.")
    # AMBIGUITY REFUSES rather than picking the later cutoff (corrected 2026-09-22).
    #
    # The original rule took max(lock_date) and called it conservative because a later cutoff admits
    # strictly fewer isolates. That reasoning is incomplete: two manifests can only BOTH verify if they
    # pin identical hashes, i.e. they describe the SAME decoder while claiming DIFFERENT cutoffs. That
    # is a protocol contradiction, not a tie -- and silently taking the later one would let an arbitrary
    # cutoff DISCARD already-accrued evidence that happens to be unfavourable, while looking cautious.
    # "Admits fewer isolates" is not the same as "is honest about which isolates it may admit".
    #
    # Identical lock_dates are NOT ambiguous: those manifests make the same commitment, so any of them
    # answers the question and the earliest path is returned for determinism.
    cutoffs = {str(m.get("lock_date", "")) for _, m in verified}
    if len(cutoffs) > 1:
        raise NoActiveLock(
            "AMBIGUOUS LOCK: several committed manifests verify against the live frozen surface but "
            f"claim DIFFERENT cutoffs {sorted(cutoffs)} ({[p.name for p, _ in verified]}). They pin the "
            "same decoder, so they cannot both be the eligibility rule. Retire one (its manifest must "
            "stop verifying) or record an explicit supersession before sweeping.")
    return min(verified, key=lambda pm: pm[0].name)


def active_lock_date(repo: Path = _REPO) -> str:
    """The eligibility cutoff of the lock that describes the live surface."""
    return str(resolve_active_lock(repo)[1]["lock_date"])
