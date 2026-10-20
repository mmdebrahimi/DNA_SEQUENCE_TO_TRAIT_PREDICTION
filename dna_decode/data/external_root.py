"""Find the project's off-volume data root when its drive letter has moved.

WHY THIS EXISTS (measured 2026-10-19). 429 occurrences across 197 files name a
literal `D:/...` path, 389 of them `D:/dna_decode_cache`. The letter is an
accident of which slot the removable disk was in when those lines were written. On
2026-10-19 the user attached that disk and it enumerated as **E:**, while
**D:** was a card-reader slot with no media -- so every one of those 389 paths
pointed at nothing, and the only symptom was datasets reporting as absent.

WHAT THIS DOES *NOT* DO, AND THAT IS THE DESIGN. It does NOT silently substitute a
discovered root for the configured one. On the same host, `E:` really did carry a
directory named `dna_decode_cache` -- nested under `PythonProjects/` -- holding only
`amr_portal` and `mpra`, and **none of CACHE_REQUIRED_SUBDIRS**. Auto-resolving to it
would have pointed 389 call sites at a cache missing the thing they wanted, turning
a loud "root absent" into a quiet "dataset absent". A wrong root that resolves is
worse than a right root that is missing, because only the second one is visible.

So the contract is: `resolve_external_root` PREFERS the configured path and returns
a discovered one only when the caller explicitly allows it; `probe` reports every
candidate together with what it actually CONTAINS, so a human (or a caller with a
required-subdir list) decides. `describe` renders that for a person.
"""
from __future__ import annotations

import string
from dataclasses import dataclass, field
from pathlib import Path

# Directory names that identify this project's off-volume root.
SIGNATURE_DIRS = ("dna_decode_cache", "dna_decode_data")

# Relative prefixes to look under on each drive. "" is the drive root; the others
# are MEASURED layouts, not guesses -- on this host the cache sits under
# `PythonProjects/`, which is why a root-only scan found nothing.
SEARCH_PREFIXES = ("", "PythonProjects", "pythonprojects")

# What a usable genome cache must contain. Named here so a caller can ask
# "is this candidate the one I need?" instead of trusting the directory name.
CACHE_REQUIRED_SUBDIRS = ("refseq",)


@dataclass(frozen=True)
class Candidate:
    """One directory that LOOKS like an off-volume root, plus what it holds."""

    path: Path
    signature: str
    entries: tuple[str, ...] = field(default=())

    def has(self, required: tuple[str, ...]) -> bool:
        have = {e.lower() for e in self.entries}
        return all(r.lower() in have for r in required)

    def missing(self, required: tuple[str, ...]) -> tuple[str, ...]:
        have = {e.lower() for e in self.entries}
        return tuple(r for r in required if r.lower() not in have)


def candidate_drive_roots() -> list[Path]:
    """Every drive root that currently resolves, cheapest possible check.

    A drive with no media (a card reader) does not resolve, which is exactly how
    the real `D:` behaved on the host that motivated this module.
    """
    out: list[Path] = []
    for letter in string.ascii_uppercase:
        root = Path(f"{letter}:/")
        try:
            if root.exists():
                out.append(root)
        except OSError:
            continue
    return out


def probe(
    *,
    drives: list[Path] | None = None,
    signatures: tuple[str, ...] = SIGNATURE_DIRS,
    prefixes: tuple[str, ...] = SEARCH_PREFIXES,
) -> list[Candidate]:
    """Every signature directory found, with a listing of its immediate entries.

    The listing is the point: a directory named `dna_decode_cache` is a NAME match,
    not evidence that it is the right one.
    """
    found: list[Candidate] = []
    seen: set[str] = set()
    for drive in drives if drives is not None else candidate_drive_roots():
        for prefix in prefixes:
            base = drive / prefix if prefix else drive
            for sig in signatures:
                cand = base / sig
                key = str(cand).lower()
                if key in seen:
                    continue
                try:
                    if not cand.is_dir():
                        continue
                    entries = tuple(sorted(c.name for c in cand.iterdir()))
                except OSError:
                    continue
                seen.add(key)
                found.append(Candidate(path=cand, signature=sig, entries=entries))
    return found


def resolve_external_root(
    configured: Path | str,
    *,
    allow_discovered: bool = False,
    required: tuple[str, ...] = (),
    candidates: list[Candidate] | None = None,
) -> tuple[Path, str]:
    """`(path, status)` for the off-volume root to use.

    `status` is one of:
      `configured`            -- the configured path exists; always preferred
      `discovered`            -- configured absent, a candidate satisfied `required`
      `absent-candidates-unfit`-- candidates exist but none satisfied `required`
      `absent`                -- configured absent and nothing was found

    With `allow_discovered=False` (the default) the configured path is returned
    unchanged even when a candidate exists, so behaviour is bit-for-bit what it was
    before this module existed. Substitution is opt-in because a same-named root
    with different CONTENTS is the failure this module was written about.
    """
    conf = Path(configured)
    if conf.is_dir():
        return conf, "configured"
    cands = probe() if candidates is None else candidates
    if not cands:
        return conf, "absent"
    fit = [c for c in cands if c.has(required)] if required else list(cands)
    if not fit:
        return conf, "absent-candidates-unfit"
    if not allow_discovered:
        return conf, "absent-candidates-unfit" if required and not fit else "absent"
    return fit[0].path, "discovered"


def describe(configured: Path | str, *, required: tuple[str, ...] = (),
             candidates: list[Candidate] | None = None) -> str:
    """A human-readable report. Never substitutes; only tells you what is there.

    `candidates` is an INJECTION SEAM, added 2026-10-20. Every other entry point here already
    had one -- `probe` takes `drives=`, `resolve_external_root` takes `candidates=` -- and
    `describe` was the only one that reached the real host unconditionally. The consequence
    was not hypothetical: `test_describe_names_the_missing_requirement` asserted the
    unfit-candidate wording against whatever volumes happened to be attached, so it PASSED
    while one carried a name-matching cache that was missing a required sub-directory, and
    FAILED the moment that volume detached -- a unit test reporting a fact about the
    operator's hardware. Pass candidates explicitly to make a test hermetic; leave it None
    for the live report, which is unchanged.
    """
    conf = Path(configured)
    lines = [f"configured off-volume root : {conf}  [{'EXISTS' if conf.is_dir() else 'ABSENT'}]"]
    cands = probe() if candidates is None else list(candidates)
    if not cands:
        lines.append("discovered candidates    : none on any resolving drive")
        return chr(10).join(lines)
    lines.append(f"discovered candidates    : {len(cands)}")
    for c in cands:
        miss = c.missing(required) if required else ()
        verdict = "UNFIT" if miss else ("fit" if required else "name match")
        lines.append(f"  {c.path}  [{verdict}]")
        lines.append(f"      holds   : {', '.join(c.entries) if c.entries else '(empty)'}")
        if miss:
            lines.append(f"      MISSING : {', '.join(miss)}  <- do not point the project here")
    return chr(10).join(lines)


if __name__ == "__main__":  # pragma: no cover - operator convenience
    print(describe("D:/dna_decode_cache", required=CACHE_REQUIRED_SUBDIRS))
