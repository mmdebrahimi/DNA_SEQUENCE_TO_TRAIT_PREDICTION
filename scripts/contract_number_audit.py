"""Check every number a cell contract CITES against the artifact it cites it FROM.

WHY
---
This repo's most frequently-recurring defect class is not wrong code, it is a claim that DRIFTED from the
artifact it rests on. Measured instances, all found by comparing a claim to its own evidence rather than by
any test: `project_status.py` quoting a figure its own data layer forbids; the report card understating the
validated surface by 4 cells; `--help` advertising two HLA screens the repo had already demoted; ResFinder's
"165 beta-lactamase genes"; and this week CYP4F2's contract citing a 54/54 concordance that was a different
axis computed on a smaller parsed set, plus a "~79% EAS" allele frequency whose real value is 0.2128.

`tests/test_cell_registry.py` already pins that every cited artifact EXISTS. Nothing checks whether the
NUMBERS in the prose are in the artifact. That is the gap this closes.

WHAT IT IS, HONESTLY
--------------------
A TRIAGE FUNNEL, not an oracle -- the same posture as scripts/kaggle_staleness_auditor.py, whose field
precision is about 1 in 6. A cited number can legitimately be absent from its artifact: it may be derived
(a ratio, a sum, a percentage of a raw count), rounded differently, or carried over from a superseded run
that the prose correctly describes as superseded. So a miss is CANDIDATE_DRIFT requiring adjudication, and
this script NEVER edits anything. Adjudication is the product.

Scope: only contracts that name a `wiki/...` artifact, and only decimal numbers (a bare integer in prose is
usually a count whose denominator lives elsewhere, which this cannot check without guessing).
"""
from __future__ import annotations

import bisect
import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

PROSE_FIELDS = ("claim", "validation_slice", "label_provenance", "demotion_rule", "claim_status")
# Contracts cite artifacts in TWO styles, and both must resolve or the audit invents drift. The second is
# shell-style brace expansion, used throughout this repo:
#     wiki/cyp4f2_getrm_concordance_2026-09-24.{md,json}
#     wiki/forward_inverse_{roundtrip,sweep,deployable}_2026-07-1{6,7}.{md,json}
# Missing it produced a false CANDIDATE_DRIFT on the inverse cell's "+53.0%", whose source
# (margin_vs_better_baseline = 0.5295) sits in a roundtrip artifact the regex never loaded.
ARTIFACT_RE = re.compile(r"wiki/[A-Za-z0-9_][A-Za-z0-9_./{},-]*")
BRACE_RE = re.compile(r"\{([^{}]*)\}")


def _expand_braces(path: str) -> list[str]:
    """Expand one or more {a,b,c} groups into every concrete path."""
    out = [path]
    while any(BRACE_RE.search(p) for p in out):
        nxt = []
        for p in out:
            m = BRACE_RE.search(p)
            if not m:
                nxt.append(p)
                continue
            for alt in m.group(1).split(","):
                nxt.append(p[:m.start()] + alt.strip() + p[m.end():])
        out = nxt
    return out
# A decimal immediately preceded by a letter/underscore is part of an IDENTIFIER, not a measurement --
# canFam3.1, GCA_000005845.2, v1.2. Capturing those produced a false CANDIDATE_DRIFT on the dog
# morphology cell's "3.1", which is the assembly name canFam3.1.
# The trailing guard rejects a DOTTED CONTINUATION (`1.2.3`, a version string) but NOT a sentence-final
# period. The original `(?![\d.])` did both, so any number ending a sentence was invisible to the audit --
# 5 of them, found only because a citation displaced the period in front of `AUROC 0.580` and the extracted
# count moved by one. `(?!\.\d)` keeps the version-string rejection that guard existed for.
DECIMAL_RE = re.compile(r"(?<![\w.])(\d{1,4}\.\d{1,4})(?![\d])(?!\.\d)")

# Numbers that are not measurements: dates, tool versions, guideline years, p-value exponents.
NON_MEASUREMENT = re.compile(r"^(19|20)\d{2}\.")

# Numbers ADJUDICATED once as legitimately absent from any artifact, with the reason. Deliberately a small
# explicit table keyed by (cell_id, number): a blanket "is it anywhere in the source tree" check was tried
# and is vacuous (it accepted 12/12 randomly-generated numbers), so each exemption is named or it flags.
ADJUDICATED_BENIGN: dict[tuple[str, str], str] = {
    ("typing:Escherichia_coli:serotype", "0.60"):
        "not a measurement: the project's own source-diversity BAR, defined in code as "
        "scripts/source_diverse_validate.MAX_SOURCE_SHARE. A threshold has no artifact to live in.",
}


# --- provenance kinds -------------------------------------------------------------------------------
# "Unverifiable" was doing the work of four different statements at once, which made 34 numbers read as
# 34 defects. They are not one class: some are real measurements whose artifact simply was not cited,
# some are external reference values no artifact of OURS could ever contain, some are pinned by a test
# rather than by a wiki file, and some are not measurements at all -- a DOI prefix and two kb lengths the
# decimal extractor cannot tell from a ratio.
#
# The two structural families are PATTERN-typed rather than listed, because a per-number list is the
# hand-enumerated-exclusion trap this repo has hit five times. The two declarative families ARE tables,
# because "this number is pinned by that test" is a claim about provenance that cannot be derived from
# the prose and must be named by a human, with a reason, exactly like ADJUDICATED_BENIGN. They are kept
# SEPARATE from ADJUDICATED_BENIGN so its `<= 5` pin stays untouched and each table stays small.
KIND_ARTIFACT = "artifact"
KIND_STRUCTURAL = "structural-non-measurement"
KIND_EXTERNAL = "external-reference"
KIND_ENFORCED = "enforced-by-test"
KIND_THRESHOLD = "threshold"
KIND_SUPERSEDED = "superseded-value"
KIND_DERIVED = "derived"
KIND_DRIFT = "candidate-drift"
KIND_UNVERIFIABLE = "unverifiable"

# A DOI prefix: the registrant half of `10.5281/zenodo.14065540`. Only structural when the very next
# character is the `/` that makes it a DOI -- `10.5281` alone would be an ordinary number.
_DOI_TOKEN_RE = re.compile(r"^10\.\d{4,9}$")
# A physical length: `5.1-kb`, `4.6 kb`, `11-bp`. At most ONE separator, so a real measurement merely in
# the same sentence as the word "kb" is never swept up -- that over-broad direction would silently
# exempt a genuine number, which is the failure this whole audit exists to catch.
_LENGTH_SUFFIX_RE = re.compile(r"^[-\s]?(kb|bp)\b", re.IGNORECASE)

EXTERNAL_REFERENCE: dict[tuple[str, str], str] = {
    ("pgx:human:nudt15", "9.5"):
        "not our measurement: the *3 allele frequency in EAS populations, a published population-genetics "
        "reference value. No artifact of ours can or should contain it; citing one would fabricate "
        "provenance for a number we did not measure.",
}

ENFORCED_BY_TEST: dict[tuple[str, str], str] = {
    ("typing:Escherichia_coli:pathotype", "0.833"):
        "pinned by an EXACT-equality assert in tests/test_pathotype_expec_recall.py (EXPEC_RECALL_CAP "
        "== 10/12). Enforced by the suite on every run, which is a STRONGER guarantee than a wiki file; "
        "the cell has no artifact by design and inventing one would be worse than naming the test.",
    ("typing:Escherichia_coli:pathotype", "1.0"):
        "pinned by the same EXACT-equality asserts (confident-supported precision == 1.0 and EPEC "
        "recall == 1.0) in tests/test_pathotype_expec_recall.py.",
}

# A number the prose deliberately records as WRONG. Distinguishing this from `unverifiable` matters: the
# contract carries it precisely so nobody quotes it, so "we cannot verify it" is the wrong description --
# and filing it as `enforced-by-test` beside the value that SUPERSEDED it would be inaccurate.
SUPERSEDED_VALUE: dict[tuple[str, str], str] = {
    ("typing:Escherichia_coli:pathotype", "0.917"):
        "deliberately recorded as the REJECTED over-rescue figure (the flat-K=1 rule's 11/12), kept in "
        "the prose so it is never quoted as the cell's recall. Superseded by 0.833 on the committed "
        "decision 'a clean 0.833 beats an overfit 0.917'.",
}

# A number the prose COMPUTES from counts it also states. This script's own posture has always said a
# cited number may legitimately be DERIVED -- a ratio, a sum, a percentage of a raw count -- but there was
# no kind for it, so a derived value could only pass by being coincidentally present, or fail as drift.
# Boundary-aware matching is what made this necessary: it correctly stopped accepting `0.862` because
# `0.86` sits inside some artifact's `0.867`, which was never evidence of anything.
#
# THE REASON MUST CARRY THE ARITHMETIC. "It's derived" with no numerator and denominator is an unaudited
# exemption; with them, a reader re-does the division in their head and the claim is checkable without any
# artifact at all. Each entry below states the fraction AND where its inputs are readable.
DERIVED_VALUE: dict[tuple[str, str], str] = {
    ("typing:bacteriophage:phage", "0.862"):
        "derived, and the prose states the fraction it comes from: '25/29 called = 0.862' "
        "(25/29 = 0.8621). Both counts are in the same clause, so the quotient is checkable by "
        "arithmetic; no artifact needs to store it.",
    ("typing:bacteriophage:phage", "0.291"):
        "derived, stated in the prose as '25/86=0.291' (25/86 = 0.2907). Same clause carries both counts.",
    ("typing:Salmonella:salmserovar", "0.900"):
        "derived coverage: the prose states 200 isolates and 'abstention 10.0%', so covered = 180/200 = "
        "0.900. Both inputs are in the same passage; the artifact stores the counts, not the rate.",
}

# Bars and cut-offs, which have no artifact to live in by definition. Same role as ADJUDICATED_BENIGN,
# held separately so that table's `<= 5` pin is untouched.
DECLARED_THRESHOLD: dict[tuple[str, str], str] = {
    ("typing:Escherichia_coli:pathotype", "0.80"):
        "a BAR, not a measurement: the per-gene coverage floor each axis must clear in the cross-axis "
        "support rule (>=1 iron-acquisition AND >=1 capsule/serum gene, each >=0.80).",
    ("typing:klebsiella:kleb", "0.90"):
        "a BAR, not a measurement: the greedy-representative clonality-correction threshold "
        "(greedy-rep @0.90) at which the leave-one-out was run.",
    ("typing:klebsiella:kleb", "0.10"):
        "not a measurement of this cell: the PRIOR NULL baseline (a 0.10 prior) that the reported lift "
        "is measured against.",
}


def _structural_kind(blob: str, tok: str) -> str | None:
    """Is every occurrence of `tok` in this prose a DOI prefix or a physical length?

    Requires ALL occurrences to be structural. A token used once as `5.1-kb` and once as a real ratio
    must NOT be exempted -- typing it structural would remove a genuine measurement from the audit
    entirely, which is a check quietly stopping checking.
    """
    spans = [m for m in DECIMAL_RE.finditer(blob) if m.group(1) == tok]
    if not spans:
        return None
    kinds = set()
    for m in spans:
        rest = blob[m.end():]
        if _DOI_TOKEN_RE.match(tok) and rest.startswith("/"):
            kinds.add(KIND_STRUCTURAL)
        elif _LENGTH_SUFFIX_RE.match(rest):
            kinds.add(KIND_STRUCTURAL)
        else:
            return None
    return KIND_STRUCTURAL if kinds else None


def _declared_kind(cell_id: str, tok: str) -> str | None:
    """Kinds a human has NAMED for this exact (cell, number), each with a mandatory reason."""
    key = (cell_id, tok)
    if key in ADJUDICATED_BENIGN or key in DECLARED_THRESHOLD:
        return KIND_THRESHOLD
    if key in ENFORCED_BY_TEST:
        return KIND_ENFORCED
    if key in EXTERNAL_REFERENCE:
        return KIND_EXTERNAL
    if key in SUPERSEDED_VALUE:
        return KIND_SUPERSEDED
    if key in DERIVED_VALUE:
        return KIND_DERIVED
    return None


def _declared_reason(cell_id: str, tok: str) -> str:
    """The mandatory reason behind a declared kind. A declaration without one is an unaudited exemption."""
    key = (cell_id, tok)
    for table in (ADJUDICATED_BENIGN, DECLARED_THRESHOLD, ENFORCED_BY_TEST,
                  EXTERNAL_REFERENCE, SUPERSEDED_VALUE, DERIVED_VALUE):
        if key in table:
            return table[key]
    return ""


def _classify_kind(cell_id: str, blob: str, tok: str, *, matched_in_artifact: bool,
                   cited: bool) -> str:
    """Precedence is deliberate: a pattern-typed non-measurement and a human declaration both outrank
    an artifact match, because a DOI prefix that happens to appear in some JSON is still not a
    measurement. `candidate-drift` is kept distinct from `unverifiable` -- cited-and-absent is a finding,
    uncited is merely unchecked, and collapsing them would hide real drift inside a coverage statistic.
    """
    structural = _structural_kind(blob, tok)
    if structural:
        return structural
    declared = _declared_kind(cell_id, tok)
    if declared:
        return declared
    if matched_in_artifact:
        return KIND_ARTIFACT
    return KIND_DRIFT if cited else KIND_UNVERIFIABLE


def _cited_artifacts(blob: str) -> list[str]:
    out, seen = [], set()
    for raw in ARTIFACT_RE.findall(blob):
        for p in _expand_braces(raw.rstrip(".,;:)")):
            p = p.rstrip(".,;:)_")
            # a contract may cite a stem without its extension (…_2026-09-24) -> try both sidecars
            for cand in ([p] if p.endswith((".json", ".md")) else [p + ".json", p + ".md", p]):
                if cand not in seen:
                    seen.add(cand)
                    out.append(cand)
    return out


def _load_artifact_text(rel: str) -> str | None:
    p = REPO / rel
    if not p.exists() or p.is_dir():
        return None
    try:
        raw = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    if p.suffix == ".json":
        try:
            # re-dump so nested numeric values are searchable as plain text
            return raw + "\n" + json.dumps(json.loads(raw))
        except json.JSONDecodeError:
            return raw
    return raw


# Discrimination bands, as RATES against the full pool. BOTH are asserted, not derived -- say so rather
# than implying a calibration that does not exist.
#
# The HIGH bar is NOT `rate == 0`, and that is a measured correction rather than a relaxation. Under the
# old 10-artifact sample, 0-of-10 was the best OBSERVABLE outcome and several cells reached it. Against
# the full ~600-trial pool, no cell in the registry reaches exactly zero -- not even
# `typing:Salmonella:salmserovar`, which carries THIRTY measured figures and still has 1 of 119 five-file
# haystacks containing all of them. A grade no cell can earn is not a grade, and collapsing everything
# below 50% into MODERATE would file resfinder (1 of 602) beside phage (246 of 602). So HIGH is the
# conventional 5%: the citation would have been falsified by more than 19 of every 20 alternatives.
# The per-cell `decoy_match_rate` is reported regardless, so a reader is never limited to the band.
_HIGH_RATE = 0.05
_LOW_RATE = 0.50


def _grade(decoy_full: int, n_decoys: int, nums: list[str]) -> str:
    if not n_decoys or not nums:
        return "n/a"
    rate = decoy_full / n_decoys
    if rate < _HIGH_RATE:
        return "HIGH"
    if rate >= _LOW_RATE:
        return "LOW"
    return "MODERATE"


def _cached_artifact(p: str, cache: dict) -> tuple[str | None, list[float]]:
    if p not in cache:
        t = _load_artifact_text(p)
        cache[p] = (t, _artifact_numbers(t) if t else [])
    return cache[p]


def _decoy_full_match_count(nums: list[str], own: set[str], pool: list[str], cache: dict,
                            group_size: int = 1) -> tuple[int, int]:
    """How many UNRELATED artifact haystacks also contain every one of this cell's numbers.

    A citation only verifies a number if the number could have failed to match. This repo already
    learned that the hard way one level up: the audit once reached a clean 0/130 by tracing numbers to
    the package source, which accepted 12 of 12 RANDOMLY GENERATED numbers. The same vacuity applies to
    a cited artifact -- a cell whose figures are few and round (`0.3`, `0.94`, `14.0`) matches a numeric
    dense JSON about anything. So the decoys are the control: 0 full matches means the citation is
    doing real work; a high match rate means the cell is CHECKED but not meaningfully verified.

    DETERMINISTIC FULL-POOL SCAN, and the reason it is not a sample
    --------------------------------------------------------------
    This was a 10-artifact random sample with a fixed seed, and the sample was UNDERPOWERED to the point
    of being wrong in both directions. Measured against the full pool: `typing:Klebsiella:ktype` is 308
    of 602 (51%) where the sample said 5 of 10, and `essentiality:any:essentiality` matches 46 unrelated
    artifacts where its sample said 0 of 10 -- which graded it HIGH, the strongest possible verdict, for
    a set of seven numbers that a yeast benchmark, two TB baselines and a pigment run each fully contain.
    A grade is the thing that says whether a citation verifies anything, so an underpowered grade is
    worse than no grade. There is no seed here because there is no sampling: it scans the whole pool.

    SIZE-MATCHED HAYSTACK, and the grouping rule
    --------------------------------------------
    The second defect in the sampled version: a cell citing N artifacts was scored against a haystack of
    N concatenated files while each decoy was a SINGLE file. Bigger haystack, more numbers, easier full
    match -- so multi-citation cells earned discrimination for free (a 33-file bundle scored 6/6 against
    single-file decoys). Each decoy trial therefore concatenates `group_size` pool artifacts, taken as
    CONSECUTIVE NON-OVERLAPPING groups over the sorted pool -- deterministic, no sampling. `group_size`
    is the count of the cell's cited artifacts that are IN THE POOL (json), because the pool is what the
    decoys are drawn from; a trailing partial group is dropped so every trial is exactly that size.
    This grouping is a design choice with no measured optimum -- a different grouping could give a
    different rate for a multi-citation cell -- so it is reported per cell rather than left implicit.

    Reported alongside drift, never folded into it -- low discrimination is weak evidence, not a defect.
    """
    decoys = [p for p in pool if p not in own]
    n = max(1, group_size)
    if len(decoys) < n or not nums:
        return 0, 0
    full, trials = 0, 0
    for i in range(0, len(decoys) - n + 1, n):
        # Keep each member's numbers SEPARATE rather than concatenating them. A rounding match against the
        # union is exactly a match against SOME member, so the union is unnecessary -- and building it
        # copied the whole cached float list on every one of ~600 trials per cell, which is what kept this
        # at minutes instead of seconds. Same for the text: only join when the group is actually >1.
        texts, num_lists = [], []
        for p in decoys[i:i + n]:
            t, a = _cached_artifact(p, cache)
            if t is None:
                continue
            texts.append(t)
            num_lists.append(a)
        if not texts:
            continue
        trials += 1
        hay = texts[0] if len(texts) == 1 else "\n".join(texts)
        if all(any(_bounded_in(v, hay) for v in _number_variants(x))
               or any(_matches_by_rounding(x, an) is not None for an in num_lists)
               for x in nums):
            full += 1
    return full, trials


def _artifact_numbers(text: str) -> list[float]:
    """Every numeric token in the artifact, as floats, SORTED, for rounding-aware comparison.

    Sorted so `_matches_by_rounding` can bisect three ranges instead of scanning tens of thousands of
    tokens per candidate. Nothing depends on document order -- the only consumer is the tolerance check.
    """
    out = []
    for m in re.findall(r"-?\d+(?:\.\d+)?", text):
        try:
            out.append(float(m))
        except ValueError:
            pass
    out.sort()
    return out


def _matches_by_rounding(tok: str, art_nums: list[float]) -> float | None:
    """Is the cited number a ROUNDING of some number in the artifact (or vice versa)?

    Load-bearing: prose legitimately rounds. This audit's first run flagged SLCO1B1's cited 0.398 as drift
    when the artifact holds 0.3977 -- a false positive produced by string matching alone, which would have
    sent a reader to adjudicate a number that is correct. Compare numerically at the cited precision, and
    also allow the percent/fraction pairing.
    """
    try:
        cited = float(tok)
    except ValueError:
        return None
    dec = len(tok.split(".")[1]) if "." in tok else 0
    # Tolerance, NOT round(): round(52.95, 1) is 52.9 because 52.95 has no exact binary representation,
    # so a cited "53.0" would spuriously fail against an artifact's 0.5295. Half a unit at the cited
    # precision is what "rounds to" means, and it is immune to float-repr edge cases.
    tol = 0.5 * (10.0 ** -dec) + 1e-9
    # `cand in [cited-tol, cited+tol]` for cand in (a, a*100, a/100) is equivalent to `a` falling in one of
    # three ranges. On a SORTED art_nums that is three bisects instead of a full scan -- which matters
    # because whole-number matching sends far more numbers down this fallback than substring matching did
    # (it took the audit from 14 s to 7m46s before this).
    lo, hi = cited - tol, cited + tol
    for a, b in ((lo, hi), (lo / 100.0, hi / 100.0), (lo * 100.0, hi * 100.0)):
        i = bisect.bisect_left(art_nums, a)
        if i < len(art_nums) and art_nums[i] <= b:
            return art_nums[i]
    return None


_SOURCE_BLOB: str | None = None


def _package_source_blob() -> str:
    """Concatenated package + scripts source, for tracing a cited THRESHOLD to the code that defines it."""
    global _SOURCE_BLOB
    if _SOURCE_BLOB is None:
        parts = []
        for d in ("dna_decode", "scripts"):
            for p in sorted((REPO / d).rglob("*.py")):
                try:
                    parts.append(p.read_text(encoding="utf-8", errors="replace"))
                except OSError:
                    pass
        _SOURCE_BLOB = "\n".join(parts)
    return _SOURCE_BLOB


def _in_package_source(tok: str) -> bool:
    return any(v in _package_source_blob() for v in _number_variants(tok))


def _bounded_in(v: str, hay: str) -> bool:
    """Is `v` present in `hay` as a WHOLE number rather than a fragment of a longer one?

    THE DEFECT THIS FIXES. The prose side has always been boundary-guarded (`DECIMAL_RE` carries
    `(?<![\\w.])` and `(?![\\d.])`), but the artifact side was a raw `v in hay` substring test. So a short
    token matched INSIDE a longer number and the audit accepted it as provenance:
        `2.1`    matched inside `72.1`   (a coverage metric, in a pneumo artifact)
        `0.51`   matched inside `0.5189` (a different protein's Spearman, for the forward cell)
        `0.9`    matched inside `0.9217` (a different accuracy, for salmserovar)
    Measured exposure when this was found: 4 of 157 then-clean numbers rested on nothing else.

    A comparison that is strict on one side and loose on the other reports agreement it has not earned,
    so this is applied to BOTH the cited haystack AND the decoy haystacks. Holding a citation to a
    stricter bar than its own control would make the control measure a different question than the one
    it exists to control for.

    Implemented with `str.find` rather than the equivalent regex `(?<![\\d.])v(?![\\d])`, and that choice is
    MEASURED, not stylistic: profiling the decoy scan showed 7 regex searches consuming 0.814 of 1.031
    seconds -- 0.116 s each -- because a lookbehind forces a full scan of a multi-megabyte artifact. The
    find loop is C-level and does O(1) work per occurrence.
    """
    n = len(v)
    end = len(hay)
    start = 0
    while True:
        i = hay.find(v, start)
        if i < 0:
            return False
        # no digit or dot immediately before (so `2.1` does not match inside `72.1` or `0.2.1`),
        # no digit immediately after (so `0.51` does not match inside `0.5189`)
        before_ok = i == 0 or not (hay[i - 1].isdigit() or hay[i - 1] == ".")
        j = i + n
        after_ok = j >= end or not hay[j].isdigit()
        if before_ok and after_ok:
            return True
        start = i + 1


def _number_variants(tok: str) -> list[str]:
    """A cited 0.83 may live in the artifact as 0.8300, 0.83, or 83.0 (percent). Accept those forms."""
    v = {tok, tok.rstrip("0").rstrip(".") if "." in tok else tok}
    try:
        f = float(tok)
    except ValueError:
        return sorted(v)
    for nd in (1, 2, 3, 4):
        v.add(f"{f:.{nd}f}")
    if f < 1:                      # fraction cited, percent stored (or vice versa)
        v.add(f"{f * 100:.1f}")
        v.add(f"{f * 100:.0f}")
        v.add(f"{f * 100:g}")
    if f > 1:
        v.add(f"{f / 100:.4f}".rstrip("0"))
    return sorted(x for x in v if x)


def audit(decoy_pool: list[str] | None = None) -> dict:
    """`decoy_pool` PINS the decoy corpus to an explicit repo-relative file list.

    Default None keeps the live `wiki/*.json` glob, so `main()` and the reported rates continue to
    measure discrimination against the CURRENT corpus -- that is where a corpus effect belongs.
    The override exists for the anchor test, which pins EXACT match counts: against a live glob those
    counts move whenever any memo is added, and NON-MONOTONICALLY, because the pool is consumed in
    groups so new files RE-PARTITION it (grouped cells fell, ungrouped rose). A matcher-regression
    guard that fires on ordinary artifact accrual is measuring the wrong thing.
    """
    sys.path.insert(0, str(REPO))
    from dna_decode.data.cell_registry import cells

    rows = []
    # Cells carrying candidate numbers that cite NO artifact. `scope` below has always said the audit
    # covers "cells whose prose names a wiki/ artifact", so the exclusion is disclosed -- but its SIZE
    # was not, and that is what tells a reader whether the audited set is near-complete or a fraction of
    # the surface. It is a fraction: the skipped set is LARGER than the checked one. A cell citing only
    # its SCRIPT (typing:bacteria:mlst names scripts/mlst_serotype_purity.py, never the wiki artifact
    # holding the numbers) drops out entirely, so its measured purity figures are unverifiable here.
    # Reported, never counted as drift -- an uncited number is unchecked, not wrong.
    unverifiable = []
    fully_typed: list[dict] = []
    # Per-kind tally across EVERY extracted number, cited cells included. Asserted at the end to sum to
    # n_extracted, so a typing bug cannot silently drop a number out of the audit.
    kind_counts: dict[str, int] = {}
    n_extracted = 0
    # The pool EXCLUDES this script's own dated outputs, and that is not tidiness -- it is the difference
    # between a control and a circular one. `main()` writes wiki/contract_number_audit_<date>.json into the
    # very directory the pool is globbed from, and that artifact RECORDS THE NUMBER TOKENS it checked
    # (`cells_unverifiable[].numbers`, `matched_only_after_rounding`, `candidate_drift`). So every accrued
    # run is a decoy guaranteed to contain the numbers under test. Caught by measurement, not by reading:
    # `finder:any:forward` moved 74/300 -> 75/300 between two consecutive runs, the second run having
    # inherited the first run's artifact. Left in, the grades would also drift downward every run day.
    if decoy_pool is None:
        _decoy_pool = sorted(p.relative_to(REPO).as_posix()
                             for p in (REPO / "wiki").glob("*.json")
                             if not p.name.startswith("contract_number_audit_"))
    else:
        # A pinned pool must honour the SAME circularity exclusion, or a snapshot taken on a day the
        # audit had already run would smuggle its own output back in as a guaranteed-matching decoy.
        _decoy_pool = sorted(f for f in decoy_pool
                             if not f.rsplit("/", 1)[-1].startswith("contract_number_audit_"))
    _pool_set = set(_decoy_pool)
    # One shared cache across every cell: the full-pool scan is dominated by reading the corpus once
    # (~15 MB of json), not by the per-cell comparison, so sharing it is what makes the control affordable.
    _art_cache: dict = {}
    for c in cells():
        blob = " ".join(str(getattr(c, f, "") or "") for f in PROSE_FIELDS)
        arts = _cited_artifacts(blob)
        if not arts:
            tokens = [t for t in dict.fromkeys(DECIMAL_RE.findall(blob)) if not NON_MEASUREMENT.match(t)]
            kinds = {t: _classify_kind(c.cell_id, blob, t, matched_in_artifact=False, cited=False)
                     for t in tokens}
            for t, k in kinds.items():
                kind_counts[k] = kind_counts.get(k, 0) + 1
            n_extracted += len(tokens)
            bare = [t for t, k in kinds.items() if k == KIND_UNVERIFIABLE]
            typed = {t: k for t, k in kinds.items() if k != KIND_UNVERIFIABLE}
            if bare:
                unverifiable.append({"cell_id": c.cell_id, "n_numbers": len(bare),
                                     "numbers": bare, "why": "prose cites no wiki/ artifact",
                                     # what the cell's OTHER numbers turned out to be, so a shrinking
                                     # residual is auditable rather than merely smaller
                                     "kinds": kinds, "n_typed_not_unverifiable": len(typed)})
            elif typed:
                # every number named -> this cell leaves the residual, but it must not leave the REPORT,
                # or the audit would look like it had fewer numbers rather than better-described ones.
                fully_typed.append({"cell_id": c.cell_id, "kinds": kinds})
            continue
        texts = {a: _load_artifact_text(a) for a in arts}
        resolved = {a: t for a, t in texts.items() if t is not None}
        if not resolved:
            rows.append({"cell_id": c.cell_id, "status": "NO_ARTIFACT_RESOLVED",
                         "cited": arts, "checked": 0, "found": 0, "candidate_drift": []})
            continue
        haystack = "\n".join(resolved.values())
        nums = [t for t in dict.fromkeys(DECIMAL_RE.findall(blob)) if not NON_MEASUREMENT.match(t)]
        art_nums = _artifact_numbers(haystack)
        found, missing, rounded, in_code = [], [], {}, []
        for tok in nums:
            if any(_bounded_in(v, haystack) for v in _number_variants(tok)):
                found.append(tok)
                continue
            src = _matches_by_rounding(tok, art_nums)
            if src is not None:
                found.append(tok)
                rounded[tok] = src          # matched only after rounding; recorded for transparency
                continue
            # A cited number may be a THRESHOLD / BAR rather than a measurement, correctly absent from
            # every artifact. Those are handled by the SMALL explicit ADJUDICATED_BENIGN table above --
            # NOT by searching the package source, which was tried and is VACUOUS: variant matching of a
            # 2-4 digit token against a multi-megabyte source blob accepted 17/17 tokens including 12
            # randomly-generated ones, turning real signal into a clean-looking zero.
            # Any HUMAN-DECLARED kind means this number is legitimately absent from the artifact -- a bar,
            # an external reference value, a test-pinned figure, a deliberately-recorded rejected value, or
            # a quotient derived from counts the prose states. Each declaration carries a mandatory reason,
            # so this is an adjudication ON THE RECORD, not a silent exemption. Only a number with NO
            # declared kind and no artifact match is candidate drift.
            declared = _declared_kind(c.cell_id, tok)
            if declared:
                found.append(tok)
                in_code.append({"number": tok, "kind": declared,
                                "why": _declared_reason(c.cell_id, tok)})
            else:
                missing.append(tok)
        found_set = set(found)
        kinds = {t: _classify_kind(c.cell_id, blob, t,
                                   matched_in_artifact=t in found_set, cited=True) for t in nums}
        for k in kinds.values():
            kind_counts[k] = kind_counts.get(k, 0) + 1
        n_extracted += len(nums)
        # Group size = how many of THIS cell's resolved citations live in the decoy pool, so the decoy
        # haystack is the same size as the cited one. See _decoy_full_match_count.
        group_size = max(1, sum(1 for a in resolved if a in _pool_set))
        decoy_full, n_decoys = _decoy_full_match_count(nums, set(resolved), _decoy_pool, _art_cache,
                                                       group_size=group_size)
        rows.append({
            "cell_id": c.cell_id,
            "status": "CANDIDATE_DRIFT" if missing else "ALL_CITED_NUMBERS_PRESENT",
            "cited": sorted(resolved), "unresolved_citations": sorted(set(arts) - set(resolved)),
            "checked": len(nums), "found": len(found), "candidate_drift": missing,
            "number_kinds": kinds,
            "matched_only_after_rounding": rounded,
            "adjudicated_benign": in_code,
            # Does the citation do any work? See _decoy_full_match_count. `n_decoys` now carries the
            # number of full-pool trials (hundreds), not a sample size of 10; the field names are kept so
            # existing consumers and the cells_low_discrimination summary keep working.
            "decoy_full_match": decoy_full, "n_decoys": n_decoys,
            "decoy_match_rate": round(decoy_full / n_decoys, 4) if n_decoys else None,
            "decoy_group_size": group_size,
            "discrimination": _grade(decoy_full, n_decoys, nums),
        })

    checked = sum(r["checked"] for r in rows)
    drift = sum(len(r["candidate_drift"]) for r in rows)
    # Nothing may vanish through a typing bug. This is the invariant that makes a SHRINKING unverifiable
    # count trustworthy: the residual got smaller because numbers were named, not because they were lost.
    if sum(kind_counts.values()) != n_extracted:
        raise AssertionError(f"per-kind counts {sum(kind_counts.values())} != extracted {n_extracted}")
    return {
        "schema": "contract-number-audit-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "posture": "TRIAGE FUNNEL, not an oracle -- a cited number may legitimately be derived, rounded, or "
                   "describe a superseded run. Every hit needs adjudication; this script edits nothing.",
        "scope": "cells whose prose names a wiki/ artifact; decimal numbers only",
        "n_cells_audited": len(rows),
        "n_numbers_checked": checked,
        "n_candidate_drift": drift,
        "hit_rate": round(drift / checked, 4) if checked else None,
        "verdict": "NOTHING_TO_ADJUDICATE" if not drift else "ADJUDICATION_REQUIRED",
        # Coverage, so a clean verdict cannot be read as "every contract number is verified". These do
        # NOT enter n_numbers_checked / n_candidate_drift / verdict -- unchecked is not drift.
        "n_cells_low_discrimination": sum(1 for r in rows if r["discrimination"] == "LOW"),
        "cells_low_discrimination": sorted(r["cell_id"] for r in rows if r["discrimination"] == "LOW"),
        # Provenance kinds. "Unverifiable" used to mean four things at once, which made 34 numbers read
        # as 34 defects; it now means ONLY the residual that has no provenance of any kind.
        "n_numbers_extracted": n_extracted,
        "n_by_kind": dict(sorted(kind_counts.items())),
        "kind_definitions": {
            KIND_ARTIFACT: "present in a wiki/ artifact the prose cites",
            KIND_STRUCTURAL: "pattern-typed non-measurement: a DOI prefix, or a kb/bp length",
            KIND_EXTERNAL: "a published external reference value, not our measurement",
            KIND_ENFORCED: "pinned by an exact-equality assert in the test suite, not by an artifact",
            KIND_THRESHOLD: "a bar or cut-off, which has no artifact to live in",
            KIND_SUPERSEDED: "recorded deliberately as the REJECTED value, so it is never quoted",
            KIND_DERIVED: "computed from counts the prose itself states; the reason carries the arithmetic",
            KIND_DRIFT: "cited but absent from the cited artifact -- the thing this audit looks for",
            KIND_UNVERIFIABLE: "no provenance of any kind: an uncited measurement",
        },
        "n_cells_all_numbers_typed_no_citation": len(fully_typed),
        "cells_all_numbers_typed_no_citation": sorted(fully_typed, key=lambda u: u["cell_id"]),
        "n_cells_with_numbers_but_no_citation": len(unverifiable),
        "n_numbers_unverifiable": sum(u["n_numbers"] for u in unverifiable),
        "cells_unverifiable": sorted(unverifiable, key=lambda u: u["cell_id"]),
        "cells": rows,
    }


def main(argv=None) -> int:
    rep = audit()
    if rep["n_numbers_checked"] == 0:
        print("REFUSED: no cited numbers were checkable -- a clean report over zero checks is not a result.",
              file=sys.stderr)
        return 3
    stamp = rep["analysis_date"]
    (REPO / "wiki" / f"contract_number_audit_{stamp}.json").write_text(json.dumps(rep, indent=2),
                                                                      encoding="utf-8")
    for r in rep["cells"]:
        flag = "" if r["status"] == "ALL_CITED_NUMBERS_PRESENT" else f"  <- {r['status']}"
        disc = "" if r["discrimination"] in ("HIGH", "n/a") else \
               f"  [discrimination={r['discrimination']}: {r['decoy_full_match']}/{r['n_decoys']} " \
               f"({r['decoy_match_rate'] * 100:.0f}%) unrelated haystacks of " \
               f"{r['decoy_group_size']} also match all its numbers]"
        print(f"{r['cell_id']:42} checked={r['checked']:3} found={r['found']:3} "
              f"drift={len(r['candidate_drift'])}{flag}{disc}")
        if r["candidate_drift"]:
            print(f"    cited but absent from {', '.join(r['cited'])}: "
                  f"{', '.join(r['candidate_drift'])}")
    print(f"\n{rep['verdict']}  {rep['n_candidate_drift']}/{rep['n_numbers_checked']} numbers need "
          f"adjudication across {rep['n_cells_audited']} cells")
    print(f"\nPROVENANCE KINDS over all {rep['n_numbers_extracted']} extracted numbers:")
    for k, n in rep["n_by_kind"].items():
        print(f"    {k:28} {n:4}   {rep['kind_definitions'][k]}")
    if rep["cells_all_numbers_typed_no_citation"]:
        print("  every number named (no citation needed): "
              + ", ".join(u["cell_id"] for u in rep["cells_all_numbers_typed_no_citation"]))
    if rep["cells_unverifiable"]:
        print(f"COVERAGE: {rep['n_numbers_unverifiable']} number(s) in "
              f"{rep['n_cells_with_numbers_but_no_citation']} further cell(s) are UNVERIFIABLE here — "
              f"they cite no wiki/ artifact, so the clean verdict above does not cover them:")
        for u in rep["cells_unverifiable"]:
            print(f"    {u['cell_id']:42} {u['n_numbers']:3} number(s): "
                  f"{', '.join(u['numbers'][:8])}{' …' if u['n_numbers'] > 8 else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
