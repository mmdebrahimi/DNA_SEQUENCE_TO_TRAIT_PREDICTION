"""A field present on SOME rungs and absent on others reads as a global statement. Refuse that.

**The measured failure this exists to stop (2026-10-13).** `wiki/mpra_substrate_census_2026-10-12.json`
describes four substrate rungs. The field `size_lever` is present on **exactly one** of them --
`peak_tiles`, where it reads `"EXHAUSTED -- 46,713 is the whole file"`. The three rungs that actually
have headroom carry no such field at all. A reader looking for a size statement finds exactly one, and
it says exhausted; so a rung-scoped fact was quoted as a project-level wall, and 321,123 unused
fragments in the same cache directory were reported as "no rows available".

**An absent field is not a neutral absence.** It reads as "nothing to say here", which is a silent
default -- the same class as a filter that removes nothing being reported as an active control. The fix
is not a better reader: it is refusing a census whose fields are asymmetric across its own rungs.

**What was tried and rejected, so it is not re-proposed.** Per-file orphan detection in
`scripts/data_inventory.py` would NOT have caught this. Its `Dataset.orphan` is sound, but the fragment
files are *referenced in code* (`dna_decode/glm/genomewide.py:48-49`) and so are not orphans by any
definition. The signal is not "unreferenced"; it is "under-used relative to its size" -- 30,000 of
321,123 rows, 9.3%.

**Ownership (OT1).** Code owns the symmetry verdict, the utilization banding and the ceiling refusal --
pure, deterministic, no I/O, recomputed from a dict the caller supplies. The MODEL/author owns what the
rungs are and what their fields say. So: *"a code-owned symmetry verdict over an author-declared
census"* -- never "a census cannot be asymmetric". Nothing here blocks a write; it reports.
"""
from __future__ import annotations

from dataclasses import dataclass, field

#: Fields that MUST be present on every rung or absent from every rung. A hand-declared set, and that
#: is a named residual: a NEW field added to one rung without being listed here is invisible to the
#: check (declared, not detected). `test_symmetry_required_cannot_be_emptied` pins that the set is
#: non-empty so the rule cannot be silently disabled.
SYMMETRY_REQUIRED: tuple[str, ...] = (
    "size_lever",
    "best_measured_arm",
    "ceiling_for_a_sequence_model",
)

#: >= this fraction of available rows actually used -> the size lever really is spent on that rung.
EXHAUSTED_AT = 0.95
#: < this fraction -> loudly under-used. Between the two -> `partial`.
UNDER_USED_BELOW = 0.50


@dataclass(frozen=True)
class SymmetryReport:
    verdict: str
    detail: str
    asymmetric_fields: dict[str, list[str]] = field(default_factory=dict)
    n_rungs: int = 0

    @property
    def ok(self) -> bool:
        return self.verdict == "ok"


@dataclass(frozen=True)
class UtilizationVerdict:
    rung: str
    verdict: str
    rows_available: int | None = None
    largest_n_measured: int | None = None
    fraction: float | None = None

    @property
    def exhausted(self) -> bool:
        """ONLY an explicitly measured >=EXHAUSTED_AT counts. `unknown` must never satisfy this --
        an unmeasured rung is not evidence that its data is spent, and treating it as such is the one
        failure direction this module exists to prevent."""
        return self.verdict == "exhausted"


def _rungs(census: dict) -> dict[str, dict]:
    r = census.get("rungs")
    return r if isinstance(r, dict) else {}


def check_symmetry(census: dict, required: tuple[str, ...] = SYMMETRY_REQUIRED) -> SymmetryReport:
    """Ordered ladder, first match wins. `indeterminate` is never silently a pass.

    `indeterminate:single-rung` is deliberate: symmetry is UNDEFINED over one rung, and returning `ok`
    there would mean a census grows its second rung asymmetrically and the check that was "passing"
    the whole time starts failing for reasons the author cannot connect to their edit.
    """
    rungs = _rungs(census)
    if not rungs:
        return SymmetryReport("indeterminate:no-rungs", "census carries no 'rungs' mapping", {}, 0)

    # A rung whose value is not a mapping used to be skipped by every per-field comprehension below,
    # so the census returned `ok` AND the rung vanished from `audit()`'s table -- this module's own
    # stated failure class ("a rung absent reads as nothing to say here") reproduced inside the
    # validator built to prevent it. `n_rungs` counted it while nothing scored it. Named ahead of the
    # ladder because an unscoreable rung makes every downstream verdict a statement about a subset.
    malformed = sorted(n for n, r in rungs.items() if not isinstance(r, dict))
    if malformed:
        return SymmetryReport(
            f"malformed:{malformed[0]}",
            f"rung(s) {malformed} are not mappings and cannot be scored; every verdict below would "
            f"describe only the {len(rungs) - len(malformed)} scoreable of {len(rungs)} rungs",
            {}, len(rungs))

    if len(rungs) == 1:
        return SymmetryReport(
            "indeterminate:single-rung",
            f"symmetry is undefined over one rung ({next(iter(rungs))}); refusing rather than passing",
            {}, 1)

    asym: dict[str, list[str]] = {}
    for f in required:
        present = [n for n, r in rungs.items() if isinstance(r, dict) and f in r]
        missing = [n for n, r in rungs.items() if isinstance(r, dict) and f not in r]
        if present and missing:
            asym[f] = sorted(missing)

    if asym:
        worst = sorted(asym)[0]
        parts = "; ".join(f"{f} present on {len(rungs) - len(m)}/{len(rungs)}, missing from {m}"
                          for f, m in sorted(asym.items()))
        return SymmetryReport(
            f"asymmetric:{worst}",
            f"a field on SOME rungs reads as a statement about ALL of them -- {parts}",
            asym, len(rungs))
    return SymmetryReport("ok", f"all {len(required)} required fields are symmetric across "
                                f"{len(rungs)} rungs", {}, len(rungs))


#: Row-count keys in PRECEDENCE order: post-filter count first (what a model actually sees), then the
#: raw file count. **The census is NOT uniform here, and a short list is a silent `unknown`.** Measured
#: 2026-10-13: a first version listed only the tile and fragment keys, so `sigma70_designed` reported
#: `unknown:no row count` while recording 10,898 rows under THREE different names
#: (`designed_rows`, `usable_joined`, `n_distinct_sequences`). The missing key did not read as "the
#: list is incomplete" -- it read as "that rung has no data", which is the exact error class this
#: module exists to prevent, reproduced inside the fix for it.
ROW_COUNT_KEYS: tuple[str, ...] = (
    "rows_after_published_filters",   # peak tiles: post-filter
    "usable_joined",                  # sigma70 designed: post-join
    "n",                              # fragments
    "designed_rows",                  # sigma70 designed: as designed
    "n_distinct_sequences",
    "rows_in_file",                   # raw, last resort
)


def _is_count(v: object) -> bool:
    """A positive integer row count, and **`bool` is REFUSED**.

    `bool` is an `int` subclass in Python, so `isinstance(True, int)` is True and a census carrying
    `{"n": true}` read as ONE available row -- which then satisfied `>= EXHAUSTED_AT` and reported
    **`exhausted`** off a type error. That is the worst fail direction this module has: a malformed
    value claiming the data is spent. `type(v) is int` excludes bool without excluding anything real.
    """
    return type(v) is int and v > 0


def _rows_available(rung: dict) -> int | None:
    """Row count under whichever key this rung uses, by `ROW_COUNT_KEYS` precedence."""
    for k in ROW_COUNT_KEYS:
        v = rung.get(k)
        if _is_count(v):
            return v
    return None


def utilization(name: str, rung: dict) -> UtilizationVerdict:
    """How much of this rung's data any measurement has actually used."""
    avail = _rows_available(rung)
    used = rung.get("largest_n_measured")
    if avail is None or not _is_count(used):
        why = "no row count" if avail is None else "no largest_n_measured"
        return UtilizationVerdict(name, f"unknown:{why}", avail,
                                  used if _is_count(used) else None, None)
    if used > avail:
        # SURFACED, not silently rounded to `exhausted`. Using more rows than exist is a census
        # arithmetic error, and its most plausible cause is a stale `largest_n_measured` left behind
        # when a row count is corrected DOWNWARD -- which is exactly the edit this census just
        # underwent. Reporting `exhausted` here would hide the error behind the reassuring answer.
        return UtilizationVerdict(name, f"inconsistent:used>{avail}", avail, used,
                                  round(used / avail, 4))
    frac = used / avail
    if frac >= EXHAUSTED_AT:
        v = "exhausted"
    elif frac < UNDER_USED_BELOW:
        v = f"under-used:{frac * 100:.1f}%"
    else:
        v = f"partial:{frac * 100:.1f}%"
    return UtilizationVerdict(name, v, avail, used, round(frac, 4))


def ceiling_guard(name: str, rung: dict) -> list[str]:
    """Refuse single_replicate_agreement presented as a ceiling.

    The census already states this rule in prose (`ceiling_rule`), and a prose rule enforced by nothing
    is exactly what got violated this session: the fragment rung's 0.4583 single-replicate agreement
    was quoted as its ceiling, understating real headroom by +0.3345 against the true 0.7928. When the
    readout is the MEAN of k replicates its reliability is higher than one replicate's, so the ceiling
    is the Spearman-Brown step-up square-rooted -- never the raw agreement.
    """
    out: list[str] = []
    sra = rung.get("single_replicate_agreement")
    ceil = rung.get("ceiling_for_a_sequence_model")
    if isinstance(sra, (int, float)) and isinstance(ceil, (int, float)) and abs(ceil - sra) > 1e-9:
        out.append(
            f"{name}: do NOT quote single_replicate_agreement {sra} as a ceiling; this rung's "
            f"ceiling_for_a_sequence_model is {ceil} (Spearman-Brown step-up, square-rooted)")
    return out


def audit(census: dict) -> dict:
    """The whole read-only verdict set, shaped for a renderer or an artifact."""
    rungs = _rungs(census)
    sym = check_symmetry(census)
    # A non-dict rung gets a ROW in the table rather than being filtered out of it. Omitting it was
    # the defect: the renderer showed a clean 3-rung table over a 4-rung census and nothing said so.
    utils = []
    for n, r in sorted(rungs.items()):
        utils.append(utilization(n, r) if isinstance(r, dict)
                     else UtilizationVerdict(n, f"malformed:not-a-mapping({type(r).__name__})"))
    guards: list[str] = []
    for n, r in sorted(rungs.items()):
        if isinstance(r, dict):
            guards.extend(ceiling_guard(n, r))
    return {
        "symmetry": {"verdict": sym.verdict, "detail": sym.detail,
                     "asymmetric_fields": sym.asymmetric_fields, "n_rungs": sym.n_rungs},
        "utilization": [{"rung": u.rung, "verdict": u.verdict, "rows_available": u.rows_available,
                         "largest_n_measured": u.largest_n_measured, "fraction": u.fraction}
                        for u in utils],
        "ceiling_warnings": guards,
        "required_fields": list(SYMMETRY_REQUIRED),
        "thresholds": {"exhausted_at": EXHAUSTED_AT, "under_used_below": UNDER_USED_BELOW,
                       "are": "ASSERTED, not derived -- the only measured anchors are 9.3% "
                              "(fragments, clearly under-used) and ~94% (tiles, clearly spent); any "
                              "cut between those separates them equally well"},
    }
