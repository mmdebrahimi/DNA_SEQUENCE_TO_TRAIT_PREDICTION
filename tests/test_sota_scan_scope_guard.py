"""A zero-shot likelihood result must never be presented as replicating the SUPERVISED encoder gate.

This repo has already paid for this exact conflation. CLAUDE.md records it verbatim: *"HEADLINE
RE-SCOPED 2026-09-30, and the unscoped version cost three months: it read 'A learned variant scorer
does NOT fill...', which a skimming session takes as closing the LEARNED direction. The body below is
entirely about ZERO-SHOT likelihood."* And the repo holds a measured case where the two readouts give
OPPOSITE verdicts on the same data -- supervised 0.8142 vs zero-shot 0.4750, a PASS where zero-shot
FAILED.

On 2026-10-22 it happened again, in `wiki/ai_genomics_sota_scan_2026-10-22.md`: two external 2026
benchmarks that score by LOG-LIKELIHOOD (the Mechanistic Invariance Test's CSS/SCR are LL-derived;
Nullsettes defines a detection as a mutant LL significantly below its nonmutant original) were
reported as "reproducing" the encoder gate -- which is a FITTED conv encoder judged by per-seed
deltas and a paired CI, and whose own artifact contains the word "likelihood" zero times.

WHAT IS GUARDED IS A RELATIONSHIP, NOT A NUMBER. A value-only check cannot catch this class: every
figure involved is correctly transcribed. The defect is the verb connecting them. So this scans for
the *claim shape* -- an external zero-shot benchmark presented as reproducing or replicating our
result -- the same way `tests/test_essentiality_overclaim.py` guards the 0.911 prose relationship
rather than the numeral.

Scope: `wiki/*.md`. Prose guards guard prose; the JSON sidecar carries its own structural fields
(`not_a_replication_of_supervised_gate`, per-result `readout`) which are asserted at write time.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WIKI = ROOT / "wiki"

#: The artifact whose own result is the supervised side of the comparison. Its premise is pinned
#: below rather than trusted.
GATE_ARTIFACT = WIKI / "glm_encoder_gate_result_2026-10-08.md"

#: Names of the external ZERO-SHOT benchmarks. Any of these in a segment makes it in-scope.
_EXTERNAL = re.compile(
    r"mechanistic invariance|nullsettes|2604\.06549|2506\.10271", re.I
)
#: Ways a segment can refer to OUR supervised result.
_OURS = re.compile(r"encoder[- ]gate|\+?0\.0122|GLM-arm negative", re.I)
#: The verbs that assert identity of result rather than agreement of direction.
#:
#: A NEGATED form is not a violation and must not be flagged -- "not a replication of it" is exactly
#: the corrected phrasing this guard exists to permit, and "un-replicated studies" is a true remark
#: about the external work's own status. The first version of this guard fired on both, which is how
#: the negation handling got written: a guard that cannot express the correct sentence is useless.
_CLAIMS_REPLICATION = re.compile(
    r"(?<!un-)(?<!un)(?<!not )(?<!not a )(?<!never )(?<!no )\b(?:reproduc|replicat)\w*", re.I
)
_NEGATED_REPLICATION = re.compile(
    r"(?:un-?|not\s+(?:a\s+)?|never\s+|no\s+)(?:reproduc|replicat)\w*", re.I
)
#: Any of these in the same segment shows the readout distinction is being made.
_READOUT_MARKER = re.compile(r"zero[- ]shot|log[- ]?likelihood|likelihood probe|frozen", re.I)
#: The specific compression that flattened three different experiments into one claim.
_COMPRESSION = re.compile(r"three teams,? three benchmarks,? one mechanism", re.I)


def _segments(text: str) -> list[str]:
    """The unit a human actually reads, so the unit to judge.

    Blank-line separated, then split again on list bullets. Without the second split a whole bullet
    list is one segment, so a correct remark in one bullet silently satisfies -- or trips -- a rule
    about a different bullet. That is not hypothetical: it produced this guard's first false
    positive, where a bullet calling the external studies "un-replicated" sat in the same block as a
    bullet naming the encoder gate.
    """
    out: list[str] = []
    for block in re.split(r"\n\s*\n", text):
        parts = re.split(r"\n(?=\s*[-*]\s)", block)
        out.extend(p for p in parts if p.strip())
    return out


def _claims_replication(seg: str) -> bool:
    """A replication assertion that is not negated anywhere in the segment."""
    if not _CLAIMS_REPLICATION.search(seg):
        return False
    # Count only matches that a negation does not already cover.
    hits = len(_CLAIMS_REPLICATION.findall(seg))
    negs = len(_NEGATED_REPLICATION.findall(seg))
    return hits > negs


def violations(text: str) -> list[str]:
    """Every way a passage can present a zero-shot result as replicating the supervised gate."""
    bad: list[str] = []
    for seg in _segments(text):
        names_external = bool(_EXTERNAL.search(seg))
        names_ours = bool(_OURS.search(seg))
        claims_repl = _claims_replication(seg)
        has_marker = bool(_READOUT_MARKER.search(seg))

        # Rule 1: a replication claim tying the external work to our result, unqualified.
        if claims_repl and names_ours and not has_marker:
            bad.append("a reproduce/replicate claim attaches to our result with no readout marker")
        # Rule 2: both sides named in one segment with no readout distinction drawn.
        if names_external and names_ours and not has_marker:
            bad.append("external benchmark and our result in one segment with no readout marker")
        # Rule 3: the compression itself.
        if _COMPRESSION.search(seg) and not has_marker:
            bad.append("'three teams, three benchmarks, one mechanism' with no readout marker")
    return bad


# --- the live surfaces --------------------------------------------------------------------------

def _wiki_md_citing_both() -> list[Path]:
    out = []
    for p in sorted(WIKI.glob("*.md")):
        try:
            t = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if _EXTERNAL.search(t) and _OURS.search(t):
            out.append(p)
    return out


def test_at_least_one_wiki_artifact_is_actually_in_scope():
    """If the glob matched nothing this whole module would pass vacuously forever."""
    found = _wiki_md_citing_both()
    assert found, "no wiki/*.md cites both an external benchmark and our gate -- guard is vacuous"


@pytest.mark.parametrize("path", _wiki_md_citing_both() or [GATE_ARTIFACT], ids=lambda p: p.name)
def test_no_wiki_artifact_presents_zero_shot_as_replicating_the_supervised_gate(path: Path):
    v = violations(path.read_text(encoding="utf-8", errors="replace"))
    assert v == [], f"{path.name}: {sorted(set(v))}"


# --- non-vacuity, against the text that actually shipped ----------------------------------------

#: Verbatim from the committed pre-fix artifact. Recovered with
#: `git show cb1b2a05:wiki/ai_genomics_sota_scan_2026-10-22.md` -- kept ONLY to prove this guard
#: fails on the real thing. A guard that flagged nothing would pass forever.
PRE_FIX_STRINGS = (
    "## F1 — Two independent 2026 results reproduce this project's encoder-gate negative",

    "Stated plainly because the opposite would matter more: this scan found **no external result "
    "that\noverturns** the GLM-arm negatives. It found two that independently reproduce them, which "
    "is a different\nand more useful thing (see F1).",

    "**Why these three agree.** All three say the same thing by different routes: **the signal these\n"
    "representations carry on regulatory sequence is compositional (GC/AT), and position/mechanism is\n"
    "largely absent.** Ours measured it as \"a conv encoder does not beat one GC number\"; the MIT "
    "measured it\nas \"apparent sensitivity is AT-content correlation\"; Nullsettes measured it as "
    "\"detection tracks\nlikelihood, not function.\" Three teams, three benchmarks, one mechanism — "
    "and it is the mechanism that\nmakes a GC baseline so hard to beat.",
)


@pytest.mark.parametrize("pre", PRE_FIX_STRINGS, ids=("f1_heading", "c2_corollary", "why_they_agree"))
def test_the_guard_is_non_vacuous_it_rejects_the_pre_fix_text(pre: str):
    assert violations(pre), f"guard failed to flag shipped over-claim: {pre[:70]!r}"


def test_the_corrected_text_passes_so_the_guard_is_not_merely_strict():
    """Non-vacuity cuts both ways: a guard nothing can satisfy is also useless."""
    ok = (
        "F1 — Two independent 2026 results point the same way as our encoder-gate negative, on a "
        "different readout. Both are ZERO-SHOT LOG-LIKELIHOOD probes; the gate is a fitted model."
    )
    assert violations(ok) == [], violations(ok)


# --- the premise, pinned rather than asserted from memory ---------------------------------------

def test_the_gate_artifact_is_not_a_likelihood_probe():
    """The whole rationale rests on this. If a future rewrite made the gate a likelihood probe the
    guard would be guarding nothing, so the rationale is checked rather than trusted."""
    if not GATE_ARTIFACT.exists():
        pytest.skip(f"{GATE_ARTIFACT.name} absent")
    t = GATE_ARTIFACT.read_text(encoding="utf-8", errors="replace").lower()
    assert t.count("likelihood") == 0, (
        "the encoder-gate artifact now mentions likelihood; re-read it before trusting this guard"
    )


def test_pre_fix_strings_are_really_in_the_pre_fix_commit():
    """Pins that PRE_FIX_STRINGS is history, not invention -- the strings must exist at cb1b2a05."""
    try:
        blob = subprocess.run(
            ["git", "show", "cb1b2a05:wiki/ai_genomics_sota_scan_2026-10-22.md"],
            cwd=ROOT, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=60,
        )
    except (OSError, subprocess.SubprocessError):  # pragma: no cover
        pytest.skip("git unavailable")
    if blob.returncode != 0:
        pytest.skip("commit cb1b2a05 not reachable in this checkout")
    text = blob.stdout
    missing = [p[:50] for p in PRE_FIX_STRINGS if p.split("\n")[0][:40] not in text]
    assert not missing, f"PRE_FIX_STRINGS not found at cb1b2a05: {missing}"
