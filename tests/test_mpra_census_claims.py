"""The census artifact's OWN claims, which nothing checked.

`tests/test_substrate_symmetry.py` guards two things about the live census: that it is symmetric, and
that every rung's utilization is knowable. Neither touches the claims the file makes about the rest
of the repo -- and that file makes several load-bearing ones:

  - `corrections_2026_10_13[].enforced_by` names a module AND a test file. Rename the test file and
    the artifact silently asserts an enforcement that no longer exists -- the "retired instrument
    leaves live claims" shape.
  - `rungs.sigma70_designed.canonical_loader` names ONE loader as canonical, recorded immediately
    after a second loader was written beside it and reproduced the exact naive join the canonical one
    documents, fixes and test-pins. A stale symbol name here points the next session at the wrong
    file, which is how that incident happened.
  - `ceiling_rule` states in prose that ceilings are PER-RUNG and must not be inherited, and that
    `single_replicate_agreement` must never be quoted as one. `substrate_symmetry.py`'s own docstring
    says a prose rule enforced by nothing is exactly what got violated -- so the guard's non-vacuity
    AGAINST THE REAL FILE is what makes the prose a rule rather than a hope.

Every test skips (never silently passes) when the artifact is absent.
"""
from __future__ import annotations

import importlib
import json
import re
from pathlib import Path

import pytest

from dna_decode.eval.substrate_symmetry import audit, ceiling_guard

ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT / "wiki" / "mpra_substrate_census_2026-10-12.json"

pytestmark = pytest.mark.skipif(
    not CENSUS.exists(), reason="census artifact absent (gitignored / fresh clone)")

#: Repo-relative file citations: `wiki/x.json`, `tests/x.py`, `dna_decode/x/y.py`, `scripts/x.py`.
_PATH = re.compile(r"(?:wiki|tests|scripts|dna_decode)/[A-Za-z0-9_./-]+\.(?:py|json|md|txt)")
#: Dotted module-or-symbol references, e.g. `dna_decode.glm.expression.load_pairs`.
_DOTTED = re.compile(r"dna_decode(?:\.[a-z_][a-z0-9_]*)+")


def _census() -> dict:
    return json.loads(CENSUS.read_text(encoding="utf-8"))


def _strings(o):
    """Every string anywhere in the artifact -- claims live in prose fields, not just in a schema."""
    if isinstance(o, dict):
        for v in o.values():
            yield from _strings(v)
    elif isinstance(o, list):
        for v in o:
            yield from _strings(v)
    elif isinstance(o, str):
        yield o


def _resolve_dotted(ref: str):
    """Import the longest importable module prefix, then walk the remainder as attributes."""
    parts = ref.split(".")
    for cut in range(len(parts), 0, -1):
        try:
            mod = importlib.import_module(".".join(parts[:cut]))
        except ImportError:
            continue
        obj = mod
        for attr in parts[cut:]:
            obj = getattr(obj, attr)  # AttributeError here IS the failure we want surfaced
        return obj
    raise ImportError(ref)


def test_every_file_the_census_cites_actually_exists():
    """A citation is provenance only if it resolves. Scanned from the artifact rather than
    hand-listed, so a citation added to a new rung is covered without editing this test."""
    cited = sorted({m for s in _strings(_census()) for m in _PATH.findall(s)})
    assert cited, "the scan found no citations at all -- the regex has drifted from the artifact"
    missing = [c for c in cited if not (ROOT / c).exists()]
    assert not missing, f"the census cites files that do not exist: {missing}"


def test_every_loader_and_helper_the_census_names_is_importable():
    """NON-VACUITY MATTERS HERE: the census's own correction records that a duplicate loader was
    written beside `dna_decode/glm/expression.py` and reproduced the naive join that file fixes. The
    census answers "which one is canonical", so a dotted name that no longer resolves re-opens
    exactly that ambiguity. Resolved through the attribute, not just the module, because
    `load_pairs`/`leave_element_out` being gone from a module that still imports is the live risk."""
    refs = sorted({m for s in _strings(_census()) for m in _DOTTED.findall(s)})
    assert len(refs) >= 3, f"expected the loader claims to be present, found {refs}"
    for ref in refs:
        try:
            assert _resolve_dotted(ref) is not None, ref
        except (ImportError, AttributeError) as exc:
            pytest.fail(f"the census names {ref} but it does not resolve: {exc!r}")


def test_the_enforcement_claim_names_a_test_file_that_still_exercises_the_check():
    """A path that resolves is not an enforcement. The named test file must actually reference
    `check_symmetry`, or the claim is false while every path in it is valid."""
    cors = [c for c in _census()["corrections_2026_10_13"] if "enforced_by" in c]
    assert cors, "the symmetry correction no longer declares what enforces it"
    for c in cors:
        cited = _PATH.findall(c["enforced_by"])
        assert cited, f"enforced_by names nothing checkable: {c['enforced_by']}"
        for path in cited:
            p = ROOT / path
            assert p.exists(), f"enforced_by cites a missing {path}"
            if p.suffix == ".py":
                assert "check_symmetry" in p.read_text(encoding="utf-8"), (
                    f"{path} is cited as the enforcement but never calls check_symmetry")


def test_the_ceiling_guard_FIRES_on_the_real_census():
    """Non-vacuity against real data. Three live rungs carry both a single-replicate agreement and a
    stepped-up ceiling, and the guard must warn on every one -- a guard that fires on synthetic input
    but is silent on the shipped file is the documented 'a filter that removes nothing is not a
    control' shape."""
    census = _census()
    expect = sorted(n for n, r in census["rungs"].items()
                    if isinstance(r.get("single_replicate_agreement"), (int, float))
                    and isinstance(r.get("ceiling_for_a_sequence_model"), (int, float))
                    and r["single_replicate_agreement"] != r["ceiling_for_a_sequence_model"])
    assert expect, "no live rung carries both numbers -- the guard has nothing to prove here"
    warnings = audit(census)["ceiling_warnings"]
    assert len(warnings) == len(expect)
    for rung in expect:
        hit = [w for w in warnings if w.startswith(f"{rung}:")]
        assert hit, f"the guard is silent on {rung}, which carries both numbers"
        sra = census["rungs"][rung]["single_replicate_agreement"]
        assert str(sra) in hit[0], "the warning must quote the number not to use"


def test_no_rung_INHERITS_another_rungs_ceiling():
    """The other half of the census's `ceiling_rule`, previously prose-only: ceilings are PER-RUNG and
    0.9410 belongs to the peak tiles alone. Two rungs sharing one non-null ceiling is the signature of
    an inherited value, and the recorded incident (a fragment rung's headroom understated by ~0.33)
    came from reading one rung's number on another."""
    seen: dict[float, str] = {}
    for name, rung in sorted(_census()["rungs"].items()):
        c = rung.get("ceiling_for_a_sequence_model")
        if not isinstance(c, (int, float)):
            continue  # an explicit null with a stated reason is a statement, not an inheritance
        assert c not in seen, f"{name} carries the same ceiling {c} as {seen[c]} -- inherited?"
        seen[c] = name
    assert len(seen) >= 3, "too few numeric ceilings for this guard to discriminate"


def test_a_rung_declaring_no_ceiling_also_declares_WHY():
    """`None` must be a statement, not a gap. The sigma70 rung has no computable ceiling (its two
    published tables are one monotone transform apart), and that reason is the difference between
    'not computable' and 'nobody filled this in' -- the absence-reads-as-nothing-to-say class."""
    for name, rung in _census()["rungs"].items():
        if rung.get("ceiling_for_a_sequence_model") is None:
            assert rung.get("ceiling_reason"), f"{name} declares a null ceiling with no reason"
            assert ceiling_guard(name, rung) == [], f"{name}: guard must be silent on a null pair"
