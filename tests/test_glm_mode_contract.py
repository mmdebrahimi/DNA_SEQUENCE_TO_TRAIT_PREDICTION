"""The (reference, conditioning) pairing contract — Step 1 of the sampling-sweep plan.

WHY (2026-10-19/20): a sweep ran `--natural-set promoter` with NO `--conditional`, which
scores an UNCONDITIONAL sample (whole-genome prior, mostly coding) against ALIGNED promoter
windows while the Markov null is fitted on the promoter FIT half. The null holds the target
distribution and the model does not, so the comparison is close to tautological — it produced
distinguishability 0.8141 against a null of 0.3876.

The repo already documented this in TWO places before the run (`corpus.py:168-175` and the
`--conditional` help) and nothing enforced it. These tests pin the enforcement.
"""
from __future__ import annotations

import importlib.util

import pytest

from dna_decode.glm.knob_contract import (
    VALID_COMPARISON,
    ComparisonModeError,
    validate_comparison_mode,
)

# ----------------------------------------------------------- the five combinations
def test_the_ONE_valid_pairing_passes():
    """promoter + conditional is the meaningful task: continue into the promoter."""
    validate_comparison_mode("promoter", True)


def test_promoter_without_conditional_is_REFUSED():
    with pytest.raises(ComparisonModeError, match="without --conditional is refused"):
        validate_comparison_mode("promoter", False)


def test_genomic_without_conditional_is_REFUSED_on_MEASURED_grounds():
    with pytest.raises(ComparisonModeError, match="MEASURED grounds"):
        validate_comparison_mode("genomic", False)


def test_genomic_with_conditional_is_REFUSED():
    with pytest.raises(ComparisonModeError, match="no shared offset structure"):
        validate_comparison_mode("genomic", True)


@pytest.mark.parametrize("conditional", [True, False])
def test_both_is_REFUSED_either_way(conditional):
    """The runner prompts from `kinds[0]` only, so the second set would be cross-referenced."""
    with pytest.raises(ComparisonModeError, match="kinds\\[0\\]"):
        validate_comparison_mode("both", conditional)


def test_an_unknown_reference_is_REFUSED_rather_than_defaulted():
    with pytest.raises(ComparisonModeError, match="unknown natural_set"):
        validate_comparison_mode("promotor", True)   # typo, not a silent pass


# ----------------------------------------------------------- non-vacuity of the messages
def test_the_promoter_refusal_TEACHES_rather_than_just_blocking():
    """A refusal that does not say why sends the reader to re-derive it.

    The message must carry the recorded defect, so the next person learns the reason
    instead of reaching for a flag to silence the guard.
    """
    with pytest.raises(ComparisonModeError) as e:
        validate_comparison_mode("promoter", False)
    msg = str(e.value)
    assert "corpus.py:168-175" in msg, "the refusal must cite the repo's own documentation"
    assert "0.8141" in msg and "0.3876" in msg, "the measured consequence must be named"
    assert "--conditional" in msg, "the fix must be stated"


def test_the_genomic_refusal_carries_the_BELOW_NOISE_FLOOR_measurement():
    """The genomic reference is refused on evidence, not taste: null 0.5086 vs 0.5255."""
    with pytest.raises(ComparisonModeError) as e:
        validate_comparison_mode("genomic", False)
    msg = str(e.value)
    assert "0.5086" in msg and "0.5255" in msg
    assert "no discriminating power" in msg


# ----------------------------------------------------------- the declared constant
# `VALID_COMPARISON` is exported, so a reader can take it as the contract -- which makes it
# a second statement of the same rule, free to drift from the function that enforces it.
_KNOWN_SETS = ("promoter", "genomic", "both")


def test_EXACTLY_ONE_pairing_passes_and_it_is_the_DECLARED_constant():
    """Pin the constant to behaviour by enumerating the whole documented cross-product.

    The five per-pairing tests above each check one cell; none of them would notice a
    SIXTH pairing being let through, nor `VALID_COMPARISON` being edited on its own. Set
    equality catches both directions: relaxing a refusal, and renaming the survivor.
    """
    passed = [
        (s, c)
        for s in _KNOWN_SETS
        for c in (True, False)
        if _accepts(s, c)
    ]
    assert passed == [VALID_COMPARISON], passed


def _accepts(natural_set, conditional) -> bool:
    try:
        validate_comparison_mode(natural_set, conditional)
    except ComparisonModeError:
        return False
    return True


# ----------------------------------------------------------- ordering in the real runner
def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "_glm_falsifier", str(__import__("pathlib").Path(__file__).resolve().parent.parent
                             / "scripts" / "glm_generator_falsifier.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_runner_REFUSES_before_building_a_corpus_or_loading_a_model():
    """Ordering is the whole point: a mis-specified run must cost milliseconds.

    If the guard sat after the corpus build this would raise something else (or spend
    ~30 min generating), so the exception TYPE is the ordering evidence.
    """
    mod = _load_runner()
    with pytest.raises(ComparisonModeError):
        mod.main(["--natural-set", "promoter"])          # no --conditional


def test_the_runner_also_refuses_both():
    mod = _load_runner()
    with pytest.raises(ComparisonModeError):
        mod.main(["--natural-set", "both", "--conditional"])


# ----------------------------------------------------------- Step 2: the invalidation memo
# The plan specified this test strategy but listed no test file for Step 2, so it lives here
# beside the Wave-0 sibling rather than in a file of its own. Deviation recorded in the
# execution log.
def test_the_invalidation_memo_names_BOTH_defects_and_cites_the_committed_baseline():
    """A memo that asserts the conclusion without its evidence is the thing to avoid.

    Both defects are independently sufficient to void the run, so both must be named, and
    the replacement baseline must be a real committed artifact rather than a recollection.
    """
    import json
    import pathlib

    root = pathlib.Path(__file__).resolve().parent.parent
    memo = root / "wiki" / "glm_sweep_invalidation_2026-10-20.md"
    assert memo.is_file(), "Step 2 memo missing"
    t = memo.read_text(encoding="utf-8")

    # defect 1: the inert knob, with its measured signature
    assert "INERT" in t and "6/6 byte-identical" in t
    # defect 2: the mode mismatch, with the repo's own prior documentation
    assert "conditional: False" in t and "corpus.py:168-175" in t
    # the replacement baseline must be named AND must actually exist
    cited = "wiki/glm_generator_falsifier_REAL_2026-10-07.json"
    assert cited in t
    art = root / cited
    assert art.is_file(), f"memo cites a missing artifact: {cited}"
    d = json.loads(art.read_text(encoding="utf-8"))
    assert d.get("conditional") is True, "the cited baseline must be the CONDITIONAL run"
    # and the memo must not present the recorded premise as a discovery
    assert "orientation error #3" in t


# ------------------------------------------- ordering, proven and reference-independent
# The two tests above take the exception TYPE as the ordering evidence, which leaves two
# gaps. (1) They need the real MG1655 reference: `data/` is gitignored and that reference is
# untracked, so on a fresh clone -- or on this host the morning after the worktree lost its
# gitignored files -- `main()` returns 2 at the reference check and never reaches the guard,
# and the tests FAIL rather than skip. (2) Nothing shows the guard precedes the corpus build
# specifically; the same error would surface from anywhere in `main()`.
class _CorpusReached(Exception):
    """Raised from the first expensive step the guard must come before."""


def _runner_with_stub_reference(tmp_path, monkeypatch):
    """A runner whose reference gate passes on a stub and whose corpus build is a tripwire."""
    mod = _load_runner()
    ref = tmp_path / "GCF_000005845.2"
    ref.mkdir()
    (ref / "genome.fna").write_text(">stub\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(mod, "REFSEQ", ref)

    def _tripwire(*_a, **_k):
        raise _CorpusReached("build_natural ran")

    monkeypatch.setattr(mod, "build_natural", _tripwire)
    return mod


def test_the_guard_fires_BEFORE_the_corpus_build_without_needing_the_reference(
    tmp_path, monkeypatch
):
    mod = _runner_with_stub_reference(tmp_path, monkeypatch)
    with pytest.raises(ComparisonModeError):
        mod.main(["--natural-set", "promoter"])          # no --conditional


def test_and_that_ordering_proof_is_NOT_vacuous(tmp_path, monkeypatch):
    """The companion that makes the test above mean something.

    If the tripwire were unreachable -- wrong stub path, an earlier return, a renamed
    helper -- the test above would pass with no guard in the file at all. This asserts the
    VALID pairing DOES reach the corpus build, so the guard is the only difference between
    the two runs.
    """
    mod = _runner_with_stub_reference(tmp_path, monkeypatch)
    with pytest.raises(_CorpusReached):
        mod.main(["--natural-set", "promoter", "--conditional"])
