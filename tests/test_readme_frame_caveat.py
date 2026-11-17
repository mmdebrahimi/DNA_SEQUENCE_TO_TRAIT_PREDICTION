"""Guards for the README's scoring-frame passage.

Two of these are verdict-independent and land whatever the gate returned.

The README used to attribute the unmeasured human check to "detached volume + an uninstalled
dependency" when `openpyxl` 3.1.5 was installed and `scripts/fetch_lentimpra.py` takes
`--dest` -- a false statement on the most-read surface this work touches, and an under-claim,
which this project treats as the same class of trust-surface falsehood as an over-claim.

And the 20-25x headline sat 13 lines above its own caveat, so a reader who stopped at the
headline got the uncaveated inference. Adjacency is now pinned per paragraph.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

README = Path("README.md")
RESOLUTION = "wiki/glm_scoring_frame_resolution.md"

CAVEAT_WORDS = ("frame", "confound", "caveat", "unresolved", "partially", "indeterminate")


@pytest.fixture(scope="module")
def text():
    return README.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def paragraphs(text):
    return [p for p in text.split("\n\n") if p.strip()]


# ----------------------------------------------------------------------------------------
# verdict-independent 1: the false dependency attribution
# ----------------------------------------------------------------------------------------

def test_the_readme_never_blames_an_uninstalled_dependency_for_the_human_check(paragraphs):
    """openpyxl 3.1.5 is installed and the fetcher takes --dest; the check was merely unrun.

    A plain absence check CANNOT work here and the first version of this guard failed for that
    reason: the correction has to QUOTE the false claim in order to retract it. So the test is
    that the phrase never appears as an assertion -- wherever it occurs, the same paragraph
    must also say it was false.
    """
    offenders = []
    for p in paragraphs:
        if "uninstalled dependency" in p.lower():
            low = p.lower()
            retracted = ("installed all along" in low or "neither blocker" in low
                         or "was false" in low or "simply unrun" in low)
            if not retracted:
                offenders.append(p[:220])
    assert not offenders, (
        "the README asserts an uninstalled dependency without retracting it in the same "
        "paragraph. openpyxl 3.1.5 is installed and scripts/fetch_lentimpra.py accepts "
        "--dest:\n" + "\n---\n".join(offenders)
    )


def test_the_correction_paragraph_actually_names_what_was_wrong(paragraphs):
    """Non-vacuity for the guard above: the retraction language must really be present."""
    hit = [p for p in paragraphs if "uninstalled dependency" in p.lower()]
    assert hit, "the phrase is absent entirely, so the guard above is vacuous"
    joined = " ".join(hit).lower()
    assert "installed all along" in joined
    assert "unrun" in joined


def test_openpyxl_really_is_importable_so_the_correction_is_true():
    """A non-vacuity check on the correction itself, not on the prose."""
    import importlib

    mod = importlib.import_module("openpyxl")
    assert mod is not None
    assert hasattr(mod, "__version__")


def test_the_fetcher_really_accepts_dest_so_the_correction_is_true():
    src = Path("scripts/fetch_lentimpra.py").read_text(encoding="utf-8")
    assert '"--dest"' in src


# ----------------------------------------------------------------------------------------
# verdict-independent 2: adjacency
# ----------------------------------------------------------------------------------------

def test_every_paragraph_stating_the_ratio_carries_a_caveat_in_the_same_paragraph(paragraphs):
    """A reader who stops at the headline must not get the uncaveated inference."""
    offenders = []
    for p in paragraphs:
        if re.search(r"20\s*[-–]\s*25x", p):
            if not any(w in p.lower() for w in CAVEAT_WORDS):
                offenders.append(p[:220])
    assert not offenders, (
        "a paragraph states the 20-25x ratio with no caveat word in the SAME paragraph:\n"
        + "\n---\n".join(offenders)
    )


def test_the_headline_paragraph_specifically_warns_against_quoting_the_ratio(paragraphs):
    hit = [p for p in paragraphs if "20–25x more on human regulatory DNA" in p]
    assert hit, "the headline paragraph was not found; this guard would be vacuous"
    joined = " ".join(hit).lower()
    assert "do not quote that" in joined or "frame caveat" in joined


def test_the_substrate_was_holding_the_method_back_claim_is_qualified_in_place(paragraphs):
    """That sentence is the inference the frame work questions, so it must not stand bare."""
    hit = [p for p in paragraphs if "holding the method back" in p]
    assert hit, "the claim was not found; this guard would be vacuous"
    assert any(w in " ".join(hit).lower() for w in CAVEAT_WORDS)


# ----------------------------------------------------------------------------------------
# verdict-dependent: the interim wording is gone and the artifact is cited
# ----------------------------------------------------------------------------------------

def test_the_interim_unresolved_wording_is_replaced_now_that_a_run_has_happened(text):
    assert "Read the ratio as unresolved, not as 20-25x." not in text, (
        "the interim placeholder is still present even though the measurements have run"
    )


def test_the_readme_cites_the_resolution_artifact(text):
    assert RESOLUTION in text, f"the README does not cite {RESOLUTION}"
    assert Path(RESOLUTION).exists(), "the README cites an artifact that does not exist"


def test_the_readme_reports_the_gate_as_ungraded_rather_than_quoting_its_margin(text):
    """The gate did not grade, so no margin from it may appear as a result."""
    assert "INDETERMINATE_NULL_NOT_CLEAN" in text
    # the ungraded active-frame conv margin must not be quoted anywhere in the README
    for forbidden in ("0.2190", "+0.2190"):
        assert forbidden not in text, (
            f"{forbidden} is the UNGRADED active-frame conv margin and must not appear as a "
            f"result while the gate reads INDETERMINATE"
        )


def test_the_curve_numbers_the_readme_does_quote_carry_their_null(text):
    """A trend without its band is the retracted statistic's mistake in a new costume."""
    if "+0.9879" in text:
        para = next(p for p in text.split("\n\n") if "+0.9879" in p)
        assert "0.4303" in para, "the trend is quoted without its attenuation-null maximum"


def test_the_headline_paragraph_qualifies_0_941_rather_than_presenting_it_as_single(paragraphs):
    """Scoped to the HEADLINE paragraph, which is where the claim is actually made.

    The first version of this guard matched ANY paragraph containing both `0.941` and the
    ratio, and so flagged the caveat paragraph -- where `0.941` appears as "none was borrowed
    from the bacterial 0.941", a correctly-scoped mention about the absent HUMAN ceiling. That
    was a guard defect, not a prose defect.
    """
    hit = [p for p in paragraphs if "20–25x more on human regulatory DNA" in p]
    assert hit, "the headline paragraph was not found; this guard would be vacuous"
    for p in hit:
        if "0.941" in p:
            assert "pooled average" in p.lower() or "0.50" in p, (
                "the headline paragraph cites 0.941 without the per-decile range that "
                "qualifies it"
            )


def test_wherever_0_941_is_called_a_ceiling_the_pooled_average_is_stated_somewhere(text):
    """The weaker document-level claim, kept because the per-paragraph one is scoped."""
    assert "pooled average" in text.lower(), (
        "the README never states that 0.941 is a pooled average over deciles, which is the "
        "one correction that changes how every bacterial number should be read"
    )
