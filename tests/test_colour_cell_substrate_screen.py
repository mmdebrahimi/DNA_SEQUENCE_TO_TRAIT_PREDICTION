"""The colour-cell substrate screen: pure classifier + the trust-surface corrections it forced.

WHY THIS EXISTS. The animal colour/plumage family reached 19 CLI cells, all KNOWLEDGE_BASELINE, before
anyone asked whether they COULD be validated. The screen (scripts/colour_cell_substrate_screen.py) derives
that per-cell from the committed catalogs and found two walls: 40 of 65 loci (62%) record no causal variant
at all, and 14 of the 25 that DO are indel/structural -- off-panel for any SNP array.

These tests pin the classifier (a text heuristic, so it needs anchoring) + the corrections that landed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from colour_cell_substrate_screen import (  # noqa: E402
    _CATALOG_GAPS, classify_variant, collect, self_check, snv_panel_scorable, summarise,
)

SCREEN_MD = ROOT / "wiki" / "colour_cell_substrate_screen_2026-08-26.md"
NEG_MAP = ROOT / "wiki" / "negative_results_map_2026-06-13.md"


# ------------------------------------------------------------------------ the classifier

@pytest.mark.parametrize("text,expect", [
    ("MC1R p.Arg306Ter = recessive red `e`", "SNV"),
    ("TYR c.604C>G p.His202Asp, recessive", "SNV"),
    ("MLPH c.-22G>A `d`", "SNV"),                      # negative-offset promoter coordinate
    ("CBD103 beta-defensin c.67_69delGGT", "INDEL"),
    ("ASIP non-agouti black c.181_184delTTCA", "INDEL"),
    ("exon-2 frameshift", "INDEL"),
    ("ASIP A^y/a^t SINE insertion + coding", "STRUCTURAL"),
    ("a 190kb duplication", "STRUCTURAL"),
    ("rabbit C (TYR): C full > chinchilla > Himalayan > c albino", "UNRECORDED"),
    ("mouse a locus (ASIP): A agouti > at > a non-agouti", "UNRECORDED"),
])
def test_classify_variant_cases(text, expect):
    assert classify_variant(text) == expect


def test_structural_beats_indel_because_a_sine_insertion_also_matches_ins():
    """ORDER IS LOAD-BEARING. 'SINE insertion' contains 'insertion'; if INDEL were tested first, every
    structural variant would be mis-filed as an indel and the STRUCTURAL count would read zero."""
    assert classify_variant("ASIP SINE insertion") == "STRUCTURAL"
    assert classify_variant("a 190kb duplication of ASIP") == "STRUCTURAL"


def test_a_frameshift_that_also_cites_a_point_coordinate_is_still_an_indel():
    """Real provenance strings mix notations in one sentence; the coarser class must win."""
    assert classify_variant("MLPH c.667_668insC frameshift (see also c.-22G>A)") == "INDEL"


def test_snv_panel_scorability_is_tri_state_not_boolean():
    """UNRECORDED must be None, never False -- 'we did not write the variant down' is NOT the same claim
    as 'a SNP panel cannot carry it', and collapsing them would fabricate evidence about the substrate."""
    assert snv_panel_scorable("SNV") is True
    assert snv_panel_scorable("INDEL") is False
    assert snv_panel_scorable("STRUCTURAL") is False
    assert snv_panel_scorable("UNRECORDED") is None


def test_empty_source_is_unrecorded_not_a_crash():
    assert classify_variant("") == "UNRECORDED"


# ------------------------------------------------------------------------ the summary verdicts

def _rows(*classes):
    return [{"locus": str(i), "variant_class": c} for i, c in enumerate(classes)]


@pytest.mark.parametrize("classes,verdict", [
    (("UNRECORDED", "UNRECORDED"), "UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED"),
    (("SNV", "SNV"), "FULLY_SNV_TRACTABLE"),
    (("SNV", "UNRECORDED"), "SNV_TRACTABLE_WHERE_RECORDED"),
    (("INDEL", "STRUCTURAL"), "NO_LOCUS_SNV_TRACTABLE"),
    (("SNV", "INDEL"), "PARTIALLY_SNV_TRACTABLE"),
])
def test_summarise_verdicts(classes, verdict):
    assert summarise(_rows(*classes))["verdict"] == verdict


def test_blocked_count_excludes_unrecorded():
    """A cell with no recorded variants has ZERO blocked loci -- it is unscreenable, not blocked. Counting
    unrecorded as blocked would overstate the substrate wall, which is the whole error this screen exists
    to avoid making."""
    s = summarise(_rows("UNRECORDED", "UNRECORDED", "SNV"))
    assert s["n_snv_panel_blocked"] == 0


# ------------------------------------------------------------------------ real catalogs

def test_the_screen_runs_on_every_committed_colour_catalog():
    """No network, no D: — the catalogs are committed Python."""
    data = collect()
    assert len(data) >= 19, f"expected >=19 colour cells, got {len(data)}: {sorted(data)}"
    assert "dog" in data and "rabbit" in data


def test_self_check_against_the_dog_catalog_passes():
    """The classifier is a text heuristic, so it is anchored on the ONE cell with measured ground truth.
    This check EARNED ITS KEEP on the first run: it flagged dog A as UNRECORDED against an expectation of
    STRUCTURAL, and the classifier was RIGHT -- the expectation had encoded the literature, not the
    catalog (the dog ASIP entry never names the SINE)."""
    assert self_check(collect()) == []


def test_the_dog_asip_catalog_gap_is_recorded_rather_than_papered_over():
    """The gap between what the measured artifact knows and what the catalog records is itself a finding:
    even the most-developed colour cell omits one of its five causal variants."""
    assert "dog/A" in _CATALOG_GAPS
    from dna_decode.pigment import dog_coat
    src = getattr(dog_coat.LOCI["A"], "source", "")
    assert "SINE" not in src, "the catalog now records the SINE — retire the gap entry deliberately"


def test_the_two_headline_counts_are_reproducible_from_the_catalogs():
    """Pins the memo's numbers to the code. If a catalog gains a causal variant these MOVE — update the
    memo in the same commit rather than loosening the test."""
    data = collect()
    tot = {}
    for rows in data.values():
        for r in rows:
            tot[r["variant_class"]] = tot.get(r["variant_class"], 0) + 1
    assert sum(tot.values()) == 65, f"loci total moved: {tot}"
    assert tot.get("UNRECORDED") == 40
    assert tot.get("INDEL", 0) + tot.get("STRUCTURAL", 0) == 14


# ------------------------------------------------------------------------ trust-surface corrections

def test_coatcolor_reports_its_measured_result_instead_of_calling_it_pending():
    """UNDER-CLAIM regression. The contract framed a run that HAPPENED (2026-07-30) as 'the v0.1 measured
    tier', and never reported black 0.994. Under-claiming is as much a trust-surface falsehood as
    over-claiming."""
    from dna_decode.cli import TRAITS
    v = TRAITS["coatcolor"]["validation"]
    assert "the v0.1 measured tier" not in v
    assert "0.994" in v and "UNSCORABLE ON THAT SUBSTRATE" in v
    assert (ROOT / "wiki" / "dog_coat_darwins_ark_measured_2026-07-30.md").exists()


def test_the_two_decoder_side_gates_are_in_the_negative_results_map():
    """G1-G8 all gate the LABEL; a curated-catalog cell can fail before a label is even relevant."""
    t = NEG_MAP.read_text(encoding="utf-8", errors="replace")
    assert "| G9 |" in t and "| G10 |" in t
    assert "## The 10 rejection gates" in t, "the gate-count heading must track the table"
    assert "screen it against G1–G10" in t, "the how-to-use line must reference every gate"


def _flat(p: Path) -> str:
    """Markdown prose is hard-wrapped, so a phrase spanning a line break defeats a naive `in` check.
    Collapse whitespace before matching -- reflowing the MEMO to suit a test would be backwards."""
    return " ".join(p.read_text(encoding="utf-8", errors="replace").split())


@pytest.mark.skipif(not SCREEN_MD.exists(), reason="screen memo absent")
def test_the_memo_does_not_claim_unrecorded_loci_are_evidence_about_substrate():
    t = _flat(SCREEN_MD)
    assert "statement about the catalog, never evidence about the substrate" in t
    assert "unvalidatable as written" in t.lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
