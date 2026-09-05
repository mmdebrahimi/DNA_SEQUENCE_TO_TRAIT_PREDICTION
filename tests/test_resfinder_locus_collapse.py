"""The ResFinder caller reports one gene per LOCUS, not one per matching catalog allele.

THE BUG THIS PINS. beta-lactamase variants differ by one to three point mutations, so a single blaTEM
locus matches ~180 catalog TEM alleles above the 90% identity bar. The pre-fix caller keyed its output
on the ALLELE name, so every one of them cleared independently and was reported as a separately present
gene -- including ESBLs (blaTEM-52, blaTEM-12) in a genome carrying only the narrow-spectrum blaTEM-1.
That is a wrong clinical reading of a genome, not a cosmetic count.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.resfinder.runner import (  # noqa: E402
    LOCUS_OVERLAP_FRACTION, _overlap_fraction, cluster_alleles_by_locus,
)
from resfinder_locus_collapse_validate import normalize_symbol  # noqa: E402

ARTIFACT = ROOT / "wiki" / "resfinder_locus_collapse_2026-09-05.json"


def hit(pid: float, cov: float, contig: str = "c1", s: int = 100, e: int = 900) -> dict:
    return {"called": True, "percent_identity": pid, "percent_coverage": cov,
            "contig": contig, "sstart": s, "send": e}


# --- the collapse ---------------------------------------------------------------------------------

def test_many_variants_of_one_locus_collapse_to_one_call():
    """The actual bug: 180 TEM variants at one position are ONE gene, not 180."""
    called = [(f"blaTEM-{i}_1_ACC", hit(99.0 + (i % 10) / 100, 100.0)) for i in range(1, 181)]
    called.append(("blaTEM-1B_1_AY458016", hit(100.0, 100.0)))
    clusters = cluster_alleles_by_locus(called)
    assert len(clusters) == 1
    assert len(clusters[0]) == 181


def test_the_winner_is_the_highest_identity_allele_not_the_first_seen():
    """Identity-primary. Within a locus every variant sits at ~100% coverage, so a coverage-first
    tiebreak is decided by dict order rather than by sequence -- the serotype/salmserovar defect."""
    called = [("blaTEM-52B_1_X", hit(99.7, 100.0)), ("blaTEM-1B_1_Y", hit(100.0, 100.0)),
              ("blaTEM-12_1_Z", hit(99.3, 100.0))]
    cluster = cluster_alleles_by_locus(called)[0]
    winner = max(cluster, key=lambda kv: (kv[1]["percent_identity"], kv[1]["percent_coverage"]))
    assert winner[0] == "blaTEM-1B_1_Y"


def test_the_esbl_overcall_is_gone():
    """A genome carrying only narrow-spectrum blaTEM-1 must not be reported as carrying blaTEM-52."""
    called = [("blaTEM-1B_1_Y", hit(100.0, 100.0)), ("blaTEM-52B_1_X", hit(99.7, 100.0)),
              ("blaTEM-12_1_Z", hit(99.3, 100.0)), ("blaTEM-10_1_W", hit(99.8, 100.0))]
    reported = {max(c, key=lambda kv: (kv[1]["percent_identity"], kv[1]["percent_coverage"]))[0]
                for c in cluster_alleles_by_locus(called)}
    assert reported == {"blaTEM-1B_1_Y"}


# --- what must NOT collapse -----------------------------------------------------------------------

def test_different_loci_stay_separate():
    """blaOXA-1 and blaOXA-48 share a name prefix, are functionally unrelated (narrow-spectrum vs
    carbapenemase), and can genuinely co-occur. Grouping by NAME would have merged them; position
    keeps them apart."""
    called = [("blaOXA-1_1_A", hit(100.0, 100.0, contig="c1", s=100, e=900)),
              ("blaOXA-48_1_B", hit(100.0, 100.0, contig="c2", s=100, e=900))]
    assert len(cluster_alleles_by_locus(called)) == 2


def test_a_tandem_array_does_not_chain_into_one_call():
    """Single-linkage would merge adjacent copies through a chain of pairwise overlaps. This project
    has seen a real 7-copy tandem blaTEM array, so the clustering is greedy-representative."""
    called = [(f"blaTEM-1_{i}_A", hit(100.0 - i * 0.01, 100.0, contig="c1",
                                      s=100 + i * 900, e=900 + i * 900)) for i in range(7)]
    assert len(cluster_alleles_by_locus(called)) == 7


def test_an_allele_with_no_position_becomes_its_own_cluster():
    """What cannot be placed must not be silently merged into someone else's locus."""
    called = [("a_1_A", hit(100.0, 100.0)),
              ("b_1_B", {"called": True, "percent_identity": 99.0, "percent_coverage": 100.0,
                         "contig": None, "sstart": None, "send": None})]
    assert len(cluster_alleles_by_locus(called)) == 2


# --- the overlap primitive ------------------------------------------------------------------------

def test_overlap_is_zero_across_contigs():
    assert _overlap_fraction(("c1", 1, 100), ("c2", 1, 100)) == 0.0


def test_overlap_is_measured_against_the_shorter_interval():
    """A short hit fully inside a long one is the same locus, which a union-based fraction would miss."""
    assert _overlap_fraction(("c1", 1, 1000), ("c1", 400, 600)) == pytest.approx(1.0)


def test_minus_strand_intervals_are_normalized():
    """A minus-strand HSP arrives as sstart>send; unnormalized it yields an empty interval."""
    called = [("a_1_A", hit(100.0, 100.0, s=900, e=100)),
              ("b_1_B", hit(99.0, 100.0, s=100, e=900))]
    assert len(cluster_alleles_by_locus(called)) == 1


def test_the_overlap_bar_is_a_named_constant_not_a_literal():
    assert 0.0 < LOCUS_OVERLAP_FRACTION <= 1.0


# --- symbol normalization -------------------------------------------------------------------------

def test_trailing_variant_letter_is_stripped_but_a_bare_gene_is_untouched():
    assert normalize_symbol("blaTEM-1B") == "blaTEM-1"
    assert normalize_symbol("blaTEM-1") == "blaTEM-1"
    assert normalize_symbol("aph(3'')-Ib") == "aph(3'')-Ib"   # no digits before the letter
    assert normalize_symbol("aadA1") == "aadA1"


# --- the committed artifact -----------------------------------------------------------------------

@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("locus-collapse artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_artifact_reports_improvement_in_every_class(art):
    assert art["verdict"] == "LOCUS_COLLAPSE_IMPROVES_AGREEMENT_WITH_AN_INDEPENDENT_CALLER"
    for cls, s in art["summary"].items():
        assert s["jaccard_vs_amrfinder_normalized_new"] > s["jaccard_vs_amrfinder_normalized_old"], cls


def test_artifact_is_non_vacuous(art):
    """Both rules agree on the empty set if blastn never called anything."""
    assert sum(s["total_alleles_called"] for s in art["summary"].values()) > 0
    assert sum(s["n_genomes"] for s in art["summary"].values()) > 0


def test_artifact_says_the_comparator_is_a_tool_not_a_wetlab_label(art):
    """Agreement with an independent implementation is not correctness; the tier does not move."""
    joined = " ".join(art["honest_limits"]).lower()
    assert "not a wet-lab label" in joined or "not correctness" in joined
    assert "faithful_to_tool" in joined
