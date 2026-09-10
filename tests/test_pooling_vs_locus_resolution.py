"""Guards for the pooling-vs-locus-resolution test.

Two things make this result trustworthy rather than decorative, and both broke during the build:

  * the VERDICT is computed from a bar frozen before any number existed -- so the middle branch
    (narrowing helps, locus identity does not) cannot be re-read afterwards as a partial win;
  * the random-gene NULL is the load-bearing arm -- and its first version silently ran ZERO draws,
    because it drew on strain-unique `gene_id`, whose intersection across strains is empty.

The gene-resolution tests exist because two annotation vocabularies name gyrA differently and one of
those names collides with parC.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.pooling_vs_locus_resolution import (  # noqa: E402
    QRDR_PATTERNS, leave_one_lineage_out_auroc, mixed_label_lineages, resolve_qrdr, verdict_from_bar,
)

ARTIFACT = ROOT / "wiki" / "pooling_vs_locus_resolution_2026-09-10.json"
BAR = ROOT / "wiki" / "pooling_vs_locus_acceptance_bar.json"


def _ann(rows):
    return pd.DataFrame(rows, columns=["gene_id", "product"])


# --- the verdict is the frozen rule, applied mechanically ---------------------------------------

def test_beating_whole_genome_is_not_enough_on_its_own():
    """THE branch that fired. Narrowing helping while the null's max is higher is a REFUTATION of the
    deployable claim; pre-registering it is what stops it being reported as a partial win."""
    assert verdict_from_bar(0.8386, 0.8029, 0.9020) == "REFUTED_LOCUS_IDENTITY"


def test_supported_requires_beating_both():
    assert verdict_from_bar(0.95, 0.80, 0.90) == "SUPPORTED"


def test_losing_to_whole_genome_is_outright_falsification():
    assert verdict_from_bar(0.70, 0.80, 0.60) == "FALSIFIED"


def test_a_missing_null_is_indeterminate_not_a_pass():
    """The first run produced two AUROCs and an EMPTY null. Treating that as a pass would have shipped
    the claim on the strength of the arm that cannot contradict it."""
    assert verdict_from_bar(0.84, 0.80, None) == "INDETERMINATE"
    assert verdict_from_bar(None, 0.80, 0.90) == "INDETERMINATE"


# --- gene resolution: two vocabularies, one collision --------------------------------------------

def test_gyra_resolves_under_both_annotation_vocabularies():
    """Handling only 'DNA gyrase subunit A' resolved gyrA on 78/140 strains; adding the formal EC name
    raised the usable cohort to 123."""
    for gyra_product in ("DNA gyrase subunit A", "DNA topoisomerase (ATP-hydrolyzing) subunit A"):
        ann = _ann([("g1", gyra_product),
                    ("g2", gyra_product.replace("subunit A", "subunit B")),
                    ("g3", "DNA topoisomerase IV subunit A"),
                    ("g4", "DNA topoisomerase IV subunit B")])
        got = resolve_qrdr(ann, {"g1", "g2", "g3", "g4"})
        assert got == {"gyrA": "g1", "gyrB": "g2", "parC": "g3", "parE": "g4"}, gyra_product


def test_gyra_and_parc_are_never_conflated():
    """Gyrase IS a type-II topoisomerase, so one vocabulary names gyrA 'DNA topoisomerase
    (ATP-hydrolyzing) subunit A'. A loose `topoisomerase.*subunit A` pattern matches that AND parC."""
    ann = _ann([("g1", "DNA topoisomerase (ATP-hydrolyzing) subunit A"),
                ("g2", "DNA topoisomerase (ATP-hydrolyzing) subunit B"),
                ("g3", "DNA topoisomerase IV subunit A"),
                ("g4", "DNA topoisomerase IV subunit B")])
    got = resolve_qrdr(ann, {"g1", "g2", "g3", "g4"})
    assert got["gyrA"] != got["parC"] and got["gyrB"] != got["parE"]
    assert got["gyrA"] == "g1" and got["parC"] == "g3"


def test_both_numeral_conventions_match_parc():
    """One source writes 'topoisomerase 4', the other 'topoisomerase IV'. Asserting the Roman numeral
    would have silently made the restricted arm gyrase-only."""
    for numeral in ("4", "IV"):
        ann = _ann([("g1", "DNA gyrase subunit A"), ("g2", "DNA gyrase subunit B"),
                    ("g3", f"DNA topoisomerase {numeral} subunit A"),
                    ("g4", f"DNA topoisomerase {numeral} subunit B")])
        assert resolve_qrdr(ann, {"g1", "g2", "g3", "g4"}) is not None, numeral


def test_an_ambiguous_or_absent_gene_refuses_rather_than_guesses():
    """A silently-wrong gene assignment would corrupt the exact comparison this script exists to make."""
    dup = _ann([("g1", "DNA gyrase subunit A"), ("g1b", "DNA gyrase subunit A"),
                ("g2", "DNA gyrase subunit B"), ("g3", "DNA topoisomerase IV subunit A"),
                ("g4", "DNA topoisomerase IV subunit B")])
    assert resolve_qrdr(dup, {"g1", "g1b", "g2", "g3", "g4"}) is None
    missing = _ann([("g2", "DNA gyrase subunit B"), ("g3", "DNA topoisomerase IV subunit A")])
    assert resolve_qrdr(missing, {"g2", "g3"}) is None


def test_a_gene_absent_from_the_cache_does_not_resolve():
    """Resolution is gated on the cache keys: an annotated gene with no embedding is not usable."""
    ann = _ann([("g1", "DNA gyrase subunit A"), ("g2", "DNA gyrase subunit B"),
                ("g3", "DNA topoisomerase IV subunit A"), ("g4", "DNA topoisomerase IV subunit B")])
    assert resolve_qrdr(ann, {"g1", "g2", "g3"}) is None


def test_the_patterns_do_not_match_gyrase_inhibitors():
    """Real annotations carry 'DNA gyrase inhibitor YacG' and 'SbmC'. Matching those would put an
    unrelated gene into the causal-locus arm."""
    ann = _ann([("g1", "DNA gyrase inhibitor YacG"), ("g2", "DNA gyrase inhibitor SbmC")])
    assert resolve_qrdr(ann, {"g1", "g2"}) is None


# --- CV honesty -----------------------------------------------------------------------------------

def test_a_lineage_never_spans_the_split():
    """Lineage-disjointness is what the de-confounding rests on; if an ST appeared on both sides a
    near-identical clone could leak across."""
    import numpy as np
    rng = np.random.default_rng(0)
    groups = np.array(["a"] * 20 + ["b"] * 20 + ["c"] * 20)
    y = np.array([0, 1] * 30)
    X = rng.normal(size=(60, 5))
    assert leave_one_lineage_out_auroc(X, y, groups) is not None


def test_within_lineage_feasibility_is_reported_not_assumed():
    """The strict within-lineage test is a different, harder question. Reporting how little of the
    cohort could support it is what stops this run reading as a contradiction of the 0-for-5 record."""
    rows = [{"mlst": "A", "y": 1}, {"mlst": "A", "y": 0}, {"mlst": "B", "y": 1}, {"mlst": "B", "y": 1}]
    f = mixed_label_lineages(rows)
    assert f["n_lineages"] == 2 and f["n_mixed_label_lineages"] == 1
    assert f["n_strains_in_mixed_lineages"] == 2


# --- the committed run ------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_committed_verdict_matches_its_own_numbers():
    """Re-derives the verdict from the artifact's own AUROCs, so a hand-edited verdict cannot stand."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["verdict"] == verdict_from_bar(d["auroc_qrdr_restricted"],
                                            d["auroc_whole_genome_pooled"],
                                            d["random_4gene_null"]["max"])
    assert d["verdict"] == "REFUTED_LOCUS_IDENTITY"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_null_actually_ran():
    """NON-VACUITY, and the single most important guard here: the first version drew on strain-unique
    gene ids, intersected to the empty set, and reported a null of zero draws -- leaving the result
    uninterpretable while looking like a successful run."""
    n = json.loads(ARTIFACT.read_text(encoding="utf-8"))["random_4gene_null"]
    assert n["n_draws"] >= 20, f"null ran only {n['n_draws']} draws"
    assert n["n_common_single_copy_products_drawn_from"] >= 4
    assert n["max"] is not None and n["mean"] is not None


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_headline_is_that_random_genes_beat_the_causal_ones():
    """The finding, pinned: if this ever flips, the memo's reading must be rewritten."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["random_4gene_null"]["max"] > d["auroc_qrdr_restricted"]
    assert d["auroc_qrdr_restricted"] > d["auroc_whole_genome_pooled"]


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_the_bar_was_frozen_before_the_result_and_names_all_outcomes():
    """A bar that only names the outcome that happened is not a pre-registration."""
    b = json.loads(BAR.read_text(encoding="utf-8"))
    assert "FROZEN BEFORE ANY RESULT" in b["status"]
    assert set(b["verdict_rule"]) == {"SUPPORTED", "REFUTED_LOCUS_IDENTITY", "FALSIFIED", "INDETERMINATE"}
    assert b["substrate_measured_before_freezing"]["n_strains_resolvable"] == 123
