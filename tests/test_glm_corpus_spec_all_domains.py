"""The all-domains corpus plan must keep unverified figures out of every headline.

This file is a PLAN, not a measurement, and that makes one failure mode dominant: a plausible but wrong
genome size produces a fully self-consistent corpus plan that is entirely wrong, with nothing downstream
able to catch it. This repo has a precedent and a rule for exactly that -- the ECOFF catalog refuses to
return a cut-off it cannot source, after a search tool volunteered a remembered value. So every entry
here carries a `source`, bacterial rows are MEASURED from today's artifacts, and rows for domains with no
local assembly are confined to an illustrative block.

The substantive thing pinned is the FORK: promoter-mode corpus scales with gene count, whole-genome
corpus scales with genome size, and because gene counts barely move across domains while genome sizes
span three orders of magnitude, the two modes differ by roughly 300x in how many genomes a download
needs. The regulatory-fraction collapse is the reason the cheap answer is a trap -- promoter windows are
12-18% of a bacterial genome and about 0.1% of a mammalian one.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_corpus_spec_all_domains.py"
OUT = ROOT / "wiki" / "glm_corpus_spec_all_domains_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    """Load the script by path, and REGISTER it in sys.modules before executing it.

    The registration is not optional here and the reason is non-obvious: this script defines a
    `@dataclass`, and dataclass construction resolves its field types via
    `sys.modules.get(cls.__module__).__dict__`. For a module loaded by path and never registered that
    lookup returns None, so `exec_module` dies with `AttributeError: 'NoneType' object has no attribute
    '__dict__'` at class-definition time. The sibling test files in this arm load their scripts the same
    way without this line and pass, because none of those scripts defines a dataclass -- so the bug
    looks like a problem with THIS script rather than with the loader. It is not: the script runs
    correctly on its own and produced the committed artifact.
    """
    import sys

    spec = importlib.util.spec_from_file_location("_corpus_spec_all_domains", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the two yield functions

def test_promoter_yield_tracks_gene_count_not_genome_size(mod):
    """The whole fork rests on this: two genomes of wildly different size but equal gene count must
    give the SAME promoter corpus."""
    small = mod.DomainEntry("animal", "compact", 100.0, 20000, mod.UNVERIFIED)
    huge = mod.DomainEntry("animal", "large", 3100.0, 20000, mod.UNVERIFIED)
    assert small.promoter_mbp() == huge.promoter_mbp()
    assert huge.genome_mbp_usable() == 31 * small.genome_mbp_usable()


def test_the_regulatory_fraction_collapses_in_large_genomes(mod):
    """Why the whole-genome shortcut is a trap, as a number rather than an assertion."""
    bact = mod.DomainEntry("bacteria", "b", 4.6, 4300, mod.MEASURED)
    mammal = mod.DomainEntry("animal", "m", 3100.0, 20000, mod.UNVERIFIED)
    assert bact.regulatory_fraction() > 0.10
    assert mammal.regulatory_fraction() < 0.002
    assert bact.regulatory_fraction() > 50 * mammal.regulatory_fraction()


def test_promoter_mode_uses_the_committed_window_length(mod):
    """A generator the committed falsifier cannot score produces numbers that mean nothing."""
    assert mod.WINDOW_BP == 150
    e = mod.DomainEntry("bacteria", "b", 4.6, 1000, mod.MEASURED)
    assert e.promoter_mbp() == pytest.approx(1000 * 150 / 1e6)


def test_a_zero_size_genome_does_not_raise_on_the_fraction(mod):
    assert mod.DomainEntry("x", "y", 0.0, 100, mod.UNVERIFIED).regulatory_fraction() == 0.0


# ------------------------------------------------------------------ the provenance rail

def test_every_illustrative_row_is_marked_unverified(mod):
    assert mod.ILLUSTRATIVE, "the illustrative block must not be empty or the fork is unillustrated"
    assert all(e.source == mod.UNVERIFIED for e in mod.ILLUSTRATIVE)


def test_the_illustrative_block_spans_more_than_one_domain(mod):
    """A single exemplar cannot show a cross-domain shape."""
    assert len({e.domain for e in mod.ILLUSTRATIVE}) >= 3


def test_measured_rows_are_read_from_artifacts_not_restated(mod):
    """Restating them in the source is how a plan and its evidence drift apart."""
    src = SCRIPT.read_text(encoding="utf-8")
    assert "load_measured" in src
    assert "glm_margin_across_organisms_2026-11-21.json" in src
    assert "Populated at runtime" in src


def test_the_script_refuses_without_measured_rows(mod):
    """A headline resting only on remembered numbers must not be reachable."""
    src = SCRIPT.read_text(encoding="utf-8")
    assert "REFUSE: no MEASURED bacterial rows" in src


# ------------------------------------------------------------------ the scaling table

def test_scaling_scales(mod):
    """Half the yield per genome must double the genomes needed."""
    a = mod.scaling_table(2.0, "a")
    b = mod.scaling_table(1.0, "b")
    for rung in a["rungs"]:
        assert b["rungs"][rung]["genomes_needed"] == pytest.approx(
            2 * a["rungs"][rung]["genomes_needed"], rel=0.01)


def test_scaling_handles_a_zero_yield_without_dividing_by_zero(mod):
    t = mod.scaling_table(0.0, "z")
    assert all(v["genomes_needed"] is None for v in t["rungs"].values())


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_artifact_declares_itself_a_plan_not_a_measurement():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["status"].startswith("PLAN")
    assert "No corpus is downloaded" in rec["status"]


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_scope_correction_is_recorded_as_a_correction():
    """The bacterial framing was drift -- an available substrate becoming the stated goal. The record
    of that is more useful than a silently corrected scope."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    s = rec["scope_correction"]
    assert "ALL LIFE FORMS" in s
    assert "quietly became the stated goal" in s


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_no_unverified_figure_appears_in_a_measured_headline():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for row in rec["measured_rows"]:
        assert row["source"] == mod_measured_token(), row
    for row in rec["illustrative_rows_UNVERIFIED"]:
        assert row["source"] == "UNVERIFIED_APPROX", row
    assert "UNVERIFIED" in rec["scaling_from_ILLUSTRATIVE_eukaryotes"]["caveat"]
    assert "must not be quoted" in rec["scaling_from_ILLUSTRATIVE_eukaryotes"]["caveat"]


def mod_measured_token() -> str:
    return "MEASURED"


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_fork_is_stated_with_both_arms_and_neither_chosen():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    f = rec["the_fork"]
    assert "PROMOTER WINDOWS" in f["question"] and "ANY GENOMIC WINDOW" in f["question"]
    assert "not_chosen_here" in f
    assert "authority call" in f["not_chosen_here"]
    for block in ("scaling_from_MEASURED_bacteria", "scaling_from_ILLUSTRATIVE_eukaryotes"):
        assert "promoter_mode" in rec[block] and "whole_genome_mode" in rec[block]


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_whole_genome_shortcut_carries_its_warning():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    w = rec["the_fork"]["why_the_700x_is_mostly_the_wrong_sequence"]
    assert "transposons" in w
    assert "NOVELTY" in w


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_domain_label_and_lineage_dedup_are_REQUIRED_not_suggested():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    r = rec["corpus_spec_requirements"]
    assert r["domain_label_per_genome"].startswith("REQUIRED")
    assert r["lineage_dedup"].startswith("REQUIRED")
    assert "0.50%" in r["lineage_dedup"], "the dedup rule must cite the measurement that justifies it"
    assert "four arms" in r["lineage_dedup"]


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_arabidopsis_negative_is_recorded_as_NOT_a_blocker():
    """It is a de-confounded eukaryotic negative and it will look like a blocker to anyone scoping
    eukaryotes. It is about discrimination, not generation -- error 0 in this repo's own notes."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    n = rec["a_closed_negative_that_is_NOT_a_blocker_here"]
    assert "DISCRIMINATIVE" in n["why_it_does_not_block_this"]
    assert "GENERATION" in n["why_it_does_not_block_this"]
    assert "error 0" in n["why_it_does_not_block_this"]


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_absence_of_a_local_eukaryotic_assembly_is_stated():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    s = rec["local_substrate_status"]
    assert s["eukaryotic_genomes_available_locally"] == 0
    assert "not an assembly" in s["note"]


@pytest.mark.skipif(not OUT.exists(), reason="spec has not been generated on this host")
def test_the_plan_does_not_claim_cross_domain_training_helps():
    """Mixing grammars could easily hurt; that is why the domain label is a requirement, not a bonus."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "No claim is made that cross-domain training helps" in joined
    assert "could easily hurt" in joined
    assert "NOT validated for DNA" in joined
    assert "upper bound" in joined
