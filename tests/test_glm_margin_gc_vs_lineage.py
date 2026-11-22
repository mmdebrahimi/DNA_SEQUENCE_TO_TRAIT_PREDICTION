"""The GC-versus-lineage probe must prove its own precondition before it reports a GC effect.

This probe exists to close a limit declared earlier the same day: the cross-organism margin correlated
with GC at r = -0.941 and had to be called uninterpretable, because the two low-GC genomes were both
*Campylobacter* and the three high-GC ones all *Enterobacteriaceae*. GC and phylogeny were the same
variable, so the effective sample size was two groups rather than five points.

That makes the collinearity guard the load-bearing part, and it has a specific non-vacuity requirement:
**it must reject the data that motivated building it.** A guard that passed the four-species case would
be decoration. `test_the_guard_rejects_the_confounded_case_that_motivated_it` pins exactly that.

The second subtlety is the within-group reading. Centring each phylum's margin and GC removes anything
lineage can explain, so a relationship that survives is about GC for its own sake. But with one genome
per phylum there is nothing left after centring, and an undefined statistic reported as zero would turn
a sampling design into a finding. That case returns a reason rather than a number.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_margin_gc_vs_lineage.py"
OUT = ROOT / "wiki" / "glm_margin_gc_vs_lineage_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_gcvlin", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the collinearity guard

def test_the_guard_rejects_the_confounded_case_that_motivated_it(mod):
    """THE non-vacuity test. The four-species shape -- two low-GC Campylobacter, three high-GC
    Enterobacteriaceae -- must come back NOT decoupled, or the guard is decoration."""
    gcs = [0.266, 0.272, 0.453, 0.457, 0.505]
    phyla = ["Campylobacterota", "Campylobacterota",
             "Pseudomonadota", "Pseudomonadota", "Pseudomonadota"]
    out = mod.collinearity(gcs, phyla)
    assert out["decoupled"] is False, out
    assert abs(out["corr_gc_vs_phylum_mean"]) >= mod.COLLINEARITY_CEILING


def test_the_guard_accepts_a_genuinely_decoupled_corpus(mod):
    """Control: if GC varies WITHIN phyla as much as between them, the question is answerable."""
    gcs = [0.29, 0.44, 0.67, 0.31, 0.46, 0.69, 0.27, 0.43, 0.72]
    phyla = ["A", "A", "A", "B", "B", "B", "C", "C", "C"]
    out = mod.collinearity(gcs, phyla)
    assert out["decoupled"] is True, out
    assert abs(out["corr_gc_vs_phylum_mean"]) < mod.COLLINEARITY_CEILING


def test_one_genome_per_phylum_is_maximally_confounded(mod):
    """Each phylum mean IS that genome's GC, so phylum predicts GC perfectly by construction.

    This is the shape the expanded corpus actually has, and reporting it honestly is the point: more
    species does not by itself decouple the two variables.
    """
    gcs = [0.27, 0.40, 0.55, 0.70]
    phyla = ["A", "B", "C", "D"]
    out = mod.collinearity(gcs, phyla)
    assert out["corr_gc_vs_phylum_mean"] == pytest.approx(1.0, abs=1e-6)
    assert out["decoupled"] is False


def test_the_guard_refuses_rather_than_dividing_by_zero(mod):
    out = mod.collinearity([0.4, 0.4, 0.4], ["A", "A", "A"])
    assert out["decoupled"] is False
    assert "too little variation" in out["reason"]


# ------------------------------------------------------------------ the within-group reading

def test_a_pure_lineage_effect_collapses_after_centring(mod):
    """If margin is set by phylum and GC merely tracks phylum, centring must remove the relationship."""
    groups = ["A", "A", "A", "B", "B", "B"]
    gcs = [0.30, 0.31, 0.29, 0.70, 0.71, 0.69]
    margins = [0.10, 0.10, 0.10, 0.30, 0.30, 0.30]   # set purely by group
    out = mod.within_group(margins, gcs, groups)
    assert out["corr"] is None or abs(out["corr"]) < 0.2, out


def test_a_real_within_group_gc_effect_survives_centring(mod):
    """Control: if margin tracks GC INSIDE each group, the relationship must survive."""
    groups = ["A", "A", "A", "B", "B", "B"]
    gcs = [0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    margins = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
    out = mod.within_group(margins, gcs, groups)
    assert out["corr"] is not None and out["corr"] > 0.9, out


def test_one_member_per_group_returns_a_reason_not_a_zero(mod):
    """An undefined statistic reported as zero turns a sampling design into a finding."""
    out = mod.within_group([0.1, 0.2, 0.3], [0.3, 0.5, 0.7], ["A", "B", "C"])
    assert out["corr"] is None
    assert "NOT a null result" in out["reason"]
    assert out["n_groups_with_replicates"] == 0


def test_groups_with_replicates_are_counted(mod):
    out = mod.within_group([0.1, 0.2, 0.3, 0.4], [0.3, 0.4, 0.5, 0.6], ["A", "A", "B", "C"])
    assert out["n_groups_with_replicates"] == 1


# ------------------------------------------------------------------ geometry and thresholds

def test_the_geometry_matches_the_committed_kernel(mod):
    """The margin must be the SAME quantity measured elsewhere, or it is not comparable."""
    assert mod.LENGTH == 150 and mod.LIMIT == 1200 and mod.N_GENERATED == 120
    assert mod.MARKOV_ORDER == 3 and mod.KMER_K == 4 and mod.POSITIONAL_BINS == 6


def test_positional_is_measured_alongside_kmer_and_both(mod):
    """Positional alone carries the signal; omitting it would miss the mode that matters."""
    assert set(mod.MODES) == {"kmer", "positional", "both"}


def test_the_claim_thresholds_are_stated(mod):
    assert mod.MIN_SPECIES_FOR_A_CLAIM >= 5
    assert mod.MIN_PHYLA_FOR_A_CLAIM >= 3
    assert 0 < mod.COLLINEARITY_CEILING < 1


def test_all_three_verdicts_are_reachable_in_the_source(mod):
    src = SCRIPT.read_text(encoding="utf-8")
    for v in ("STILL_CONFOUNDED__NO_GC_CLAIM_POSSIBLE",
              "DECOUPLED_BUT_UNDERPOWERED",
              "DECOUPLED_AND_POWERED__GC_READING_IS_INTERPRETABLE"):
        assert v in src, v


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_artifact_states_the_precondition_before_any_gc_reading():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert "collinearity_precondition" in rec
    assert "decoupled" in rec["collinearity_precondition"]
    assert rec["verdict"].startswith(("STILL_CONFOUNDED", "DECOUPLED"))


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_a_gc_claim_is_only_made_when_the_precondition_holds():
    """The guard must actually gate the verdict, not sit beside it."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    if not rec["collinearity_precondition"]["decoupled"]:
        assert rec["verdict"] == "STILL_CONFOUNDED__NO_GC_CLAIM_POSSIBLE"


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_what_this_closes_names_the_earlier_uninterpretable_number():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert "-0.941" in rec["what_this_closes"] or "0.941" in rec["what_this_closes"]
    assert "effective n" in rec["what_this_closes"]


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_undefined_within_group_case_is_declared_in_limits():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "undefined" in joined
    assert "not as zero" in joined or "rather than as zero" in joined
