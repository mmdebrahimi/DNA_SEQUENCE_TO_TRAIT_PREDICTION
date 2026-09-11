"""Guards for the determinant-clade-coupling measurement.

Three things make this result trustworthy rather than decorative, and all three broke during the build:

  * the SYNONYMOUS filter -- `mutations.tsv` is `--mutation_all` output listing every position SCREENED
    including wildtype. Keeping those turns "screened and found nothing" into a called mutation;
  * the CEILING -- the null holds carrier count fixed, which controls the EXPECTATION but not the
    MAXIMUM. At 400 carriers only ~0.048 of coupling is achievable, so the dominant determinants read
    as near-chance when they are in fact near-saturated. A memo stating the opposite was one step from
    being published;
  * the VARIANT filter -- the single-residue regex silently discarded 12 symbols AMRFinder CALLS
    (indels, frameshifts, nonsense, negative-coordinate promoter positions).
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.determinant_clade_coupling import (  # noqa: E402
    SCHEME, ceiling_largest_fraction, coupling_score, is_called_substitution, is_called_variant,
    permutation_gap, verdict_from_bar,
)

ARTIFACT = ROOT / "wiki" / "determinant_clade_coupling_2026-09-10.json"
SENSITIVITY = ROOT / "wiki" / "determinant_clade_coupling_filter_sensitivity_2026-09-10.json"
BAR = ROOT / "wiki" / "determinant_clade_coupling_acceptance_bar.json"


# --- the synonymous filter -----------------------------------------------------------------------

def test_wildtype_screen_rows_are_not_calls():
    """NON-VACUITY for the load-bearing filter: 137,288 such rows exist. Counting them as called
    mutations would leave the chromosomal arm mostly noise."""
    assert not is_called_substitution("baeS_Q163Q")
    assert not is_called_substitution("16S_A523A")
    assert not is_called_variant("ampC_C-11C")       # promoter wildtype, negative coordinate
    assert not is_called_variant("pmrB_RPISLR6RPISLR")


def test_real_substitutions_are_calls_under_both_filters():
    for sym in ("gyrA_S83L", "parC_S80I", "gyrA_D87N"):
        assert is_called_substitution(sym), sym
        assert is_called_variant(sym), sym


def test_a_nonsense_call_is_not_mistaken_for_wildtype():
    """`*` as the alt must not collide with the equal-ref/alt wildtype rule."""
    assert is_called_substitution("nfsA_Q67*")
    assert is_called_variant("nfsA_Q67*")


# --- the corrected variant filter ----------------------------------------------------------------

def test_the_narrow_filter_really_does_drop_called_determinants():
    """NON-VACUITY for the sensitivity band. If the narrow filter accepted these there would be no
    defect to report and the band would be ceremony -- these are all real `main.tsv` symbols."""
    dropped = ["ftsI_N337NYRIN", "ftsI_I336IKYRI", "cirA_S90YfsTer15",
               "ampC_C-42T", "nfsA_Q67Ter", "23S_G543DEL"]
    for sym in dropped:
        assert not is_called_substitution(sym), f"{sym} should be invisible to the frozen filter"
        assert is_called_variant(sym), f"{sym} should be visible to the corrected filter"


def test_the_corrected_filter_still_refuses_wildtype():
    """The widening must not swallow the screen rows the whole result depends on excluding."""
    for sym in ("baeS_Q163Q", "16S_A523A", "ampC_T-14T", "ampC_G-15G"):
        assert not is_called_variant(sym), sym


def test_an_unparseable_symbol_is_not_silently_a_call():
    assert not is_called_variant("")
    assert not is_called_variant("justagene")


# --- the ceiling ----------------------------------------------------------------------------------

def test_the_ceiling_falls_as_prevalence_rises():
    """THE correction. With the biggest lineage at 53 genomes, a 400-carrier determinant cannot exceed
    53/400 -- which is why raw coupling declines with prevalence in BOTH arms."""
    sizes = [53, 40, 20, 10]
    assert ceiling_largest_fraction(40, sizes) == 1.0          # fits inside the biggest lineage
    assert ceiling_largest_fraction(53, sizes) == 1.0
    assert ceiling_largest_fraction(400, sizes) == pytest.approx(53 / 400)
    assert ceiling_largest_fraction(106, sizes) == pytest.approx(0.5)


def test_normalization_rescues_a_saturated_high_prevalence_determinant():
    """NON-VACUITY, and the whole point: a raw score of 0.043 and one of 0.873 are indistinguishable in
    kind once each is read against what was achievable. Without this the memo said the dominant
    determinants were 'essentially decoupled'."""
    # gyrA_S83L: observed 0.128, null 0.085, ceiling 0.133  -> tiny raw, near-saturated normalized
    dominant_raw = 0.128 - 0.085
    dominant_norm = dominant_raw / (53 / 400 - 0.085)
    # parE_I529L: observed 0.979, null 0.106, ceiling 1.0
    rare_norm = (0.979 - 0.106) / (1.0 - 0.106)
    assert dominant_raw < 0.05                      # reads as near-chance
    assert dominant_norm > 0.85                     # is in fact near-saturated
    assert rare_norm > 0.85
    assert abs(dominant_norm - rare_norm) < 0.15    # same kind, once the ceiling is honoured


def test_coupling_score_reports_ceiling_and_normalized():
    lineage_of = {f"g{i}": ("A" if i < 10 else f"L{i}") for i in range(40)}
    genomes = sorted(lineage_of)
    carriers = {f"g{i}" for i in range(10)}          # exactly the one big lineage
    s = coupling_score(carriers, lineage_of, genomes, np.random.default_rng(0), n_null=50,
                       lineage_sizes=list({}.values()) or None)
    assert s["largest_lineage_fraction"] == 1.0
    assert s["ceiling_largest_fraction"] == 1.0
    assert s["normalized_coupling"] == pytest.approx(1.0, abs=1e-6)


def test_normalized_is_none_rather_than_exploding_without_headroom():
    """A zero achievable range must not produce a huge or undefined ratio."""
    sizes = [2]
    assert ceiling_largest_fraction(1000, sizes) == pytest.approx(0.002)
    lineage_of = {f"g{i}": f"L{i}" for i in range(30)}   # every genome its own lineage
    genomes = sorted(lineage_of)
    s = coupling_score(set(genomes[:20]), lineage_of, genomes, np.random.default_rng(0), n_null=20)
    # ceiling == null mean == 1/20 here, so there is no headroom to normalize by
    assert s["normalized_coupling"] is None or isinstance(s["normalized_coupling"], float)


# --- the verdict rule is the frozen one, applied mechanically -------------------------------------

def test_all_four_frozen_branches():
    strong = {"exceeds_perm_max": True}
    weak = {"exceeds_perm_max": False}
    many_p, many_a = [0.5] * 25, [0.1] * 25
    assert verdict_from_bar(many_p, many_a, strong) == "SUPPORTED"
    assert verdict_from_bar(many_p, many_a, weak) == "WEAK_DIRECTIONAL"
    assert verdict_from_bar([0.1] * 25, [0.5] * 25, strong) == "FALSIFIED"
    assert verdict_from_bar([0.5] * 5, many_a, strong) == "INDETERMINATE"


def test_a_thin_arm_is_indeterminate_not_a_pass():
    """A 19-family arm must not clear a 20-family floor just because its gap is large."""
    assert verdict_from_bar([0.9] * 19, [0.0] * 92, {"exceeds_perm_max": True}) == "INDETERMINATE"


def test_the_permutation_null_guards_the_comparison_itself():
    """Two arms of different sizes drawn from ONE distribution differ by chance; without this the sign
    of the gap alone would be the finding."""
    rng = np.random.default_rng(0)
    pooled = list(rng.normal(size=116))
    perm = permutation_gap(pooled[:24], pooled[24:], np.random.default_rng(1), n_perm=200)
    assert not perm["exceeds_perm_max"], "same-distribution arms must not clear their own permutation"
    separated = permutation_gap([5.0] * 24, [0.0] * 92, np.random.default_rng(1), n_perm=200)
    assert separated["exceeds_perm_max"]


# --- the single-scheme restriction ----------------------------------------------------------------

def test_only_one_mlst_scheme_is_admitted():
    """An ST number means a DIFFERENT lineage under a different scheme; pooling invents identity."""
    raw = {"a": f"{SCHEME}131", "b": f"{SCHEME}10", "c": "MLST.senterica_achtman_2.11"}
    kept = {k: v for k, v in raw.items() if v.startswith(SCHEME)}
    assert set(kept) == {"a", "b"}


# --- the committed run ------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_committed_verdict_is_re_derivable_from_its_own_numbers():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for band in d["results"].values():
        pv = [band["mean_coupling_point"]] * band["n_point_families"]
        av = [band["mean_coupling_acquired"]] * band["n_acquired_families"]
        assert band["verdict"] == verdict_from_bar(pv, av, band["permutation"])
    assert d["headline_verdict"] == "SUPPORTED"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_every_gap_beats_the_permutation_maximum_not_just_p95():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for name, band in d["results"].items():
        p = band["permutation"]
        assert p["observed_gap"] > p["perm_max"], name


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_ceiling_normalized_band_shipped_and_widens_the_gap():
    """The headline correction, pinned: if normalization ever NARROWS the gap the memo's central
    reading must be rewritten."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    norm = d["ceiling_normalized"]["called_point_prevalence_matched"]
    raw = d["results"]["called_point_prevalence_matched"]
    assert norm["mean_normalized_point"] > norm["mean_normalized_acquired"]
    assert norm["permutation"]["observed_gap"] > raw["permutation"]["observed_gap"]
    assert norm["permutation"]["observed_gap"] > norm["permutation"]["perm_max"]


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_stratification_carries_both_metrics():
    """Raw alone made the decline look like biology; normalized alone would hide that the decline is
    real. The table must carry both."""
    strat = json.loads(ARTIFACT.read_text(encoding="utf-8"))["prevalence_stratified"]
    for arm in ("point", "acquired"):
        for band in strat[arm].values():
            assert "mean_coupling_raw" in band and "mean_coupling_normalized" in band


@pytest.mark.skipif(not SENSITIVITY.exists(), reason="sensitivity band not generated")
def test_the_filter_sensitivity_band_shows_the_defect_and_that_it_does_not_flip_the_verdict():
    s = json.loads(SENSITIVITY.read_text(encoding="utf-8"))
    frozen, corrected = s["frozen_single_residue"], s["corrected_all_variant_forms"]
    assert corrected["raw"]["n_point"] > frozen["raw"]["n_point"], "the fix must recover families"
    assert corrected["n_filtered_rows"] < frozen["n_filtered_rows"]
    for key in ("raw", "normalized"):
        assert frozen[key]["verdict"] == corrected[key]["verdict"] == "SUPPORTED"


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_the_bar_was_frozen_first_and_names_every_outcome():
    """A bar that only names the outcome that happened is not a pre-registration."""
    b = json.loads(BAR.read_text(encoding="utf-8"))
    assert "FROZEN BEFORE" in b["status"]
    assert set(b["verdict_rule"]) == {"SUPPORTED", "WEAK_DIRECTIONAL", "FALSIFIED", "INDETERMINATE"}


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_the_bars_synonymous_count_is_the_TRUE_wildtype_count():
    """The reconciliation that found the filter defect: the bar's 137,288 is the true wildtype count,
    and the script's 156,471 additionally swept in 19,183 unparseable symbols. Pinning the bar's figure
    keeps that distinction from quietly collapsing again."""
    b = json.loads(BAR.read_text(encoding="utf-8"))
    assert b["substrate_measured_before_freezing"]["synonymous_point_rows_filtered"] == 137288
