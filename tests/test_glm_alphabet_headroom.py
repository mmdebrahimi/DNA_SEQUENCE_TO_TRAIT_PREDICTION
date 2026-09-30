"""Guards for the GLM headroom head-to-head (2026-09-30).

The real run needs gitignored Stanford data + the ESM masked-marginal cache; these pin the PURE bar logic
and the anti-laundering property that makes the comparison meaningful.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "glm_alphabet_headroom", REPO / "scripts" / "glm_alphabet_headroom.py")
G = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(G)

ART = REPO / "wiki" / "glm_alphabet_headroom_2026-09-30.json"


def _art():
    if not ART.exists():
        pytest.skip("artifact absent (gitignored data or ESM cache missing at build time)")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_bar_is_the_one_the_zero_shot_arm_FAILED():
    """THE ANTI-LAUNDERING PIN. Judging a new method by a fresh bar is how a recorded negative becomes a
    positive. The bar is copied verbatim from the committed 2026-07-09 artifact, and the number that
    failed it must STILL fail it -- if this ever passes, the bar drifted and the comparison is void."""
    z = G.PREREGISTERED["zero_shot_result_being_compared"]
    assert z["passed"] is False
    assert G.verdict(z["esm_mean_damage_auroc"], 0.4514, 0.4974) != "PASS"
    assert G.PREREGISTERED["min_auroc"] == 0.65
    assert G.PREREGISTERED["max_null_auroc"] == 0.55
    assert "hiv_esm_vs_catalog_2026-07-09" in G.PREREGISTERED["source"]


def test_the_committed_bar_matches_the_artifact_it_cites():
    """The bar claims to be verbatim from another artifact. Verify that, rather than trusting the claim."""
    src = REPO / "wiki" / "hiv_esm_vs_catalog_2026-07-09.json"
    if not src.exists():
        pytest.skip("cited 2026-07-09 artifact absent")
    d = json.loads(src.read_text(encoding="utf-8"))
    assert d["pass_bar"] == G.PREREGISTERED["rule_verbatim"], (
        "the bar drifted from the artifact it cites")
    # The ONLY licensed difference is the subject noun. My first version PARAPHRASED the bar
    # ('esm' -> 'score') and this test caught it -- a paraphrased bar cannot be verified against source.
    for thresh in ("0.65", "0.55", "burden"):
        assert thresh in G.PREREGISTERED["rule_verbatim"]
        assert thresh in G.PREREGISTERED["rule_as_applied"]
    assert (G.PREREGISTERED["rule_verbatim"].replace("esm", "score")
            == G.PREREGISTERED["rule_as_applied"]), "the generalisation changed more than the subject noun"
    assert d["passed"] is False
    assert d["esm_mean_damage_auroc"] == G.PREREGISTERED[
        "zero_shot_result_being_compared"]["esm_mean_damage_auroc"]


def test_every_bar_clause_can_independently_fail():
    """A bar whose clauses cannot each bite is decoration."""
    assert G.verdict(0.81, 0.45, 0.50) == "PASS"
    assert G.verdict(0.60, 0.45, 0.50) == "FAIL_BELOW_MIN_AUROC"
    assert G.verdict(0.70, 0.75, 0.50) == "FAIL_NOT_BETTER_THAN_MUTATION_BURDEN"
    assert G.verdict(0.70, 0.45, 0.60) == "FAIL_NULL_TOO_HIGH"
    assert G.verdict(0.65, 0.64, 0.5499) == "PASS", "the boundary must be inclusive as written"


def test_self_check_passes():
    G.self_check()


def test_the_pipeline_anchor_reproduced_the_committed_zero_shot_number():
    """The comparison is only trustworthy if our subset construction matches the 2026-07-09 run's. The
    anchor is that check; a DISAGREES anchor means nothing else in the artifact can be read."""
    d = _art()
    assert "REPRODUCES" in d["pipeline_anchor"], f"anchor did not reproduce: {d['pipeline_anchor']}"


def test_the_strictest_split_is_the_reported_headline_and_is_grouped_on_the_full_set():
    """An ungrouped k-fold leaks (1.35x patient duplication). The full set must be scored by
    leave-one-STUDY-out, and the headline must be that number -- not the flattering ungrouped one."""
    d = _art()
    full = next((c for c in d["cells"] if c["dataset"] == "full"), None)
    if full is None:
        pytest.skip("full dataset cell absent")
    assert full["strictest_split"] == "leave_one_study_out"
    sup = full["supervised_genotype_tokens_auroc"]
    assert full["strictest_supervised_auroc"] == min(sup.values()), \
        "the headline must be the STRICTEST (lowest) grouped number, never the most flattering"


def test_the_supervised_arm_clears_the_bar_the_zero_shot_arm_failed():
    """The finding itself. Pinned so a regression in the scan surfaces as a failure, not a quiet reversal."""
    d = _art()
    for c in d["cells"]:
        assert c["verdict_zero_shot"].startswith("FAIL"), "zero-shot must still fail its own bar"
        assert c["verdict_supervised"] == "PASS"
        assert c["delta_supervised_minus_zero_shot"] > 0.2


def test_the_blind_spot_is_catalog_negative_so_it_cannot_rediscover_the_catalog():
    """The whole claim rests on this: if any isolate carried a catalogued major DRM, the model could be
    re-deriving the catalog instead of finding headroom beyond it. Verified against the DEPLOYED catalog."""
    import dna_decode.data.hiv_amr as H
    d = _art()
    majors = H.NNRTI_RT_MAJOR_DRMS
    for c in d["cells"]:
        for t in c["top_resistance_features_full_fit"]:
            tok = t["token"]                                     # '<pos><aa>', e.g. '179D'
            pos = int("".join(ch for ch in tok if ch.isdigit()))
            aa = tok[-1]
            wt = H._RT_WT.get(pos)
            if wt is not None:
                assert f"{wt}{pos}{aa}" not in majors, \
                    f"{tok} IS a catalogued major DRM -- the subset is not catalog-negative"


def test_honest_scope_is_recorded_and_names_the_linear_floor():
    """A pass here is a FLOOR on the supervised family, not evidence that context/attention helps. If that
    caveat ever drops out of the artifact the result reads as far stronger than it is."""
    d = _art()
    scope = d["honest_scope"].lower()
    assert "weakest" in scope and "floor" in scope
    assert "in-distribution" in scope
