"""Pin the TMP-SMX external scorer: strata math, artifact schema/branding, and the frozen-leak guard."""
from scripts.tmp_smx_external_validate import build_artifact, score_cohort


def _fixture():
    """Synthetic cohort: sul+dfr high-R, sul-only S, dfr-only S, neither S (the expected separation)."""
    geno = {
        "a": ["sul1", "dfrA17"], "b": ["sul2", "dfrA1"], "c": ["sul1", "dfrA12"],  # sul+dfr
        "d": ["sul2"], "e": ["sul1"],                                              # sul-only
        "f": ["dfrA1"],                                                            # dfr-only
        "g": ["blaCTX-M-15"], "h": [],                                             # neither
    }
    mics = {"a": 32, "b": 16, "c": 8, "d": 1, "e": 0.5, "f": 1, "g": 0.5, "h": 1}
    return geno, mics


def test_strata_and_reproduction():
    cells = score_cohort(*_fixture())
    s = cells["strata"]
    assert s["sul+dfr"]["n"] == 3 and s["sul+dfr"]["r_rate"] == 1.0
    assert s["sul-only"]["r_rate"] == 0.0
    assert s["neither"]["r_rate"] == 0.0
    assert cells["strata_reproduced"] is True


def test_binary_metric_primary():
    cells = score_cohort(*_fixture())
    b = cells["binary"]
    # a,b,c are sul+dfr with MIC>=4 (R) and the rule predicts R -> 3 TP; S strains predicted S
    assert b["tp"] == 3 and b["fp"] == 0
    assert b["sens"] == 1.0 and b["spec"] == 1.0


def test_strata_fail_marks_indeterminate():
    # sul-only is high-R (rule would mis-call) -> reproduction fails -> INDETERMINATE headline
    geno = {"a": ["sul1", "dfrA1"], "d": ["sul2"], "e": ["sul1"]}
    mics = {"a": 32, "d": 32, "e": 16}  # sul-only strains are R here
    cells = score_cohort(geno, mics)
    assert cells["strata"]["sul-only"]["r_rate"] >= 0.5
    assert cells["strata_reproduced"] is False
    art = build_artifact("synthetic", cells, "runX", "indep", "leak")
    assert art["headline"] == "INDETERMINATE"


def test_artifact_schema_and_branding():
    art = build_artifact("sci234", score_cohort(*_fixture()), "runX", "indep note", "leak note")
    assert art["_schema"] == "external-validation-v1"
    assert art["drug"] == "trimethoprim-sulfamethoxazole"
    assert art["run_id"] == "runX"
    assert art["rule_status"] == "EXPERIMENTAL_SCORED"
    assert art["rule_scope"] == "scorer_local"
    assert art["not_in_shipped_surface"] is True
    assert art["headline"] == "SCORED"
    assert "strict" in art and "binary" in art and "strata" in art


def test_binary_intermediate_mic_dropped_from_binary_and_strata():
    # MIC in (2, 4) is the _binary None gap: such a strain must NOT count toward the
    # binary pairs OR the strata R/S tallies (only its tier may reach strict/relaxed).
    geno = {"x": ["sul1", "dfrA17"]}
    cells = score_cohort(geno, {"x": 3})          # 2 < 3 < 4 -> _binary returns None
    assert cells["binary"]["n_scored"] == 0
    assert cells["strata"]["sul+dfr"]["n"] == 0    # not tallied into the stratum either
    assert cells["strata"]["sul+dfr"]["r_rate"] is None


def test_dfr_only_stratum_assigned():
    # dfr-only routing: a strain with only a dfr gene lands in the dfr-only stratum (not neither).
    cells = score_cohort({"x": ["dfrA1"]}, {"x": 1})  # MIC 1 -> S in binary
    assert cells["strata"]["dfr-only"]["n"] == 1 and cells["strata"]["dfr-only"]["S"] == 1
    assert cells["strata"]["sul-only"]["n"] == 0


def test_reproduction_holds_when_sul_only_stratum_empty():
    # so (sul-only r_rate) is None when no sul-only strain exists -> the `so is None`
    # branch must still allow reproduction when sul+dfr is the max-R stratum.
    geno = {"a": ["sul1", "dfrA17"], "g": ["blaCTX-M-15"]}
    cells = score_cohort(geno, {"a": 32, "g": 0.5})   # sul+dfr R, neither S; no sul-only
    assert cells["strata"]["sul-only"]["r_rate"] is None
    assert cells["strata_reproduced"] is True


def test_frozen_leak_guard():
    """TMP-SMX must NOT leak into the frozen deployed surface."""
    from dna_decode.data.mic_tiers import supported_drugs
    from dna_decode.eval.amr_rules import DRUG_RULE
    assert "trimethoprim-sulfamethoxazole" not in supported_drugs()
    assert "trimethoprim-sulfamethoxazole" not in DRUG_RULE
    assert "cotrimoxazole" not in supported_drugs()
