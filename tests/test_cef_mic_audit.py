from dna_decode.data.cohort import CandidateStrain, StrainCohort
from scripts.cef_mic_audit import _tier_detail
from dna_decode.data.mic_tiers import breakpoints_for, classify_tier


def test_cef_high_r_tier():
    breakpoints = breakpoints_for("ceftriaxone")
    assert classify_tier([16.0], {"R"}, breakpoints) == "HIGH_R"


def test_cef_high_s_tier():
    breakpoints = breakpoints_for("ceftriaxone")
    assert classify_tier([0.125], {"S"}, breakpoints) == "HIGH_S"


def test_tier_detail_contains_calls():
    breakpoints = breakpoints_for("ceftriaxone")
    detail = _tier_detail([4.0, 8.0], breakpoints)
    assert detail["median_mic"] == 6.0
    assert detail["clsi_call"] == "R"
    assert detail["eucast_call"] == "R"


def test_cohort_ids_prefers_saved_pool_ids():
    from scripts.cef_mic_audit import _cohort_ids_for_drug

    cohort = StrainCohort(
        strains=[
            CandidateStrain("a", ast_labels={"ceftriaxone": 1}),
            CandidateStrain("b", ast_labels={"ceftriaxone": 0}),
        ],
        per_drug_strain_ids={"ceftriaxone": ["a"]},
        three_drug_intersection=[],
    )

    ids, scope = _cohort_ids_for_drug(cohort, "ceftriaxone")
    assert scope == "pool"
    assert ids == {"a"}


def test_cohort_ids_can_force_all_labeled():
    from scripts.cef_mic_audit import _cohort_ids_for_drug

    cohort = StrainCohort(
        strains=[
            CandidateStrain("a", ast_labels={"ceftriaxone": 1}),
            CandidateStrain("b", ast_labels={"ceftriaxone": 0}),
        ],
        per_drug_strain_ids={"ceftriaxone": ["a"]},
        three_drug_intersection=[],
    )

    ids, scope = _cohort_ids_for_drug(cohort, "ceftriaxone", all_labeled=True)
    assert scope == "labeled"
    assert ids == {"a", "b"}
