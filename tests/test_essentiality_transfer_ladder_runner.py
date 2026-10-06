"""Guards for the ladder runner — the reconcile gate above all, since everything else depends on it."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import essentiality_transfer_ladder as run  # noqa: E402

from dna_decode.eval import transfer_ladder as tl  # noqa: E402

CACHE = Path("D:/dna_decode_cache/essentiality")
_HAVE = all((CACHE / f).exists() for f in (
    "CEGv2.txt", "NEGv1.txt", "Homo_sapiens.gene_info.gz",
    "goodall_TableS1_essential.xlsx", "ecoli_k12_feature_table.txt.gz"))


def _ok():
    return [{"key": "ecoli", "scored": True, "coverage_lift": 0.3125,
             "auroc_SECONDARY_not_cross_rung_comparable": 0.6952,
             "n_essential": 351, "n_nonessential": 3432},
            {"key": "human", "scored": True, "coverage_lift": 0.1633,
             "auroc_SECONDARY_not_cross_rung_comparable": 0.5805,
             "n_essential": 681, "n_nonessential": 899}]


# --------------------------------------------------------------------------------------------------
# the reconcile gate — live, not decorative
# --------------------------------------------------------------------------------------------------
def test_matching_values_reconcile():
    out = run.reconcile_or_raise(_ok())
    assert out["ecoli"]["coverage_lift"] == 0.3125
    assert out["human"]["auroc"] == 0.5805


@pytest.mark.parametrize("perturb", [
    {"coverage_lift": 0.3126},                                  # 4th dp on the primary
    {"auroc_SECONDARY_not_cross_rung_comparable": 0.6962},      # 3rd dp on the secondary
    {"n_essential": 350},
    {"n_nonessential": 3431},
    {"scored": False},
], ids=["lift", "auroc", "n_ess", "n_non", "unscored"])
def test_the_gate_RAISES_on_every_perturbation(perturb):
    """Non-vacuity. A gate that cannot fail is decoration — and this one already earned its place by
    catching a real divergence during execution (a first-wins vs last-wins duplicate-symbol rule that
    shifted coverage_lift to 0.3122 and AUROC to 0.6951)."""
    recs = [dict(_ok()[0], **perturb), _ok()[1]]
    with pytest.raises(run.ReconcileMismatch):
        run.reconcile_or_raise(recs)


def test_a_missing_rung_refuses_rather_than_skipping():
    with pytest.raises(run.ReconcileMismatch, match="did not score"):
        run.reconcile_or_raise([_ok()[0]])


# --------------------------------------------------------------------------------------------------
# walls
# --------------------------------------------------------------------------------------------------
def test_a_walled_rung_carries_no_metric_and_names_its_route():
    w = run.walled_record(tl.RUNGS_BY_KEY["saureus"], "no W0 record for ntml_nebraska",
                          "WALL_SCHEMA_UNVERIFIED")
    assert w["scored"] is False
    assert w["wall"] == "WALL_SCHEMA_UNVERIFIED"
    assert "ntml" in w["route_tried"]
    for forbidden in ("coverage_lift", "coverage_essential",
                      "auroc_SECONDARY_not_cross_rung_comparable"):
        assert forbidden not in w


# --------------------------------------------------------------------------------------------------
# null + thresholds
# --------------------------------------------------------------------------------------------------
def test_the_permutation_null_centres_on_zero():
    n = run.permutation_null([0] * 900 + [1] * 100, 500, draws=200)
    assert abs(n["null_mean"]) < 0.03, n
    assert n["null_p95"] <= n["null_max"]


def test_the_null_is_deterministic_under_its_seed():
    a = run.permutation_null([0] * 90 + [1] * 10, 50, draws=50)
    b = run.permutation_null([0] * 90 + [1] * 10, 50, draws=50)
    assert a == b


def test_the_frozen_thresholds_ship_with_their_derivation():
    t, d = run.FROZEN_THRESHOLDS, run.THRESHOLD_DERIVATION
    assert set(t) == {"plateau_tol", "cliff_drop", "min_rungs"}
    assert t["cliff_drop"] == 2 * t["plateau_tol"], "the recorded derivation says 2x the noise floor"
    assert d["measured_nulls"]["ecoli"]["null_max"] <= t["plateau_tol"]
    assert d["observed_two_rung_spread"] > t["cliff_drop"], "a cliff must be reachable"
    assert d["observed_two_rung_spread"] > t["plateau_tol"], "a plateau must not be trivially satisfied"


def test_the_bar_is_labelled_a_consistency_lock_not_an_endpoint_test():
    """It is derived from the two numbers that motivated the question, so it cannot also test them.
    Claiming 'frozen before the numbers' here would be pre-registration theatre."""
    lab = run.THRESHOLD_DERIVATION["HONEST_LABEL"]
    assert "CONSISTENCY LOCK" in lab.upper()
    assert "not a frozen-before-the-numbers" in lab


# --------------------------------------------------------------------------------------------------
# self-check + the real run
# --------------------------------------------------------------------------------------------------
def test_self_check_passes_offline():
    assert run._self_check() == 0


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_the_real_run_reconciles_reports_all_five_rungs_and_returns_INDETERMINATE():
    """The pre-committed outcome for the only-two-rungs case, asserted on the real data."""
    records = run.build_records()
    assert {r["key"] for r in records} == {r.key for r in tl.RUNGS}, "every rung must be reported"
    run.reconcile_or_raise(records)                      # raises if the committed numbers moved

    scored = [r for r in records if r.get("scored")]
    walled = [r for r in records if not r.get("scored")]
    assert {r["key"] for r in scored} == {"ecoli", "human"}
    assert all(r["wall"] == "WALL_SCHEMA_UNVERIFIED" for r in walled), walled

    v = tl.classify_ladder(records, **run.FROZEN_THRESHOLDS)
    assert v.verdict == tl.INDETERMINATE_INSUFFICIENT_RUNGS
    assert v.n_scored == 2


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_both_scored_rungs_clear_their_own_null_by_a_wide_margin():
    records = [r for r in run.build_records() if r.get("scored")]
    for r in records:
        assert r["coverage_lift"] > r["null_max"] * 2.5, r
        assert abs(r["null_mean"]) < 0.01, r


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_the_adjusted_lift_RISES_on_both_real_rungs():
    """Floor genes are misses, so removing them raises the lift. Reported, never clamped."""
    for r in (x for x in run.build_records() if x.get("scored")):
        assert r["coverage_lift_adjusted"] > r["coverage_lift"], r["key"]


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_the_duplicate_symbol_rule_is_what_reproduces_the_published_numbers():
    """Pins the divergence the reconcile gate caught: `first` must NOT reproduce, `last` must."""
    from dna_decode.essentiality import annotation_join as aj
    from dna_decode.essentiality.core_decoder import score_gene
    from essentiality_e3_learned import load_labels

    p = str(CACHE / "ecoli_k12_feature_table.txt.gz")
    lab = load_labels()

    def lift(rule):
        ann = aj.load_annotation(p, "ncbi_feature_table", on_duplicate=rule)
        e = [score_gene(g, ann[g]).core_score for g, v in lab.items() if v == 1 and g in ann]
        n = [score_gene(g, ann[g]).core_score for g, v in lab.items() if v == 0 and g in ann]
        return round(tl.coverage_lift(e, n), 4)

    assert lift("last") == 0.3125, "the default must reproduce the committed number"
    assert lift("first") == 0.3122, "and the alternative must measurably differ, else the rule is moot"


def test_an_invalid_duplicate_rule_raises():
    from dna_decode.essentiality import annotation_join as aj
    with pytest.raises(aj.AnnotationError, match="on_duplicate"):
        aj.load_annotation("x", "ncbi_feature_table", on_duplicate="middle")


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_no_rung_RECORD_carries_a_bare_auroc_key():
    """A field name is the cheapest disclosure that a number is not cross-rung comparable.

    Checked on the RECORDS rather than by counting source strings: the first version of this test
    asserted an arbitrary occurrence count (<=4) on the file text, which is brittle and said nothing
    about the data. `RECONCILE` and the reconciled summary legitimately use a short 'auroc' key -- what
    must never happen is a per-rung record exposing one.
    """
    for r in run.build_records():
        if r.get("scored"):
            assert "auroc" not in r, r["key"]
            assert "auroc_SECONDARY_not_cross_rung_comparable" in r, r["key"]


@pytest.mark.skipif(not _HAVE, reason="D: essentiality cache absent")
def test_the_emitted_artifact_is_parseable_and_carries_the_derivation():
    hits = sorted(ROOT.glob("wiki/essentiality_transfer_ladder_*.json"))
    if not hits:
        pytest.skip("ladder artifact not emitted on this host")
    art = json.loads(hits[-1].read_text(encoding="utf-8"))
    assert art["record"] == "essentiality-transfer-ladder-v1"
    assert art["frozen_thresholds"] == run.FROZEN_THRESHOLDS
    assert "HONEST_LABEL" in art["threshold_derivation"]
    assert len(art["rungs"]) == len(tl.RUNGS)
    assert art["verdict"]["verdict"].startswith("INDETERMINATE")
    # the eukaryote tie must be visible in the artifact, not just in code
    tied = [c for c in art["distance_classes"] if len(c["rungs"]) > 1]
    assert tied and sorted(tied[0]["rungs"]) == ["human", "scerevisiae"]
