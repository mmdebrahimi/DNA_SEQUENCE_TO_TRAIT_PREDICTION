"""Offline tests for the conv ⊕ max_tx rank-average hybrid arm.

The arm's own hypothesis is made falsifiable here, on synthetic data, BEFORE it meets the real table:
an orthogonal second modality must lift, and a near-chance one must not. If those two did not hold on
data we construct, a real-data result either way would be uninterpretable.
"""
from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import glm_tile_hybrid_arm as arm  # noqa: E402
from glm_tile_hybrid_arm import (  # noqa: E402
    PREREGISTERED,
    as_table,
    combine,
    orthogonality,
    render_md,
    shuffled,
    verdict,
)
from glm_tile_predictions import TilePredictionError, table_payload  # noqa: E402

from dna_decode.forward.variant_effect import rank_average_hybrid  # noqa: E402
from dna_decode.glm.genomewide import spearman  # noqa: E402


# --------------------------------------------------------------------------------------------------
# the combine plumbing
# --------------------------------------------------------------------------------------------------
def test_as_table_keys_are_positions_as_strings():
    assert as_table([0.5, 0.1, 0.9]) == {"0": 0.5, "1": 0.1, "2": 0.9}


def test_combine_returns_original_index_order():
    """`rank_average_hybrid` returns a dict over sorted STRING keys ('10' < '2'), so the caller must
    re-index by position or every vector past index 9 silently reorders."""
    n = 12
    preds = {"a": [float(i) for i in range(n)], "b": [float(i) for i in range(n)]}
    out = combine(preds, ("a", "b"))
    assert len(out) == n
    assert out == sorted(out), "two identical ascending scorers must combine to an ascending vector"


def test_combine_is_monotone_when_both_modalities_agree():
    preds = {"a": [1.0, 2.0, 3.0, 4.0], "b": [10.0, 20.0, 30.0, 40.0]}
    assert combine(preds, ("a", "b")) == sorted(combine(preds, ("a", "b")))


def test_rank_average_hybrid_contract_is_respected():
    """Its own refusals: fewer than 2 tables, and an empty key intersection."""
    with pytest.raises(ValueError, match="needs >=2"):
        rank_average_hybrid([{"0": 1.0}])
    with pytest.raises(ValueError, match="no variant"):
        rank_average_hybrid([{"0": 1.0}, {"1": 1.0}])


def test_combine_survives_the_two_digit_key_boundary():
    """Element-wise at n=11, the first length where a string sort disagrees with a numeric one.

    `rank_average_hybrid` keys on sorted STRINGS, so '10' sorts between '1' and '2'. At n<=10 the two
    orders coincide and the bug is invisible; n=11 is the smallest length that exposes it. The
    expected vector is written out rather than checked for monotonicity, so a reindex that merely
    happened to stay sorted could not pass.
    """
    n = 11
    v = [float(i) for i in range(n)]
    assert combine({"a": v, "b": v}, ("a", "b")) == [i / (n - 1) for i in range(n)]


def test_combine_preserves_the_input_ORDER_not_merely_sortedness():
    """The decisive form: a SHUFFLED vector combined with itself must rank-correlate 1.0 with itself.

    The existing ascending-input test passes even if positions are permuted into sorted order, because
    its correct answer IS sorted. With a shuffled input any mis-indexing shows up as rho < 1.
    """
    v = [float(i) for i in range(25)]
    random.Random(4).shuffle(v)
    out = combine({"a": v, "b": v}, ("a", "b"))
    assert spearman(out, v) == pytest.approx(1.0, abs=1e-9)


def test_combine_of_three_modalities_also_holds_the_order():
    """The `hybrid_conv_maxtx_gc` combo takes three tables, not two."""
    v = [float(i) for i in range(30)]
    random.Random(5).shuffle(v)
    out = combine({"a": v, "b": v, "c": v}, ("a", "b", "c"))
    assert len(out) == 30
    assert spearman(out, v) == pytest.approx(1.0, abs=1e-9)


def test_combine_on_a_length_mismatch_raises_one_way_and_TRUNCATES_the_other():
    """Recorded as it IS, not as it should be — and why that is safe here.

    Positional keys intersect to the shorter vector. When the FIRST name is the longer one, the
    re-index walks past the intersection and raises. When it is the shorter one, `n` is the short
    length and the result is silently truncated. The protection against the second case is upstream,
    not here: `table_payload` and `load_table` both refuse a ragged table, so a mismatch cannot reach
    `combine` from a real seed file. Pinned so the asymmetry is visible rather than assumed away.
    """
    with pytest.raises(KeyError):
        combine({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0]}, ("a", "b"))
    assert len(combine({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0]}, ("b", "a"))) == 2


def test_shuffled_preserves_the_multiset_but_not_the_order():
    v = [float(i) for i in range(50)]
    s = shuffled(v, seed=0)
    assert sorted(s) == sorted(v)
    assert s != v


def test_shuffled_is_deterministic_per_seed():
    v = [float(i) for i in range(30)]
    assert shuffled(v, 3) == shuffled(v, 3)
    assert shuffled(v, 3) != shuffled(v, 4)


# --------------------------------------------------------------------------------------------------
# the arm's HYPOTHESIS, both directions, on constructed data
# --------------------------------------------------------------------------------------------------
def _truth(n: int = 400) -> list[float]:
    return [float(i) for i in range(n)]


def test_an_ORTHOGONAL_second_modality_LIFTS():
    """The positive direction. Two noisy views of the truth with INDEPENDENT noise must combine
    better than either alone -- this is the published mechanism the arm is testing for."""
    rng = random.Random(0)
    truth = _truth()
    a = [t + rng.gauss(0, 120) for t in truth]
    b = [t + rng.gauss(0, 120) for t in truth]
    sa, sb = spearman(a, truth), spearman(b, truth)
    sh = spearman(combine({"a": a, "b": b}, ("a", "b")), truth)
    assert sh > max(sa, sb), f"hybrid {sh:.4f} did not beat max(a={sa:.4f}, b={sb:.4f})"


def test_a_NEAR_CHANCE_second_modality_REDUCES_the_signal():
    """The negative direction, and the measured failure mode (ESM2+ProSST 0.204 vs ESM2 0.352).

    Without this test a real-data negative could not be distinguished from a broken combine.
    """
    rng = random.Random(1)
    truth = _truth()
    good = [t + rng.gauss(0, 60) for t in truth]
    junk = [rng.gauss(0, 1) for _ in truth]
    sg = spearman(good, truth)
    sh = spearman(combine({"a": good, "b": junk}, ("a", "b")), truth)
    assert sh < sg, f"hybrid {sh:.4f} should be WORSE than the good modality alone {sg:.4f}"


def test_orthogonality_separates_redundant_from_independent():
    rng = random.Random(2)
    truth = _truth()
    a = [t + rng.gauss(0, 50) for t in truth]
    assert orthogonality(a, a) == pytest.approx(1.0, abs=1e-9)
    b = [rng.gauss(0, 1) for _ in truth]
    assert abs(orthogonality(a, b)) < 0.25


# --------------------------------------------------------------------------------------------------
# the frozen bar
# --------------------------------------------------------------------------------------------------
def _p(median: float, wins: int, n: int = 10) -> dict:
    return {"median": median, "wins": wins, "n": n, "ci95": [median - 0.01, median + 0.01]}


def test_verdict_hybrid_adds():
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 9),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "HYBRID_ADDS"


def test_verdict_weak_when_below_margin():
    arms = {"_paired_hybrid_vs_comparator": _p(0.005, 9),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "WEAK"


def test_verdict_weak_when_below_seed_count():
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 5),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "WEAK"


def test_verdict_no_gain_on_a_negative_margin():
    arms = {"_paired_hybrid_vs_comparator": _p(-0.004, 2),
            "_paired_control_vs_comparator": _p(0.000, 4)}
    assert verdict(arms)[0] == "NO_GAIN"


def test_verdict_REFUSES_when_the_control_also_lifts():
    """A positive that the control reproduces is an artifact of rank-averaging, not evidence.

    Ordered BEFORE the HYBRID_ADDS branch on purpose: otherwise a run where both lift would print a
    headline and bury its own control.
    """
    arms = {"_paired_hybrid_vs_comparator": _p(0.03, 9),
            "_paired_control_vs_comparator": _p(0.03, 9)}
    v, detail = verdict(arms)
    assert v == "INDETERMINATE_CONTROL_ALSO_LIFTS"
    assert "artifact" in detail


def test_verdict_REFUSES_on_a_drifted_table():
    arms = {"_reconcile_problems": ["gc seed 0 = 0.2000, not 0.3 +-0.01"],
            "_paired_hybrid_vs_comparator": _p(0.05, 10),
            "_paired_control_vs_comparator": _p(0.0, 3)}
    assert verdict(arms)[0] == "INDETERMINATE_RECONCILE_FAILED"


def test_verdict_REFUSES_on_a_missing_scorer():
    assert verdict({"_missing_scorer": ["conv"]})[0] == "INDETERMINATE_MISSING_SCORER"


def test_refusals_outrank_a_headline():
    """A drifted table AND a huge apparent gain must still refuse."""
    arms = {"_missing_scorer": ["conv"], "_reconcile_problems": ["x"],
            "_paired_hybrid_vs_comparator": _p(0.99, 10),
            "_paired_control_vs_comparator": _p(0.0, 0)}
    assert verdict(arms)[0].startswith("INDETERMINATE")


def test_frozen_bar_mirrors_the_decomposition_arm():
    """Same margin and seed count, so 'does more information help?' is comparable across the two."""
    assert PREREGISTERED["margin"] == 0.02
    assert PREREGISTERED["min_seeds_winning"] == 8
    assert PREREGISTERED["comparator"] == "max_tx"


def test_the_precondition_is_reported_not_used_as_a_gate():
    """Predicting a combine will fail is not a substitute for measuring it (R6)."""
    assert PREREGISTERED["precondition"]["measured_first"] is True
    assert "not_a_gate" in PREREGISTERED["precondition"]


# --------------------------------------------------------------------------------------------------
# render_md -- the human-facing surface. A caveat carried only in the JSON is not a disclosure.
# --------------------------------------------------------------------------------------------------
def _art(**over) -> dict:
    art = {
        "date": "2026-10-14",
        "verdict": "WEAK",
        "verdict_detail": "hybrid +0.0054 over max_tx, 7/10 seeds",
        "orthogonality": {"median": 0.6421, "min": 0.5900, "max": 0.7100,
                          "per_seed": [0.59, 0.64, 0.71]},
        "arms": {
            "max_tx": {"median": 0.3968, "n": 2, "per_seed": [0.39, 0.40]},
            "conv": {"median": 0.3234, "n": 2, "per_seed": [0.32, 0.33]},
            "hybrid_conv_maxtx": {"median": 0.4022, "n": 2, "per_seed": [0.40, 0.41]},
            "_paired_hybrid_vs_comparator": _p(0.0054, 7),
            "_paired_control_vs_comparator": _p(0.0001, 4),
        },
        "honest_limits": ["it RANKS and does not dose", "one rung, one combine rule"],
    }
    art.update(over)
    return art


def test_render_md_carries_the_verdict_and_its_detail():
    md = render_md(_art())
    assert "**VERDICT: `WEAK`**" in md
    assert "hybrid +0.0054 over max_tx" in md


def test_render_md_reports_the_orthogonality_precondition_with_its_SPREAD():
    """The median alone would hide a seed where the two scorers were redundant."""
    md = render_md(_art())
    assert "0.6421" in md and "0.5900" in md and "0.7100" in md


def test_render_md_renders_one_row_per_PUBLIC_arm_and_never_a_private_key():
    """Non-vacuous: the `_paired_*` dicts carry no `per_seed`, so a leaked private key would raise
    rather than merely render an odd row. That this call SUCCEEDS is the assertion."""
    md = render_md(_art())
    assert "| `max_tx` |" in md and "| `conv` |" in md and "| `hybrid_conv_maxtx` |" in md
    assert "_paired_hybrid_vs_comparator" not in md.split("## Paired margins")[0]


def test_render_md_omits_a_paired_line_whose_arm_did_not_run():
    """`hybrid+gc` exists only when `gc` is in the table; its absence must not print an empty row."""
    md = render_md(_art())
    assert "hybrid vs max_tx" in md
    assert "CONTROL (shuffled conv) vs max_tx" in md
    assert "hybrid+gc vs max_tx" not in md


def test_render_md_includes_the_three_way_line_when_that_arm_ran():
    a = _art()
    a["arms"]["_paired_hybrid3_vs_comparator"] = _p(0.0070, 8)
    assert "hybrid+gc vs max_tx" in render_md(a)


def test_render_md_always_renders_the_control_line():
    """A positive the control reproduces is an artifact; the control must be visible beside it."""
    assert "CONTROL (shuffled conv)" in render_md(_art(verdict="HYBRID_ADDS"))


def test_render_md_carries_every_honest_limit():
    md = render_md(_art())
    for lim in _art()["honest_limits"]:
        assert lim in md


# --------------------------------------------------------------------------------------------------
# main() -- the refusal and write paths, on synthetic seed tables. No D:, no model, no network.
# --------------------------------------------------------------------------------------------------
def _noisy(n: int, sd: float, seed: int) -> list[float]:
    r = random.Random(seed)
    return [float(i) + r.gauss(0, sd) for i in range(n)]


def _write_seed(d: Path, seed: int, preds: dict[str, list[float]],
                sp: dict[str, float], *, smoke: bool = False, n: int = 40) -> None:
    """One synthetic seed table on disk, built through the script's own validating payload builder."""
    d.mkdir(parents=True, exist_ok=True)
    payload = table_payload(seed=seed, held=["peakA"], tile_indices=list(range(n)),
                            measured=[float(i) for i in range(n)], predictions=preds,
                            spearmans=sp, smoke=smoke)
    with gzip.open(d / f"seed-{seed:02d}.json.gz", "wt", encoding="utf-8") as fh:
        json.dump(payload, fh)


def _full_seed(d: Path, seed: int, *, smoke: bool = False) -> None:
    """A table carrying all three scorers, with seed 0's recomputed Spearmans ON the published
    anchors so the reconcile gate is clean."""
    _write_seed(d, seed,
                {"gc": _noisy(40, 8, 10 + seed), "max_tx": _noisy(40, 6, 20 + seed),
                 "conv": _noisy(40, 7, 30 + seed)},
                {"gc": 0.3000 if seed == 0 else 0.3100,
                 "max_tx": 0.3968 if seed == 0 else 0.4000, "conv": 0.3234},
                smoke=smoke)


def test_main_REFUSES_an_empty_input_directory(tmp_path, monkeypatch):
    """Nothing to combine is not a NO_GAIN. It must name the producer to run first."""
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path)])
    with pytest.raises(SystemExit, match="no seed tables"):
        arm.main()


def test_main_REFUSES_a_table_without_the_scorers_it_combines(tmp_path, monkeypatch, capsys):
    """A gc-only table cannot answer this arm's question; main must exit non-zero and write nothing."""
    _write_seed(tmp_path, 0, {"gc": _noisy(40, 5, 1)}, {"gc": 0.3000})
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1"])
    assert arm.main() == 1
    out = capsys.readouterr().out
    assert "INDETERMINATE_MISSING_SCORER" in out
    assert not list(tmp_path.glob("*.json")), "a refusal must not leave an artifact behind"


def test_main_REFUSES_a_SMOKE_table_at_the_consumer_boundary(tmp_path, monkeypatch):
    """The laundering guard where it matters: a smoke table skipped the reconcile gate, so the entry
    point itself must refuse rather than relying on a caller to have checked."""
    _full_seed(tmp_path, 0, smoke=True)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1"])
    with pytest.raises(TilePredictionError, match="SMOKE"):
        arm.main()


def test_main_writes_both_the_json_and_the_markdown_sidecar(tmp_path, monkeypatch):
    """End to end on two synthetic seeds. The MD sidecar is the only surface a human reads, so a run
    that wrote JSON alone would be a silent half-report."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    assert arm.main() == 0
    assert out.is_file() and out.with_suffix(".md").is_file()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["schema"] == "glm-tile-hybrid-arm-v1"
    assert art["n_seeds"] == 2
    assert art["preregistered"]["margin"] == PREREGISTERED["margin"]
    assert art["arm_name"] == "promcalc-plus-conv-rank-average-hybrid-150bp"
    assert art["verdict"] in PREREGISTERED["verdicts"]
    assert "**VERDICT:" in out.with_suffix(".md").read_text(encoding="utf-8")


def test_main_cannot_reach_HYBRID_ADDS_on_fewer_seeds_than_the_bar_needs(tmp_path, monkeypatch):
    """The bar wants 8 winning seeds of 10. A 2-seed run must therefore not print the headline, even
    when both seeds win — the same under-scoped-run trap the selection script carries a branch for."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["verdict"] != "HYBRID_ADDS"
    assert art["arms"]["_paired_hybrid_vs_comparator"]["n"] == 2


def test_main_records_the_orthogonality_and_the_control_arm(tmp_path, monkeypatch):
    """Both are preconditions for reading the headline, so neither may be absent from the artifact."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    for s in (0, 1):
        _full_seed(tables, s)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "2",
                                      "--out", str(out)])
    arm.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert set(art["orthogonality"]) >= {"per_seed", "median", "min", "max"}
    assert len(art["orthogonality"]["per_seed"]) == 2
    assert "control_shuffled_conv" in art["arms"]
    assert "_paired_control_vs_comparator" in art["arms"]
