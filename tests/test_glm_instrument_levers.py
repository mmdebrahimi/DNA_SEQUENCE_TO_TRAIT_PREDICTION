"""The two instrument levers, and the one result that explains the whole kmer-versus-both story.

Two findings are pinned here because each reverses something the arm had been assuming.

1. **`positional` alone matches `both`.** The committed grid scores only `kmer` and `both`, and `both` is
   the concatenation of the two blocks, so the positional block was always the natural explanation for
   `both` beating `kmer` -- and was never run on its own. Seed-paired at the committed order,
   `both - positional` is -0.0364 with an interval containing zero, while `positional - kmer` is +0.1411
   and resolvable. The 256 k-mer features add nothing resolvable to the 96 positional ones.

2. **The null's ORDER sets the kmer margin and leaves positional alone.** A Markov chain's sufficient
   statistic IS the k-mer conditional distribution, so an order-3 chain saturates a 4-mer discriminator by
   construction: kmer's margin falls from +0.3259 at order 1 to +0.0329 at order 3. A positional
   discriminator measures WHERE composition sits, which no order of Markov chain models, so its margin is
   flat across orders (range 0.0344, nothing resolvably different from the committed order). That is the
   mechanism behind the kmer mode's small margin, and it means the small margin is the CORRECT answer
   rather than a defect of the mode.

`positional_features`' own docstring already argued point 2 from first principles. Nobody had measured the
order dependence that tests it, which is why these tests pin the measurement and not the argument.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_instrument_levers.py"
OUT = ROOT / "wiki" / "glm_instrument_levers_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_levers", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the feature-count clamp

def test_the_positional_block_respects_the_k_clamp_inside_the_builder(mod):
    """`positional_features` is called with k=min(k,2), so the block is bins*16 and NOT bins*256.

    Getting this wrong would misreport every p-against-n figure in the artifact, and the clamp is easy to
    miss because it lives inside the builder rather than at the call site.
    """
    assert mod.feature_count("kmer") == 256
    assert mod.feature_count("positional") == 6 * 16 == 96
    assert mod.feature_count("both") == 256 + 96 == 352


def test_the_clamp_matches_what_build_features_actually_produces(mod):
    """Non-vacuity: compare against the real builder rather than against my arithmetic."""
    from dna_decode.glm.falsifier import build_features

    seqs = ["ACGT" * 40, "TGCA" * 40, "AACC" * 40]
    for mode in ("kmer", "positional", "both"):
        rows = build_features(seqs, mode=mode, k=mod.KMER_K, bins=mod.POSITIONAL_BINS)
        assert len(rows[0]) == mod.feature_count(mode), mode


# ------------------------------------------------------------------ the paired comparison

def test_the_mode_comparison_is_seed_paired(mod):
    """Two modes share the corpus draw and split at each seed, so differencing within a seed removes
    the draw-to-draw variation that dominates either arm alone. An unpaired comparison at n=10 would
    hide a real difference behind corpus noise."""
    a = [0.20, 0.25, 0.15, 0.22, 0.18, 0.24, 0.16, 0.21, 0.19, 0.23]
    b = [v - 0.05 for v in a]
    out = mod.paired_diff(a, b)
    assert out["mean_diff"] == pytest.approx(0.05, abs=1e-6)
    assert out["sd_diff"] == pytest.approx(0.0, abs=1e-9)
    assert out["resolvably_nonzero"] is True
    assert out["n_positive"] == 10


def test_a_zero_variance_difference_is_resolvable_not_unresolvable(mod):
    """Exactly-zero spread, built so no floating-point residue can rescue it.

    The version of this test above subtracts a constant, which leaves se tiny-but-nonzero rather than
    exactly zero -- and it therefore PASSED against a `bool(se and ...)` guard that reported sd == 0 as
    unresolvable, the exact inversion of the truth. Integer-exact inputs make the degenerate branch
    actually execute.
    """
    out = mod.paired_diff([0.10] * 10, [0.20] * 10)
    assert out["sd_diff"] == 0.0
    assert out["mean_diff"] == pytest.approx(-0.10, abs=1e-9)
    assert out["resolvably_nonzero"] is True
    out_tied = mod.paired_diff([0.10] * 10, [0.10] * 10)
    assert out_tied["sd_diff"] == 0.0 and out_tied["mean_diff"] == 0.0
    assert out_tied["resolvably_nonzero"] is False, "a zero difference must NOT read as resolvable"


def test_a_paired_difference_of_zero_is_not_resolvable(mod):
    a = [0.20, 0.25, 0.15, 0.22, 0.18, 0.24, 0.16, 0.21, 0.19, 0.23]
    b = [0.21, 0.24, 0.16, 0.21, 0.19, 0.23, 0.17, 0.20, 0.20, 0.22]
    out = mod.paired_diff(a, b)
    assert out["resolvably_nonzero"] is False


def test_paired_diff_detects_a_negative_difference_too(mod):
    """One-sided logic here would call a resolvable NEGATIVE difference 'not resolvable', which is
    exactly the case that matters for `both - positional`."""
    a = [0.10] * 10
    b = [0.20] * 10
    out = mod.paired_diff(a, b)
    assert out["mean_diff"] < 0 and out["resolvably_nonzero"] is True


def test_resolve_uses_the_standard_error(mod):
    m = [0.17, 0.18, 0.16, 0.175, 0.165, 0.17, 0.18, 0.17, 0.17, 0.175]
    out = mod.resolve(m)
    assert out["stderr"] == pytest.approx(out["sd"] / math.sqrt(10), abs=1e-4)


# ------------------------------------------------------------------ the design

def test_the_sweep_covers_the_committed_cell_and_the_unmeasured_mode(mod):
    assert mod.COMMITTED_ORDER == 3
    assert "positional" in mod.MODES, "the mode the committed grid never scored"
    assert mod.COMMITTED_ORDER in mod.ORDERS
    assert len(mod.ORDERS) >= 4, "one order on each side of the committed value is not a sweep"


def test_the_floor_is_computed_once_per_mode_and_seed(mod):
    """The floor is real-against-real and never involves the null, so it CANNOT depend on the null's
    order. Recomputing it per order would add noise and invite a reader to think order moves the floor."""
    src = SCRIPT.read_text(encoding="utf-8")
    assert "floors: dict[tuple[str, int], float]" in src
    assert "order-independent" in src or "order_independent" in src


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_floor_is_identical_across_orders_within_a_mode():
    """The structural guarantee, checked rather than trusted."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "positional", "both"):
        floors = {rec["cells"][f"{mode}_order{o}"]["floor_mean"] for o in rec["geometry"]["seeds"][:1]
                  for o in (1, 2, 3, 4, 5)}
        assert len(floors) == 1, (mode, floors)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_positional_alone_matches_both_and_beats_kmer():
    """Finding 1, pinned. The k-mer block is 256 of 352 features and adds nothing resolvable."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    q1 = rec["q1_feature_blocks"]
    assert q1["both_minus_positional"]["resolvably_nonzero"] is False
    assert q1["positional_minus_kmer"]["resolvably_nonzero"] is True
    assert q1["positional_minus_kmer"]["mean_diff"] > 0
    assert q1["verdict"].startswith("POSITIONAL_ALONE_MATCHES_BOTH")


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_kmer_margin_collapses_with_null_order_and_positional_does_not():
    """Finding 2, pinned -- the mechanism behind the whole kmer-versus-both story.

    A Markov chain's sufficient statistic IS the k-mer conditional distribution, so raising its order
    saturates a k-mer discriminator by construction while leaving a positional one untouched.
    """
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    q2 = rec["q2_null_order"]
    assert q2["kmer"]["range"] > 0.2, q2["kmer"]["range"]
    assert q2["positional"]["range"] < 0.06, q2["positional"]["range"]
    assert q2["kmer"]["range"] > 5 * q2["positional"]["range"]
    # and positional is not merely flatter -- no order is resolvably different from the committed one
    assert q2["positional"]["n_orders_resolvably_different_from_committed"] == 0
    assert q2["kmer"]["n_orders_resolvably_different_from_committed"] >= 2


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_kmer_margin_is_large_at_a_weak_null():
    """Rules out 'the kmer discriminator is simply blind': against an order-1 chain it reads +0.33.

    So the small margin at order 3 is a property of the NULL matching the features, not of the mode
    lacking power -- which is what makes 'retire the mode' the wrong conclusion.
    """
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["cells"]["kmer_order1"]["margin"]["mean"] > 0.3
    assert rec["cells"]["kmer_order1"]["margin"]["resolvable_above_zero"] is True


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_overfitting_alternative_for_high_order_nulls_is_declared():
    """A high-order chain fitted on 600 windows approaches memorising them, so part of any order effect
    is overfitting the fit half rather than modelling the substrate. The sweep does not separate those."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "memorising" in joined or "memorizing" in joined
    assert "does NOT separate" in joined or "not separate" in joined


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_new_cells_are_not_presented_as_a_reread_of_the_committed_grid():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    if "note_on_the_committed_grid" in rec:
        assert "NEW cells" in rec["note_on_the_committed_grid"]


# ------------------------------------------------------------------ the structural cross-check

def test_the_structural_prediction_is_registered_with_its_falsifier(mod):
    """The mechanism claim is strong, so the attempt to break it must name what breaking looks like."""
    p = mod.STRUCTURAL_PREDICTION
    assert p["registered_before_the_second_substrate_was_run"] is True
    assert "sufficient statistic" in p["claim"]
    assert "falsified_if" in p and len(p["falsified_if"]) > 40
    src = SCRIPT.read_text(encoding="utf-8")
    assert "ORDER_DEPENDENCE_DID_NOT_REPLICATE__MECHANISM_CLAIM_UNSUPPORTED" in src, (
        "the failing branch must be reachable, or the cross-check can only ever confirm")


def test_the_second_substrate_is_chosen_for_being_different(mod):
    """A second genome that resembled the first would not test whether the pattern is structural."""
    assert "Campylobacter" in mod.SECOND_LABEL
    assert "0.266" in mod.SECOND_LABEL and "0.453" in mod.REF_LABEL


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_order_dependence_replicates_in_a_very_different_genome():
    """THE mechanism test. If a Markov chain's sufficient statistic IS the k-mer distribution, the
    saturation is a property of the two objects and must reproduce at GC 0.266 as it does at 0.453."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    cc = rec["structural_cross_check_second_substrate"]
    assert cc is not None, "the cross-check did not run; the mechanism claim is then unsupported"
    assert cc["prediction_held"] is True
    assert cc["verdict"] == "ORDER_DEPENDENCE_IS_STRUCTURAL__REPLICATES_ACROSS_SUBSTRATES"
    prof = cc["profile"]
    assert prof["kmer"]["range"] > 5 * prof["positional"]["range"]
    assert prof["positional"]["n_orders_resolvably_different_from_committed"] == 0
    assert prof["kmer"]["n_orders_resolvably_different_from_committed"] >= 2


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_both_substrates_show_the_same_shape_not_the_same_numbers():
    """The claim is about SHAPE. Demanding equal margins across two organisms would be the wrong test
    -- the cross-organism probe already measured that Campylobacter carries a LARGER margin."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    primary = rec["q2_null_order"]
    second = rec["structural_cross_check_second_substrate"]["profile"]
    for mode in ("kmer", "both"):
        assert primary[mode]["range"] > 0.15 and second[mode]["range"] > 0.15, mode
    assert primary["positional"]["range"] < 0.06 and second["positional"]["range"] < 0.06
    # and the margins are genuinely NOT equal, so this is not accidentally a sameness test
    assert primary["kmer"]["range"] != second["kmer"]["range"]


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_two_substrate_limit_is_declared():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "do not characterise how it varies" in joined or "not characterise" in joined
    assert "ONE Campylobacter" in joined or "single draft" in joined
