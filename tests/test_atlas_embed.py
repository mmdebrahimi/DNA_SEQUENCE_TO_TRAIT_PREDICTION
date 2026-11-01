"""Guards for the local PromoterAtlas embedding harness.

The load-bearing tests are the two REFUSALS -- ragged input and a truncating width -- because both
failure modes produce a plausible number rather than an error. Everything that needs the real
checkpoint skips when the clone is absent, so this file is green in CI without the 686 MB clone.

**Mask-awareness is tested as a MEASURED property, not as an invariance.** It would be natural to
assert that appending zero columns leaves the pooled vector for the real region unchanged -- that is
what masked pooling is supposed to buy. On this encoder it is FALSE and known to be false: attention
lets the all-zero columns reach the real positions, so the pooled vector moves. The arm measured max
|delta| 1.232 on real tiles. The test therefore asserts the perturbation is NON-ZERO, which is the
fact that demoted the padded-200 geometry to a secondary cell. Asserting the invariance would pin a
belief the project has already refuted.
"""

from __future__ import annotations

import numpy as np
import pytest

from dna_decode.glm.atlas_embed import (
    EXPECTED_PARAMS,
    NUC,
    SAFE_GLOBAL_PREFIXES,
    AtlasEmbedError,
    DEFAULT_PA_REPO,
    audit_checkpoint,
    load_atlas,
    masked_mean_pool,
    ridge_spearman,
)

CLONE = DEFAULT_PA_REPO / "trained_weights/base_model/promoteratlas-base.pt"
SRC = DEFAULT_PA_REPO / "src/promoter_atlas/models/dna_transformer.py"
needs_clone = pytest.mark.skipif(
    not (CLONE.is_file() and SRC.is_file()), reason="PromoterAtlas clone absent on this host"
)


# ---------------------------------------------------------------- offline contract


def test_nucleotide_map_matches_the_pretraining_channel_order():
    """A permuted channel map would embed every sequence as a different sequence and still return a
    (N,128) matrix, so this is pinned rather than trusted."""
    assert NUC == {"A": 0, "C": 1, "G": 2, "T": 3}


def test_expected_param_count_is_the_arms_pinned_value():
    assert EXPECTED_PARAMS == 1_847_428


def test_allowlist_does_not_admit_arbitrary_callables():
    for prefix in SAFE_GLOBAL_PREFIXES:
        assert not prefix.startswith(("os", "subprocess", "builtins.eval", "builtins.exec"))
    assert not any(p in ("builtins", "") for p in SAFE_GLOBAL_PREFIXES)


def test_absent_clone_raises_a_named_error_not_an_import_error(tmp_path):
    with pytest.raises(AtlasEmbedError) as ei:
        load_atlas(tmp_path / "no-such-clone")
    assert "absent" in str(ei.value)


def test_audit_refuses_a_non_zip_file(tmp_path):
    p = tmp_path / "not-a-checkpoint.pt"
    p.write_bytes(b"not a zip")
    with pytest.raises(Exception):
        audit_checkpoint(p)


# ---------------------------------------------------------------- refusals (no model needed)


class _Stub:
    """Stands in for the encoder so the input refusals are testable with no checkpoint.

    It is only ever reached AFTER the geometry checks, so a test that gets a result from it has
    proven the refusal did not fire.
    """

    def __call__(self, X, return_latent=False):  # pragma: no cover - refusals fire first
        import torch

        n, _, w = X.shape
        return None, None, torch.zeros((n, 128, w))


def test_ragged_input_is_REFUSED_not_padded():
    """The decisive refusal. Padding a ragged cohort would silently produce a plausible number, and
    on this encoder padding is measured to perturb the embedding."""
    with pytest.raises(AtlasEmbedError) as ei:
        masked_mean_pool(_Stub(), ["ACGT", "ACG"])
    assert "refuses to pad" in str(ei.value)


def test_width_shorter_than_the_sequence_is_REFUSED_not_truncated():
    with pytest.raises(AtlasEmbedError) as ei:
        masked_mean_pool(_Stub(), ["ACGTACGT"], width=4)
    assert "truncate" in str(ei.value)


def test_empty_input_raises_rather_than_returning_an_empty_matrix():
    with pytest.raises(AtlasEmbedError):
        masked_mean_pool(_Stub(), [])


def test_native_geometry_is_the_default():
    """`width=None` must mean len(seq); a default of 200 would silently run the demoted cell."""
    out = masked_mean_pool(_Stub(), ["ACGTACGT"])
    assert out.shape == (1, 128)


# ---------------------------------------------------------------- ridge protocol


def test_ridge_spearman_mirrors_fit_predict_ridge_on_the_same_input():
    """The incumbents' protocol is StandardScaler -> Ridge(alpha=1.0) -> fit on log1p -> score
    against RAW. If these two drift apart, every local number becomes incomparable with the
    committed tables, so they are compared directly on one feature matrix."""
    from dna_decode.glm.expression import Pair, fit_predict_ridge

    rng = np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), 20)) for _ in range(80)]
    expr = [float(s.count("G") + s.count("C")) + float(rng.normal(0, 0.3)) for s in seqs]
    pairs = [
        Pair(name="p%d" % i, seq=s, expression=max(e, 0.0))
        for i, (s, e) in enumerate(zip(seqs, expr))
    ]
    tr, te = pairs[:60], pairs[60:]

    ref = fit_predict_ridge(tr, te, feature="gc")

    # The matrix is built with the EXACT function fit_predict_ridge dispatches to. Using
    # genomewide.gc_features instead would pass today (the two are behaviourally identical) and
    # would silently stop testing parity the moment either drifts.
    from dna_decode.glm.expression import FEATURES

    Xtr = np.asarray(FEATURES["gc"]([p.seq for p in tr]), dtype=float)
    Xte = np.asarray(FEATURES["gc"]([p.seq for p in te]), dtype=float)
    got, _ = ridge_spearman(
        Xtr, [p.expression for p in tr], Xte, [p.expression for p in te]
    )
    # fit_predict_ridge returns round(spearman, 4) (expression.py:374) while ridge_spearman returns
    # full precision, so parity is asserted by EXACT equality at the reference's own precision --
    # stronger than widening a tolerance, which would also pass a real protocol divergence.
    assert round(got, 4) == ref["spearman"], (got, ref["spearman"])


def test_shuffle_seed_changes_the_result_so_the_null_is_not_vacuous():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(120, 5))
    y = list(np.abs(X[:, 0]) * 3 + 1)
    real, _ = ridge_spearman(X[:90], y[:90], X[90:], y[90:])
    null, _ = ridge_spearman(X[:90], y[:90], X[90:], y[90:], shuffle_seed=0)
    assert real != null


# ---------------------------------------------------------------- real checkpoint


@needs_clone
def test_real_checkpoint_loads_with_the_pinned_param_count():
    model = load_atlas()
    assert sum(p.numel() for p in model.parameters()) == EXPECTED_PARAMS


@needs_clone
def test_real_checkpoint_audit_finds_only_allowlisted_globals():
    globs = audit_checkpoint(CLONE)
    assert globs, "an audit that finds NO globals has not parsed the pickle"
    for g in globs:
        assert any(g == p or g.startswith(p) for p in SAFE_GLOBAL_PREFIXES), g


@needs_clone
def test_real_embedding_is_deterministic_across_two_passes():
    model = load_atlas()
    probe = ["ACGT" * 37 + "AA"] * 4
    a = masked_mean_pool(model, probe)
    b = masked_mean_pool(model, probe)
    assert float(np.abs(a - b).max()) == 0.0


@needs_clone
def test_padding_PERTURBS_the_pooled_embedding_which_is_why_200_was_demoted():
    """Measured, not assumed: masked pooling does NOT undo zero-padding on an attention encoder."""
    model = load_atlas()
    s = ["ACGT" * 37 + "AA"]
    native = masked_mean_pool(model, s)
    padded = masked_mean_pool(model, s, width=200)
    assert float(np.abs(native - padded).max()) > 0.1


@needs_clone
def test_real_embedding_is_not_degenerate():
    model = load_atlas()
    rng = np.random.default_rng(7)
    seqs = ["".join(rng.choice(list("ACGT"), 150)) for _ in range(24)]
    X = masked_mean_pool(model, seqs)
    assert X.shape == (24, 128)
    assert int((X.var(axis=0) < 1e-8).sum()) < 128, "every dimension constant: encoder is dead"
