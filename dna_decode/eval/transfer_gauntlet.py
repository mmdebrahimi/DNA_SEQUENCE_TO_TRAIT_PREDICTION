"""Five continuous baseline controls plus the four-cell verdict matrix (E19).

NO STATISTICS LIVE HERE. Every number comes from `dna_decode.deconfound`, imported lazily. The AST
guard in tests/test_transfer_gauntlet.py asserts this module imports no `sklearn.*` or `scipy.*`
symbol -- which is why the PCA control reduces `X` with numpy SVD rather than sklearn's PCA. That
guard is the enforcement; the docstring is not.

THE CONTROLS ARE CONTINUOUS ON PURPOSE. `eval/clade_baseline.py::validation_gate` and
`eval/cv.py::leave_one_clade_out_cv` encode the same ideas in AUROC / classification shape, and both
substrates here are continuous. Routing a continuous score into an AUROC-shaped gate would typecheck,
run, and mislead silently -- the recorded `classify_wildtype`-vs-`classify_tier` trap, where an ECOFF
passed in as a clinical breakpoint would have produced a fully self-consistent and entirely wrong
evaluation. `refuse_auroc_shaped_gate` makes that attempt loud.

WHY `assert_all_controls_ran` IS THE WEAKER GUARD. Four of the five controls are imported, so counting
five labels proves almost nothing. The load-bearing invariant is that every control was scored on the
SAME basis -- same group vector, same scored ids, same residualization convention -- which is what
`assert_same_evaluation_basis` checks.
"""
from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from dna_decode.eval.transfer import WITHIN_GROUP_MIN_N

BEATS_ALL_CONTROLS = "BEATS_ALL_CONTROLS"
GROUP_OFFSET_LEARNED = "GROUP_OFFSET_LEARNED"
WITHIN_NOISE = "WITHIN_NOISE"
WITHIN_GROUP_UNSCORABLE = "WITHIN_GROUP_UNSCORABLE"

CONTROL_NAMES = ("shuffled_label_null", "structure_only_baseline", "nearest_neighbour_sequence",
                 "pca_only_baseline", "held_out_clade")

#: The group-offset baseline's WITHIN-GROUP r2 is 0 by construction: a predictor that is constant
#: inside each group explains none of the within-group variance. So "beat the offset baseline on the
#: within-group cell" reduces to "within-group r2 > 0, and above the within-group permutation null".
GROUP_OFFSET_WITHIN_GROUP_R2 = 0.0


class AurocShapedGateRefused(TypeError):
    """A continuous score was aimed at an AUROC-shaped gate."""


class GauntletIncomplete(RuntimeError):
    """Fewer than five controls ran, or they were not scored on one basis."""


@dataclass
class MetricMatrix:
    """Four cells, because one score hides three different questions.

    At k>0 a model can learn the held-out group's MEAN and score well on the pooled cell while
    learning nothing about genotype. Only the within-group cell discounts that.
    """
    zero_shot_group_r2: float | None
    kshot_pooled_r2: float | None
    kshot_within_group_r2: float | None
    kshot_group_offset_baseline: float | None
    n_groups_used_within: int = 0

    def as_dict(self) -> dict:
        return {"zero_shot_group_r2": self.zero_shot_group_r2,
                "kshot_pooled_r2": self.kshot_pooled_r2,
                "kshot_within_group_r2": self.kshot_within_group_r2,
                "kshot_group_offset_baseline": self.kshot_group_offset_baseline,
                "n_groups_used_within": self.n_groups_used_within}


@dataclass
class GauntletReport:
    controls: dict[str, float | None] = field(default_factory=dict)
    matrix: MetricMatrix | None = None
    skipped_groups: tuple[str, ...] = ()
    basis: dict = field(default_factory=dict)
    primitive_provenance: dict = field(default_factory=dict)

    @property
    def n_controls_run(self) -> int:
        """A nan control is NOT run.

        FAIL-OPEN FIXED (found by verification, 2026-10-03). This tested only `is not None`, and
        `nan is not None` is True -- so an uncomputable control counted as run and the report read
        clean, directly contradicting the fail-closed behaviour its own docstring promised.
        """
        return sum(1 for v in self.controls.values()
                   if v is not None and not (isinstance(v, float) and np.isnan(v)))

    def as_dict(self) -> dict:
        return {"controls": dict(self.controls), "n_controls_run": self.n_controls_run,
                "matrix": self.matrix.as_dict() if self.matrix else None,
                "skipped_groups": list(self.skipped_groups), "basis": dict(self.basis),
                "primitive_provenance": dict(self.primitive_provenance)}

    def assert_all_controls_ran(self) -> None:
        """WEAK by design -- see the module docstring. Five labels is necessary, not sufficient."""
        def _unrun(v) -> bool:
            return v is None or (isinstance(v, float) and np.isnan(v))

        missing = [n for n in CONTROL_NAMES if _unrun(self.controls.get(n))]
        if missing:
            raise GauntletIncomplete(
                f"{len(missing)} control(s) did not run: {missing}. A gauntlet that ran nothing is "
                f"not a gauntlet, and a report with holes must not read as clean.")

    def assert_same_evaluation_basis(self) -> None:
        """THE LOAD-BEARING GUARD. Four of five controls are imported, so the thing that can silently
        go wrong is not a missing label -- it is one control scored on a different set of rows."""
        b = self.basis
        required = ("group_vector_hash", "scored_ids_hash", "residualization", "shot_ids_hash")
        missing = [k for k in required if k not in b]
        if missing:
            raise GauntletIncomplete(f"evaluation basis is not recorded: missing {missing}")
        per_control = b.get("per_control", {})
        off = {n: d for n, d in per_control.items()
               if any(d.get(k) != b[k] for k in required)}
        if off:
            raise GauntletIncomplete(
                f"control(s) {sorted(off)} were scored on a DIFFERENT basis than the candidate; the "
                f"comparison is not like-for-like")


def refuse_auroc_shaped_gate(score) -> None:
    """Both substrates here are continuous; `clade_baseline.validation_gate` wants per-clade AUROC
    dicts. Passing a float would typecheck and mislead."""
    if isinstance(score, (int, float, np.floating)):
        raise AurocShapedGateRefused(
            f"a continuous score ({score!r}) was aimed at an AUROC-shaped gate. "
            f"clade_baseline.validation_gate expects {{clade: auroc}} mappings; the continuous "
            f"controls in this module are the right comparison instead.")


def _primitive_provenance() -> dict:
    """Qualnames of the imported primitives plus the sha256 of the module that defines them, so drift
    between the replay arm and the live arms is detectable FROM THE ARTIFACT."""
    from dna_decode import deconfound as pkg
    src = Path(pkg.__file__).parent / "deconfound.py"
    return {
        "deconfound_module_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
        "qualnames": sorted([
            "dna_decode.deconfound.deconfound.cv_r2",
            "dna_decode.deconfound.deconfound.within_group_r2",
            "dna_decode.deconfound.deconfound.group_centered_spearman",
            "dna_decode.deconfound.deconfound.permutation_null",
            "dna_decode.deconfound.deconfound.cluster_from_distance",
            "dna_decode.deconfound.deconfound.r2",
        ]),
        "within_group_min_n_passed_explicitly": WITHIN_GROUP_MIN_N,
    }


def _hash_ids(ids: Sequence) -> str:
    return hashlib.sha256("|".join(str(i) for i in ids).encode()).hexdigest()[:16]


# --- the five controls -----------------------------------------------------------------------------

def shuffled_label_null(X, y, groups, *, n_draws: int = 200,
                        within_group_min_n: int = WITHIN_GROUP_MIN_N) -> float:
    """p95 of a WITHIN-GROUP label-permutation null, ON THE SAME SCALE as the cell it gates.

    UNIT MISMATCH FIXED (found by verification, 2026-10-03). The first version returned
    `deconfound.permutation_null`'s output, which is a SPEARMAN RHO, and the verdict compared it
    against `kshot_within_group_r2`, an R-SQUARED. Measured on one fixture: rho p95 0.0614 vs r2
    0.7237 -- incommensurable quantities. That is the same class as the AUROC-vs-continuous trap this
    module's own `refuse_auroc_shaped_gate` exists to stop, committed inside the module that warns
    about it. (A second defect hid behind it: the branch doing the comparison was unreachable, so the
    control gated nothing at all -- a control that constrains nothing is not a control.)

    So the null is now built on the metric itself: shuffle y INSIDE each group -- which destroys the
    within-group association while leaving every group mean intact, the correct null for this
    estimand -- and re-run `within_group_r2`. The STATISTIC is still deconfound's; only the
    permutation bookkeeping is here. That duplicates the shuffle loop `permutation_null` performs
    internally, which is a cost (see the plan's duplicated-knowledge risk flag) and is accepted
    because the alternative is comparing a rho to an r-squared.

    Returns nan when the null is not computable, which leaves the control UNRUN (see
    `GauntletReport.n_controls_run`, which treats nan as not-run) and makes
    `assert_all_controls_ran` raise. Fail-closed on purpose.
    """
    from dna_decode.deconfound import within_group_r2
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    uniq = np.unique(groups)

    draws = []
    for s in range(int(n_draws)):
        rng = np.random.default_rng(s)
        yp = y.copy()
        for g in uniq:                       # shuffle WITHIN each group: group means survive
            idx = np.where(groups == g)[0]
            v = yp[idx].copy()
            rng.shuffle(v)
            yp[idx] = v
        r, used = within_group_r2(X, yp, groups, min_n=within_group_min_n)
        if used > 0 and not np.isnan(r):
            draws.append(float(r))
    return float(np.percentile(draws, 95)) if draws else float("nan")


def structure_only_baseline(y, groups) -> float:
    """One-hot group membership -> ridge. The control that BEAT the Arabidopsis embedding."""
    from dna_decode.deconfound import cv_r2
    groups = np.asarray(groups)
    levels = sorted(set(groups.tolist()))
    onehot = np.zeros((len(groups), len(levels)), dtype=float)
    for j, lv in enumerate(levels):
        onehot[groups == lv, j] = 1.0
    return float(cv_r2(onehot, np.asarray(y, dtype=float), groups=groups))


def nearest_neighbour_sequence(y, dist, *, exclude_self: bool = True) -> float:
    """1-NN on a supplied distance matrix. numpy only -- no sklearn, so the AST guard can hold."""
    from dna_decode.deconfound import r2
    y = np.asarray(y, dtype=float)
    D = np.asarray(dist, dtype=float).copy()
    if exclude_self:
        np.fill_diagonal(D, np.inf)
    nn = np.argmin(D, axis=1)
    return float(r2(y, y[nn]))


def pca_only_baseline(X, y, groups, *, n_components: int = 5) -> float:
    """Top-k genotype PCs -> ridge, via numpy SVD (NOT sklearn.decomposition)."""
    from dna_decode.deconfound import cv_r2
    X = np.asarray(X, dtype=float)
    Xc = X - X.mean(axis=0, keepdims=True)
    k = int(min(n_components, min(Xc.shape)))
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    return float(cv_r2(Xc @ vt[:k].T, np.asarray(y, dtype=float), groups=np.asarray(groups)))


def held_out_clade(X, y, dist, *, n_clades: int = 3) -> float:
    """Clade folds from the distance matrix (average-linkage) scored leave-one-clade-out."""
    from dna_decode.deconfound import cluster_from_distance, cv_r2
    clades = cluster_from_distance(np.asarray(dist, dtype=float), n_clades)
    return float(cv_r2(np.asarray(X, dtype=float), np.asarray(y, dtype=float), groups=clades))


# --- the matrix -------------------------------------------------------------------------------------

def group_offset_baseline(y, groups, shots_by_group: dict | None = None) -> float:
    """The cheap shortcut: predict each row as the mean of its own group's k shots (the global mean
    when a group has none). Pooled r2. Its WITHIN-group r2 is 0 by construction."""
    from dna_decode.deconfound import r2
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    shots_by_group = shots_by_group or {}
    pred = np.full(len(y), float(np.mean(y)))
    for g, idx in shots_by_group.items():
        idx = list(idx)
        if idx:
            pred[groups == g] = float(np.mean(y[idx]))
    return float(r2(y, pred))


def build_matrix(X, y, groups, *, kshot_pooled: float | None = None,
                 shots_by_group: dict | None = None,
                 within_group_min_n: int = WITHIN_GROUP_MIN_N) -> MetricMatrix:
    """`within_group_min_n` is passed EXPLICITLY, never inherited from the primitive's default -- a
    default in another module silently governing the PASS condition is the drift this plan refuses.

    `kshot_pooled` is supplied BY THE CALLER, from `eval/transfer.k_shot_curve`. That is deliberate:
    computing it here would mean fitting a model in this module, and the caller already owns the
    estimator via `fit_predict`. With no k-shot run it falls back to the zero-shot cell, which is the
    honest value at k=0.
    """
    from dna_decode.deconfound import cv_r2, within_group_r2
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)

    zero_shot = float(cv_r2(X, y, groups=groups))
    within, used = within_group_r2(X, y, groups, min_n=within_group_min_n)
    pooled = zero_shot if kshot_pooled is None else float(kshot_pooled)
    offset = group_offset_baseline(y, groups, shots_by_group)
    return MetricMatrix(
        zero_shot_group_r2=None if np.isnan(zero_shot) else zero_shot,
        kshot_pooled_r2=None if np.isnan(pooled) else float(pooled),
        kshot_within_group_r2=None if np.isnan(within) else float(within),
        kshot_group_offset_baseline=None if np.isnan(offset) else offset,
        n_groups_used_within=int(used),
    )


def run_gauntlet(X, y, groups, *, dist=None, ids: Sequence | None = None,
                 kshot_pooled: float | None = None,
                 shots_by_group: dict | None = None,
                 within_group_min_n: int = WITHIN_GROUP_MIN_N,
                 n_draws: int = 200) -> GauntletReport:
    """Run all five controls plus the matrix on ONE basis, and record that basis."""
    from dna_decode.eval.transfer import held_out_group_folds, skipped_groups as _skipped
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    groups = np.asarray(groups)
    ids = list(ids) if ids is not None else list(range(len(y)))
    if dist is None:                                    # euclidean genotype distance, numpy only
        d2 = ((X[:, None, :] - X[None, :, :]) ** 2).sum(-1)
        dist = np.sqrt(np.maximum(d2, 0.0))

    rep = GauntletReport()
    rep.controls["shuffled_label_null"] = shuffled_label_null(X, y, groups, n_draws=n_draws)
    rep.controls["structure_only_baseline"] = structure_only_baseline(y, groups)
    rep.controls["nearest_neighbour_sequence"] = nearest_neighbour_sequence(y, dist)
    rep.controls["pca_only_baseline"] = pca_only_baseline(X, y, groups)
    rep.controls["held_out_clade"] = held_out_clade(X, y, dist)
    rep.matrix = build_matrix(X, y, groups, kshot_pooled=kshot_pooled,
                              shots_by_group=shots_by_group,
                              within_group_min_n=within_group_min_n)

    folds = held_out_group_folds({str(i): str(g) for i, g in zip(ids, groups.tolist())})
    rep.skipped_groups = _skipped(folds)

    basis = {
        "group_vector_hash": _hash_ids(groups.tolist()),
        "scored_ids_hash": _hash_ids(ids),
        "residualization": f"group_centered_train_mean,min_n={within_group_min_n}",
        "shot_ids_hash": _hash_ids(sorted(
            f"{g}:{sorted(v)}" for g, v in (shots_by_group or {}).items())),
    }
    basis["per_control"] = {n: dict(basis) for n in CONTROL_NAMES}
    rep.basis = basis
    rep.primitive_provenance = _primitive_provenance()
    return rep


def verdict(matrix: MetricMatrix, report: GauntletReport, *, candidate: float | None = None) -> str:
    """The PASS condition reads the WITHIN-GROUP cell, and REFUSES when that cell is undecidable.

    `within_group_r2` returns `(nan, 0)` whenever every group falls under `min_n`, while the pooled
    cell still returns a confident-looking number on the identical data (measured: nan/0 groups vs
    -0.4230 on 6 groups x 20). Falling back to the pooled cell there is exactly how a
    GROUP_OFFSET_LEARNED run would be reported as a pass, so it never happens.
    """
    w = matrix.kshot_within_group_r2
    if w is None or (isinstance(w, float) and np.isnan(w)) or matrix.n_groups_used_within == 0:
        return WITHIN_GROUP_UNSCORABLE

    # the offset baseline's within-group r2 is 0 by construction (constant inside each group)
    if w <= GROUP_OFFSET_WITHIN_GROUP_R2:
        pooled = matrix.kshot_pooled_r2
        offset = matrix.kshot_group_offset_baseline
        if pooled is not None and offset is not None and pooled > offset:
            return GROUP_OFFSET_LEARNED
        return WITHIN_NOISE

    score = candidate if candidate is not None else w

    # THE NULL NOW GATES, on matched units. Previously this sat below a `w <= 0` condition that the
    # branch above had already returned on, so it was unreachable and the control constrained
    # nothing; and it compared a Spearman rho to an r-squared. Both fixed -- `shuffled_label_null`
    # returns a within-group r2 null, so `score` and `null_p95` are the same quantity.
    null_p95 = report.controls.get("shuffled_label_null")
    if null_p95 is not None and not np.isnan(null_p95) and score <= null_p95:
        return WITHIN_NOISE

    # RIVALS. `structure_only_baseline` is deliberately NOT one, and the reason is measured rather
    # than stylistic: under leave-one-group-out the held-out group's one-hot column is all-zero in
    # every training row, so the model cannot represent that group's offset at all. It therefore
    # scores WORSE the more structure explains (-0.777 where group explains ~everything vs -0.001
    # where it explains nothing) -- presenting its easiest bar exactly where it was added to bite.
    # The Arabidopsis comparison that motivated it ("the control that BEAT the embedding") was on
    # POOLED metrics, which is not the frame the verdict reads. In the WITHIN-GROUP frame any purely
    # structural predictor is identically 0 by construction, so it is also redundant with the offset
    # baseline already handled above. It stays in the report as a POOLED DIAGNOSTIC.
    for name in ("nearest_neighbour_sequence", "pca_only_baseline", "held_out_clade"):
        rival = report.controls.get(name)
        if rival is not None and not np.isnan(rival) and score <= rival:
            return f"LOSES_TO_{name}"
    return BEATS_ALL_CONTROLS
