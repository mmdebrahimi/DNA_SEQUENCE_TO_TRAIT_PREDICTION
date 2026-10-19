"""Genotype-derived population structure: IBS distance, a derived clustering threshold, and
the two guards that make a grouping honest (F3 Step 3).

Pure, offline, no Docker, no network. Mash is deliberately NOT used here: these are dense
genotype matrices, and IBS over the markers is the DIRECT measurement where a Mash sketch
would be a proxy for it.

The threshold is DERIVED from each panel's own distance distribution rather than swept over a
grid. That is a reaction to a measured defect in this repo's Mash work: `pick_best_threshold`
minimised `variance_ratio = intra/inter`, which tends to 0 as clusters tend to singletons, so
MINIMIZE was monotonically biased toward over-splitting -- extending its grid downward picked
0.002 and produced 101 clades with 56.5% singletons. Its qualifying criteria guarded only
against a DOMINANT cluster, never against over-splitting. `panel_groups` enforces BOTH.

MEASURED CAVEAT -- A PERCENTILE THRESHOLD ALWAYS FINDS SOMETHING, AND THE GUARDS CANNOT SEE IT.
On an i.i.d. binary panel with NO structure at all, the 10th-percentile cut still produced 19
well-shaped clusters over 50 individuals (largest 14%), clearing both guards. By construction
~`percentile`% of pairs fall below any quantile whether or not any structure exists, so a
clean-looking grouping is NOT evidence that a panel is structured. Two consequences:

  1. `pc1_variance_fraction` -- which has no free parameter -- is the PRIMARY structure
     statistic. Cluster shape is corroboration, never the existence claim.
  2. This is a FEATURE for the F3 negative control, not a bug to patch out. On a bi-parental
     cross ancestry is randomized by construction, so an arbitrary grouping is precisely the
     thing the negative control must exercise; refusing on a gapless distribution would make
     that control unrunnable. The calibration's verdict is therefore DIFFERENTIAL (does the
     de-confounder respond MORE where structure is measured to exist) rather than absolute.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
from dna_decode.eval.phylogeny import DistanceMatrix

# --- asserted constants (NOT derived -- stated so a reader is not misled) -------------------
# A single fixed percentile is WRONG IN BOTH DIRECTIONS, measured: at 10.0 an i.i.d. panel gets
# 19 manufactured clusters while a genuinely structured two-block panel (pc1 0.894) is REFUSED
# for 84% singletons, because 0.02 sits below typical within-block distances. So the threshold
# is SWEPT over this ladder of the panel's own distance percentiles and selected by the guards.
# This is the fix CLAUDE.md prescribes for the Mash threshold defect -- "add an over-split
# guard, THEN extend the grid" -- in that order, which is what keeps the sweep from inheriting
# the monotone bias toward singletons that sank `pick_best_threshold`.
IBS_PERCENTILE_LADDER = (1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0)
# Retained as the ladder's historical starting rung; `panel_groups` no longer defaults to it.
IBS_THRESHOLD_PERCENTILE = 10.0
# A cluster holding more than this much of the panel makes grouped CV leave-a-large-fraction-out.
MAX_LARGEST_CLUSTER_FRACTION = 0.60
# More singletons than this means the "grouping" barely groups -- it approximates random CV
# while carrying the authority of a de-confounded number. This is the guard the Mash work lacked.
MAX_SINGLETON_FRACTION = 0.50


class PanelStructureError(ValueError):
    """A panel whose structure cannot support an honest grouping.

    Raised rather than returned: a grouping that silently degrades toward random folds would
    report a de-confounded number that was never de-confounded.
    """


@dataclass
class PanelGroups:
    ids: list[str]
    groups: list[int]                 # aligned to `ids`
    n_individuals: int
    n_clusters: int
    largest_cluster_fraction: float
    singleton_fraction: float
    threshold: float
    percentile: float                 # the ladder rung the guards selected
    ibs_path: str                     # "binary_exact" | "general_loop"
    missing_fraction: float
    sweep: list = None                # every rung + why the rejected ones failed

    def as_dict(self) -> dict:
        return {
            "n_individuals": self.n_individuals,
            "n_clusters": self.n_clusters,
            "largest_cluster_fraction": round(self.largest_cluster_fraction, 4),
            "singleton_fraction": round(self.singleton_fraction, 4),
            "threshold": round(self.threshold, 6),
            "percentile": self.percentile,
            "ibs_path": self.ibs_path,
            "missing_fraction": round(self.missing_fraction, 4),
            "sweep": [
                {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                 if k != "refusal"}
                for r in (self.sweep or [])
            ],
        }


def _as_matrix(X) -> np.ndarray:
    A = np.asarray(X, dtype=float)
    if A.ndim != 2:
        raise PanelStructureError(f"genotype matrix must be 2-D, got shape {A.shape}")
    if A.shape[0] < 2:
        raise PanelStructureError(f"need >=2 individuals, got {A.shape[0]}")
    if A.shape[1] < 1:
        raise PanelStructureError("genotype matrix has no markers")
    return A


def _binary_recode(A: np.ndarray) -> np.ndarray | None:
    """Return a {0,1} recoding when the observed values admit one, else None.

    Covers the two codings real marker files here actually use: {0,1} and {-1,+1} (the Bloom
    segregant panel ships B/R parent-of-origin calls). A dosage coding {0,1,2} is NOT binary
    and falls through to the exact general path.
    """
    obs = A[np.isfinite(A)]
    if obs.size == 0:
        return None
    vals = np.unique(obs)
    if vals.size > 2:
        return None
    if vals.size == 1:
        return np.where(np.isfinite(A), 0.0, np.nan)
    lo, hi = float(vals[0]), float(vals[1])
    return np.where(np.isfinite(A), (A - lo) / (hi - lo), np.nan)


@dataclass
class IBSResult:
    """The distance matrix plus which code path produced it (auditable in an artifact)."""

    matrix: DistanceMatrix
    path: str                         # "binary_exact" | "general_loop"


def ibs_distance_matrix(X, ids) -> DistanceMatrix:
    """Pairwise identity-by-state distance; see `ibs_distance_result` for the audit field."""
    return ibs_distance_result(X, ids).matrix


def ibs_distance_result(X, ids) -> IBSResult:
    """Pairwise identity-by-state distance: the fraction of JOINTLY OBSERVED markers at which
    two individuals differ.

    Missing-aware: NaN markers are excluded pairwise, so a pair sharing no observed marker
    yields a non-finite distance, which `derive_threshold` refuses rather than silently
    treating as distance 0 (identical).
    """
    A = _as_matrix(X)
    ids = [str(i) for i in ids]
    if len(ids) != A.shape[0]:
        raise PanelStructureError(f"{len(ids)} ids for {A.shape[0]} rows")
    if len(set(ids)) != len(ids):
        raise PanelStructureError("duplicate ids: the clusterer keys on id, so a duplicate "
                                  "would silently drop an individual")

    B = _binary_recode(A)
    if B is not None:
        # Exact for a binary coding via |a-b| = a + b - 2ab, which turns the whole O(n^2 m)
        # pairwise computation into four matmuls (BLAS). Not an approximation.
        M = np.isfinite(B).astype(np.float64)
        Z = np.where(np.isfinite(B), B, 0.0)
        num = (Z @ M.T) + (M @ Z.T) - 2.0 * (Z @ Z.T)
        den = M @ M.T
        with np.errstate(invalid="ignore", divide="ignore"):
            D = num / den
        D[den == 0] = np.nan
        path = "binary_exact"
    else:
        n = A.shape[0]
        rng = float(np.nanmax(A) - np.nanmin(A))
        if not np.isfinite(rng) or rng == 0:
            raise PanelStructureError("genotype matrix has no variation (zero value range)")
        D = np.full((n, n), np.nan, dtype=np.float64)
        for i in range(n):
            with np.errstate(invalid="ignore"):
                D[i] = np.nanmean(np.abs(A[i] - A), axis=1) / rng
        path = "general_loop"

    D = (D + D.T) / 2.0          # enforce exact symmetry against float asymmetry
    np.fill_diagonal(D, 0.0)
    return IBSResult(matrix=DistanceMatrix(strain_ids=ids, matrix=D), path=path)


def derive_threshold(distance_matrix: DistanceMatrix,
                     percentile: float = IBS_THRESHOLD_PERCENTILE) -> float:
    """Threshold from the panel's OWN off-diagonal distance distribution.

    Refuses on a non-finite distance: a pair with no jointly observed marker is unknown, not
    close, and letting NaN through would have the clusterer compare NaN <= threshold (False)
    and quietly treat an unmeasurable pair as distant.
    """
    m = np.asarray(distance_matrix.matrix, dtype=float)
    n = m.shape[0]
    off = m[~np.eye(n, dtype=bool)]
    if off.size == 0:
        raise PanelStructureError("no off-diagonal distances")
    n_bad = int(np.sum(~np.isfinite(off)))
    if n_bad:
        raise PanelStructureError(
            f"{n_bad} of {off.size} off-diagonal distances are non-finite (pairs sharing no "
            "observed marker); an unmeasurable pair is not a close pair"
        )
    return float(np.percentile(off, percentile))


def pc1_variance_fraction(X) -> float:
    """Leading-eigenvalue share of the genotype covariance -- the parameter-free structure
    statistic. High on a structured panel, low on one whose ancestry is randomized.

    Missing values are mean-imputed per marker (the standard PCA convention here); a marker
    observed nowhere is dropped.
    """
    A = _as_matrix(X).copy()
    with warnings.catch_warnings():
        # an entirely-missing marker is expected and dropped two lines below; the "Mean of
        # empty slice" warning is noise, and a noisy suite trains readers to discount noise
        warnings.simplefilter("ignore", RuntimeWarning)
        col_mean = np.nanmean(np.where(np.isfinite(A), A, np.nan), axis=0)
    keep = np.isfinite(col_mean)
    if not keep.any():
        raise PanelStructureError("every marker is entirely missing")
    A = A[:, keep]
    col_mean = col_mean[keep]
    bad = ~np.isfinite(A)
    if bad.any():
        A[bad] = np.take(col_mean, np.where(bad)[1])
    Ac = A - A.mean(axis=0, keepdims=True)
    if not np.any(np.abs(Ac) > 0):
        raise PanelStructureError("genotype matrix has no variation after centering")
    s = np.linalg.svd(Ac, full_matrices=False, compute_uv=False)
    var = s ** 2
    total = float(var.sum())
    if total <= 0:
        raise PanelStructureError("zero total genotypic variance")
    return float(var[0] / total)


def cluster_shape(groups) -> tuple[int, float, float]:
    """(n_clusters, largest_cluster_fraction, singleton_fraction) -- pure, no thresholds.

    Report cluster STRUCTURE, not just count: this repo's TB work found 43 "discordant
    clusters" that between them hid 97% of the isolates, so a count alone is uninformative.
    """
    g = list(groups)
    if not g:
        raise PanelStructureError("empty grouping")
    sizes: dict[int, int] = {}
    for x in g:
        sizes[x] = sizes.get(x, 0) + 1
    n_clusters = len(sizes)
    largest = max(sizes.values()) / len(g)
    singletons = sum(1 for v in sizes.values() if v == 1) / n_clusters
    return n_clusters, largest, singletons


def assert_usable_shape(n_clusters: int, largest: float, singletons: float,
                        threshold: float, max_largest: float | None = None) -> None:
    """REFUSE a grouping that fails in either direction. Pure, so both guards are directly
    testable without having to induce them through a random panel.

    `max_largest` tightens the dominant-cluster bound. A CONSUMER needs this: a grouping can
    clear the panel-level 60% bound and still be unusable for k-fold CV, where a group above
    n/k makes one fold larger than its share. Passing 1/n_splits here surfaces that at
    GROUPING time instead of letting `genomic_prediction._resolve_folds` raise later, which
    would look like a plumbing failure rather than a property of the panel.
    """
    limit = MAX_LARGEST_CLUSTER_FRACTION if max_largest is None else float(max_largest)
    if largest > limit:
        raise PanelStructureError(
            f"largest cluster holds {largest:.1%} of the panel (limit "
            f"{limit:.0%}) at threshold {threshold:.6g}; grouped CV "
            f"would become leave-{largest:.0%}-out"
        )
    if singletons > MAX_SINGLETON_FRACTION:
        raise PanelStructureError(
            f"{singletons:.1%} of {n_clusters} clusters are singletons (limit "
            f"{MAX_SINGLETON_FRACTION:.0%}) at threshold {threshold:.6g}; this grouping "
            "approximates random CV while carrying the authority of a de-confounded number"
        )


def sweep_thresholds(distance_matrix: DistanceMatrix, ids,
                     ladder=IBS_PERCENTILE_LADDER,
                     max_largest: float | None = None) -> list[dict]:
    """Cluster shape at every rung of the percentile ladder. Pure; no selection, no refusal."""
    ids_s = [str(i) for i in ids]
    rows = []
    for p in ladder:
        t = derive_threshold(distance_matrix, p)
        assigned = greedy_representative_clusters_from_matrix(distance_matrix, t)
        groups = [assigned[i] for i in ids_s]
        n_c, largest, singles = cluster_shape(groups)
        try:
            assert_usable_shape(n_c, largest, singles, t, max_largest)
            passes, why = True, None
        except PanelStructureError as e:
            passes, why = False, str(e)
        rows.append({"percentile": p, "threshold": t, "n_clusters": n_c,
                     "largest_cluster_fraction": largest, "singleton_fraction": singles,
                     "passes": passes, "refusal": why})
    return rows


def select_threshold(distance_matrix: DistanceMatrix, ids,
                     ladder=IBS_PERCENTILE_LADDER,
                     max_largest: float | None = None) -> tuple[float, float, list[dict]]:
    """Pick the SMALLEST ladder rung whose cluster shape clears BOTH guards.

    PRE-COMMITTED, and its bias is named rather than hidden: the smallest passing rung
    UNDER-merges, so the grouping sits closer to random CV than a coarser one would. For the
    F3 differential criterion that biases AGAINST detecting a positive-control response --
    against the instrument -- so a response that survives this choice is not an artifact of it.
    (The opposite rule, picking the coarsest passing rung, would flatter the instrument.)

    Returns (percentile, threshold, sweep). Refuses when no rung passes, with the whole sweep
    in the message so the failure is diagnosable rather than a bare "no threshold found".
    """
    sweep = sweep_thresholds(distance_matrix, ids, ladder, max_largest)
    for row in sweep:
        if row["passes"]:
            return float(row["percentile"]), float(row["threshold"]), sweep
    detail = "; ".join(
        f"p{r['percentile']:g}: {r['n_clusters']}cl largest={r['largest_cluster_fraction']:.2f} "
        f"single={r['singleton_fraction']:.2f}" for r in sweep
    )
    raise PanelStructureError(
        f"no percentile on the ladder yields a usable grouping -- {detail}. The panel has no "
        "scale at which clusters are neither dominant nor mostly singletons."
    )


def panel_groups(X, ids, *, percentile: float | None = None,
                 ladder=IBS_PERCENTILE_LADDER,
                 max_largest: float | None = None) -> PanelGroups:
    """Derive a grouping and REFUSE if it is unusable in either direction.

    `percentile=None` (default) SWEEPS the ladder and selects by the guards; an explicit value
    forces that single rung (used to drive each guard deterministically in tests).

    Reuses `clonality.greedy_representative_clusters_from_matrix`, which is
    chaining-resistant. The alternative in this repo (`phylogeny.cluster_by_ani`, single
    linkage) is RECORDED as defective: it chains A~B~C into one cluster and produced a 57.8%
    largest cluster where greedy-representative gave 31.3% on the same data.
    """
    res = ibs_distance_result(X, ids)
    if percentile is None:
        chosen_p, threshold, sweep = select_threshold(res.matrix, ids, ladder, max_largest)
    else:
        chosen_p, threshold, sweep = percentile, derive_threshold(res.matrix, percentile), []
    assigned = greedy_representative_clusters_from_matrix(res.matrix, threshold)

    ids_s = [str(i) for i in ids]
    n = len(ids_s)
    groups = [assigned[i] for i in ids_s]
    n_clusters, largest, singletons = cluster_shape(groups)

    A = _as_matrix(X)
    missing = float(np.mean(~np.isfinite(A)))

    assert_usable_shape(n_clusters, largest, singletons, threshold, max_largest)

    return PanelGroups(
        ids=ids_s, groups=groups, n_individuals=n, n_clusters=n_clusters,
        largest_cluster_fraction=largest, singleton_fraction=singletons,
        threshold=threshold, percentile=chosen_p, ibs_path=res.path,
        missing_fraction=missing, sweep=sweep,
    )


def structure_census(X, ids, *, percentile: float | None = None) -> dict:
    """Structure statistics WITHOUT the refusals -- for a census that must describe every panel,
    including the unusable ones, rather than stopping at the first.

    `pc1_variance_fraction` is the primary statistic (no free parameter); the cluster shape is
    corroboration and carries `groups_refused` with the reason when the guards would refuse.
    """
    out: dict = {"pc1_variance_fraction": None, "pc1_error": None,
                 "groups": None, "groups_refused": None}
    try:
        out["pc1_variance_fraction"] = pc1_variance_fraction(X)
    except Exception as e:                        # noqa: BLE001 - census must not stop
        out["pc1_error"] = f"{type(e).__name__}: {e}"
    try:
        out["groups"] = panel_groups(X, ids, percentile=percentile).as_dict()
    except Exception as e:                        # noqa: BLE001
        out["groups_refused"] = f"{type(e).__name__}: {e}"
    return out
