"""Mash-cluster orchestration for the cipro N=147 cohort -- PASS-path artifact.

Triggered by `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` Step P
sub-step 1 only when the falsifier verdict is PASS. Codex on the Precision 7780
runs this (Docker Desktop + Mash image required); Claude on the GTX 860M laptop
drafts + tests the pure-logic helpers locally.

Threshold sweep (per brainstorm B5):
  sweep [0.02, 0.03, 0.04, 0.05, 0.07, 0.10]
  pick = the threshold satisfying ALL of:
    - n_clades >= 3
    - max-clade fraction < 0.60 (no single clade dominates)
    - intra-clade-vs-inter-clade variance ratio MINIMIZED (lower = clades tight + well separated).
      NOTE: this docstring previously said "maximized"; the implementation and its tests MINIMIZE.
      Corrected 2026-07-09.
  fallback = 0.05 (matches plan default, flagged in JSON sidecar)

KNOWN DEFECTS (documented 2026-07-09, deliberately NOT fixed here -- the selection rule feeds downstream
contracts, so the fix is surfaced for ratification. See
`wiki/cipro_mash_clades_n147_threshold_diagnostic_2026-07-09.md`):
  1. 0.02 is the GRID FLOOR and, on the real N=147 cohort, the ONLY qualifying rung -- so the chosen
     threshold is a boundary artifact. It yields a 57.8%-largest clade, which makes leave-one-clade-out
     CV effectively leave-58%-out. (0.015 would beat it at the 4th decimal if merely added.)
  2. `variance_ratio` -> 0 as clusters -> singletons, so MINIMIZING it is monotonically biased toward
     over-splitting. Extending the grid downward WITHOUT an over-split guard (max singleton fraction /
     min clade size) selects 0.002 -> 101 clades, 56.5% singletons. The coarse grid is load-bearing by
     accident; the qualifying criteria guard only against a dominant clade, never against over-splitting.
  3. `cluster_by_ani` is single-linkage and CHAINS (at 0.015 it reports a 57.8% clade where
     `clonality.greedy_representative_clusters_from_matrix` reports 31.3%). The greedy-representative
     method is the chaining-resistant one the frozen lineage layer uses.
  Also: a missing FASTA is WARN-and-skipped, so an incomplete refseq cache silently clusters a subset
  (43 of 146 accessions were cached before the 2026-07-09 run fetched the rest).

Outputs:
  wiki/cipro_mash_clades_n147_<DATE>.json  (machine-readable: clade per strain + scoring summary)
  wiki/cipro_mash_clades_n147_<DATE>.md    (narrative: clade-size distribution + per-clade R/S balance)
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path

import numpy as np


CANDIDATE_THRESHOLDS: tuple[float, ...] = (0.02, 0.03, 0.04, 0.05, 0.07, 0.10)
FALLBACK_THRESHOLD: float = 0.05
MIN_CLADES: int = 3
MAX_SINGLE_CLADE_FRACTION: float = 0.60


@dataclass(frozen=True)
class ThresholdScore:
    """Pure scoring record for one threshold candidate."""

    threshold: float
    n_clades: int
    max_clade_fraction: float
    variance_ratio: float  # intra-clade variance / inter-clade variance (smaller = tighter)
    satisfies_min_clades: bool
    satisfies_max_fraction: bool

    @property
    def fully_satisfied(self) -> bool:
        return self.satisfies_min_clades and self.satisfies_max_fraction


def score_threshold(
    distance_matrix: np.ndarray,
    cluster_assignments: dict[str, int],
    strain_ids: list[str],
    threshold: float,
) -> ThresholdScore:
    """Score a clustering result against the B5 criteria.

    variance_ratio = mean(intra-clade pairwise distance) / mean(inter-clade pairwise distance)
    Lower ratio = clades are tighter + better-separated. If only 1 clade -> ratio = NaN.
    """
    n = len(strain_ids)
    if n == 0:
        return ThresholdScore(
            threshold=threshold,
            n_clades=0,
            max_clade_fraction=0.0,
            variance_ratio=float("nan"),
            satisfies_min_clades=False,
            satisfies_max_fraction=False,
        )

    cluster_sizes = Counter(cluster_assignments.values())
    n_clades = len(cluster_sizes)
    max_clade_size = max(cluster_sizes.values())
    max_clade_fraction = max_clade_size / n

    # Compute variance ratio
    intra_dists: list[float] = []
    inter_dists: list[float] = []
    sid_to_clade = cluster_assignments
    for i in range(n):
        for j in range(i + 1, n):
            d = float(distance_matrix[i, j])
            if sid_to_clade[strain_ids[i]] == sid_to_clade[strain_ids[j]]:
                intra_dists.append(d)
            else:
                inter_dists.append(d)
    if not intra_dists or not inter_dists:
        variance_ratio = float("nan")
    else:
        intra_mean = float(np.mean(intra_dists))
        inter_mean = float(np.mean(inter_dists))
        variance_ratio = intra_mean / inter_mean if inter_mean > 0 else float("nan")

    return ThresholdScore(
        threshold=threshold,
        n_clades=n_clades,
        max_clade_fraction=max_clade_fraction,
        variance_ratio=variance_ratio,
        satisfies_min_clades=n_clades >= MIN_CLADES,
        satisfies_max_fraction=max_clade_fraction < MAX_SINGLE_CLADE_FRACTION,
    )


def pick_best_threshold(
    distance_matrix: np.ndarray,
    strain_ids: list[str],
    cluster_at_threshold_fn,
    candidates: tuple[float, ...] = CANDIDATE_THRESHOLDS,
    fallback: float = FALLBACK_THRESHOLD,
) -> tuple[float, list[ThresholdScore], bool]:
    """Pick the threshold satisfying all B5 criteria with the LOWEST variance_ratio.

    Args:
      distance_matrix: N x N float matrix of pairwise distances.
      strain_ids: row/column ordering for the matrix.
      cluster_at_threshold_fn: callable(threshold) -> dict[strain_id, cluster_id].
        Injected so this function is pure + testable.
      candidates: thresholds to evaluate.
      fallback: threshold returned if none of the candidates satisfy all criteria.

    Returns:
      (chosen_threshold, per-candidate-scores, fellback_to_default)
    """
    scores: list[ThresholdScore] = []
    for t in candidates:
        assignments = cluster_at_threshold_fn(t)
        scores.append(score_threshold(distance_matrix, assignments, strain_ids, t))

    qualifying = [s for s in scores if s.fully_satisfied and not np.isnan(s.variance_ratio)]
    if not qualifying:
        return fallback, scores, True

    # Pick the threshold with the LOWEST variance_ratio among qualifying candidates.
    best = min(qualifying, key=lambda s: s.variance_ratio)
    return best.threshold, scores, False


def per_clade_label_balance(
    cluster_assignments: dict[str, int],
    labels_by_strain: dict[str, int | None],
) -> dict[int, dict[str, int]]:
    """Per-clade counts of R / S / unknown labels.

    Returns:
      clade_id -> {"R": count, "S": count, "unknown": count, "n": total}.
    """
    out: dict[int, dict[str, int]] = {}
    for sid, clade in cluster_assignments.items():
        slot = out.setdefault(clade, {"R": 0, "S": 0, "unknown": 0, "n": 0})
        label = labels_by_strain.get(sid)
        if label == 1:
            slot["R"] += 1
        elif label == 0:
            slot["S"] += 1
        else:
            slot["unknown"] += 1
        slot["n"] += 1
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Mash-cluster N=147 cipro cohort with threshold sweep.")
    parser.add_argument("--cohort", type=Path, required=True,
                        help="Cohort parquet (e.g., data/processed/stage2_n150_cipro_cohort.parquet)")
    parser.add_argument("--refseq-cache", type=Path, required=True,
                        help="RefSeq cache root (FASTA per accession lives at <root>/<accession>/genome.fna)")
    parser.add_argument("--drug", default="ciprofloxacin",
                        help="Drug for per-clade R/S balance reporting")
    parser.add_argument("--use-docker", action="store_true", default=True,
                        help="Route Mash through Docker (required on Windows hosts without native binary)")
    parser.add_argument("--output-prefix", type=Path,
                        default=Path(f"wiki/cipro_mash_clades_n147_{_date.today().isoformat()}"))
    args = parser.parse_args(argv)

    from dna_decode.data.cohort import load_cohort
    from dna_decode.data.refseq import fasta_path
    from dna_decode.eval.phylogeny import (
        cluster_by_ani,
        compute_mash_distances,
    )

    cohort = load_cohort(args.cohort)
    drug_lower = args.drug.lower()

    strain_genomes: dict[str, Path] = {}
    labels_by_strain: dict[str, int | None] = {}
    for s in cohort.strains:
        fna = fasta_path(s.assembly_accession, args.refseq_cache)
        if not fna.exists():
            print(f"[mash-cluster] WARN no FASTA for {s.strain_id} ({s.assembly_accession}); skipping", file=sys.stderr)
            continue
        strain_genomes[s.strain_id] = fna
        labels_by_strain[s.strain_id] = s.ast_labels.get(drug_lower)

    if len(strain_genomes) < MIN_CLADES:
        print(f"[mash-cluster] FAIL only {len(strain_genomes)} strains with FASTA; need >= {MIN_CLADES}", file=sys.stderr)
        return 2

    print(f"[mash-cluster] computing Mash distances for {len(strain_genomes)} strains "
          f"(docker={args.use_docker})")
    dm = compute_mash_distances(strain_genomes, use_docker=args.use_docker)

    chosen, all_scores, fellback = pick_best_threshold(
        dm.matrix,
        dm.strain_ids,
        cluster_at_threshold_fn=lambda t: cluster_by_ani(dm, threshold=t),
    )
    print(f"[mash-cluster] chose threshold={chosen:.3f} "
          f"(fellback_to_default={fellback})")

    final_assignments = cluster_by_ani(dm, threshold=chosen)
    clade_balance = per_clade_label_balance(final_assignments, labels_by_strain)

    # ---- JSON sidecar ----
    out_json = args.output_prefix.with_suffix(".json")
    out_md = args.output_prefix.with_suffix(".md")
    out_json.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "run_date": _date.today().isoformat(),
        "cohort_path": str(args.cohort),
        "drug": args.drug,
        "n_strains": len(strain_genomes),
        "chosen_threshold": chosen,
        "fellback_to_default": fellback,
        "min_clades_required": MIN_CLADES,
        "max_single_clade_fraction": MAX_SINGLE_CLADE_FRACTION,
        "candidate_threshold_scores": [
            {
                "threshold": s.threshold,
                "n_clades": s.n_clades,
                "max_clade_fraction": s.max_clade_fraction,
                "variance_ratio": s.variance_ratio,
                "satisfies_min_clades": s.satisfies_min_clades,
                "satisfies_max_fraction": s.satisfies_max_fraction,
                "fully_satisfied": s.fully_satisfied,
            }
            for s in all_scores
        ],
        "clade_assignments": final_assignments,  # strain_id -> clade_id
        "per_clade_label_balance": clade_balance,
    }
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[mash-cluster] wrote {out_json}")

    # ---- MD sidecar ----
    lines = [
        f"# Cipro Mash-clades N={len(strain_genomes)} cohort ({_date.today().isoformat()})",
        "",
        f"**Drug:** {args.drug}",
        f"**Chosen threshold:** {chosen:.3f}" + (" (fallback)" if fellback else ""),
        f"**Number of clades:** {len(clade_balance)}",
        "",
        "## Threshold-sweep scoring",
        "",
        "| threshold | n_clades | max_clade_frac | variance_ratio | min_clades_ok | max_frac_ok | qualified |",
        "|---:|---:|---:|---:|:---:|:---:|:---:|",
    ]
    for s in all_scores:
        lines.append(
            f"| {s.threshold:.3f} | {s.n_clades} | {s.max_clade_fraction:.3f} | "
            f"{s.variance_ratio:.4f} | "
            f"{'Y' if s.satisfies_min_clades else 'n'} | "
            f"{'Y' if s.satisfies_max_fraction else 'n'} | "
            f"{'Y' if s.fully_satisfied else 'n'} |"
        )
    lines += ["", "## Per-clade R/S balance", "", "| clade_id | n | R | S | unknown | R_fraction |", "|---:|---:|---:|---:|---:|---:|"]
    for cid, c in sorted(clade_balance.items()):
        r_frac = c["R"] / c["n"] if c["n"] else 0.0
        lines.append(f"| {cid} | {c['n']} | {c['R']} | {c['S']} | {c['unknown']} | {r_frac:.3f} |")
    lines += ["", "Generated by `scripts/mash_cluster_n147.py`."]
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"[mash-cluster] wrote {out_md}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
