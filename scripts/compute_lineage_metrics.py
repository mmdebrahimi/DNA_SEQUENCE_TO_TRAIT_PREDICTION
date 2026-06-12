"""Per-cohort lineage-metrics recompute for the disclosure layer.

For each provenance-disjoint cohort under data/raw/*_provdisjoint_*, this script
emits the honest companion to the report card's raw-isolate sens/spec:

  - raw_N + raw sens/spec (RECONCILED against the cohort's committed artifact — M4)
  - effective_lineage_N @ {0.001, 0.005} per class (greedy-representative clustering)
  - cluster-weighted sens/spec, each with a Wilson CI + effective-N (C3)
  - a graded lineage annotation (buckets, not a binary >=20 tier — M3)
  - DISCORDANT-lineage count

M1 genome-completeness gate: every selected accession must have a present, non-empty,
valid-FASTA genome before Mash runs. Missing genomes are fetched (restartable); a
cohort that still can't complete is marked `partial` + `n_genomes_missing` and emits
NO lineage tier (a wrong lineage number is worse than an honest "incomplete").

The Mash call (Docker) lives only in `main()`; all the math (clustering aggregation,
reconciliation, cell assembly, upsert) is in pure helpers so it is unit-testable
offline against on-disk fixtures with no network and no Docker.

Output: wiki/provdisjoint_lineage_metrics.json (schema provdisjoint-lineage-metrics-v1),
keyed by canonical_cell_key, checkpointed per-cohort (Docker-wedge-safe), idempotent.

Usage: .venv/Scripts/python.exe scripts/compute_lineage_metrics.py [--cohort-dir data/raw/<one>]
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.cell_key import canonical_cell_key
from dna_decode.data.refseq import download_genome, fasta_path
from dna_decode.eval.amr_rules import call_resistance
from dna_decode.eval.clonality import (
    cluster_weighted_confusion,
    effective_lineage_n,
    wilson_ci,
)
from scripts.independent_cohort_validate import _conf
from scripts.organism_drug_validate import _run_dir

SCHEMA = "provdisjoint-lineage-metrics-v1"
SIDECAR = Path("wiki/provdisjoint_lineage_metrics.json")
THRESHOLDS: tuple[float, ...] = (0.001, 0.005)
# Grade is taken at the COARSER threshold (more merging -> fewer, more-conservative
# lineages) so the headline annotation is the honest lower bound on diversity.
GRADE_THRESHOLD = 0.005


class ReconcileMismatch(Exception):
    """Recomputed raw metrics disagree with the cohort's committed artifact (M4)."""


# --------------------------------------------------------------------------- #
# cohort <-> artifact resolution
# --------------------------------------------------------------------------- #
def parse_cohort_dir(name: str) -> tuple[str, str]:
    """('klebsiella_provdisjoint_ciprofloxacin') -> ('klebsiella', 'ciprofloxacin')."""
    if "_provdisjoint_" not in name:
        raise ValueError(f"not a provdisjoint cohort dir: {name!r}")
    slug, drug = name.split("_provdisjoint_", 1)
    return slug, drug


def find_artifact(slug: str, drug: str, wiki: Path = Path("wiki")) -> Path | None:
    """Locate the committed provenance-disjoint validation JSON for a cohort.

    Filenames are provenance_disjoint_validation_{slug}_{drug[:5]}_{date}.json;
    returns the latest by name when several dates exist."""
    matches = sorted(glob.glob(str(wiki / f"provenance_disjoint_validation_{slug}_{drug[:5]}_*.json")))
    return Path(matches[-1]) if matches else None


def read_selected(selected: Path) -> dict[str, int]:
    """Parse selected.tsv -> {accession: 1 if R else 0}."""
    out: dict[str, int] = {}
    for ln in selected.read_text(encoding="utf-8").splitlines():
        if "\t" in ln:
            acc, rs = ln.split("\t", 1)
            out[acc.strip()] = 1 if rs.strip() == "R" else 0
    return out


# --------------------------------------------------------------------------- #
# M1 genome completeness gate
# --------------------------------------------------------------------------- #
def _fasta_ok(p: Path) -> bool:
    """Present + non-empty + first non-blank line is a FASTA header (>= 1 contig)."""
    if not p.exists() or p.stat().st_size == 0:
        return False
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.strip():
                    return line.startswith(">")
    except OSError:
        return False
    return False


def ensure_cohort_genomes(
    selected: dict[str, int], refseq_root: Path, *, fetch: bool = True
) -> tuple[dict[str, Path], list[str]]:
    """Gate (and optionally fetch) genomes for every selected accession.

    Returns (present {acc: fasta_path}, missing [acc...]). download_genome is
    skip-complete + restartable via its `.complete` sentinel, so re-runs are cheap.
    """
    present: dict[str, Path] = {}
    missing: list[str] = []
    for acc in selected:
        fp = fasta_path(acc, refseq_root)
        if not _fasta_ok(fp) and fetch:
            try:
                download_genome(acc, refseq_root)
                fp = fasta_path(acc, refseq_root)
            except Exception:  # noqa: BLE001 — any fetch failure -> mark missing (M1)
                pass
        if _fasta_ok(fp):
            present[acc] = fp
        else:
            missing.append(acc)
    return present, missing


# --------------------------------------------------------------------------- #
# M4 raw reconciliation (raises on mismatch)
# --------------------------------------------------------------------------- #
def reconcile_raw_metrics(
    selected: dict[str, int],
    own_runs: Path,
    reuse_glob: str,
    drug: str,
    registry_organism: str,
    artifact_metrics: dict,
) -> tuple[dict, dict[str, str]]:
    """Recompute raw sens/spec with the cohort's ORIGINAL rule-path args and assert
    it reconciles with the committed artifact before any weighted number is trusted.

    Returns (raw_conf, preds_by_acc). Raises ReconcileMismatch if the recomputed
    confusion disagrees with `artifact_metrics`.
    """
    applied: list[tuple[str, int]] = []
    preds: dict[str, str] = {}
    for acc, y in selected.items():
        rd = _run_dir(acc, own_runs, reuse_glob)
        if rd is None:
            continue
        pred = call_resistance(rd / "main.tsv", drug, organism=registry_organism)["prediction"]
        preds[acc] = pred
        applied.append((pred, y))
    raw = _conf(applied)
    for k in ("tp", "fp", "tn", "fn", "sens", "spec", "n_scored"):
        if artifact_metrics.get(k) != raw.get(k):
            raise ReconcileMismatch(
                f"{k}: artifact={artifact_metrics.get(k)} recomputed={raw.get(k)} "
                f"(refusing to trust the cluster-weighted metric)"
            )
    return raw, preds


# --------------------------------------------------------------------------- #
# pure metric assembly (no Docker / network)
# --------------------------------------------------------------------------- #
def graded_lineage_bucket(n_eff: int) -> str:
    """Graded lineage-diversity annotation (M3) — buckets, not a binary tier."""
    if n_eff >= 15:
        return "moderate (>=15 effective lineages)"
    if n_eff >= 8:
        return "limited (8-14 effective lineages)"
    if n_eff >= 3:
        return "scarce (3-7 effective lineages)"
    return "clonal (<3 effective lineages)"


def _weighted_block(preds: dict[str, str], labels: dict[str, int], clusters: dict[str, int]) -> dict:
    """cluster_weighted_confusion + inline Wilson CIs (C3 — never a point without CI)."""
    w = cluster_weighted_confusion(preds, labels, clusters)
    sens_ci = wilson_ci(w["tp"], w["tp"] + w["fn"])
    spec_ci = wilson_ci(w["tn"], w["tn"] + w["fp"])
    return {
        "sens": w["sens"], "sens_ci": list(sens_ci), "sens_eff_n": w["tp"] + w["fn"],
        "spec": w["spec"], "spec_ci": list(spec_ci), "spec_eff_n": w["tn"] + w["fp"],
        "tp": w["tp"], "fp": w["fp"], "tn": w["tn"], "fn": w["fn"],
        "n_clusters_R": w["n_clusters_R"], "n_clusters_S": w["n_clusters_S"],
        "n_discordant": w["n_discordant"], "n_cluster_abstain": w["n_cluster_abstain"],
    }


def build_threshold_results(
    preds: dict[str, str],
    labels: dict[str, int],
    clusters_by_threshold: dict[float, dict[str, int]],
) -> dict[str, dict]:
    """Per-threshold {effective_lineage_N per class + weighted block}. Pure."""
    out: dict[str, dict] = {}
    for t, clusters in clusters_by_threshold.items():
        out[str(t)] = {
            "effective_lineage_N_R": effective_lineage_n(clusters, labels, "R"),
            "effective_lineage_N_S": effective_lineage_n(clusters, labels, "S"),
            "cluster_weighted": _weighted_block(preds, labels, clusters),
        }
    return out


def build_cell(
    *,
    organism: str,
    drug: str,
    cohort: str,
    raw: dict,
    raw_reconciled: bool,
    partial: bool,
    n_genomes_missing: int,
    threshold_results: dict[str, dict] | None,
) -> dict:
    """Assemble one lineage-metrics cell. A partial/unreconciled cohort emits raw +
    flags but NO lineage tier (threshold_results None, lineage_tier_emitted False)."""
    cell = {
        "organism": organism, "drug": drug, "cohort": cohort,
        "raw_N": raw.get("n_scored"),
        "raw_sens": raw.get("sens"), "raw_spec": raw.get("spec"),
        "raw_tp": raw.get("tp"), "raw_fp": raw.get("fp"),
        "raw_tn": raw.get("tn"), "raw_fn": raw.get("fn"),
        "raw_reconciled": raw_reconciled,
        "partial": partial,
        "n_genomes_missing": n_genomes_missing,
    }
    if partial or not raw_reconciled or not threshold_results:
        cell["thresholds"] = {}
        cell["lineage_grade"] = None
        cell["lineage_tier_emitted"] = False
        return cell
    cell["thresholds"] = threshold_results
    grade_n = threshold_results.get(str(GRADE_THRESHOLD), {}).get("effective_lineage_N_R", 0)
    cell["lineage_grade"] = graded_lineage_bucket(grade_n)
    cell["lineage_tier_emitted"] = True
    return cell


# --------------------------------------------------------------------------- #
# idempotent sidecar upsert
# --------------------------------------------------------------------------- #
def load_sidecar(path: Path = SIDECAR) -> dict:
    if path.exists():
        try:
            d = json.loads(path.read_text(encoding="utf-8"))
            d.setdefault("cells", [])
            return d
        except Exception:  # noqa: BLE001
            pass
    return {"_schema": SCHEMA, "date": _date.today().isoformat(), "cells": []}


def upsert_cell(sidecar: dict, cell: dict) -> dict:
    """Replace the cell with the same canonical (organism, drug); else append."""
    key = canonical_cell_key(cell["organism"], cell["drug"])
    cells = [c for c in sidecar.get("cells", [])
             if canonical_cell_key(c["organism"], c["drug"]) != key]
    cells.append(cell)
    cells.sort(key=lambda c: canonical_cell_key(c["organism"], c["drug"]))
    sidecar["cells"] = cells
    return sidecar


def write_sidecar(sidecar: dict, path: Path = SIDECAR) -> None:
    sidecar["date"] = _date.today().isoformat()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(sidecar, indent=2), encoding="utf-8")
    tmp.replace(path)


# --------------------------------------------------------------------------- #
# orchestration (Docker Mash here only)
# --------------------------------------------------------------------------- #
def process_cohort(cohort_dir: Path, *, use_docker: bool = True) -> dict:
    """Full per-cohort pipeline: reconcile -> genome gate -> Mash cluster -> cell."""
    from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
    from dna_decode.eval.phylogeny import compute_mash_distances

    slug, drug = parse_cohort_dir(cohort_dir.name)
    artifact_path = find_artifact(slug, drug)
    if artifact_path is None:
        raise FileNotFoundError(f"no committed artifact for {cohort_dir.name}")
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    organism = artifact["organism"]
    reg_org = artifact.get("registry_organism") or organism
    a_drug = artifact["drug"]
    a_metrics = artifact.get("metrics", {})

    selected = read_selected(cohort_dir / "selected.tsv")
    own_runs = cohort_dir / "amrfinder_runs"
    reuse_glob = f"data/raw/{slug}_*/amrfinder_runs"

    try:
        raw, preds = reconcile_raw_metrics(selected, own_runs, reuse_glob, a_drug, reg_org, a_metrics)
    except ReconcileMismatch as e:
        print(f"  [{cohort_dir.name}] RECONCILE FAILED: {e}", file=sys.stderr)
        return build_cell(organism=organism, drug=a_drug, cohort=cohort_dir.name,
                          raw=a_metrics, raw_reconciled=False, partial=False,
                          n_genomes_missing=0, threshold_results=None)

    present, missing = ensure_cohort_genomes(selected, cohort_dir / "refseq")
    if missing:
        print(f"  [{cohort_dir.name}] PARTIAL: {len(missing)}/{len(selected)} genomes missing "
              f"-> no lineage tier", file=sys.stderr)
        return build_cell(organism=organism, drug=a_drug, cohort=cohort_dir.name,
                          raw=raw, raw_reconciled=True, partial=True,
                          n_genomes_missing=len(missing), threshold_results=None)

    # ONE Mash sketch+dist per cohort; cluster at every threshold from the same matrix
    # (avoids redundant Docker invocations — each container spin-up is a churn-corruption risk on this host).
    dm = compute_mash_distances(present, use_docker=use_docker)
    clusters_by_threshold = {
        t: greedy_representative_clusters_from_matrix(dm, t) for t in THRESHOLDS
    }
    tr = build_threshold_results(preds, selected, clusters_by_threshold)
    return build_cell(organism=organism, drug=a_drug, cohort=cohort_dir.name,
                      raw=raw, raw_reconciled=True, partial=False,
                      n_genomes_missing=0, threshold_results=tr)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-dir", type=Path, default=None,
                    help="process a single cohort dir (default: all data/raw/*_provdisjoint_*)")
    ap.add_argument("--no-docker", action="store_true", help="use a native mash binary instead of Docker")
    a = ap.parse_args(argv)

    if a.cohort_dir:
        cohort_dirs = [a.cohort_dir]
    else:
        cohort_dirs = [Path(p) for p in sorted(glob.glob("data/raw/*_provdisjoint_*")) if Path(p).is_dir()]
    if not cohort_dirs:
        print("no provdisjoint cohort dirs found", file=sys.stderr)
        return 1

    sidecar = load_sidecar()
    for cd in cohort_dirs:
        print(f"[lineage] {cd.name} ...", flush=True)
        try:
            cell = process_cohort(cd, use_docker=not a.no_docker)
        except Exception as e:  # noqa: BLE001 — checkpoint discipline: one bad cohort doesn't lose the rest
            print(f"  [{cd.name}] ERROR: {type(e).__name__}: {e}", file=sys.stderr)
            continue
        sidecar = upsert_cell(sidecar, cell)
        write_sidecar(sidecar)  # checkpoint per-cohort (Docker-wedge-safe)
        tier = cell.get("lineage_grade") or ("partial" if cell["partial"] else "no-tier")
        print(f"  [{cd.name}] raw_N={cell['raw_N']} tier={tier}", flush=True)

    print(f"wrote {SIDECAR} ({len(sidecar['cells'])} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
