"""External-validation roll-up — render the external cells + inline clonality.

Globs the SEPARATE `wiki/external_validation_*.json` namespace (NOT the frozen
provenance_disjoint_validation_* surface) and renders
`wiki/external_validation_report_card.{md,json}`. The FROZEN
`build_validation_report_card.py` + `compute_lineage_metrics.py` are NOT touched
or invoked — this is the Fix-C parallel roll-up.

Per cell:
  - raw STRICT + RELAXED acc/sens/spec (from the artifact)
  - cluster-weighted STRICT sens/spec + Wilson CI + effective-lineage-N, computed
    INLINE by reusing clonality.py {greedy_representative_clusters,
    cluster_weighted_confusion, wilson_ci, effective_lineage_n} on the cohort
    genomes (Mash via the existing Docker path). On a non-Docker host the lineage
    block degrades to `{"status": "unavailable"}` with a reason — raw stays.

The Mash call lives only in `compute_lineage_for_cohort`/`main`; all aggregation,
CI, and rendering is pure + offline-testable.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.eval.clonality import (
    cluster_weighted_confusion,
    effective_lineage_n,
    wilson_ci,
)

FROZEN_REPORT_CARD = "decoder_validation_report_card"      # must NOT be our output
OUR_REPORT_CARD = "external_validation_report_card"


# --------------------------------------------------------------------------- #
# Pure: artifact loading + lineage block + rendering
# --------------------------------------------------------------------------- #
def load_external_artifacts(wiki_dir: str | Path) -> list[dict]:
    """Load every external-validation-v1 artifact (excludes the roll-up itself)."""
    out: list[dict] = []
    for f in sorted(glob.glob(str(Path(wiki_dir) / "external_validation_*.json"))):
        if Path(f).stem == OUR_REPORT_CARD:
            continue
        try:
            d = json.loads(Path(f).read_text(encoding="utf-8"))
        except (ValueError, OSError):
            continue
        if d.get("_schema") == "external-validation-v1":
            d["_path"] = f
            out.append(d)
    return out


def cluster_weighted_with_ci(preds: dict[str, str], labels: dict[str, object],
                             clusters: dict[str, int]) -> dict:
    """clonality.cluster_weighted_confusion + Wilson CIs + effective-lineage-N.

    sens CI on (tp, tp+fn); spec CI on (tn, tn+fp). Reuses the already-tested
    clonality math — this adds only the CI/effective-N wrapper the report card needs.
    """
    conf = cluster_weighted_confusion(preds, labels, clusters)
    tp, fp, tn, fn = conf["tp"], conf["fp"], conf["tn"], conf["fn"]
    conf["sens_ci"] = wilson_ci(tp, tp + fn)
    conf["spec_ci"] = wilson_ci(tn, tn + fp)
    conf["effective_lineage_n_R"] = effective_lineage_n(clusters, labels, "R")
    conf["effective_lineage_n_S"] = effective_lineage_n(clusters, labels, "S")
    return conf


def build_cell(artifact: dict, lineage: dict | None) -> dict:
    """One report-card cell from an artifact + an optional clonality block."""
    return {
        "cohort": artifact.get("cohort"),
        "organism": artifact.get("organism"),
        "drug": artifact.get("drug"),
        "evidence_tier": artifact.get("evidence_tier"),
        "strict": artifact.get("strict", {}),
        "relaxed": artifact.get("relaxed", {}),
        "lineage": lineage or {"status": "unavailable", "reason": "no clonality computed"},
        "independence_tier": artifact.get("independence_tier"),
    }


def render_json(cells: list[dict]) -> dict:
    return {
        "_schema": "external-validation-report-card-v1",
        "date": _date.today().isoformat(),
        "n_cells": len(cells),
        "note": ("external clinical re-validation of the frozen decoder; strict-tier is the "
                 "primary metric, relaxed secondary; raw sens/spec is clonality-inflated — see "
                 "the cluster-weighted block. Separate from the frozen decoder report card."),
        "cells": cells,
    }


def _fmt_ci(ci) -> str:
    if not ci or not isinstance(ci, (list, tuple)) or len(ci) != 2:
        return "—"
    return f"[{ci[0]}, {ci[1]}]"


def render_md(cells: list[dict]) -> str:
    lines = [
        f"# External-validation report card — {_date.today().isoformat()}",
        "",
        "External clinical re-validation of the FROZEN decoder on independent measured-MIC "
        "cohorts (different country / lab / AST method than the US-NCBI-PD tuning provenance). "
        "**Strict tier (HIGH_R/HIGH_S) is the primary metric**; relaxed (+DECISIVE) is secondary. "
        "Raw sens/spec is clonality-inflated — the cluster-weighted block (one vote per lineage, "
        "Wilson CI) is the honest companion. This is SEPARATE from the frozen decoder report card.",
        "",
        "| cohort | drug | strict sens | strict spec | strict n | lineage-wt sens (CI) | lineage-wt spec (CI) | eff-N R/S |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for c in cells:
        s = c["strict"]
        lin = c["lineage"]
        if lin.get("status") == "unavailable":
            lw_sens = lw_spec = effn = "n/a"
        else:
            lw_sens = f"{lin.get('sens')} {_fmt_ci(lin.get('sens_ci'))}"
            lw_spec = f"{lin.get('spec')} {_fmt_ci(lin.get('spec_ci'))}"
            effn = f"{lin.get('effective_lineage_n_R')}/{lin.get('effective_lineage_n_S')}"
        lines.append(
            f"| {c['cohort']} | {c['drug']} | {s.get('sens')} | {s.get('spec')} | "
            f"{s.get('n_scored')} | {lw_sens} | {lw_spec} | {effn} |"
        )
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Docker Mash composition (not unit-tested) + main
# --------------------------------------------------------------------------- #
def compute_lineage_for_cohort(artifact: dict, *, threshold: float = 0.001,
                               use_docker: bool = True) -> dict:
    """Compute the cluster-weighted STRICT block for one cohort artifact.

    Reads predictions_strict.json + the cohort refseq FASTAs, runs greedy-representative
    Mash clustering, and applies cluster_weighted_with_ci. Returns {"status":"unavailable",
    "reason":...} if Docker/Mash or the genomes are not available.
    """
    cohort, drug = artifact.get("cohort"), artifact.get("drug")
    base = Path(f"data/raw/{cohort}_extval_{drug}")
    pred_path = base / "predictions_strict.json"
    if not pred_path.exists():
        return {"status": "unavailable", "reason": f"no {pred_path}"}
    records = json.loads(pred_path.read_text(encoding="utf-8"))
    # keep only HIGH_R/HIGH_S-scored strains that produced an R/S call
    records = [r for r in records if str(r.get("prediction", "")).upper() in ("R", "S")]
    if not records:
        return {"status": "unavailable", "reason": "no R/S strict predictions"}
    try:
        from dna_decode.data import refseq
        from dna_decode.eval.clonality import greedy_representative_clusters
        gcache = base / "refseq"
        genomes = {r["gca"]: Path(refseq.fasta_path(r["gca"], gcache)) for r in records}
        genomes = {g: p for g, p in genomes.items() if p.exists()}
        if len(genomes) < 2:
            return {"status": "unavailable", "reason": f"only {len(genomes)} genome FASTAs on disk"}
        clusters = greedy_representative_clusters(genomes, threshold, use_docker=use_docker)
        preds = {r["gca"]: r["prediction"] for r in records if r["gca"] in clusters}
        labels = {r["gca"]: r["label"] for r in records if r["gca"] in clusters}
        block = cluster_weighted_with_ci(preds, labels, clusters)
        block["status"] = "ok"
        block["threshold"] = threshold
        return block
    except Exception as e:  # noqa: BLE001 — Mash/Docker absent or failed -> degrade, keep raw
        return {"status": "unavailable", "reason": f"{type(e).__name__}: {e}"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--wiki-dir", default="wiki")
    ap.add_argument("--no-clonality", action="store_true", help="skip the Docker Mash lineage block")
    ap.add_argument("--threshold", type=float, default=0.001)
    a = ap.parse_args()

    artifacts = load_external_artifacts(a.wiki_dir)
    cells = []
    for art in artifacts:
        lineage = None if a.no_clonality else compute_lineage_for_cohort(art, threshold=a.threshold)
        cells.append(build_cell(art, lineage))

    out_json = Path(a.wiki_dir) / f"{OUR_REPORT_CARD}.json"
    out_md = Path(a.wiki_dir) / f"{OUR_REPORT_CARD}.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(render_json(cells), indent=2), encoding="utf-8")
    out_md.write_text(render_md(cells), encoding="utf-8")
    print(f"wrote {out_md} + {out_json} ({len(cells)} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
