"""BacDive carbon-utilization substrate feasibility census (EP-6 entry gate).

The go/no-go gate BEFORE any NT-cache populate on the carbon-utilization
substrate. Mirrors `scripts/bvbrc_strict_mic_4drug_census.py`: a layered census
that answers "which carbon source, if any, has a de-confoundable E. coli
utilizer/non-utilizer cohort large enough to test whether NT embeddings beat
gene-content on a sampling-independent label with no curated catalog?"

Per-carbon-source pipeline:
  Stage 1: distinct E. coli strains with a binary utilization label
  Stage 2: both classes present + minority fraction >= --min-minority-frac
  Stage 3: strains with a downloadable assembly_accession (NT-cacheable)
  Stage 4: de-confound gate (within-MLST utilizer/non-utilizer contrast) when an
           MLST sidecar is supplied — reuses dna_decode.eval.cohort_deconfound

A carbon source is FEASIBLE iff it clears Stage 2 + Stage 3 floors AND (if MLST
provided) the de-confound gate returns DE_CONFOUNDED (promotable). Without MLST
the verdict is FEASIBLE_PENDING_DECONFOUND (counts pass; gate not yet runnable).

The headline picks the carbon source where a learned decoder has the best shot:
both-class balance + downloadable genomes + within-lineage contrast. Li et al.
2023 (PMC10729968) warn that easy carbon sources are already near-ceiling for
gene-content RF out-of-clade, so the ranking also surfaces minority_fraction
(harder/balanced sources are where embeddings could add value).

Output: wiki/bacdive_carbon_util_feasibility_<date>.{md,json}.
Exit 0 if >=1 carbon source is FEASIBLE / FEASIBLE_PENDING_DECONFOUND; 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path
from typing import Optional

from dna_decode.data.bacdive import (
    census_carbon_sources,
    get_binary_labels,
    load_bacdive_carbon,
)
from dna_decode.eval.cohort_deconfound import (
    DE_CONFOUNDED,
    confound_report,
    render_report,
)

DEFAULT_MIN_STRAINS = 100        # ≥100 de-confoundable floor (intent-contract bar)
DEFAULT_MIN_MINORITY_FRAC = 0.15  # at least 15% in the minority class


def _load_mlst_sidecar(path: Optional[Path]) -> dict[str, str]:
    """Optional strain_id → MLST map (JSON or 2-col CSV/TSV). Empty if absent."""
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"MLST sidecar not found at {p}")
    if p.suffix.lower() == ".json":
        return {str(k): str(v) for k, v in json.loads(p.read_text()).items()}
    import pandas as pd
    df = pd.read_csv(p, sep=None, engine="python", dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    sid_col = next((c for c in ("strain_id", "id", "genome_id") if c in df.columns), df.columns[0])
    mlst_col = next((c for c in ("mlst", "st", "sequence_type") if c in df.columns), df.columns[1])
    return {str(r[sid_col]).strip(): str(r[mlst_col]).strip() for _, r in df.iterrows()}


def assess_carbon_source(
    table,
    carbon_source: str,
    mlst_map: dict[str, str],
    *,
    min_strains: int,
    min_minority_frac: float,
) -> dict:
    """Layered feasibility assessment for one carbon source."""
    labels = get_binary_labels(table, carbon_source)
    n = len(labels)
    pos = sum(1 for v in labels.values() if v == 1)
    neg = n - pos
    minority_frac = (min(pos, neg) / n) if n else 0.0

    sub = table[table["carbon_source"] == carbon_source.lower()]
    with_acc = set(sub[sub["assembly_accession"].str.len() > 0]["strain_id"])

    result = {
        "carbon_source": carbon_source,
        "n_strains": n,
        "n_positive": pos,
        "n_negative": neg,
        "minority_fraction": round(minority_frac, 4),
        "n_with_accession": len(with_acc),
        "deconfound_verdict": None,
        "deconfound_reason": None,
        "verdict": None,
    }

    # Stage 2: both-class + balance floor
    if n < min_strains:
        result["verdict"] = "INFEASIBLE_TOO_FEW"
        return result
    if pos == 0 or neg == 0:
        result["verdict"] = "INFEASIBLE_SINGLE_CLASS"
        return result
    if minority_frac < min_minority_frac:
        result["verdict"] = "INFEASIBLE_IMBALANCED"
        return result
    # Stage 3: downloadable genomes
    if len(with_acc) < min_strains:
        result["verdict"] = "INFEASIBLE_NO_ACCESSIONS"
        return result

    # Stage 4: de-confound gate (only when MLST is available)
    if mlst_map:
        sids = list(labels.keys())
        ys = [labels[s] for s in sids]
        lin = [mlst_map.get(s) for s in sids]
        rep = confound_report(ys, lin)
        result["deconfound_verdict"] = rep["verdict"]
        result["deconfound_reason"] = render_report(rep)
        result["verdict"] = "FEASIBLE" if rep["verdict"] == DE_CONFOUNDED else "BLOCKED_CONFOUNDED"
    else:
        result["verdict"] = "FEASIBLE_PENDING_DECONFOUND"
    return result


def run_census(
    export_path: Path,
    mlst_path: Optional[Path],
    *,
    min_strains: int,
    min_minority_frac: float,
    organism: str,
) -> dict:
    table = load_bacdive_carbon(export_path, organism=organism)
    mlst_map = _load_mlst_sidecar(mlst_path)

    overview = census_carbon_sources(table, min_strains=1)
    assessments = [
        assess_carbon_source(
            table, c.carbon_source, mlst_map,
            min_strains=min_strains, min_minority_frac=min_minority_frac,
        )
        for c in overview
    ]
    feasible = [a for a in assessments
                if a["verdict"] in ("FEASIBLE", "FEASIBLE_PENDING_DECONFOUND")]
    # rank feasible by balance (closer to 0.5 minority = harder = embedding niche),
    # then by downloadable cohort size
    feasible.sort(key=lambda a: (-a["minority_fraction"], -a["n_with_accession"]))
    return {
        "organism": organism,
        "n_carbon_sources_total": len(overview),
        "n_carbon_sources_assessed": len(assessments),
        "min_strains_floor": min_strains,
        "min_minority_frac": min_minority_frac,
        "mlst_provided": bool(mlst_map),
        "n_feasible": len(feasible),
        "feasible_ranked": feasible,
        "all_assessments": assessments,
    }


def render_md(census: dict) -> str:
    d = _date.today().isoformat()
    lines = [
        f"# BacDive carbon-utilization substrate feasibility census — {d}",
        "",
        "> EP-6 entry gate. Which E. coli carbon source has a de-confoundable",
        "> utilizer/non-utilizer cohort big enough to test NT embeddings vs gene-content",
        "> on a sampling-INDEPENDENT label with NO curated catalog?",
        "",
        f"- Organism: `{census['organism']}`",
        f"- Carbon sources in export: {census['n_carbon_sources_total']}",
        f"- Min-strains floor: {census['min_strains_floor']} · "
        f"min-minority-frac: {census['min_minority_frac']}",
        f"- MLST sidecar provided (de-confound gate runnable): {census['mlst_provided']}",
        f"- **Feasible carbon sources: {census['n_feasible']}**",
        "",
        "## Feasible (ranked: balance, then downloadable cohort size)",
        "",
        "| carbon source | N | +/- | minority frac | with-accession | de-confound | verdict |",
        "|---|---:|---|---:|---:|---|---|",
    ]
    for a in census["feasible_ranked"]:
        dv = a["deconfound_verdict"] or "—"
        lines.append(
            f"| {a['carbon_source']} | {a['n_strains']} | "
            f"{a['n_positive']}/{a['n_negative']} | {a['minority_fraction']} | "
            f"{a['n_with_accession']} | {dv} | {a['verdict']} |"
        )
    if not census["feasible_ranked"]:
        lines.append("| (none) | | | | | | |")
    lines += [
        "",
        "## All assessments",
        "",
        "| carbon source | N | +/- | with-acc | verdict |",
        "|---|---:|---|---:|---|",
    ]
    for a in census["all_assessments"]:
        lines.append(
            f"| {a['carbon_source']} | {a['n_strains']} | "
            f"{a['n_positive']}/{a['n_negative']} | {a['n_with_accession']} | {a['verdict']} |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "- `FEASIBLE` = clears count + balance + accession floors AND de-confound gate",
        "  returned DE_CONFOUNDED (within-MLST utilizer/non-utilizer contrast exists).",
        "- `FEASIBLE_PENDING_DECONFOUND` = count floors clear; MLST not supplied so the",
        "  de-confound gate could not run. Supply `--mlst` to resolve.",
        "- `BLOCKED_CONFOUNDED` = no within-lineage contrast → predicting utilization would",
        "  predict lineage, not metabolism (the Li et al. 2023 phylogeny-dominance trap).",
        "- Ranking prioritizes HIGHER minority_fraction (balanced/harder sources) because",
        "  easy carbon sources are already near-ceiling for gene-content RF out-of-clade —",
        "  the embedding niche is the hard ones.",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--export", required=True, type=Path,
                    help="BacDive carbon-utilization export (long-format CSV/TSV)")
    ap.add_argument("--mlst", type=Path, default=None,
                    help="optional strain_id→MLST sidecar (JSON or 2-col CSV/TSV)")
    ap.add_argument("--organism", default="Escherichia coli")
    ap.add_argument("--min-strains", type=int, default=DEFAULT_MIN_STRAINS)
    ap.add_argument("--min-minority-frac", type=float, default=DEFAULT_MIN_MINORITY_FRAC)
    ap.add_argument("--out-dir", type=Path, default=Path("wiki"))
    args = ap.parse_args(argv)

    census = run_census(
        args.export, args.mlst,
        min_strains=args.min_strains,
        min_minority_frac=args.min_minority_frac,
        organism=args.organism,
    )

    d = _date.today().isoformat()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    md_path = args.out_dir / f"bacdive_carbon_util_feasibility_{d}.md"
    json_path = args.out_dir / f"bacdive_carbon_util_feasibility_{d}.json"
    md_path.write_text(render_md(census), encoding="utf-8")
    json_path.write_text(json.dumps(census, indent=2), encoding="utf-8")

    # ascii-safe console summary (Windows stdout is cp1252; md file is utf-8)
    print(f"carbon sources assessed: {census['n_carbon_sources_assessed']} | "
          f"FEASIBLE: {census['n_feasible']}")
    for a in census["feasible_ranked"]:
        print(f"  + {a['carbon_source']}: N={a['n_strains']} "
              f"({a['n_positive']}/{a['n_negative']}) minority={a['minority_fraction']} "
              f"acc={a['n_with_accession']} deconfound={a['deconfound_verdict'] or 'pending'} "
              f"-> {a['verdict']}")
    print(f"Wrote {md_path} + {json_path}")
    return 0 if census["n_feasible"] >= 1 else 1


if __name__ == "__main__":
    sys.exit(main())
