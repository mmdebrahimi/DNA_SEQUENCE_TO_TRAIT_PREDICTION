"""Finer Mash-based lineage collapse of the INDEPENDENT TB AMR-Portal cohort.

Closes the named-deferred item from `wiki/tb_independent_lineage_collapsed_result_2026-07-02.md`:
the shipped independent lineage number collapses by the pinned **Napier barcode**, which is COARSE
(~67 lineages over 2,845 isolates) and therefore EXCLUDES 43-44 mixed-label "discordant" lineages
from sens/spec. A finer, distance-based collapse (Mash greedy-representative, the SAME clustering the
frozen provenance-disjoint lineage layer uses) resolves sub-lineage structure, so fewer isolates are
lost to discordance and the effective-N / Wilson CI tighten.

**What is REUSED, unchanged (frozen math):**
  - `dna_decode.eval.clonality.greedy_representative_clusters_from_matrix` (chaining-resistant greedy rep)
  - `dna_decode.eval.clonality.cluster_weighted_confusion` (one vote per same-label lineage; DISCORDANT
    lineages excluded + counted, never majority-voted)
  - `dna_decode.eval.clonality.wilson_ci` / `effective_lineage_n`

**What is NEW here:** only the Mash invocation. `phylogeny.compute_mash_distances` cannot be used at
this scale on this host for two grounded reasons:
  1. it passes every FASTA path as argv (2,845 paths ~= 71 KB) -> exceeds the Windows 32 KB
     CreateProcess limit;
  2. it `shutil.copy`s the whole cohort (12 GB) into a temp dir on C: (17 GB free) -> disk blowout.
So we invoke mash directly against a READ-ONLY bind mount of the assembly cache, using `mash sketch -l`
(file-of-file-names) + `mash triangle` (compact relaxed-Phylip lower triangle instead of an 8.1M-line
`dist` dump). Both flags verified against `mash 2.3 -h` in the pinned container.

**Sketch size:** `-s 10000`, NOT the default 1000. M. tuberculosis is monomorphic — within/between
sub-lineage Mash distances live near 1e-4..1e-3, which a 1000-hash sketch quantizes away. A coarser
sketch would silently merge sub-lineages and REINTRODUCE the very discordance this run exists to reduce.

**Honesty rails (inherited, unchanged):** the cohort IS accession-level provenance-disjoint from the
CRyPTIC build (leaked=0; biosample cross-archive overlap 0/30). The lineage-collapsed number is the
honest headline; the raw per-isolate number is clonality-inflated. Callability is unassessed (no regeno),
so the determinant calls remain a CONSERVATIVE lower bound. This script changes the CLUSTERING only --
it does not re-score any genome and does not touch the frozen AMR surface.

Run:  uv run python -m scripts.tb_independent_mash_collapse
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

MASH_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"
SKETCH_SIZE = 10000
MASH_THREADS = 4

# Swept thresholds. 0.001 + 0.005 mirror the frozen provdisjoint lineage layer; the finer rungs
# resolve TB sub-lineage structure the barcode cannot see.
THRESHOLDS: tuple[float, ...] = (0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005)
HEADLINE_THRESHOLD = 0.001  # same rung the frozen provdisjoint layer headlines

WORK = Path(os.environ.get("TB_INDEP_WORK", "D:/dna_decode_cache/tb_indep"))

# Barcode-collapsed reference numbers (wiki/tb_independent_lineage_collapsed_result_2026-07-02.md)
BARCODE_REFERENCE = {
    "rifampicin": {"sens": 0.444, "spec": 0.979, "n_clusters_R": 20, "n_clusters_S": 47, "n_discordant": 43},
    "isoniazid": {"sens": 0.321, "spec": 0.972, "n_clusters_R": 30, "n_clusters_S": 36, "n_discordant": 44},
}


def parse_phylip_lower_triangle(text: str) -> tuple[list[str], np.ndarray]:
    """Parse `mash triangle` relaxed-Phylip lower-triangular output.

    Format:
        <n>
        name1
        name2\td21
        name3\td31\td32
    Returns (names, symmetric NxN float matrix with zero diagonal).
    Names are basenames with a `.fna` suffix stripped.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise ValueError("empty mash triangle output")
    n = int(lines[0].strip())
    names: list[str] = []
    m = np.zeros((n, n), dtype=float)
    for row, ln in enumerate(lines[1 : 1 + n]):
        parts = ln.split("\t")
        raw = parts[0].strip()
        names.append(Path(raw).name.removesuffix(".fna"))
        for col, val in enumerate(parts[1:]):
            d = float(val)
            m[row, col] = d
            m[col, row] = d
    if len(names) != n:
        raise ValueError(f"expected {n} rows, parsed {len(names)}")
    return names, m


def load_results(work: Path) -> dict[str, dict]:
    """{strain_id: {rif_label, inh_label, rif_pred, inh_pred}} from the runner's checkpoint."""
    out: dict[str, dict] = {}
    with open(work / "results.jsonl", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            out[r["strain_id"]] = r
    return out


def run_mash(asm_dir: Path, strain_ids: list[str], scratch: Path,
             reuse: bool = False) -> tuple[list[str], np.ndarray]:
    """Sketch + all-vs-all triangle over `strain_ids`, mounting `asm_dir` READ-ONLY.

    `reuse=True` re-parses an existing `tri.txt` instead of re-sketching (the sketch+triangle is the
    only expensive step; threshold analysis is cheap and iterated).
    """
    from tools.docker_runner import run as drun

    scratch.mkdir(parents=True, exist_ok=True)
    tri = scratch / "tri.txt"
    if reuse and tri.exists() and tri.stat().st_size > 0:
        print(f"[mash] reusing cached triangle {tri}", flush=True)
        return parse_phylip_lower_triangle(tri.read_text(encoding="utf-8"))
    listing = "\n".join(f"/asm/{sid}.fna" for sid in strain_ids) + "\n"
    # LF explicitly: Path.write_text on Windows translates "\n" -> "\r\n", and mash then reads each
    # entry as "/asm/<id>.fna\r" and fails with "could not open" on every genome.
    (scratch / "list.txt").write_bytes(listing.encode("utf-8"))

    mounts = {str(asm_dir): "/asm", str(scratch): "/work"}

    print(f"[mash] sketch {len(strain_ids)} genomes (s={SKETCH_SIZE}, p={MASH_THREADS})", flush=True)
    drun(
        MASH_IMAGE,
        ["mash", "sketch", "-l", "/work/list.txt", "-o", "/work/sketch",
         "-s", str(SKETCH_SIZE), "-p", str(MASH_THREADS)],
        mounts=mounts, capture_output=True, check=True, timeout=14400,
    )

    print("[mash] triangle (all-vs-all lower triangle)", flush=True)
    drun(
        MASH_IMAGE,
        ["sh", "-c",
         f"mash triangle -p {MASH_THREADS} /work/sketch.msh > /work/tri.txt"],
        mounts=mounts, capture_output=True, check=True, timeout=14400,
    )

    return parse_phylip_lower_triangle((scratch / "tri.txt").read_text(encoding="utf-8"))


def collapse_at(matrix, names, results, threshold: float) -> dict:
    from dna_decode.eval.clonality import (
        cluster_weighted_confusion,
        effective_lineage_n,
        greedy_representative_clusters_from_matrix,
        wilson_ci,
    )
    from dna_decode.eval.phylogeny import DistanceMatrix

    dm = DistanceMatrix(strain_ids=names, matrix=matrix)
    clusters = greedy_representative_clusters_from_matrix(dm, threshold)

    out: dict = {"threshold": threshold, "n_clusters_total": len(set(clusters.values())), "drugs": {}}
    for drug, code in (("rifampicin", "rif"), ("isoniazid", "inh")):
        labels = {s: results[s][f"{code}_label"] for s in names
                  if results[s].get(f"{code}_label") in ("R", "S")}
        preds = {s: results[s][f"{code}_pred"] for s in labels}
        sub = {s: c for s, c in clusters.items() if s in labels}
        conf = cluster_weighted_confusion(preds, labels, sub)
        conf["sens_ci95"] = wilson_ci(conf["tp"], conf["tp"] + conf["fn"])
        conf["spec_ci95"] = wilson_ci(conf["tn"], conf["tn"] + conf["fp"])
        conf["effective_lineage_n_R"] = effective_lineage_n(sub, labels, "R")
        conf["effective_lineage_n_S"] = effective_lineage_n(sub, labels, "S")
        conf["n_isolates"] = len(labels)
        out["drugs"][drug] = conf
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", type=Path, default=WORK)
    ap.add_argument("--scratch", type=Path, default=WORK / "mash_collapse")
    ap.add_argument("--output-prefix", type=Path,
                    default=REPO / "wiki" / f"tb_independent_mash_lineage_{_date.today().isoformat()}")
    a = ap.parse_args(argv)

    results = load_results(a.work)
    asm = a.work / "asm"
    strain_ids = sorted(s for s in results if (asm / f"{s}.fna").exists())
    missing = sorted(set(results) - set(strain_ids))
    print(f"[tb-mash] {len(strain_ids)} assemblies present; {len(missing)} of "
          f"{len(results)} scored isolates lack a FASTA", flush=True)
    if len(strain_ids) < 50:
        print("[tb-mash] FAIL too few assemblies", file=sys.stderr)
        return 2

    names, matrix = run_mash(asm, strain_ids, a.scratch)
    print(f"[tb-mash] distance matrix {matrix.shape}", flush=True)

    sweep = [collapse_at(matrix, names, results, t) for t in THRESHOLDS]
    headline = next(s for s in sweep if s["threshold"] == HEADLINE_THRESHOLD)

    payload = {
        "_schema": "tb-independent-mash-lineage-v1",
        "run_date": _date.today().isoformat(),
        "cohort": "EBI AMR-Portal TB provenance-disjoint (leaked=0)",
        "n_isolates_scored": len(results),
        "n_isolates_clustered": len(strain_ids),
        "n_isolates_missing_fasta": len(missing),
        "mash": {"image": MASH_IMAGE, "sketch_size": SKETCH_SIZE, "threads": MASH_THREADS},
        "headline_threshold": HEADLINE_THRESHOLD,
        "threshold_sweep": sweep,
        "headline": headline,
        "barcode_reference": BARCODE_REFERENCE,
        "honesty": (
            "GENUINELY INDEPENDENT (out-of-CRyPTIC-build): accession-level provenance-disjoint, measured "
            "phenotype, WHO catalogue applied UNCHANGED. The lineage-collapsed number is the honest "
            "headline; the raw per-isolate number is clonality-inflated. This run REFINES the clustering "
            "(Mash greedy-representative) vs the coarser pinned Napier barcode; it re-scores nothing."),
        "scope_limits": (
            "Callability unassessed (no regeno) -> determinant calls are a CONSERVATIVE lower bound. "
            "Assembly-available subset (not prevalence-preserving) -> status stays TB_SUBSET_PLUMBING. "
            "asm5 minimap2 VCFs miss some determinants."),
    }

    out_json = a.output_prefix.with_suffix(".json")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        f"# Independent TB cohort — finer Mash lineage collapse ({payload['run_date']})",
        "",
        f"**Cohort:** {payload['cohort']} · **clustered:** {len(strain_ids)} / {len(results)} scored isolates",
        f"**Mash:** `{MASH_IMAGE}` sketch `-s {SKETCH_SIZE}` · greedy-representative (chaining-resistant)",
        "",
        "Refines the coarse Napier-barcode collapse. Frozen clustering + confusion math reused unchanged;",
        "only the lineage definition changes. Nothing is re-scored.",
        "",
        "## Threshold sweep",
        "",
        "| threshold | clusters | drug | lineage sens [95% CI] | lineage spec [95% CI] | R-lin | S-lin | discordant |",
        "|---:|---:|---|---|---|---:|---:|---:|",
    ]
    for s in sweep:
        for drug, c in s["drugs"].items():
            mark = " **" if s["threshold"] == HEADLINE_THRESHOLD else " "
            lines.append(
                f"|{mark}{s['threshold']}{mark.strip() and '**' or ''} | {s['n_clusters_total']} | {drug} | "
                f"{c['sens']} [{c['sens_ci95'][0]}–{c['sens_ci95'][1]}] | "
                f"{c['spec']} [{c['spec_ci95'][0]}–{c['spec_ci95'][1]}] | "
                f"{c['n_clusters_R']} | {c['n_clusters_S']} | {c['n_discordant']} |")

    lines += ["", f"## Headline (threshold {HEADLINE_THRESHOLD})", ""]
    for drug, c in headline["drugs"].items():
        b = BARCODE_REFERENCE[drug]
        lines += [
            f"### {drug}",
            f"- Mash-collapsed: **sens {c['sens']}** [{c['sens_ci95'][0]}–{c['sens_ci95'][1]}] · "
            f"**spec {c['spec']}** [{c['spec_ci95'][0]}–{c['spec_ci95'][1]}]",
            f"- lineages: R={c['n_clusters_R']} S={c['n_clusters_S']} discordant={c['n_discordant']} "
            f"(barcode: R={b['n_clusters_R']} S={b['n_clusters_S']} discordant={b['n_discordant']})",
            f"- barcode reference: sens {b['sens']} · spec {b['spec']}",
            "",
        ]
    lines += ["## Honesty", "", payload["honesty"], "", "## Scope limits", "", payload["scope_limits"], "",
              "Generated by `scripts/tb_independent_mash_collapse.py`."]
    a.output_prefix.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")

    print(f"[tb-mash] headline threshold {HEADLINE_THRESHOLD}: "
          f"{headline['n_clusters_total']} clusters")
    for drug, c in headline["drugs"].items():
        print(f"  {drug}: sens={c['sens']} {c['sens_ci95']} spec={c['spec']} {c['spec_ci95']} "
              f"(R-lin={c['n_clusters_R']} S-lin={c['n_clusters_S']} discordant={c['n_discordant']})")
    print(f"artifact -> {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
