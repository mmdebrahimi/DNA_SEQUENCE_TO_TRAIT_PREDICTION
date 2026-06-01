"""Laptop-local pathotype model pipeline (EP-4).

End-to-end on a single machine, no workhorse / no VirulenceFinder / no emails:

    accession CSV -> download genome (GCA) -> Bakta annotate -> parse GFF3
    -> extract CDS -> NT-embed per gene -> mean-pool per strain
    -> binary XGBoost + leave-one-strain-out CV -> AUROC + result packet

Substrate is re-derived from public accessions (data/ is gitignored), so this
script + a committed cohort CSV reproduce the run on any machine. The genome
caches + Bakta DB live on D: (C: is space-constrained).

This is the LEARNED-MODEL track (NT embeddings + classifier), distinct from the
deterministic VirulenceFinder resolver (v0). Discovered constraints baked in:
- NCBI GCA assemblies ship FASTA but no GFF3 and 0 CDS in the GBK -> annotation
  (Bakta) is a MANDATORY stage, not optional.
- Bakta appends an embedded ##FASTA block to its GFF3 that the repo parse_gff3
  cannot handle -> we strip it before parsing.
- The Bakta image entrypoint is already `bakta`, so docker args must NOT repeat it.

Smoke vs real: --max-genes subsamples CDS per strain for a fast pipeline-validation
smoke; omit (or set 0) for a full mean-pool over all CDS.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np

from dna_decode.data.refseq import download_genome
from dna_decode.data import annotations as ann
from dna_decode.models.foundation import model_factory
from tools.docker_runner import run as docker_run, DockerRunnerError

BAKTA_IMAGE = "oschwengers/bakta:v1.11.4"


def bakta_annotate(acc: str, refseq_cache: str, bakta_db: str, bakta_out: str,
                   threads: int = 4, timeout: float = 2400) -> Path:
    """Annotate one genome with Bakta (Docker). Idempotent: skips if GFF3 exists."""
    gff = Path(bakta_out) / f"{acc}.gff3"
    if gff.exists():
        return gff
    Path(bakta_out).mkdir(parents=True, exist_ok=True)
    docker_run(
        image=BAKTA_IMAGE,
        # NOTE: image entrypoint is already `bakta` -> do NOT prepend "bakta".
        args=["--db", "/db/db-light", "--output", "/out", "--prefix", acc,
              "--skip-plot", "--force", "--threads", str(threads), "/data/genome.fna"],
        mounts={
            f"{refseq_cache}/{acc}": "/data:ro",
            bakta_db: "/db:ro",
            bakta_out: "/out",
        },
        capture_output=True, check=True, timeout=timeout,
    )
    if not gff.exists():
        raise RuntimeError(f"Bakta produced no GFF3 for {acc}")
    return gff


def parse_bakta_gff3(gff: Path):
    """Strip Bakta's embedded ##FASTA block, then parse with repo parse_gff3."""
    raw = Path(gff).read_text(encoding="utf-8")
    stripped = raw.split("##FASTA")[0]
    tmp = Path(gff).with_suffix(".nofasta.gff3")
    tmp.write_text(stripped, encoding="utf-8")
    return ann.parse_gff3(str(tmp))


def strain_embedding(model, genome_fasta: Path, ann_table, max_genes: int,
                     rng: np.random.Generator) -> np.ndarray | None:
    """Mean-pool NT embeddings over (a subsample of) the strain's CDS."""
    cds = ann.extract_cds_sequences(str(genome_fasta), ann_table)
    keys = [k for k, s in cds.items() if 60 <= len(s) <= 12000]
    if not keys:
        return None
    if max_genes and len(keys) > max_genes:
        keys = list(rng.choice(np.array(keys, dtype=object), size=max_genes, replace=False))
    vecs = []
    for k in keys:
        v = model.embed(cds[k])
        if getattr(v, "ndim", 1) == 2:
            v = v.mean(axis=0)
        vecs.append(np.asarray(v, dtype=np.float32))
    return np.mean(np.stack(vecs), axis=0) if vecs else None


def loso_auroc(X: np.ndarray, y: np.ndarray):
    """Manual leave-one-strain-out CV; plain XGBoost (calibrate=False at small N)."""
    from xgboost import XGBClassifier
    from sklearn.metrics import roc_auc_score
    n = len(y)
    proba = np.zeros(n, dtype=float)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        clf = XGBClassifier(n_estimators=120, max_depth=3, learning_rate=0.1,
                            subsample=0.9, eval_metric="logloss", n_jobs=4)
        clf.fit(X[tr], y[tr])
        proba[i] = clf.predict_proba(X[i:i+1])[0, 1]
    auroc = float(roc_auc_score(y, proba)) if len(set(y.tolist())) == 2 else float("nan")
    return auroc, proba


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Laptop-local pathotype NT-classifier pipeline")
    ap.add_argument("--cohort", required=True, help="CSV with gca_accession + y columns")
    ap.add_argument("--refseq-cache", default="D:/dna_decode_cache/refseq")
    ap.add_argument("--bakta-db", default="C:/Users/Farshad/dna_decode_stage2/bakta_db")
    ap.add_argument("--bakta-out", default="D:/dna_decode_cache/bakta_out")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--max-genes", type=int, default=200,
                    help="Subsample CDS per strain (0 = all; smoke uses ~200)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True, help="Result packet path (.json; .md sibling written)")
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    rows = list(csv.DictReader(open(args.cohort, encoding="utf-8")))
    print(f"[pipeline] cohort={args.cohort} n={len(rows)} max_genes={args.max_genes} device={args.device}")

    print("[pipeline] loading NT model...")
    model = model_factory("nucleotide_transformer",
                          config_path=str(REPO / "config" / "datasources.yaml"),
                          device=args.device)

    X, y, ids, per = [], [], [], []
    for r in rows:
        acc = r["gca_accession"]
        label = int(r["y"])
        print(f"[pipeline] {acc} (y={label}) download...", flush=True)
        download_genome(acc, args.refseq_cache)
        print(f"[pipeline] {acc} bakta annotate...", flush=True)
        gff = bakta_annotate(acc, args.refseq_cache, args.bakta_db, args.bakta_out)
        atab = parse_bakta_gff3(gff)
        gfna = Path(args.refseq_cache) / acc / "genome.fna"
        print(f"[pipeline] {acc} NT embed...", flush=True)
        vec = strain_embedding(model, gfna, atab, args.max_genes, rng)
        if vec is None:
            print(f"[pipeline] WARN {acc}: no usable CDS, skipped")
            continue
        X.append(vec); y.append(label); ids.append(acc)
        per.append({"acc": acc, "y": label, "n_cds": int((atab["type"] == "CDS").sum())})

    X = np.stack(X); y = np.array(y)
    print(f"[pipeline] embedded {len(y)} strains; X={X.shape}; running LOSO...")
    auroc, proba = loso_auroc(X, y)
    print(f"[pipeline] LOSO AUROC = {auroc:.4f}")

    result = {
        "cohort": args.cohort, "n_strains": int(len(y)),
        "n_pos": int(y.sum()), "n_neg": int((y == 0).sum()),
        "max_genes_per_strain": args.max_genes, "embedding_dim": int(X.shape[1]),
        "cv": "leave_one_strain_out", "loso_auroc": auroc,
        "model": "nucleotide_transformer_v2_100m", "device": args.device,
        "per_strain": [{**p, "held_out_proba": float(proba[i])} for i, p in enumerate(per)],
        "caveat": "Pipeline-validation smoke if max_genes>0 (subsampled mean-pool). "
                  "Toxin-label (LT/ST) contrast is partly circular for ETEC -> "
                  "resolver-conformance-adjacent, not strong external validity.",
    }
    outp = Path(args.out); outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(result, indent=2), encoding="utf-8")
    md = outp.with_suffix(".md")
    md.write_text(
        f"# Pathotype Laptop Pipeline Result\n\n"
        f"- cohort: `{args.cohort}`\n- n={len(y)} ({int(y.sum())} pos / {int((y==0).sum())} neg)\n"
        f"- model: NT v2 100M | max_genes/strain: {args.max_genes} | device: {args.device}\n"
        f"- **LOSO AUROC: {auroc:.4f}**\n\n{result['caveat']}\n",
        encoding="utf-8")
    print(f"[pipeline] wrote {outp} + {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
