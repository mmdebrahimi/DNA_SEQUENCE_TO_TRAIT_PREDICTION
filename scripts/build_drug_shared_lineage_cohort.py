"""Build a shared-lineage cohort from a fetch manifest: download genomes + run AMRFinder + emit a parquet.

Restartable (skip-existing -> Docker-corruption-safe): re-running picks up where it left off. Genomes
cache on D: (disk-tight C:); AMRFinder TSVs land in data/amrfinder_runs/<accession>/. The emitted parquet
(strain_id / assembly_accession / mlst / ast_<drug>) is the substrate for scripts/functional_alphabet_probe.py.

Usage:
  DNA_DECODE_AMRFINDER_DB=C:/Users/Farshad/dna_decode_stage2/amrfinder_db \
  uv run python scripts/build_drug_shared_lineage_cohort.py \
    --manifest wiki/dna_llm_shared_lineage_manifest_tetracycline_2026-06-26.json \
    --drug tetracycline --refseq-cache D:/dna_decode_cache/refseq [--limit N]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--drug", required=True)
    ap.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    ap.add_argument("--runs-root", type=Path, default=ROOT / "data/amrfinder_runs")
    ap.add_argument("--organism", default="Escherichia")
    ap.add_argument("--limit", type=int, default=None, help="pilot: only the first N strains")
    ap.add_argument("--download-only", action="store_true",
                    help="fetch genomes to D: only (skip AMRFinder); for pre-fetching before offline")
    ap.add_argument("--out-parquet", type=Path, default=None)
    args = ap.parse_args(argv)

    from dna_decode.data.refseq import download_genome
    from scripts.drug_mechanism_audit import _run_amrfinder

    man = json.loads(args.manifest.read_text(encoding="utf-8"))
    strains = man["strains"][: args.limit] if args.limit else man["strains"]
    args.refseq_cache.mkdir(parents=True, exist_ok=True)
    args.runs_root.mkdir(parents=True, exist_ok=True)

    rows, done, skipped, failed = [], 0, 0, 0
    for i, s in enumerate(strains, 1):
        acc, mlst, lab = s["assembly_accession"], s["mlst"], 1 if s["label"] == "R" else 0
        run_dir = args.runs_root / acc
        gpath = args.refseq_cache / acc / "genome.fna"
        try:
            if not gpath.exists():
                t0 = time.time()
                download_genome(acc, args.refseq_cache)
                print(f"[{i}/{len(strains)}] {acc} ({s['label']}) download {time.time()-t0:.0f}s", flush=True)
                done += 1
            else:
                skipped += 1
            if args.download_only:
                # genome-only pre-fetch: write the cohort row regardless (AMRFinder runs offline later)
                if gpath.exists():
                    rows.append({"strain_id": acc, "assembly_accession": acc, "mlst": mlst,
                                 f"ast_{args.drug}": lab})
                continue
            if not (run_dir / "main.tsv").exists():
                run_dir.mkdir(parents=True, exist_ok=True)
                _run_amrfinder(gpath, run_dir, organism=args.organism)
            if (run_dir / "main.tsv").exists():
                rows.append({"strain_id": acc, "assembly_accession": acc, "mlst": mlst,
                             f"ast_{args.drug}": lab})
        except Exception as e:  # restartable: one bad accession doesn't kill the batch
            print(f"[{i}/{len(strains)}] {acc} FAILED: {type(e).__name__}: {e}", flush=True)
            failed += 1

    df = pd.DataFrame(rows)
    out = args.out_parquet or (ROOT / "data/processed" /
                               f"shared_lineage_{args.drug}_cohort.parquet")
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    n_R = int((df[f"ast_{args.drug}"] == 1).sum()) if len(df) else 0
    n_S = int((df[f"ast_{args.drug}"] == 0).sum()) if len(df) else 0
    print(f"COHORT: {len(df)} strains ({n_R}R/{n_S}S) | amrfinder_new={done} skipped={skipped} "
          f"failed={failed} | unique mlst={df['mlst'].nunique() if len(df) else 0} -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
