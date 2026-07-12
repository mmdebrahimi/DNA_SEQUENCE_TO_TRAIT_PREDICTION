"""Virulence axis sweep — the 2nd multi-axis Family C axis (2026-07-12), run LOCALLY overnight.

Adds the VIRULENCE axis (VirulenceFinder blastn) alongside the plasmid axis, so the multi-axis co-resistance
model can ask: does resistance / a plasmid backbone co-occur with virulence (do ESBL plasmids ride in
virulent E. coli; do ExPEC virulence genes travel with AMR)? Mirrors `scripts/plasmid_axis_sweep.py`:
reuses the committed E. coli VirulenceFinder DB + `pathotype.vf_runner.run_canonical_vf` (blastn 90/60);
scoped to E. coli/Shigella (the VF DB is E. coli). Checkpointed JSONL (restartable). Local blastn on cached
assemblies — no Databricks/Kaggle/fetch.

Run: uv run python scripts/virulence_axis_sweep.py   (writes data/processed/virulence_axis_cache.json)
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pathotype.vf_runner import run_canonical_vf  # noqa: E402

DB = REPO / "data" / "virulencefinder_db" / "virulence_ecoli.fsa"
CACHE = REPO / "data" / "processed" / "virulence_axis_cache.json"
CKPT = REPO / "data" / "processed" / "virulence_axis_checkpoint.jsonl"
ORG = "escherichia_coli_shigella"


def _norm(p): return p.replace(os.sep, "/")


def build_fasta_index():
    idx = {}
    for pat in ("D:/dna_decode_cache/refseq/**/*.fna", "data/raw/**/refseq/**/*.fna", "data/raw/**/*.fna"):
        for p in glob.glob(pat, recursive=True):
            m = re.search(r"(GC[AF]_\d+\.\d+)", _norm(p))
            if m and m.group(1) not in idx:
                idx[m.group(1)] = p
    return idx


def ecoli_accessions():
    accs = set()
    for f in glob.glob(str(REPO / "data" / "raw" / "escherichia_coli_shigella*" / "selected.tsv")):
        for ln in open(f, encoding="utf-8"):
            a = ln.split("\t")[0].strip()
            if a.startswith("GC"):
                accs.add(a)
    return sorted(accs)


def load_done():
    if not CKPT.exists():
        return {}
    return {json.loads(l)["acc"]: json.loads(l) for l in CKPT.read_text(encoding="utf-8").splitlines() if l.strip()}


def _genes(vf_result) -> list[str]:
    """CALLED virulence gene clusters only (per_cluster lists the best hit for EVERY cluster incl. non-called;
    filter on `called` = cleared the identity/coverage thresholds) — the real per-genome virulence profile."""
    pc = vf_result.get("per_cluster") or {}
    return sorted(k for k, v in pc.items() if v.get("called"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max", type=int, default=0)
    a = ap.parse_args(argv)
    if not DB.exists():
        print(f"ERROR: VirulenceFinder E. coli DB absent at {DB}", file=sys.stderr)
        return 2
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    idx = build_fasta_index()
    todo = [acc for acc in ecoli_accessions() if acc in idx]
    if a.max:
        todo = todo[:a.max]
    done = load_done()
    print(f"[vir-sweep] {len(todo)} E. coli/Shigella accessions with a cached FASTA, {len(done)} done", flush=True)
    with open(CKPT, "a", encoding="utf-8") as ck:
        for i, acc in enumerate(todo, 1):
            if acc in done:
                continue
            r = run_canonical_vf(idx[acc], DB)
            genes = _genes(r) if r.get("status") == "ok" else []
            rec = {"acc": acc, "organism": ORG, "status": r.get("status"), "virulence_genes": genes}
            ck.write(json.dumps(rec) + "\n"); ck.flush()
            done[acc] = rec
            if i % 25 == 0 or i == len(todo):
                nok = sum(1 for d in done.values() if d["status"] == "ok")
                print(f"    [{i}/{len(todo)}] {acc} status={r.get('status')} genes={genes[:5]} | {nok} ok", flush=True)
    cache = {acc: {"organism": d["organism"], "virulence_genes": d["virulence_genes"]}
             for acc, d in done.items() if d["status"] == "ok"}
    CACHE.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    n_with = sum(1 for v in cache.values() if v["virulence_genes"])
    print(f"[vir-sweep] DONE: {len(cache)} genomes ok, {n_with} carry >=1 virulence gene -> {CACHE}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
