"""Plasmid Inc-replicon axis sweep — the multi-axis Family C deepening (2026-07-11), run LOCALLY.

C-deep found the AMRFinder virulence CROSS-axis unavailable (AMR-only runs). This adds the PLASMID axis
instead — the canonical co-resistance vehicle: which resistance determinants travel on which plasmid
Inc-type. It is a blastn finder sweep (PlasmidFinder), NOT Bakta/GPU, and every cohort assembly is CACHED
locally, so it runs on this host with the installed blastn — no Databricks/Kaggle/fetch needed.

Scope: ENTEROBACTERALES (E. coli/Shigella, Klebsiella, Salmonella, Enterobacter) — plasmid Inc-typing is an
Enterobacterales concept and `data/plasmidfinder_db/enterobacteriales.fsa` is that DB; these are exactly
Family C's well-populated organisms. Per genome: `plasmid.runner.call_replicons` (blastn 95/60) -> the set
of called Inc-replicons. Checkpointed JSONL (restartable; a wedged blastn loses <=1 genome). Output cache
feeds `scripts/coresistance_multiaxis.py` (determinant/class x Inc-type co-occurrence).

Run: uv run python scripts/plasmid_axis_sweep.py   (writes data/processed/plasmid_axis_cache.json)
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
sys.path.insert(0, str(REPO / "scripts"))

from dna_decode.plasmid.runner import call_replicons  # noqa: E402
import determinant_cooccurrence as dcc  # noqa: E402  (organism_of)

ENTERO = ("escherichia_coli_shigella", "klebsiella", "salmonella", "enterobacter")
DB = REPO / "data" / "plasmidfinder_db" / "enterobacteriales.fsa"
CACHE = REPO / "data" / "processed" / "plasmid_axis_cache.json"
CKPT = REPO / "data" / "processed" / "plasmid_axis_checkpoint.jsonl"


def _norm(p): return p.replace(os.sep, "/")


def build_fasta_index():
    idx = {}
    for pat in ("D:/dna_decode_cache/refseq/**/*.fna", "data/raw/**/refseq/**/*.fna", "data/raw/**/*.fna"):
        for p in glob.glob(pat, recursive=True):
            m = re.search(r"(GC[AF]_\d+\.\d+)", _norm(p))
            if m and m.group(1) not in idx:
                idx[m.group(1)] = p
    return idx


def entero_accessions():
    """{accession: organism} for Enterobacterales cohort genomes."""
    acc_org = {}
    for f in glob.glob(str(REPO / "data" / "raw" / "*" / "selected.tsv")):
        org = dcc.organism_of(_norm(f).replace("/selected.tsv", "/amrfinder_runs/x/main.tsv"))
        if not any(org.startswith(e) for e in ENTERO):
            continue
        for ln in open(f, encoding="utf-8"):
            a = ln.split("\t")[0].strip()
            if a.startswith("GC"):
                acc_org.setdefault(a, org)
    return acc_org


def load_done():
    if not CKPT.exists():
        return {}
    return {json.loads(l)["acc"]: json.loads(l) for l in CKPT.read_text(encoding="utf-8").splitlines() if l.strip()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max", type=int, default=0, help="cap genomes (smoke)")
    a = ap.parse_args(argv)
    if not DB.exists():
        print(f"ERROR: PlasmidFinder DB absent at {DB}", file=sys.stderr)
        return 2
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    idx = build_fasta_index()
    acc_org = entero_accessions()
    todo = [(acc, org) for acc, org in acc_org.items() if acc in idx]
    if a.max:
        todo = todo[:a.max]
    done = load_done()
    print(f"[plasmid-sweep] {len(acc_org)} Enterobacterales accessions, {len(todo)} with a cached FASTA, "
          f"{len(done)} already done", flush=True)
    with open(CKPT, "a", encoding="utf-8") as ck:
        for i, (acc, org) in enumerate(todo, 1):
            if acc in done:
                continue
            r = call_replicons(idx[acc], DB)
            reps = sorted({x["replicon"] for x in r.get("replicons", [])})
            rec = {"acc": acc, "organism": org, "status": r["status"], "replicons": reps}
            ck.write(json.dumps(rec) + "\n"); ck.flush()
            done[acc] = rec
            if i % 25 == 0 or i == len(todo):
                nok = sum(1 for d in done.values() if d["status"] == "ok")
                print(f"    [{i}/{len(todo)}] {acc} ({org}) status={r['status']} reps={reps[:5]} "
                      f"| {nok} ok so far", flush=True)
    # final cache
    cache = {acc: {"organism": d["organism"], "replicons": d["replicons"]}
             for acc, d in done.items() if d["status"] == "ok"}
    CACHE.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    n_ok = len(cache)
    n_with_rep = sum(1 for v in cache.values() if v["replicons"])
    print(f"[plasmid-sweep] DONE: {n_ok} genomes called ok, {n_with_rep} carry >=1 replicon -> {CACHE}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
