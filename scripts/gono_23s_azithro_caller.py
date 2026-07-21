"""Custom 23S rRNA caller for N. gonorrhoeae azithromycin — the determinant AMRFinder does NOT report.

AMRFinderPlus (and thus NCBI-PD) does not call gono 23S rRNA point mutations (23S is a 4-copy rRNA;
AMRFinder is protein-centric), so the azithromycin cell scored DEGENERATE on the NCBI-PD external validation
(0/110 R carried any 23S marker). This closes that gap with a native-BLAST caller (NO Docker):

  per isolate: download the GCA assembly (ENA, cached) -> blastn a full E. coli 23S reference vs it ->
  walk the alignment to E. coli macrolide positions 2059 + 2611 -> read the assembled base. A2059G or
  C2611T (in ANY assembled 23S copy) -> azithromycin R.

**INHERENT CEILING (honest):** 23S is multi-copy and azithromycin resistance is often heteroplasmic (a
subset of the 4 copies mutated). The consensus assembly collapses the copies, so only high-level (all/most
copies mutated) R is reliably detectable from WGS -- a well-known literature limitation, analogous to tet(M)
detecting only high-level TRNG. Expect high specificity, ceiling-limited sensitivity.

  uv run python scripts/gono_23s_azithro_caller.py --limit 10   # probe
  uv run python scripts/gono_23s_azithro_caller.py              # full azithro cohort
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.amr_portal_score_independent import wilson_ci  # noqa: E402
from scripts.independent_cohort_validate import _conf  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
CACHE = Path("D:/dna_decode_cache/gono_asm")
REF_23S = None  # set at runtime (fetched to scratch)
POSITIONS = {2059: ("A", "G"), 2611: ("C", "T")}  # E. coli 23S macrolide sites: (WT, resistant)


def fetch_ref(scratch: Path) -> Path:
    ref = scratch / "ecoli23S_full.fa"
    if not ref.exists() or ref.stat().st_size == 0:
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=V00331.1"
               "&rettype=fasta&retmode=text")
        ref.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    return ref


def download_asm(gca: str) -> Path | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    dest = CACHE / f"{gca}.fna"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{gca}?download=true"
    try:
        data = urllib.request.urlopen(url, timeout=180).read()
        if len(data) < 1000:
            return None
        dest.write_bytes(data)
        return dest
    except Exception:  # noqa: BLE001
        return None


def _base_at(qseq: str, sseq: str, qstart: int, target: int) -> str | None:
    """Walk the aligned query/subject seqs; return the subject base aligned to query position `target`."""
    qpos = qstart - 1
    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            qpos += 1
        if qpos == target and qc != "-":
            return sc.upper()
    return None


def call_23s(asm: Path, ref: Path) -> dict:
    """blastn 23S ref vs assembly; read the base at each macrolide position across ALL aligned 23S copies."""
    out = subprocess.run(
        [BLASTN, "-query", str(ref), "-subject", str(asm),
         "-outfmt", "6 qstart qend sstart send sstrand qseq sseq", "-max_target_seqs", "20"],
        capture_output=True, text=True, timeout=120)
    calls = {p: set() for p in POSITIONS}
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 7:
            continue
        qstart, qseq, sseq = int(f[0]), f[5], f[6]
        for pos in POSITIONS:
            if qstart <= pos <= int(f[1]):
                b = _base_at(qseq, sseq, qstart, pos)
                if b and b != "-":
                    calls[pos].add(b)
    muts = []
    for pos, (wt, res) in POSITIONS.items():
        if res in calls[pos]:
            muts.append(f"23S_{wt}{pos}{res}")
    return {"prediction": "R" if muts else "S", "matched_23S": muts,
            "bases_seen": {p: sorted(calls[p]) for p in POSITIONS}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", default="data/raw/gono_ncbipd_extval/cohort.tsv")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    scratch = Path("C:/Users/Farshad/AppData/Local/Temp/claude/"
                   "C--Users-Farshad-PythonProjects-dna-decode/29b7c779-fe9b-4a00-ac93-2e6f871e69ca/scratchpad")
    scratch.mkdir(parents=True, exist_ok=True)
    ref = fetch_ref(scratch)
    rows = [r for r in csv.DictReader(open(a.cohort), delimiter="\t") if r.get("azithromycin") in ("R", "S")]
    if a.limit:  # balanced probe subset
        R = [r for r in rows if r["azithromycin"] == "R"][:a.limit // 2]
        S = [r for r in rows if r["azithromycin"] == "S"][:a.limit // 2]
        rows = R + S
    scored, n_missing = [], 0
    for r in rows:
        asm = download_asm(r["asm_acc"])
        if asm is None:
            n_missing += 1
            continue
        c = call_23s(asm, ref)
        pred, lab = c["prediction"], r["azithromycin"]
        scored.append((pred, 1 if lab == "R" else 0))
        print(f"  {r['biosample']} {r['asm_acc']}: lab={lab} pred={pred} 23S={c['matched_23S']} "
              f"bases={c['bases_seen']}", flush=True)
    conf = _conf(scored)
    n_R, n_S = conf["tp"] + conf["fn"], conf["tn"] + conf["fp"]
    print(f"\nazithromycin (23S caller): n={conf['n_scored']} ({n_R}R/{n_S}S) acc={conf['acc']} "
          f"sens={conf['sens']} spec={conf['spec']} missing={n_missing}")
    art = {"_schema": "gono-23s-azithro-v1", "date": _date.today().isoformat(),
           "caller": "native blastn E.coli-23S(V00331) vs assembly @ E.coli 2059/2611",
           "binary": conf, "n_R": n_R, "n_S": n_S, "n_missing": n_missing,
           "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
           "ceiling": "23S multi-copy heteroplasmy -> only high-level all-copy R reliably WGS-detectable"}
    Path(f"wiki/gono_23s_azithro_{_date.today().isoformat()}.json").write_text(json.dumps(art, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
