"""ExPEC-vs-EPEC k-mer baseline arm of the representation bake-off (EP-4).

The MANDATED classical control (CLAUDE.md: FM must beat classical by >=3pp) and
the confound diagnostic: ExPEC(Salipante) vs EPEC(Hazen) is study-confounded, so
a near-perfect k-mer AUROC means the contrast is trivially/batch-separable and
NT's value is unmeasurable on this substrate; a moderate k-mer AUROC leaves room
for NT lift.

Cheap: assembly FASTA only (ENA WGS-set fetch, no GCA, no Bakta, no GPU).
Non-circular labels: ExPEC=isolation-site, EPEC=DECA-curated.
"""
from __future__ import annotations
import csv, gzip, io, sys, json
from pathlib import Path
import urllib.request

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
import numpy as np
from sklearn.metrics import roc_auc_score
from dna_decode.models.classical_baselines import CONTIG_SEPARATOR
from dna_decode.eval.loso_kmer import run_kmer_xgboost_loso

F1 = REPO / "data/external/horesh2021_F1_genome_metadata.csv"
ENA_CACHE = Path("C:/Users/Farshad/PythonProjects/dna_decode/data/ena_wgs")
N_PER_CLASS = 12


import re
WGS_MASTER = re.compile(r"^[A-Z]{4}\d{8}(\.fa)?$", re.I)  # AAAA00000000(.fa) -> WGS set


def select_strains():
    rows = list(csv.DictReader(open(F1, encoding="utf-8")))
    clean = [r for r in rows if "(predicted)" not in r["Pathotype"] and r["Pathotype"] not in ("Not determined", "")]
    # Only genuine WGS-master Assembly_names are ENA-fetchable; some Salipante rows
    # carry SRA run ids (SRR...) which the prefix logic mangles -> filter them out.
    def ok(r): return bool(WGS_MASTER.match(r["Assembly_name"].strip()))
    expec = [r for r in clean if r["Pathotype"].strip().startswith("ExPEC") and r["Source"].startswith("Salipante") and ok(r)][:N_PER_CLASS]
    epec = [r for r in clean if r["Pathotype"].strip().startswith("EPEC") and r["Source"].startswith("Hazen") and ok(r)][:N_PER_CLASS]
    out = []
    for r in expec: out.append((r["ID"], r["Assembly_name"], "ExPEC", 0, r["Source"]))
    for r in epec: out.append((r["ID"], r["Assembly_name"], "EPEC", 1, r["Source"]))
    return out


def wgs_set_prefix(assembly_name: str) -> str:
    base = assembly_name.replace(".fa", "").replace(".fasta", "").strip()
    return base[:4] + "01"  # master AAAA00000000 -> contig set AAAA01


def fetch_assembly(setp: str) -> str | None:
    ENA_CACHE.mkdir(parents=True, exist_ok=True)
    cache = ENA_CACHE / f"{setp}.fna"
    if cache.exists() and cache.stat().st_size > 1000:
        return cache.read_text(encoding="utf-8")
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{setp}"  # bare form (params break it)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        raw = urllib.request.urlopen(req, timeout=120).read()
    except Exception as e:
        print(f"  FETCH FAIL {setp}: {type(e).__name__} {str(e)[:80]}")
        return None
    if raw[:2] == b"\x1f\x8b":  # gzip magic
        text = gzip.decompress(raw).decode("utf-8", "replace")
    else:
        text = raw.decode("utf-8", "replace")
    if ">" not in text:
        print(f"  NO FASTA in {setp} (len={len(text)})")
        return None
    cache.write_text(text, encoding="utf-8")
    return text


def concat_contigs(fasta_text: str) -> str:
    seqs, cur = [], []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return CONTIG_SEPARATOR.join(seqs)


def main() -> int:
    strains = select_strains()
    print(f"[kmer] selected {len(strains)} strains "
          f"({sum(1 for s in strains if s[3]==0)} ExPEC / {sum(1 for s in strains if s[3]==1)} EPEC)")
    seqs_by_strain, labels_by_strain, ids = {}, {}, []
    for sid, aname, cls, y, src in strains:
        setp = wgs_set_prefix(aname)
        print(f"[kmer] {sid} ({cls}, {src.split()[0]}) <- ENA {setp}", flush=True)
        ft = fetch_assembly(setp)
        if ft is None:
            print(f"  skip {sid}")
            continue
        seq = concat_contigs(ft)
        if len(seq) < 1_000_000:
            print(f"  skip {sid}: assembly too short ({len(seq)} bp)")
            continue
        key = f"{sid}|{setp}"
        seqs_by_strain[key] = seq
        labels_by_strain[key] = y
        ids.append(key)
    npos = sum(labels_by_strain.values()); nneg = len(ids) - npos
    print(f"[kmer] fetched {len(ids)} usable ({nneg} ExPEC / {npos} EPEC); running LOSO k-mer...")
    if npos < 3 or nneg < 3:
        print("[kmer] ABORT: too few per class after fetch"); return 1
    cv = run_kmer_xgboost_loso(seqs_by_strain, labels_by_strain, ids, drug="expec_vs_epec", k=8)
    yt, ys = cv.all_y_true(), cv.all_y_score()
    auroc = float(roc_auc_score(yt, ys)) if len(set(yt.tolist())) == 2 else float("nan")
    uniq = sorted(set(round(float(s), 4) for s in ys))
    degenerate = len(uniq) <= 2
    print(f"[kmer] LOSO AUROC = {auroc:.4f} | distinct scores={len(uniq)} degenerate={degenerate}")
    res = {"contrast": "ExPEC(Salipante) vs EPEC(Hazen)", "n": len(ids),
           "n_expec": nneg, "n_epec": npos, "representation": "kmer_k8_xgboost",
           "cv": "loso", "auroc": auroc, "distinct_scores": len(uniq), "degenerate": degenerate,
           "confound": "study-confounded (ExPEC=Salipante, EPEC=Hazen) -> high AUROC may be batch, not biology",
           "ids": ids}
    out = REPO / "research_outputs/pathotype_bakeoff_kmer_expec_epec_2026-05-31.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[kmer] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
