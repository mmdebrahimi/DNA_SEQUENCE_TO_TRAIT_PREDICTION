#!/usr/bin/env python
"""MaveDB CPU smoke — score a MaveDB protein-missense DMS with real ESM-2 on CPU (free, no GPU/D:).

A SECOND, independent variant-effect data source beyond ProteinGym: fetch a MaveDB score set (metadata +
scores CSV) via the public API, adapt its `hgvs_pro` variants + target protein sequence into the canonical
scorer's input, and run `scripts/esm_zeroshot_dms.py::score_assay` on CPU. Proves the world-model pipeline
generalizes to an independent DB, and captures MaveDB as a reusable substrate.

Target sequence handling: MaveDB targets are `protein` (used directly) or `dna` (translated to protein here,
so the numbering matches `hgvs_pro`). Only clean single missense variants (p.Xaa###Yaa, both real AAs,
wt!=mut) with a finite score are kept; the scorer's own reference-mismatch guard drops any that don't match.

HONEST SCOPE: MaveDB largely OVERLAPS ProteinGym (ProteinGym curates from MaveDB + others), so this is a
cross-DB pipeline check + an independent-ingest capability, not a fully independent benchmark. A per-assay
rho with a ~0 shuffled control is a real correctness signal.

Run:  python scripts/mavedb_cpu_smoke.py --urn urn:mavedb:00000001-a-1 --out wiki/mavedb_ube2i_35M_result.json
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from esm_zeroshot_dms import AA, score_assay  # noqa: E402  (canonical, drift-guarded scorer)

API = "https://api.mavedb.org/api/v1/score-sets/"
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
_CODON = {  # standard genetic code
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
_MISSENSE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def translate(dna):
    """DNA -> protein (standard code), stopping at the first stop codon."""
    dna = dna.upper().replace("U", "T")
    out = []
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON.get(dna[i:i + 3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def parse_hgvs_pro(h):
    """'p.Met1Ala' -> ('M',1,'A') for a clean single missense; None otherwise (syn/nonsense/multi/indel)."""
    m = _MISSENSE.match((h or "").strip())
    if not m:
        return None
    wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
    wt, mut = AA3.get(wt3), AA3.get(mut3)
    if not wt or not mut or wt == mut:      # both must be real AAs and differ (drops Ter/=, syn)
        return None
    return wt, pos, mut


def _get(url, timeout=120, retries=6):
    """GET with retry+backoff — MaveDB's on-demand /scores export intermittently returns 504."""
    import time
    last = None
    for i in range(retries):
        try:
            return urllib.request.urlopen(url, timeout=timeout).read()
        except Exception as e:                       # noqa: BLE001 (504/timeout/transient — retry)
            last = e
            print(f"  fetch attempt {i + 1}/{retries} failed ({str(e)[:50]}); retrying...")
            time.sleep(3 * (i + 1))
    raise last


def fetch_target_protein(urn):
    """Return the target PROTEIN sequence (translating if the MaveDB target is DNA)."""
    meta = json.loads(_get(API + urn, timeout=40))
    ts = (meta.get("targetGenes", [{}])[0].get("targetSequence") or {})
    seq, styp = (ts.get("sequence") or ""), (ts.get("sequenceType") or "").lower()
    return translate(seq) if styp == "dna" else seq, meta.get("title", "")


def build_dms(scores_csv_text):
    """MaveDB scores CSV text -> {'<wt><pos><mut>': score} for clean single missense with finite score."""
    out = {}
    for row in csv.DictReader(io.StringIO(scores_csv_text)):
        v = parse_hgvs_pro(row.get("hgvs_pro", ""))
        if not v:
            continue
        try:
            s = float(row.get("score", ""))
        except (ValueError, TypeError):
            continue
        if not np.isfinite(s):
            continue
        wt, pos, mut = v
        out[f"{wt}{pos}{mut}"] = s
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--urn", default="urn:mavedb:00000001-a-1")
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--maxlen", type=int, default=1022)
    ap.add_argument("--out", default="")
    ap.add_argument("--scores-csv", default="",
                    help="local MaveDB scores CSV (fallback when the /scores API is 504-degraded)")
    args = ap.parse_args(argv)

    prot, title = fetch_target_protein(args.urn)
    if args.scores_csv:
        scores_text = open(args.scores_csv, encoding="utf-8", errors="replace").read()
        print(f"  using local scores CSV: {args.scores_csv}")
    else:
        scores_text = _get(API + args.urn + "/scores", timeout=180).decode("utf-8", "replace")
    dms = build_dms(scores_text)
    print(f"urn={args.urn} title={title[:44]!r} protein_len={len(prot)} single_missense={len(dms)}")

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    device = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}
    res, why = score_assay(prot, dms, tok, model, device, mask_id, aa_ids, args.batch, args.maxlen, torch)
    if not res:
        print(f"NOT SCORED ({why}) — mism means hgvs_pro numbering != translated target"); return 1
    print(f"MaveDB {args.urn}: n={res['n']} rho={res['rho']:+.3f} shuffled={res['rho_shuf']:+.3f} "
          f"mism={res['mism']}")
    if args.out:
        json.dump({"source": "MaveDB", "urn": args.urn, "title": title, "model": args.model,
                   "protein_len": len(prot), **{k: res[k] for k in ("n", "rho", "rho_shuf", "mism")}},
                  open(args.out, "w", encoding="utf-8"), indent=2)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
