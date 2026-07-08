#!/usr/bin/env python
"""J2 CPU smoke — get a REAL ESM-2 zero-shot number on REAL ProteinGym DMS, FREE, no GPU / no D: / no Kaggle.

The existing J2 runners assume the ProteinGym data is pre-attached AND a GPU is present. This one is the
opposite: it FETCHES a handful of the smallest ProteinGym assays from the public HF mirror
(ICML2022/ProteinGym, v1.0) and runs the CANONICAL scorer (`scripts/esm_zeroshot_dms.py::score_assay`) on
CPU with a small ESM-2 — so any laptop can produce a real zero-shot |Spearman| in minutes and confirm the
whole pipeline works end-to-end on real wet-lab data (not just the synthetic unit tests).

HONEST SCOPE: the smallest assays are BIASED toward hard/atypical short constructs (viral / toxin-antitoxin),
and a small ESM-2 (35M) is the weakest model — so this MINI-median is a pipeline-validation + lower bound,
NOT the published 0.48 (which is ESM2-650M over the FULL 217-assay set). A per-assay rho with a ~0 shuffled
control IS a real correctness signal. Use `--model facebook/esm2_t33_650M_UR50D` for a fairer (slower) number.

Run:
  python scripts/j2_cpu_smoke.py --k 6 --model facebook/esm2_t12_35M_UR50D --data-dir <scratch-dir>
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from esm_zeroshot_dms import AA, score_assay  # noqa: E402  (canonical, drift-guarded scorer)

HF = "https://huggingface.co/datasets/ICML2022/ProteinGym/resolve/main/"
REF = HF + "ProteinGym_reference_file_substitutions.csv"
SUBS = HF + "ProteinGym_substitutions/"


def smallest_assays(ref_rows, k):
    """The k assays with the shortest target_seq that have >=20 single mutants (CPU-cheap; deterministic)."""
    def nm(r):
        try:
            return int(r.get("DMS_number_single_mutants", 0) or 0)
        except ValueError:
            return 0
    ok = [r for r in ref_rows if nm(r) >= 20 and (r.get("target_seq") or "")]
    return sorted(ok, key=lambda r: len(r["target_seq"]))[:k]


def _fetch(url, path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        urllib.request.urlretrieve(url, path)   # HF resolve 302-redirects; urllib follows


def load_dms_local(path):
    from esm_zeroshot_dms import load_dms
    return load_dms(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, help="scratch dir for the fetched reference + assay CSVs")
    ap.add_argument("--model", default="facebook/esm2_t12_35M_UR50D")
    ap.add_argument("--k", type=int, default=6, help="how many of the smallest assays to score")
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--maxlen", type=int, default=1022)
    args = ap.parse_args(argv)

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    os.makedirs(os.path.join(args.data_dir, "subs"), exist_ok=True)
    ref_path = os.path.join(args.data_dir, "ref.csv")
    _fetch(REF, ref_path)
    ref = list(csv.DictReader(open(ref_path, encoding="utf-8")))
    picked = smallest_assays(ref, args.k)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"model={args.model} device={device} assays={len(picked)} (smallest with >=20 single mutants)")
    tok = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model).to(device).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}

    rows = []
    for r in picked:
        fn, seq = r["DMS_filename"], r["target_seq"]
        p = os.path.join(args.data_dir, "subs", fn)
        try:
            _fetch(SUBS + fn, p)
        except Exception as e:                  # noqa: BLE001 (network is best-effort here)
            print(f"  fetch-fail {fn}: {e}")
            continue
        t0 = time.time()
        res, why = score_assay(seq, load_dms_local(p), tok, model, device, mask_id, aa_ids,
                               args.batch, args.maxlen, torch)
        if res:
            rows.append(res)
            print(f"  {r['DMS_id'][:34]:34s} {len(seq):4d}aa  rho={res['rho']:+.3f}  "
                  f"shuf={res['rho_shuf']:+.3f}  ({time.time() - t0:.0f}s)")
        else:
            print(f"  {r['DMS_id'][:34]:34s} skipped ({why})")
    if not rows:
        print("no assays scored"); return 1
    med = float(np.median([abs(x["rho"]) for x in rows]))
    meds = float(np.median([abs(x["rho_shuf"]) for x in rows]))
    print(f"\nMINI-MEDIAN |Spearman| over {len(rows)} smallest assays = {med:.3f} (shuffled {meds:.3f})")
    print("NOTE: smallest-assay + small-model subset — a pipeline-validation + lower bound, NOT the 0.48 "
          "full-benchmark number (see module docstring).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
