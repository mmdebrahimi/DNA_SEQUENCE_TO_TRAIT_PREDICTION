#!/usr/bin/env python
"""Fetch the ProteinGym v1.1 substitution DMS assays FREE, without the 11 GB download.

The full v1.1 benchmark ships as one 11 GB Zenodo zip (record 13936340), and the HF mirrors only carry the
older v1.0 87-assay set. But the DMS assay CSVs live in a 43 MB NESTED zip inside it, and Zenodo supports
HTTP range requests — so `remotezip` can pull just that inner zip (43 MB) and extract all 217 assays locally
(~1 GB uncompressed), no 11 GB download. This stages the current benchmark for free-CPU scoring
(scripts/j2_cpu_smoke.py / esm_zeroshot_dms.py); the tiny Tsuboyama domains (37–45 aa) are the best free-CPU
substrate (35M-model mini-median ~0.48 on them — see wiki/proteingym_v11_tiny_35M_result.json).

Run:  python scripts/proteingym_v11_fetch.py --out-dir <dir>        # stages ref.csv + subs/*.csv
"""
from __future__ import annotations

import argparse
import io
import os
import time
import zipfile

ZIP_URL = "https://zenodo.org/api/records/13936340/files/ProteinGym_v1.1.zip/content"
INNER = "ProteinGym_v1.1/DMS_ProteinGym_substitutions.zip"
REF = "ProteinGym_v1.1/DMS_substitutions.csv"


def fetch(out_dir):
    from remotezip import RemoteZip
    os.makedirs(os.path.join(out_dir, "subs"), exist_ok=True)
    t0 = time.time()
    with RemoteZip(ZIP_URL) as z:                       # reads central dir via range (~2 s)
        open(os.path.join(out_dir, "ref.csv"), "wb").write(z.read(REF))   # reference (0.2 MB)
        inner = z.read(INNER)                                             # 43 MB nested zip, range-read
    iz = zipfile.ZipFile(io.BytesIO(inner))
    csvs = [n for n in iz.namelist() if n.endswith(".csv")]
    for n in csvs:
        open(os.path.join(out_dir, "subs", n.rsplit("/", 1)[-1]), "wb").write(iz.read(n))
    print(f"staged {len(csvs)} v1.1 substitution assays + ref.csv to {out_dir} "
          f"(range-read {len(inner) // 1_000_000} MB in {time.time() - t0:.0f}s; NOT the 11 GB zip)")
    return len(csvs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    n = fetch(args.out_dir)
    return 0 if n else 1


if __name__ == "__main__":
    raise SystemExit(main())
