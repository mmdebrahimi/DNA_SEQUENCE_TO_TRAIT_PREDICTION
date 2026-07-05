#!/usr/bin/env python
"""capture_pgp_uk.py - index-capture the PGP-UK open individual-level human cohort.

WHY (the data gap): DNA-decode's individual-level human path is the GATED UK Biobank application
(pending). PGP-UK is a FREE, OPEN-CONSENT, NON-application-gated individual-level human multi-omics
cohort (WGS + WGBS methylation + RNA-seq + self-reported phenotype), deposited in OPEN EMBL-EBI
repositories (ENA / ArrayExpress) — NOT dbGaP/EGA. So it unblocks individual-level human work WITHOUT
waiting on a data-access committee. Verified live real-surface 2026-07-04 (ENA PRJEB17529 returns real
FTP download URLs).

This captures the INDEX (accessions + per-run FTP links + sample pointers) like the sibling
capture_clinvar/capture_1000g scripts — it does NOT bulk-download the ~2 TB raw panel (that is a
separate, explicit, D:-targeted fetch). Manifest lands on D: (heavy-artifacts rule), code stays on C:.

Open accessions (verified in the open repos, no gate):
  - ENA PRJEB17529  : WGS + WGBS (whole-genome + bisulfite methylation)
  - ENA PRJEB25139  : RNA-seq  (also ArrayExpress E-MTAB-6523)
  - ArrayExpress E-MTAB-5377 : raw 450k methylation IDATs (pointer only; not an ENA read_run)
  - phenotype + genome/methylome reports: https://www.personalgenomes.org.uk/data

Run:  python capture_pgp_uk.py
"""
from __future__ import annotations

import csv
import io
import os
import sys
import urllib.request

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT_DIR = os.environ.get("PGP_UK_OUT", "D:/dna_decode_cache/pgp_uk")
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
READ_RUN_ACCESSIONS = ["PRJEB17529", "PRJEB25139"]  # WGS/WGBS, RNA-seq (open, no gate)
POINTERS = {
    "methylation_450k_idat": "ArrayExpress E-MTAB-5377",
    "rnaseq_arrayexpress": "ArrayExpress E-MTAB-6523",
    "phenotype_and_reports": "https://www.personalgenomes.org.uk/data",
    "api": "https://www.personalgenomes.org/api",
}
FIELDS = "run_accession,sample_title,fastq_ftp,submitted_ftp,library_strategy,base_count"


def fetch_ena_index(accession, timeout=60):
    url = (f"{ENA_API}?accession={accession}&result=read_run&fields={FIELDS}"
           f"&format=tsv&limit=0")
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode-capture/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310 (fixed EBI host)
        text = r.read().decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    manifest_path = os.path.join(OUT_DIR, "pgp_uk_manifest.tsv")
    all_rows, per_acc = [], {}
    print("=== PGP-UK index capture (open, non-gated individual-level human cohort) ===")
    for acc in READ_RUN_ACCESSIONS:
        try:
            rows = fetch_ena_index(acc)
        except Exception as e:  # noqa: BLE001
            print(f"  {acc}: FETCH FAILED ({e}) — endpoint unreachable, skipped")
            per_acc[acc] = 0
            continue
        for row in rows:
            row["_project"] = acc
        all_rows.extend(rows)
        per_acc[acc] = len(rows)
        dl = next((r.get("fastq_ftp") or r.get("submitted_ftp") for r in rows
                   if (r.get("fastq_ftp") or r.get("submitted_ftp"))), "")
        first_dl = dl.split(";")[0] if dl else "(none)"
        print(f"  {acc}: {len(rows)} runs  e.g. {(rows[0]['run_accession'] if rows else '-')}  "
              f"dl={first_dl}")

    if all_rows:
        cols = ["_project", "run_accession", "sample_title", "fastq_ftp", "submitted_ftp",
                "library_strategy", "base_count"]
        with open(manifest_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(all_rows)
        print(f"\nmanifest: {manifest_path}  ({len(all_rows)} run records)")

    # write the pointer sidecar (non-ENA data + phenotype page + API)
    ptr_path = os.path.join(OUT_DIR, "pgp_uk_pointers.tsv")
    with open(ptr_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["key", "value"])
        for k, v in POINTERS.items():
            w.writerow([k, v])
    print(f"pointers: {ptr_path}  ({len(POINTERS)} entries)")

    total = sum(per_acc.values())
    ok = total > 0
    print(f"\n{'PASS' if ok else 'FAIL'}: captured {total} downloadable run records across "
          f"{len([a for a,n in per_acc.items() if n])} open ENA projects "
          f"(+ {len(POINTERS)} pointers). Individual-level human, NO application gate.")
    print("  Raw panel is ~2 TB — bulk fetch is a separate explicit D:-targeted step, not this index.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
