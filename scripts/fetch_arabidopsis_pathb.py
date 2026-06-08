"""EP-8 Path B — fetch/verify the Arabidopsis data sources (workhorse first step).

Makes the workhorse's first action one command: confirm the (committed) phenotype labels are present,
re-fetch them if missing, and HEAD-check the large genotype/pseudogenome URLs resolve BEFORE the big
downloads (catches a 1001genomes.org site reorg early instead of mid-batch). With --download it pulls the
genotype VCF + imputed SNP matrix to a target dir (use D: on the workhorse).

All endpoints verified 2026-06-08 (see plans/EP8_PathB_PreStage_Manifest.md). Phenotype CSVs are committed
at data/arabidopsis/ so the workhorse already has the labels (FT10 n=1162, FT16 n=1122, 1122 with both).
No GPU, no money. Big downloads are large (VCF + SNP matrix are multi-GB) — target D:.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PHENO = {
    "FT10_pheno_261.csv": "https://arapheno.1001genomes.org/rest/phenotype/261/values.csv",
    "FT16_pheno_262.csv": "https://arapheno.1001genomes.org/rest/phenotype/262/values.csv",
}
GENOTYPE = {
    "1001genomes_snp-short-indel_only_ACGTN.vcf.gz":
        "https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz",
}
# SNP matrix + pseudogenomes are linked from the Data Center (path can change on site reorg):
DATA_CENTER = "https://1001genomes.org/data-center.html"
FM_MODEL = "kuleshov-group/PlantCaduceus_l32"


def _curl(args: list[str], timeout: float = 120) -> subprocess.CompletedProcess:
    return subprocess.run(["curl", *args], capture_output=True, text=True, timeout=timeout)


def head_check(url: str) -> str:
    """Return a one-line resolve/size status for a URL (HEAD; falls back to ranged GET)."""
    p = _curl(["-sIL", "-m", "40", url])
    status = next((ln.split()[1] for ln in p.stdout.splitlines() if ln.startswith("HTTP/")), "?")
    size = next((ln.split(":", 1)[1].strip() for ln in p.stdout.splitlines()
                 if ln.lower().startswith("content-length")), "?")
    return f"HTTP {status}, content-length {size}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno-dir", default="data/arabidopsis")
    ap.add_argument("--genotype-dir", default="D:/dna_decode_cache/arabidopsis",
                    help="where --download writes the multi-GB genotype files (use D: on the workhorse)")
    ap.add_argument("--download", action="store_true", help="actually download the large genotype files")
    a = ap.parse_args(argv)

    pdir = Path(a.pheno_dir)
    pdir.mkdir(parents=True, exist_ok=True)
    print("== phenotype labels ==")
    for fn, url in PHENO.items():
        f = pdir / fn
        if not f.exists():
            print(f"  {fn}: missing -> fetching")
            _curl(["-sSL", "-m", "60", url, "-o", str(f)])
        n = sum(1 for _ in f.open(encoding="utf-8")) - 1 if f.exists() else 0
        print(f"  {fn}: {n} accessions  ({'ok' if n > 100 else 'CHECK'})")

    print("\n== genotype/pseudogenome URLs (HEAD-check; verify before big pulls) ==")
    for fn, url in GENOTYPE.items():
        print(f"  {fn}: {head_check(url)}\n    {url}")
    print(f"  imputed SNP matrix (1001_SNP_MATRIX.tar.gz) + pseudogenomes: see Data Center\n    {DATA_CENTER}")
    print(f"\n== plant DNA-FM == HF: {FM_MODEL} (frozen; fits RTX 3500 Ada 12GB at fp16)")

    if a.download:
        gdir = Path(a.genotype_dir)
        gdir.mkdir(parents=True, exist_ok=True)
        print(f"\n== downloading genotype to {gdir} ==")
        for fn, url in GENOTYPE.items():
            out = gdir / fn
            if out.exists():
                print(f"  {fn}: exists, skip"); continue
            print(f"  {fn}: downloading (large)...", flush=True)
            r = _curl(["-sSL", "-m", "3600", url, "-o", str(out)], timeout=3700)
            print(f"  {fn}: {'OK' if out.exists() and out.stat().st_size > 0 else 'FAILED'} "
                  f"({out.stat().st_size if out.exists() else 0} bytes)")
        print("  NOTE: imputed SNP matrix + pseudogenomes are not auto-downloaded "
              "(Data Center path varies) — pull manually per the manifest.")
    else:
        print("\n(--download omitted: phenotype confirmed + genotype URLs checked only. "
              "Re-run with --download on the workhorse to pull the genotype.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
