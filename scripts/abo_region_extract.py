"""Extract the full ABO-region genotype matrix for the J3 substrate samples — the git-transferable input.

WHY: the workhorse cannot mount external drives and git cannot carry the 21 GB OpenSNP zip. But the workhorse
does NOT need the zip — it only needs the ABO-locus genotypes for the 395 J3 substrate samples, which is tiny
and travels through the mmdebrahimi git repo. This extracts every array SNP in a generous ABO window (GRCh37)
for each substrate sample, long-format + lossless, aligned to the substrate `user_id` index.

The workhorse builds `frozen_fm` / `learned_rep` from (a) this per-sample ABO-window genotype matrix + (b) a
GRCh37 reference sequence it fetches itself (a standard public download — NOT gated by the zip). No external
drive, no large transfer.

Window: ABO gene GRCh37 chr9:136,125,788-136,150,617 + flanks -> chr9:136,120,000-136,155,000 (35 kb). The
23andMe dumps are build 37; genotypes are plus-strand-oriented per the file header.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import zipfile
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

DEFAULT_ZIP = Path("D:/dna_decode_cache/opensnp/opensnp_datadump.2017-12-08.zip")
ABO_CHROM = "9"
ABO_START = 136_120_000
ABO_END = 136_155_000


def _region_from_member(zf: zipfile.ZipFile, member: str) -> list[tuple[str, int, str]]:
    """Stream a 23andMe/AncestryDNA member -> [(rsid, pos, genotype)] for chr9 ABO-window SNPs."""
    out: list[tuple[str, int, str]] = []
    try:
        with zf.open(member) as fh:
            for bline in fh:
                line = bline.decode("utf-8", errors="replace")
                if not line or line[0] == "#" or line.startswith("rsid") or line.startswith('"rsid'):
                    continue
                parts = line.rstrip("\n").replace('"', "").replace(",", "\t").split("\t")
                if len(parts) < 4:
                    continue
                chrom = parts[1].strip()
                if chrom != ABO_CHROM or not parts[2].strip().isdigit():
                    continue
                pos = int(parts[2].strip())
                if not (ABO_START <= pos <= ABO_END):
                    continue
                if len(parts) >= 5:                        # AncestryDNA allele1/allele2
                    gt = (parts[3].strip() + parts[4].strip()).upper()
                else:
                    gt = parts[3].strip().upper()
                out.append((parts[0].strip(), pos, gt))
    except (KeyError, zipfile.BadZipFile, OSError):
        return out
    return out


def run(zip_path: Path, substrate_json: Path) -> dict:
    if not zip_path.exists():
        return {"status": "ZIP_NOT_PRESENT", "zip": str(zip_path)}
    if not substrate_json.exists():
        return {"status": "NO_SUBSTRATE", "hint": "run scripts/abo_opensnp_ingest.py first"}
    substrate = json.loads(substrate_json.read_text(encoding="utf-8"))
    samples = substrate["samples"]
    zf = zipfile.ZipFile(str(zip_path))
    rows = []
    positions = set()
    n_samples_with_snps = 0
    for s in samples:
        member = s["genotype_member"]
        region = _region_from_member(zf, member)
        if region:
            n_samples_with_snps += 1
        for rsid, pos, gt in region:
            positions.add(pos)
            rows.append((s["user_id"], rsid, ABO_CHROM, pos, gt))
    return {
        "status": "EXTRACTED" if rows else "NO_SNPS",
        "schema": "j3-abo-region-genotypes-v1", "date": _date.today().isoformat(),
        "assembly": "GRCh37", "window": f"chr{ABO_CHROM}:{ABO_START}-{ABO_END}",
        "n_samples": len(samples), "n_samples_with_region_snps": n_samples_with_snps,
        "n_distinct_positions": len(positions), "n_genotype_rows": len(rows),
        "note": ("git-transferable ABO-window genotype matrix aligned to the substrate user_id index; the "
                 "workhorse embeds this + a self-fetched GRCh37 reference. No zip transfer / no external drive."),
        "_rows": rows,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--substrate", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_substrate.json")
    ap.add_argument("--out-tsv", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_region_genotypes.tsv")
    ap.add_argument("--out-json", type=Path, default=REPO / "data" / "j3_abo" / "j3_abo_region_manifest.json")
    a = ap.parse_args(argv)
    res = run(a.zip, a.substrate)
    if res.get("status") == "EXTRACTED":
        a.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        with a.out_tsv.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["user_id", "rsid", "chrom", "pos", "genotype"])
            for row in sorted(res["_rows"], key=lambda r: (r[0], r[3])):
                w.writerow(row)
        manifest = {k: v for k, v in res.items() if k != "_rows"}
        a.out_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        print(f"[wrote {a.out_tsv} + {a.out_json}]")
    print(json.dumps({k: v for k, v in res.items() if k != "_rows"}, indent=2))
    return 0 if res.get("status") == "EXTRACTED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
