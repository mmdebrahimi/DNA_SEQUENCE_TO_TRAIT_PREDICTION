"""CRyPTIC feasibility probe — is the M. tuberculosis measured-MIC compendium a viable decoder substrate?

The acquisition-fork research (wiki/acquirable... / research_outputs/acquirable-label-sources-2026-06-16.md)
ranked CRyPTIC #1: 12,287 M. tb isolates x 13 drugs, reference UKMYC broth-microdilution MIC + binary
phenotype + per-sample VCF (variant calls vs H37Rv), FREE on the EBI FTP, disjoint by construction (no
TB tuning set exists). TB resistance is POINT-mutation-driven (rpoB/katG/gyrA/embB/pncA) — the decoder's
STRONGEST regime (cf. cipro QRDR). This probe answers, cheaply, WITHOUT building the full TB decoder:

  1. POWERING census (CSV-only, decisive): per drug, count R / S isolates at HIGH phenotype quality ->
     which drugs clear a >=20/class cohort bar (and by how much).
  2. GENOTYPE ACCESS: confirm the per-sample VCF is fetchable from the FTP + parseable + carries the
     expected resistance-region variants (no genome assembly needed — the VCF IS the determinant source).
  3. DETERMINANT-RULE PROOF-OF-CONCEPT (rifampicin): on a balanced HIGH-quality sample, predict R iff the
     VCF carries a variant in the rpoB RRDR (the canonical RIF region) and compare to the measured binary
     phenotype -> a first sens/spec, demonstrating the deterministic rule transfers to TB. This is a PoC
     PROXY for the WHO-catalogue codon-level rule, NOT the final rule.
  4. CLONALITY proxy: site distribution (UNIQUEID encodes collection site) — a coarse read on whether the
     cohort would survive the project's Mash-lineage clonality correction.

Run (needs network for steps 2-3; step 1 is offline once the CSV is local):
  uv run python scripts/cryptic_feasibility_probe.py            # full probe (census + access + RIF PoC)
  uv run python scripts/cryptic_feasibility_probe.py --census-only
  uv run python scripts/cryptic_feasibility_probe.py --rif-sample 20   # N per class for the PoC
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import sys
import urllib.request
from collections import Counter
from datetime import date as _date
from pathlib import Path

REUSE_CSV = Path("data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv")
FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"
VCF_CACHE = Path("data/raw/cryptic/vcf_cache")

DRUGS = {"AMI": "amikacin", "BDQ": "bedaquiline", "CFZ": "clofazimine", "DLM": "delamanid",
         "EMB": "ethambutol", "ETH": "ethionamide", "INH": "isoniazid", "KAN": "kanamycin",
         "LEV": "levofloxacin", "LZD": "linezolid", "MXF": "moxifloxacin", "RIF": "rifampicin",
         "RFB": "rifabutin"}

# rpoB RRDR window in H37Rv (NC_000962.3): rpoB starts 759807 (+strand); codons 426-452 (the canonical
# "81-bp RRDR", E. coli 507-533) span ~761082-761162. Pad slightly. A variant here is the dominant RIF
# determinant (S450L etc.) — a PoC proxy for the WHO catalogue, NOT codon-resolved.
RRDR_LO, RRDR_HI = 761055, 761165


def load_rows(csv_path: Path = REUSE_CSV) -> list[dict]:
    with open(csv_path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def census(rows: list[dict]) -> dict:
    """Per-drug R/S counts at HIGH phenotype quality (the powering census)."""
    out = {}
    for code, name in DRUGS.items():
        ph, q = f"{code}_BINARY_PHENOTYPE", f"{code}_PHENOTYPE_QUALITY"
        r = s = r_hi = s_hi = 0
        for row in rows:
            call = (row.get(ph) or "").strip().upper()
            qual = (row.get(q) or "").strip().upper()
            if call == "R":
                r += 1
                if qual == "HIGH":
                    r_hi += 1
            elif call == "S":
                s += 1
                if qual == "HIGH":
                    s_hi += 1
        powered = "POWERED" if (r_hi >= 20 and s_hi >= 20) else "UNDERPOWERED"
        out[name] = {"code": code, "R_all": r, "S_all": s, "R_high": r_hi, "S_high": s_hi,
                     "verdict": powered}
    return out


def _vcf_url(rel_path: str) -> str:
    return FTP_BASE + rel_path.split("reproducibility/", 1)[1] if "reproducibility/" in rel_path else FTP_BASE + rel_path.lstrip("./")


def fetch_vcf(rel_path: str, timeout: int = 60) -> bytes | None:
    """Fetch + cache a gzipped VCF. Returns decompressed bytes or None on failure."""
    VCF_CACHE.mkdir(parents=True, exist_ok=True)
    cache = VCF_CACHE / (rel_path.replace("/", "_").replace("..", "").strip("_"))
    if cache.exists():
        return cache.read_bytes()
    url = _vcf_url(rel_path)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            raw = r.read()
        data = gzip.decompress(raw)
        cache.write_bytes(data)
        return data
    except Exception as e:  # noqa: BLE001
        print(f"    [vcf fetch fail] {url}: {type(e).__name__}: {e}")
        return None


def rrdr_variant_present(vcf_bytes: bytes) -> bool:
    """True iff the VCF has an ALT record in the rpoB RRDR window (single-chromosome organism -> POS only)."""
    for line in io.StringIO(vcf_bytes.decode("utf-8", "replace")):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        try:
            pos = int(f[1])
        except ValueError:
            continue
        alt = f[4].strip()
        if RRDR_LO <= pos <= RRDR_HI and alt not in (".", "", "<NON_REF>"):
            return True
    return False


def rif_poc(rows: list[dict], n_per_class: int) -> dict:
    """Determinant-rule PoC: rpoB-RRDR variant -> predict R; compare to measured RIF binary phenotype."""
    def pick(call):
        out = []
        for row in rows:
            if (row.get("RIF_BINARY_PHENOTYPE") or "").strip().upper() == call \
               and (row.get("RIF_PHENOTYPE_QUALITY") or "").strip().upper() == "HIGH" \
               and (row.get("VCF") or "").strip():
                out.append(row)
            if len(out) >= n_per_class:
                break
        return out

    sample = pick("R") + pick("S")
    tp = fp = tn = fn = miss = 0
    for row in sample:
        vcf = fetch_vcf(row["VCF"].strip())
        if vcf is None:
            miss += 1
            continue
        pred_r = rrdr_variant_present(vcf)
        true_r = (row.get("RIF_BINARY_PHENOTYPE") or "").strip().upper() == "R"
        if pred_r and true_r:
            tp += 1
        elif pred_r and not true_r:
            fp += 1
        elif (not pred_r) and not true_r:
            tn += 1
        else:
            fn += 1
    n = tp + fp + tn + fn
    return {"n_scored": n, "vcf_fetch_miss": miss, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sens": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "spec": round(tn / (tn + fp), 3) if (tn + fp) else None,
            "rule": f"rpoB RRDR variant in NC_000962.3:{RRDR_LO}-{RRDR_HI} (PoC proxy for WHO catalogue)"}


def clonality_proxy(rows: list[dict]) -> dict:
    """Coarse clonality read: distinct collection sites (UNIQUEID 'site.NN.subj...')."""
    sites = Counter()
    for row in rows:
        uid = row.get("UNIQUEID") or ""
        if uid.startswith("site."):
            sites[uid.split(".")[1]] += 1
    return {"distinct_sites": len(sites), "top_sites": dict(sites.most_common(8))}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--census-only", action="store_true")
    ap.add_argument("--rif-sample", type=int, default=15, help="N per class for the RIF PoC")
    a = ap.parse_args()
    if not REUSE_CSV.exists():
        print(f"missing {REUSE_CSV} — download from "
              "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reuse/")
        return 2
    rows = load_rows()
    print(f"loaded {len(rows)} CRyPTIC isolates\n")

    cen = census(rows)
    print("=== POWERING census (R/S at HIGH phenotype quality) ===")
    print(f"{'drug':<14}{'R_high':>8}{'S_high':>8}{'R_all':>8}{'S_all':>8}  verdict")
    for name, c in sorted(cen.items(), key=lambda kv: -kv[1]['R_high']):
        print(f"{name:<14}{c['R_high']:>8}{c['S_high']:>8}{c['R_all']:>8}{c['S_all']:>8}  {c['verdict']}")
    clon = clonality_proxy(rows)
    print(f"\nclonality proxy: {clon['distinct_sites']} distinct collection sites; "
          f"top: {clon['top_sites']}")

    poc = None
    if not a.census_only:
        print(f"\n=== GENOTYPE ACCESS + RIF determinant-rule PoC (sample {a.rif_sample}/class) ===")
        poc = rif_poc(rows, a.rif_sample)
        print(f"RIF PoC: n={poc['n_scored']} sens={poc['sens']} spec={poc['spec']} "
              f"(TP{poc['tp']}/FP{poc['fp']}/TN{poc['tn']}/FN{poc['fn']}; vcf_miss {poc['vcf_fetch_miss']})")
        print(f"  rule: {poc['rule']}")

    art = {"_schema": "cryptic-feasibility-probe-v1", "date": _date.today().isoformat(),
           "source": "CRyPTIC reuse table 20240917 (EBI FTP)", "n_isolates": len(rows),
           "census": cen, "clonality_proxy": clon, "rif_poc": poc}
    out = Path(f"wiki/cryptic_feasibility_probe_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nartifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
