#!/usr/bin/env python
"""Multi-gene AF-corroboration — validation numbers for the new single-SNP PGx cells (v0.1).

Extends the DPYD AF-corroboration (scripts/dpyd_af_corroboration.py) to NUDT15 / UGT1A1 / CYP4F2 / ABCG2.
For each cell's actionable variant, confirm the ALT-allele frequency in the population where it is most
informative matches the CPIC/gnomAD-expected band. Same honest tier as DPYD: KNOWLEDGE_BASELINE, NOT an
independent per-sample concordance number (GeT-RM truth for these genes is a paper-supplement = external
wall). This is the machine-fetchable validation ceiling for the new cells, produced without manual curation.

ALT-allele frequencies fetched from Ensembl 1000G phase3 (the ALT = functional allele, selected explicitly
per population) 2026-07-07; `--live` re-fetches.
"""
from __future__ import annotations

import argparse
import datetime
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# gene, allele, rsid, ALT (functional), population, observed ALT AF [Ensembl 1000G, 2026-07-07],
# expected band (CPIC allele-functionality + gnomAD literature).
PGX_AF = [
    {"gene": "NUDT15", "allele": "*3", "rsid": "rs116855232", "alt": "T", "pop": "EAS",
     "af": 0.0952, "lo": 0.05, "hi": 0.15, "note": "no-function; EAS ~10% (thiopurine toxicity)"},
    {"gene": "UGT1A1", "allele": "*80", "rsid": "rs887829", "alt": "T", "pop": "EUR",
     "af": 0.2982, "lo": 0.22, "hi": 0.42, "note": "*28 LD-tag; EUR ~30% == *28 freq"},
    {"gene": "UGT1A1", "allele": "*6", "rsid": "rs4148323", "alt": "A", "pop": "EAS",
     "af": 0.1379, "lo": 0.08, "hi": 0.20, "note": "decreased; EAS-common ~14%"},
    {"gene": "CYP4F2", "allele": "*3", "rsid": "rs2108622", "alt": "T", "pop": "EUR",
     "af": 0.2903, "lo": 0.20, "hi": 0.40, "note": "reduced; EUR ~29% (warfarin dose-up)"},
    {"gene": "ABCG2", "allele": "141K", "rsid": "rs2231142", "alt": "T", "pop": "EAS",
     "af": 0.2907, "lo": 0.20, "hi": 0.40, "note": "poor-function; EAS ~29% (rosuvastatin)"},
]


def classify_af(af: float | None, lo: float, hi: float) -> str:
    """PURE: is the observed ALT-allele frequency within the CPIC/gnomAD-expected band?"""
    if af is None:
        return "NO_DATA"
    return "IN_BAND" if lo <= af <= hi else "OUT_OF_BAND"


def _fetch_alt_af(rsid: str, alt: str, pop: str) -> float | None:
    import urllib.request
    url = f"https://rest.ensembl.org/variation/homo_sapiens/{rsid}?content-type=application/json;pops=1"
    d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=30))
    hits = [p for p in d.get("populations", [])
            if f":{pop}" in p.get("population", "") and "phase_3" in p.get("population", "") and p.get("allele") == alt]
    return hits[0]["frequency"] if hits else None


def build_report(rows: list[dict]) -> dict:
    scored = [{**r, "verdict": classify_af(r["af"], r["lo"], r["hi"])} for r in rows]
    n_in = sum(1 for s in scored if s["verdict"] == "IN_BAND")
    return {
        "schema": "pgx-af-corroboration-v0", "analysis_date": datetime.date.today().isoformat(),
        "genes": sorted({r["gene"] for r in scored}),
        "method": ("1000G ALT-allele-frequency corroboration of each new cell's actionable variant vs the "
                   "CPIC/gnomAD-expected band, in the most-informative population; AFs from Ensembl 1000G phase3"),
        "n_variants": len(scored), "n_in_band": n_in,
        "verdict": "AF_CORROBORATED" if n_in == len(scored) else "AF_PARTIAL",
        "getrm_concordance_status": ("EXTERNAL_WALL for all: GeT-RM truth for NUDT15/UGT1A1/CYP4F2/ABCG2 is "
                                     "paper-supplement (not in the CYP-only ursaPGx benchmark); manual curation"),
        "honesty_tier": ("KNOWLEDGE_BASELINE AF-corroboration (like DPYD/VKORC1/SLCO1B1) — confirms each new "
                         "cell's actionable variant is real + correctly-positioned at CPIC-expected population "
                         "frequency; NOT an independent per-sample concordance number. NOT clinical."),
        "variants": scored,
    }


def render_md(rep: dict) -> str:
    L = [f"# Multi-gene PGx AF-corroboration ({rep['analysis_date']})", "",
         f"_{rep['honesty_tier']}_", "",
         f"**Verdict: {rep['verdict']}** — {rep['n_in_band']}/{rep['n_variants']} actionable variants in band "
         f"across {', '.join(rep['genes'])}.", "",
         f"- **GeT-RM concordance:** {rep['getrm_concordance_status']}", "",
         "| gene | allele | rsid | ALT | pop | ALT AF | expected band | verdict |",
         "|---|---|---|---|---|---|---|---|"]
    for v in rep["variants"]:
        L.append(f"| {v['gene']} | {v['allele']} | {v['rsid']} | {v['alt']} | {v['pop']} | {v['af']} | "
                 f"{v['lo']}-{v['hi']} | {v['verdict']} |")
    L += ["", "_AFs: Ensembl 1000G phase3 (ALT = functional allele, per population), 2026-07-07. "
          "Expected bands: CPIC allele-functionality + gnomAD. NOT a clinical tool._", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Multi-gene PGx AF-corroboration (NUDT15/UGT1A1/CYP4F2/ABCG2).")
    ap.add_argument("--live", action="store_true", help="re-fetch ALT AFs from Ensembl (else use committed)")
    ap.add_argument("--out-md", type=Path, default=REPO / "wiki" / "pgx_af_corroboration_2026-07-07.md")
    ap.add_argument("--out-json", type=Path, default=REPO / "wiki" / "pgx_af_corroboration_2026-07-07.json")
    args = ap.parse_args(argv)
    rows = [dict(r) for r in PGX_AF]
    if args.live:
        for r in rows:
            r["af"] = _fetch_alt_af(r["rsid"], r["alt"], r["pop"])
    rep = build_report(rows)
    args.out_md.write_text(render_md(rep), encoding="utf-8")
    args.out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    print(f"Multi-gene AF-corroboration: {rep['verdict']} ({rep['n_in_band']}/{rep['n_variants']} in band)")
    print(f"[-> {args.out_md}]")
    return 0 if rep["verdict"] == "AF_CORROBORATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
