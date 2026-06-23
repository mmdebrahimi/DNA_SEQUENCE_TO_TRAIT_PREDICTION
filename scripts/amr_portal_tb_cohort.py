"""Prepare the AMR Portal provenance-disjoint M. tuberculosis cohort for the INDEPENDENT TB number (rec #1).

The TB cell's rule is the WHO catalogue on RAW per-isolate VCF (`organism_rules/tb_amr`), so — unlike the
bacterial cells — the AMR Portal's AMRFinder genotype table is NOT a faithful input (it is AMRFinder's
narrower, catalogue-derived TB POINT calls; scoring our WHO rule on it tests AMRFinder, not the catalogue).
The independent TB number therefore needs RAW variants: fetch each isolate's assembly → call variants vs
H37Rv NC_000962.3 → masked VCF → `tb_amr.score_drug`. That fetch+call is the DOCKER-COMPUTE-gated step.

This script does the FREE, code-side part: emit the ready cohort — every provenance-disjoint (not in CRyPTIC,
not in our cohorts) M. tuberculosis isolate with a MEASURED rifampicin AND/OR isoniazid label, joined to its
assembly accession — as a TSV the variant-call+score step consumes. It pivots phenotype to per-isolate
{rif,inh} and stamps the leakage verdict. The actual independent number is then one Docker-host run away.

Reuses the verified leakage + SIR helpers from `amr_portal_feasibility`. Pure pivot logic is unit-testable.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.amr_portal_feasibility import (  # noqa: E402
    DEFAULT_CRYPTIC, DEFAULT_PARQUET, binarize_sir, iso_leaked, _load_leak_set,
)

TB = "Mycobacterium tuberculosis"
TB_DRUGS = {"rifampin": "rif", "rifampicin": "rif", "isoniazid": "inh"}
OUT_TSV = REPO / "data" / "raw" / "tb_goldset" / "amr_portal_tb_disjoint_cohort.tsv"


def pivot_tb(rows) -> dict:
    """rows: (biosample, sra, assembly, amr_drug, sir) -> {isolate_key: {rif,inh,biosample,sra,assembly}}.

    Conflicting repeat measurements for the same isolate+drug -> blanked (discordant); first clean kept."""
    iso: dict = defaultdict(lambda: {"rif": "", "inh": "", "biosample": "", "sra": "", "assembly": ""})
    for bs, sra, asm, amr_drug, sir in rows:
        col = TB_DRUGS.get((amr_drug or "").strip().lower())
        if col is None:
            continue
        key = bs or sra or asm
        if not key:
            continue
        rec = iso[key]
        rec["biosample"], rec["sra"], rec["assembly"] = bs, sra, asm
        s = binarize_sir(sir)
        if s in ("R", "S"):
            prior = rec[col]
            if prior and prior != s:
                rec[col] = ""                      # discordant repeat -> blank
            elif not prior:
                rec[col] = s
    return dict(iso)


def main(argv=None) -> int:
    import argparse, json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    ap.add_argument("--out", type=Path, default=OUT_TSV)
    a = ap.parse_args(argv)
    if not a.pheno.exists():
        print(f"ERROR: AMR Portal phenotype parquet not found: {a.pheno}", file=sys.stderr)
        return 2
    import pyarrow.parquet as pq
    leak = _load_leak_set(a.cryptic)
    cols = ["organism", "antibiotic_name", "BioSample_ID", "SRA_accession", "assembly_ID", "resistance_phenotype"]
    pf = pq.ParquetFile(str(a.pheno))
    raw = []
    for batch in pf.iter_batches(batch_size=200_000, columns=cols):
        d = batch.to_pydict()
        for org, ab, bs, sra, asm, rp in zip(d["organism"], d["antibiotic_name"], d["BioSample_ID"],
                                              d["SRA_accession"], d["assembly_ID"], d["resistance_phenotype"]):
            if org == TB and (ab or "").strip().lower() in TB_DRUGS:
                raw.append(((bs or "").strip(), (sra or "").strip(), (asm or "").strip(), ab, rp))
    iso = pivot_tb(raw)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    n_total = n_disjoint = n_disjoint_with_asm = 0
    stats = {"rif_RS": 0, "inh_RS": 0, "disjoint_rif_RS": 0, "disjoint_inh_RS": 0}
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["strain_id", "biosample", "sra", "assembly", "rif_label", "inh_label", "leaked"])
        for key, r in iso.items():
            if r["rif"] not in ("R", "S") and r["inh"] not in ("R", "S"):
                continue
            n_total += 1
            stats["rif_RS"] += r["rif"] in ("R", "S")
            stats["inh_RS"] += r["inh"] in ("R", "S")
            leaked = iso_leaked({r["biosample"], r["sra"], r["assembly"]}, leak)
            if not leaked:
                n_disjoint += 1
                stats["disjoint_rif_RS"] += r["rif"] in ("R", "S")
                stats["disjoint_inh_RS"] += r["inh"] in ("R", "S")
                if r["assembly"]:
                    n_disjoint_with_asm += 1
            w.writerow([key, r["biosample"], r["sra"], r["assembly"], r["rif"], r["inh"],
                        "1" if leaked else "0"])
    print(f"[amr-portal-tb-cohort] {n_total} TB isolates with a measured RIF/INH label -> {a.out}")
    print(f"  provenance-DISJOINT: {n_disjoint} ({n_disjoint_with_asm} with an assembly accession to fetch)")
    print(f"  disjoint measured: RIF {stats['disjoint_rif_RS']} | INH {stats['disjoint_inh_RS']}")
    print("  NEXT (Docker-compute-gated): fetch the assembly_ID genomes -> call variants vs H37Rv "
          "NC_000962.3 -> masked VCF -> tb_amr.score_drug. The cohort above is the ready input.")
    (REPO / "wiki" / "amr_portal_tb_cohort_stats.json").write_text(
        json.dumps({"n_total": n_total, "n_disjoint": n_disjoint,
                    "n_disjoint_with_assembly": n_disjoint_with_asm, **stats}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
