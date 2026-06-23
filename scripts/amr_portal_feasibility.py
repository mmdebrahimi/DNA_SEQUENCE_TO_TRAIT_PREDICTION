"""AMR Portal (EBI/CABBAGE) feasibility + leakage pass — is a FREE INDEPENDENT number in reach per cell?

The EBI AMR Portal is free, per-isolate, MEASURED-AST, accession-linked (`wiki/ebi_amr_portal_finding_2026-06-23.md`) — but it AGGREGATES NCBI-PD/PATRIC/NARMS, so it OVERLAPS our cohorts and CRyPTIC; it is NOT
auto-independent. This script measures, per (organism, drug), the PROVENANCE-DISJOINT subset: AMR-Portal
isolates whose BioSample / ERS / GCA accession is NOT in CRyPTIC (the TB in-distribution set) and NOT in any
of our existing cohorts (`cohort_manifest`). A disjoint subset with measured R+S is a FREE independent
validation set — the thing the TB gold-set saga (5 author/DUA/circular walls) could not get.

Leakage is ALIAS-AWARE (the `tb_goldset.assert_independent_aliased` discipline): an isolate is leaked iff ANY
of its accessions (BioSample SAMEA / SRA ERS / assembly GCA) is in the union leak set. The phenotype is the
MEASURED `resistance_phenotype` (lab AST) — non-circular by construction.

Pure logic (`binarize_sir` / `iso_leaked` / `summarize`) is unit-tested; the parquet + leakage-set load is live.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

DEFAULT_PARQUET = Path("D:/dna_decode_cache/data files donwload/amr_portal/phenotype.parquet")
DEFAULT_CRYPTIC = Path("D:/dna_decode_cache/data files donwload/CRyPTIC TB MIC compendium/"
                       "CRyPTIC_reuse_table_20240917.csv")
# AMR-Portal drug name -> our cell drug name (for the headline cells we care about).
DRUG_ALIASES = {"rifampin": "rifampicin", "rifampicin": "rifampicin", "isoniazid": "isoniazid",
                "ciprofloxacin": "ciprofloxacin", "tetracycline": "tetracycline", "gentamicin": "gentamicin",
                "meropenem": "meropenem", "ceftriaxone": "ceftriaxone",
                "trimethoprim-sulfamethoxazole": "trimethoprim-sulfamethoxazole"}
# the cells we most want to know about (organism, our-drug)
HEADLINE = {("Mycobacterium tuberculosis", "rifampicin"), ("Mycobacterium tuberculosis", "isoniazid"),
            ("Escherichia coli", "ciprofloxacin"), ("Escherichia coli", "tetracycline"),
            ("Escherichia coli", "ceftriaxone"), ("Escherichia coli", "gentamicin"),
            ("Escherichia coli", "meropenem")}


def binarize_sir(resistance_phenotype: str) -> str | None:
    """Measured SIR -> 'R'/'S'/None. resistant/non-susceptible -> R; susceptible(+SDD) -> S; intermediate/'' -> None."""
    v = (resistance_phenotype or "").strip().lower()
    if v in ("resistant", "non-susceptible"):
        return "R"
    if v in ("susceptible", "susceptible-dose dependent"):
        return "S"
    return None


def iso_leaked(aliases: set[str], leak_set: set[str]) -> bool:
    """True iff ANY of the isolate's accessions is in the union leak set (alias-aware, upper-cased)."""
    return bool({a.upper() for a in aliases if a} & leak_set)


def summarize(rows, leak_set: set[str]) -> dict:
    """rows: iterable of (organism, drug_our, biosample, sra, assembly, sir). Returns per-(org,drug) powering.

    Counts UNIQUE (organism, drug, isolate) — an isolate's leakage is computed once from its accessions."""
    seen: set = set()
    agg: dict = defaultdict(lambda: {"total": 0, "disjoint": 0, "disjoint_R": 0, "disjoint_S": 0,
                                     "leaked": 0, "intermediate_or_blank": 0})
    for organism, drug, biosample, sra, assembly, rp in rows:
        sir = binarize_sir(rp)
        iso_key = biosample or sra or assembly
        dedupe = (organism, drug, iso_key, sir)
        if dedupe in seen:
            continue
        seen.add(dedupe)
        cell = agg[(organism, drug)]
        if sir is None:
            cell["intermediate_or_blank"] += 1
            continue
        cell["total"] += 1
        if iso_leaked({biosample, sra, assembly}, leak_set):
            cell["leaked"] += 1
        else:
            cell["disjoint"] += 1
            cell["disjoint_R" if sir == "R" else "disjoint_S"] += 1
    return dict(agg)


def _load_leak_set(cryptic_csv: Path) -> set[str]:
    from dna_decode.organism_rules import tb_goldset
    from dna_decode.eval import cohort_manifest as cm
    leak = set()
    if cryptic_csv.exists():
        leak |= {a.upper() for a in tb_goldset.cryptic_accessions(cryptic_csv)}
    m = cm.build_manifest()
    leak |= {a.upper() for a in cm.prior_accessions(m, exclude_cohort="__none__")}
    return leak


def _load_rows(parquet: Path):
    import pyarrow.parquet as pq
    cols = ["organism", "antibiotic_name", "BioSample_ID", "SRA_accession", "assembly_ID",
            "resistance_phenotype"]
    pf = pq.ParquetFile(str(parquet))
    for batch in pf.iter_batches(batch_size=200_000, columns=cols):
        d = batch.to_pydict()
        for org, ab, bs, sra, asm, rp in zip(d["organism"], d["antibiotic_name"], d["BioSample_ID"],
                                              d["SRA_accession"], d["assembly_ID"], d["resistance_phenotype"]):
            drug = DRUG_ALIASES.get((ab or "").strip().lower())
            if drug is None:
                continue                       # only the drugs our cells route (keeps the report focused)
            yield (org or "", drug, (bs or "").strip(), (sra or "").strip(), (asm or "").strip(), rp)


def main(argv=None) -> int:
    import argparse, json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parquet", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--cryptic", type=Path, default=DEFAULT_CRYPTIC)
    ap.add_argument("--min-per-class", type=int, default=10)
    a = ap.parse_args(argv)
    if not a.parquet.exists():
        print(f"ERROR: AMR Portal phenotype parquet not found at {a.parquet}", file=sys.stderr)
        return 2
    leak = _load_leak_set(a.cryptic)
    print(f"[amr-portal-feasibility] leak set: {len(leak):,} accessions (CRyPTIC + our cohorts)")
    agg = summarize(_load_rows(a.parquet), leak)

    def powered(c):
        return c["disjoint_R"] >= a.min_per_class and c["disjoint_S"] >= a.min_per_class

    rows_sorted = sorted(agg.items(), key=lambda kv: -kv[1]["disjoint"])
    print(f"\n{'organism':<32} {'drug':<22} {'total':>7} {'disjoint':>8} {'dR':>5} {'dS':>5} {'leaked':>7} powered")
    headline_powered = []
    for (org, drug), c in rows_sorted:
        is_head = (org, drug) in HEADLINE
        p = powered(c)
        if is_head or c["disjoint"] >= 50:
            mark = "***" if is_head else "   "
        else:
            continue
        print(f"{mark}{org:<29} {drug:<22} {c['total']:>7} {c['disjoint']:>8} {c['disjoint_R']:>5} "
              f"{c['disjoint_S']:>5} {c['leaked']:>7} {'YES' if p else 'no'}")
        if p:
            headline_powered.append((org, drug, c))
    print(f"\nPROVENANCE-DISJOINT POWERED cells (>= {a.min_per_class} R and S, free independent number in reach): "
          f"{len(headline_powered)}")
    for org, drug, c in headline_powered:
        print(f"  {org} / {drug}: disjoint {c['disjoint']} ({c['disjoint_R']}R / {c['disjoint_S']}S)")

    out = REPO / "wiki" / "amr_portal_feasibility.json"
    out.write_text(json.dumps({f"{o}|{d}": c for (o, d), c in agg.items()}, indent=2), encoding="utf-8")
    print(f"artifact -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
