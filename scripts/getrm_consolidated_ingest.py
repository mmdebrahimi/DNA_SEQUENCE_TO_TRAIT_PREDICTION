"""Ingest the full GeT-RM Consolidated PGx/HLA table (CDC, public-domain) -> a clean long-format truth TSV.

The GeT-RM Consolidated table (363 samples x 34 genes/loci) + the 137-PGx-sample table (137 samples with
per-gene star-allele diplotypes + ENA sequencing-data URLs) are the free, authoritative CDC/Coriell PGx
truth set. We currently score a SLICE (ursaPGx common.tsv + a few consensus.tsv). This parses the official
xlsx -> `data/pgx_getrm/getrm_consolidated_truth.tsv` (long: coriell_id, gene, diplotype, run_accession),
the reusable substrate that (a) EXPANDS N for our scored genes, (b) lets us SCORE our built-but-unvalidated
cells (UGT1A1/SLCO1B1/CYP4F2), (c) FIXES CYP3A5's underpowered n=8, and (d) maps free-truth roadmap genes.

The source xlsx is a user browser-download (the CDC page bot-walls automated fetch); pass its path. The
DERIVED TSV is public-domain (US-gov work) and committable, mirroring the existing getrm_*_consensus.tsv.

    uv run python scripts/getrm_consolidated_ingest.py \
        --consolidated "C:/.../Consolidated_PGx-HLA_table_1-22-25-V4.xlsx" \
        --panel137 "C:/.../137-PGx-sample-table-8-1-25-V3.xlsx"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
_META = {"URL", "RUN_ACCESSION", "FASTQ_FTP", "SUBMITTED_FTP", "SRA_FTP"}
# our decoder genes (pgx catalogs + single-variant cells) + which are already SCORED
OUR_GENES = {"CYP2C19", "CYP2C9", "CYP2C8", "CYP2D6", "CYP2B6", "CYP3A5", "TPMT",
             "NUDT15", "UGT1A1", "DPYD", "VKORC1", "CYP4F2", "ABCG2", "SLCO1B1"}
SCORED = {"CYP2C19", "CYP2C9", "CYP2C8", "CYP2D6", "CYP2B6", "CYP3A5", "TPMT"}


def _is_gene_col(name) -> bool:
    return (isinstance(name, str) and name and "Reference" not in name and "ID" not in name
            and "GeT-RM" not in name and name.upper() == name and name.upper() not in _META
            and any(ch.isalpha() for ch in name))


def _looks_like_diplotype(v: str) -> bool:
    return "*" in v or ("/" in v and any(ch.isalnum() for ch in v))


def parse_xlsx(path: Path, data_sheet: str) -> tuple[list[dict], dict[str, int]]:
    """Return (rows, per-gene sample count). rows = {coriell_id, gene, diplotype, run_accession}."""
    import pandas as pd
    raw = pd.ExcelFile(path).parse(data_sheet, header=None)
    hdr = raw.iloc[1].tolist()
    body = raw.iloc[2:].reset_index(drop=True)
    # locate the Coriell-ID + run_accession columns
    def col_idx(pred):
        for j, h in enumerate(hdr):
            if isinstance(h, str) and pred(h):
                return j
        return None
    cid = col_idx(lambda h: "Coriell" in h)
    acc = col_idx(lambda h: h.strip().lower() == "run_accession")
    gene_cols = [(j, h) for j, h in enumerate(hdr) if _is_gene_col(h)]
    rows: list[dict] = []
    cov: dict[str, int] = {}
    for _, r in body.iterrows():
        sample = str(r.iloc[cid]).strip() if cid is not None else ""
        if not sample or sample.lower() == "nan":
            continue
        run = str(r.iloc[acc]).strip() if acc is not None else ""
        run = "" if run.lower() == "nan" else run
        for j, gene in gene_cols:
            v = str(r.iloc[j]).strip()
            if v and v.lower() != "nan" and _looks_like_diplotype(v):
                rows.append({"coriell_id": sample, "gene": gene, "diplotype": v, "run_accession": run})
                cov[gene] = cov.get(gene, 0) + 1
    return rows, cov


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="getrm_consolidated_ingest")
    ap.add_argument("--consolidated", help="Consolidated_PGx-HLA_table xlsx (363 samples)")
    ap.add_argument("--panel137", help="137-PGx-sample-table xlsx (137 samples + ENA URLs)")
    ap.add_argument("--out-tsv", default="data/pgx_getrm/getrm_consolidated_truth.tsv")
    ap.add_argument("--out-json", default="wiki/pgx_getrm_consolidated_coverage_2026-07-28.json")
    args = ap.parse_args(argv)
    if not args.consolidated and not args.panel137:
        print("error: give --consolidated and/or --panel137 xlsx path(s)", file=sys.stderr)
        return 2

    all_rows: dict[tuple, dict] = {}   # (coriell_id, gene) -> row (panel137 wins for run_accession)
    src_cov: dict[str, dict] = {}
    for label, path, sheet in (("consolidated", args.consolidated, "PGx HLA Genotypes"),
                               ("panel137", args.panel137, "137 PGx Panel")):
        if not path:
            continue
        rows, cov = parse_xlsx(Path(path), sheet)
        src_cov[label] = cov
        for r in rows:
            key = (r["coriell_id"], r["gene"])
            if key not in all_rows or (r["run_accession"] and not all_rows[key]["run_accession"]):
                all_rows[key] = r

    rows = sorted(all_rows.values(), key=lambda r: (r["gene"], r["coriell_id"]))
    out_tsv = REPO / args.out_tsv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as fh:
        fh.write("coriell_id\tgene\tdiplotype\trun_accession\n")
        for r in rows:
            fh.write(f"{r['coriell_id']}\t{r['gene']}\t{r['diplotype']}\t{r['run_accession']}\n")

    # coverage (merged, best per gene) + cross-reference vs our cells
    cov: dict[str, int] = {}
    for r in rows:
        cov[r["gene"]] = cov.get(r["gene"], 0) + 1
    ours = {g: cov[g] for g in OUR_GENES if g in cov}
    unscored_cells = {g: n for g, n in ours.items() if g not in SCORED}
    new_free = {g: n for g, n in cov.items() if g not in OUR_GENES and "HLA" not in g and n >= 20}
    hla = sorted(g for g in cov if "HLA" in g)
    summary = {
        "source": "GeT-RM Consolidated PGx/HLA + 137-PGx-sample tables (CDC, public-domain)",
        "n_rows": len(rows), "n_samples": len({r["coriell_id"] for r in rows}),
        "n_genes": len(cov), "n_with_ena_url": sum(1 for r in rows if r["run_accession"]),
        "our_genes_coverage": dict(sorted(ours.items(), key=lambda x: -x[1])),
        "built_but_unscored_cells_now_scoreable": dict(sorted(unscored_cells.items(), key=lambda x: -x[1])),
        "new_free_truth_pharmacogenes_no_decoder_yet": dict(sorted(new_free.items(), key=lambda x: -x[1])),
        "hla_loci": hla, "per_source_gene_counts": src_cov,
    }
    out_json = REPO / args.out_json
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[wrote {len(rows)} truth rows / {summary['n_samples']} samples / {summary['n_genes']} genes -> {out_tsv}]")
    print(json.dumps({k: summary[k] for k in
                      ("built_but_unscored_cells_now_scoreable",
                       "our_genes_coverage", "new_free_truth_pharmacogenes_no_decoder_yet")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
