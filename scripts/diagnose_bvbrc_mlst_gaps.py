"""Diagnose why 87% of BV-BRC cipro candidates lack MLST.

3357 cipro candidates after AST+metadata join; only 387 carry MLST. The 2970
MLST-missing strains include most cipro-R (R-ceiling is 49 after the filter).
This script answers: where did the MLST data go?

Layer-by-layer inspection:
  1. Raw BVBRC_genome (1).csv: count rows with non-empty `MLST` column.
  2. After load_bvbrc_genome_metadata parsing: count entries with non-empty mlst.
  3. Per-R/S breakdown: are MLST gaps uniformly distributed or class-skewed?
  4. Sample 10 MLST-missing R strains: what does their raw row look like?
  5. Recovery options: alternate MLST columns? PubMLST cross-reference?
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


AST_CSV = Path("C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv")
METADATA_CSV = Path("C:/Users/Farshad/Downloads/BVBRC_genome (1).csv")


def main() -> int:
    # --- Layer 1: raw metadata CSV inspection -----------------------------------
    print(f"[layer 1] raw metadata CSV inspection: {METADATA_CSV}")
    raw = pd.read_csv(METADATA_CSV, sep=None, engine="python", dtype=str, keep_default_na=False)
    raw.columns = [str(c).strip().lower().replace(" ", "_") for c in raw.columns]
    print(f"  total rows: {len(raw)}")
    print(f"  columns containing 'mlst': {[c for c in raw.columns if 'mlst' in c.lower()]}")
    if "mlst" in raw.columns:
        non_empty_mlst = (raw["mlst"].astype(str).str.strip() != "").sum()
        print(f"  rows with non-empty raw 'mlst' field: {non_empty_mlst} / {len(raw)} ({100*non_empty_mlst/len(raw):.1f}%)")
    else:
        print("  WARNING: no 'mlst' column in raw metadata!")

    # --- Layer 2: parser output ---------------------------------------------------
    print(f"\n[layer 2] after load_bvbrc_genome_metadata parsing")
    from dna_decode.data.bvbrc_genome import load_bvbrc_genome_metadata
    metadata = load_bvbrc_genome_metadata(METADATA_CSV)
    print(f"  metadata dict size: {len(metadata)}")
    with_mlst = sum(1 for m in metadata.values() if m.get("mlst") and str(m["mlst"]).strip() not in ("", "None"))
    print(f"  entries with non-empty parsed mlst: {with_mlst} / {len(metadata)} ({100*with_mlst/len(metadata):.1f}%)")

    # --- Layer 3: AST-join + R/S breakdown ---------------------------------------
    print(f"\n[layer 3] AST + metadata join, per-R/S MLST breakdown")
    from dna_decode.data.ast_data import load_bvbrc_ast
    ast = load_bvbrc_ast(AST_CSV)
    cipro_ast = ast[ast["antibiotic"].str.lower() == "ciprofloxacin"]
    print(f"  cipro AST rows: {len(cipro_ast)}")
    # Join to metadata, mark MLST presence
    joined_rows = []
    for _, row in cipro_ast.iterrows():
        sid = row["strain_id"]
        meta = metadata.get(sid, {})
        mlst = meta.get("mlst", "")
        has_mlst = bool(str(mlst).strip()) and str(mlst).strip() not in ("None",)
        has_assembly = bool(str(meta.get("assembly_accession", "")).strip())
        joined_rows.append({
            "strain_id": sid,
            "label": row["susceptibility_label"],
            "has_mlst": has_mlst,
            "has_assembly": has_assembly,
            "mlst": mlst,
        })
    j = pd.DataFrame(joined_rows)
    print(f"  cipro candidates with metadata entry: {(j['has_assembly']).sum()} of {len(j)}")
    print(f"  cipro candidates with MLST: {(j['has_mlst']).sum()} of {len(j)} ({100*(j['has_mlst']).sum()/len(j):.1f}%)")
    print(f"  MLST x label crosstab:")
    print(j.groupby(["label", "has_mlst"]).size().unstack(fill_value=0))

    # --- Layer 4: sample MLST-missing R rows -------------------------------------
    print(f"\n[layer 4] sample 10 MLST-missing cipro-R strains")
    missing_r = j[(j["has_mlst"] == False) & (j["label"] == 1)]
    print(f"  total MLST-missing cipro-R: {len(missing_r)}")
    print("  sample strain IDs:", list(missing_r["strain_id"].head(10)))
    # Look up raw metadata for these
    for sid in list(missing_r["strain_id"].head(5)):
        meta = metadata.get(sid, "(not in metadata dict)")
        print(f"    {sid}: {meta}")

    # --- Layer 5: recovery options ------------------------------------------------
    print(f"\n[layer 5] recovery options analysis")
    # Check if raw CSV has the strain IDs but parser dropped them
    if "genome_id" in raw.columns and "mlst" in raw.columns:
        raw_ids = set(raw["genome_id"].astype(str))
        missing_r_ids = set(missing_r["strain_id"].astype(str))
        in_raw = missing_r_ids & raw_ids
        print(f"  MLST-missing cipro-R also present in RAW metadata CSV: {len(in_raw)} / {len(missing_r_ids)}")
        if in_raw:
            sample = list(in_raw)[:5]
            print(f"  raw rows for sample IDs (first 5):")
            for sid in sample:
                rows = raw[raw["genome_id"] == sid]
                if len(rows):
                    r0 = rows.iloc[0]
                    mlst_raw = r0.get("mlst", "(no mlst col)")
                    print(f"    {sid}: raw mlst = {mlst_raw!r}")

    print("\n[diag] done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
