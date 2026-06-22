"""Ingest a candidate INDEPENDENT TB cohort -> a leakage-checked gold-set manifest (deliverable-b helper).

This is the bridge between "you obtained a non-CRyPTIC TB source" (per wiki/tb_goldset_howto_2026-06-22.md)
and the existing scoring arm (`organism_rules/tb_goldset` + `scripts/score_tb_independent_goldset`). It:
  1. reads a simple candidate TSV you assemble from the source (one row per isolate),
  2. EXCLUDES any isolate present in CRyPTIC (the independence gate — `tb_goldset.assert_independent`
     against `tb_goldset.cryptic_accessions`), writing a leaked-report so the drop is auditable, and
  3. writes a per-drug gold-set manifest (clean isolates only) in the JSON shape the scorer consumes.

Candidate TSV columns (tab-separated, header required):
  strain_id        a stable id for the isolate (your choice)
  masked_vcf       path to the isolate's masked VCF vs H37Rv NC_000962.3 (determinant CALLS)
  regeno_vcf       path to the regeno VCF (callability); leave blank if none
  rif_label        measured RIF DST: R / S / (blank = no RIF label)
  inh_label        measured INH DST: R / S / (blank = no INH label)
  # accession aliases for the CRyPTIC leakage check — supply AS MANY as you have (ALIAS-AWARE):
  ena_accession            (legacy single token) AND/OR
  run_accession            ERR/SRR/DRR
  sample_accession         ERS/SRS/DRS
  biosample_accession      SAMEA/SAMN/SAMD
A candidate is dropped if ANY alias is in CRyPTIC. India PRJNA1155695 is SRA-style (SRR/SAMN) while CRyPTIC
is ENA-style (ERR/SAMEA) — supplying run+sample+biosample maximizes the leakage catch.

The scored INDEPENDENT number is still gated on you SUPPLYING this TSV + the VCFs (the external/data wall);
this script makes the leakage + manifest step a single, auditable, executor-run command.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.organism_rules import tb_goldset  # noqa: E402

DEFAULT_REUSE = Path("D:/dna_decode_cache/data files donwload/CRyPTIC TB MIC compendium/"
                     "CRyPTIC_reuse_table_20240917.csv")
DEFAULT_OUT_DIR = Path("data/raw/tb_goldset")


def read_candidates(tsv_path: Path) -> list[dict]:
    with open(tsv_path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def manifest_rows(candidates: list[dict], drug_code: str) -> list[dict]:
    """Candidate rows -> manifest dicts for one drug (only rows with a clean R/S label for that drug)."""
    label_col = {"RIF": "rif_label", "INH": "inh_label"}[drug_code]
    out = []
    for c in candidates:
        lab = (c.get(label_col) or "").strip().upper()
        if lab not in ("R", "S"):
            continue
        out.append({"strain_id": str(c.get("strain_id", "")).strip(),
                    "masked_vcf": str(c.get("masked_vcf", "")).strip(),
                    "regeno_vcf": (str(c.get("regeno_vcf")).strip() or None) if c.get("regeno_vcf") else None,
                    "label": lab})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", type=Path, required=True, help="candidate TSV (see module docstring)")
    ap.add_argument("--reuse-csv", type=Path, default=DEFAULT_REUSE, help="CRyPTIC reuse table (leakage set)")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    a = ap.parse_args(argv)
    if not a.candidates.exists():
        print(f"ERROR: candidate TSV not found: {a.candidates}", file=sys.stderr)
        return 2
    if not a.reuse_csv.exists():
        print(f"ERROR: CRyPTIC reuse table not found at {a.reuse_csv} (needed for the leakage check)",
              file=sys.stderr)
        return 2

    cand = read_candidates(a.candidates)
    cryptic = tb_goldset.cryptic_accessions(a.reuse_csv)
    # ALIAS-AWARE leakage: gather every accession alias each candidate carries (any subset of these cols).
    alias_cols = ("ena_accession", "run_accession", "sample_accession", "biosample_accession")
    aliases = {str(c.get("strain_id", "")).strip():
               [c.get(col, "") for col in alias_cols if c.get(col)] for c in cand}
    rep = tb_goldset.assert_independent_aliased(aliases, cryptic)
    clean_ids = set(rep.clean)
    clean_cand = [c for c in cand if str(c.get("strain_id", "")).strip() in clean_ids]
    print(f"[tb-goldset] {rep.n_checked} candidates | independent {len(rep.clean)} | "
          f"CRyPTIC-leaked (dropped) {len(rep.leaked)} (alias-aware: run+sample+biosample)")

    a.out_dir.mkdir(parents=True, exist_ok=True)
    if rep.leaked:
        (a.out_dir / "leaked_excluded.txt").write_text("\n".join(rep.leaked) + "\n", encoding="utf-8")
        print(f"  leaked accessions -> {a.out_dir / 'leaked_excluded.txt'}")
    for code in ("RIF", "INH"):
        rows = manifest_rows(clean_cand, code)
        out = a.out_dir / f"goldset_{code.lower()}.json"
        out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        print(f"  {code}: {len(rows)} independent isolates with measured {code} DST -> {out}")
    print("Next: score with  uv run python -m scripts.score_tb_independent_goldset --drug rifampicin")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
