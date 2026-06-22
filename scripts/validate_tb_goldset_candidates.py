"""Pre-flight validator for a TB gold-set candidate TSV — run BEFORE expensive download/variant-calling.

Checks the candidate TSV (the input to `scripts/build_tb_goldset_manifest.py`) is well-formed: required
columns, unique strain_id, >=1 accession alias, clean R/S labels, >=1 usable measured label, and the
RIF/INH class balance (the 10/10/10 MVP buckets). VCF-existence is OPT-IN (`--require-vcf`) so you can
validate labels + accessions FIRST, then re-run with --require-vcf after variant-calling.

Stdlib only (no pandas) — pure `validate()` is unit-tested; `main` is the CLI wrapper. Exit 0 = OK.
Merged from the ChatGPT shortlist's validator suggestion (2026-06-22).
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

REQUIRED = ("strain_id", "masked_vcf", "regeno_vcf", "rif_label", "inh_label")
ALIAS_COLS = ("ena_accession", "run_accession", "sample_accession", "biosample_accession")
VALID_LABELS = {"R", "S", "", "NA"}


def validate(rows: list[dict], header: list[str], *, require_vcf: bool = False) -> tuple[bool, list[str], dict]:
    """Return (ok, errors, summary). summary has the (rif,inh) bucket counts + usable count."""
    errors: list[str] = []
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        errors.append(f"missing required columns: {missing}")
    if not any(c in header for c in ALIAS_COLS):
        errors.append(f"no accession-alias column present (need >=1 of {list(ALIAS_COLS)})")
    if errors:   # structural failure -> can't check rows meaningfully
        return False, errors, {}

    sids = [(r.get("strain_id") or "").strip() for r in rows]
    dup_sid = [s for s, n in Counter(sids).items() if s and n > 1]
    if dup_sid:
        errors.append(f"duplicate strain_id: {dup_sid[:10]}")

    for col in ("rif_label", "inh_label"):
        bad = sorted({(r.get(col) or "").strip().upper() for r in rows} - VALID_LABELS)
        if bad:
            errors.append(f"bad labels in {col}: {bad} (only R/S/blank/NA allowed)")

    def _lab(r, c):
        return (r.get(c) or "").strip().upper()
    usable = [r for r in rows if _lab(r, "rif_label") in ("R", "S") or _lab(r, "inh_label") in ("R", "S")]
    if not usable:
        errors.append("no usable RIF/INH measured labels (every row blank/NA)")

    if require_vcf:
        missing_vcf = [(r.get("masked_vcf") or "").strip() for r in rows
                       if (r.get("masked_vcf") or "").strip() and not Path((r.get("masked_vcf")).strip()).exists()]
        if missing_vcf:
            errors.append(f"{len(missing_vcf)} masked_vcf paths do not exist (run with --require-vcf only "
                          f"after variant-calling)")

    buckets = Counter((_lab(r, "rif_label") or "-", _lab(r, "inh_label") or "-") for r in rows)
    summary = {"n_rows": len(rows), "n_usable": len(usable),
               "rif_inh_buckets": {f"RIF={k[0]}/INH={k[1]}": v for k, v in sorted(buckets.items())}}
    return (not errors), errors, summary


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("candidates", type=Path)
    ap.add_argument("--require-vcf", action="store_true", help="also assert each masked_vcf path exists")
    a = ap.parse_args(argv)
    if not a.candidates.exists():
        print(f"ERROR: candidate TSV not found: {a.candidates}", file=sys.stderr)
        return 2
    with open(a.candidates, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows = list(rd)
        header = rd.fieldnames or []
    ok, errors, summary = validate(rows, header, require_vcf=a.require_vcf)
    if not ok:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(f"OK: {summary['n_rows']} rows | usable measured-label rows: {summary['n_usable']}")
    print("RIF/INH balance buckets (target ~10 each of SS / SR / RR):")
    for k, v in summary["rif_inh_buckets"].items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
