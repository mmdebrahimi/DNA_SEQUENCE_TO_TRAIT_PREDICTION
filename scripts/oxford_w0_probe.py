"""W0 empirical probe — pin the Oxford MIC-table schema + crosswalk feasibility BEFORE
building the ingester against it (the "build against real data" gate).

Given the project accession + a MIC-table path, emits a `wiki/oxford_w0_probe_<date>.json`
audit: row cardinality, unique-key counts per candidate ENA field (run/sample_alias/
secondary/native BioSample), duplicate MIC rows per isolate/drug, operator/censoring
distribution per drug, and the MIC-key -> BioSample resolution rate. This output PINS
the ingester's column map (Step 5) + the crosswalk's key candidates (Step 4).

PURE summarize helpers are split from the fetch/read so they unit-test offline.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.external_mic_labels import canonical_drug, parse_mic_value

ENA_CANDIDATE_FIELDS = ("run_accession", "sample_accession", "sample_alias", "secondary_sample_accession")


# --------------------------------------------------------------------------- #
# Pure summarizers (no network / no file IO)
# --------------------------------------------------------------------------- #
def candidate_key_cardinality(records: list[dict], fields=ENA_CANDIDATE_FIELDS) -> dict:
    """Per ENA field: how many distinct non-empty values across the records."""
    return {f: len({r[f] for r in records if r.get(f)}) for f in fields}


def mic_table_summary(rows: list[dict], key_col: str, drug_cols: dict[str, str]) -> dict:
    """Summarize a parsed MIC table.

    `drug_cols` maps file-column -> drug-alias. Returns row count, unique key count,
    duplicate (isolate, drug) row counts, and per-canonical-drug operator/censoring
    distribution.
    """
    keys = [r.get(key_col, "").strip() for r in rows if r.get(key_col, "").strip()]
    seen_pairs: dict[tuple[str, str], int] = {}
    censoring: dict[str, dict[str, int]] = {}
    for r in rows:
        k = r.get(key_col, "").strip()
        if not k:
            continue
        for col, alias in drug_cols.items():
            drug = canonical_drug(alias)
            if drug is None:
                continue
            cell = r.get(col)
            mv = parse_mic_value(cell)
            if mv is None:
                continue
            seen_pairs[(k, drug)] = seen_pairs.get((k, drug), 0) + 1
            d = censoring.setdefault(drug, {"=": 0, ">": 0, ">=": 0, "<": 0, "<=": 0})
            d[mv.operator] = d.get(mv.operator, 0) + 1
    dup_isolate_drug = {f"{k}|{drug}": n for (k, drug), n in seen_pairs.items() if n > 1}
    return {
        "n_rows": len(rows),
        "n_unique_keys": len(set(keys)),
        "n_isolate_drug_pairs": len(seen_pairs),
        "duplicate_isolate_drug": dup_isolate_drug,
        "censoring_by_drug": censoring,
    }


def resolution_summary(mic_keys, records: list[dict], fields=ENA_CANDIDATE_FIELDS) -> dict:
    """How many MIC keys match SOME ENA candidate-field value (resolution feasibility)."""
    value_set = {r[f] for r in records for f in fields if r.get(f)}
    mic_keys = sorted({k for k in mic_keys if k})
    resolved = [k for k in mic_keys if k in value_set]
    n = len(mic_keys)
    return {
        "n_mic_keys": n,
        "n_resolved": len(resolved),
        "resolution_rate": round(len(resolved) / n, 4) if n else 0.0,
        "unresolved_sample": [k for k in mic_keys if k not in value_set][:10],
    }


# --------------------------------------------------------------------------- #
# IO + main
# --------------------------------------------------------------------------- #
def read_table(path: str | Path) -> list[dict]:
    """Read a CSV/TSV into a list of dict rows (delimiter sniffed; tab/comma fallback)."""
    text = Path(path).read_text(encoding="utf-8")
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        delim = dialect.delimiter
    except csv.Error:
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    return list(csv.DictReader(text.splitlines(), delimiter=delim))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", required=True)
    ap.add_argument("--mic-table", required=True)
    ap.add_argument("--key-col", required=True, help="MIC-table column holding the isolate key")
    ap.add_argument("--drug-col", action="append", default=[], metavar="COL=ALIAS",
                    help="map a file column to a drug alias (repeatable), e.g. CIP=ciprofloxacin")
    ap.add_argument("--wiki-dir", default="wiki")
    a = ap.parse_args()
    drug_cols = dict(kv.split("=", 1) for kv in a.drug_col)

    rows = read_table(a.mic_table)
    from dna_decode.eval.biosample_resolver import BioSampleResolver
    resolver = BioSampleResolver()
    records = resolver.read_run_records_for_project(a.project)

    mic_keys = [r.get(a.key_col, "").strip() for r in rows]
    audit = {
        "_schema": "oxford-w0-probe-v1",
        "date": _date.today().isoformat(),
        "project": a.project,
        "mic_table": str(a.mic_table),
        "key_col": a.key_col,
        "drug_cols": drug_cols,
        "ena_candidate_cardinality": candidate_key_cardinality(records),
        "mic_table_summary": mic_table_summary(rows, a.key_col, drug_cols),
        "resolution": resolution_summary(mic_keys, records),
    }
    out = Path(a.wiki_dir) / f"oxford_w0_probe_{_date.today().isoformat()}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(f"W0 probe -> {out}; resolution_rate={audit['resolution']['resolution_rate']} "
          f"({audit['resolution']['n_resolved']}/{audit['resolution']['n_mic_keys']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
