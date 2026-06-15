"""Config-driven MIC-table ingester.

Turns a raw cohort MIC table (CSV/TSV, schema pinned by the W0 probe) into the
per-isolate structure the label builder consumes. NO hardcoded column names — the
caller supplies the isolate key column + a {file-column -> drug-alias} map (+ optional
categorical-call columns). Drug aliases are normalized via `canonical_drug` (non-pilot
drugs skipped); MIC cells are parsed into `MicValue` (operator preserved, Step 2);
columns mapped to nothing are surfaced in an `unmapped_columns` report (loud, not silent).
"""
from __future__ import annotations

import csv
from pathlib import Path

from dna_decode.data.external_mic_labels import MicValue, canonical_drug, parse_mic_value


def _read_table(path: str | Path) -> list[dict]:
    text = Path(path).read_text(encoding="utf-8")
    sample = text[:4096]
    try:
        delim = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    return list(csv.DictReader(text.splitlines(), delimiter=delim))


def ingest_mic_table(path: str | Path, *, key_col: str, drug_cols: dict[str, str],
                     call_cols: dict[str, str] | None = None) -> dict:
    """Ingest a MIC table -> {"data": {native_key: {drug: {"mics":[MicValue], "calls":set}}},
    "unmapped_columns": [...]}.

    `drug_cols` / `call_cols` map file-columns to drug aliases. Raises ValueError if
    `key_col` is absent from the table header.
    """
    rows = _read_table(path)
    call_cols = call_cols or {}
    if rows and key_col not in rows[0]:
        raise ValueError(f"key_col {key_col!r} not in MIC-table columns {sorted(rows[0].keys())}")
    header = set(rows[0].keys()) if rows else set()
    mapped = {key_col, *drug_cols.keys(), *call_cols.keys()}
    unmapped = sorted(header - mapped)

    data: dict[str, dict[str, dict]] = {}
    for r in rows:
        k = (r.get(key_col) or "").strip()
        if not k:
            continue
        per = data.setdefault(k, {})
        for col, alias in drug_cols.items():
            drug = canonical_drug(alias)
            if drug is None:
                continue
            slot = per.setdefault(drug, {"mics": [], "calls": set()})
            mv = parse_mic_value(r.get(col))
            if mv is not None:
                slot["mics"].append(mv)
        for col, alias in call_cols.items():
            drug = canonical_drug(alias)
            if drug is None:
                continue
            call = (r.get(col) or "").strip().upper()
            if call:
                per.setdefault(drug, {"mics": [], "calls": set()})["calls"].add(call)
    return {"data": data, "unmapped_columns": unmapped}
