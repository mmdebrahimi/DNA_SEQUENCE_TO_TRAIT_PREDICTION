"""Label-emission driver: MIC table -> BioSample-keyed labels + the cohort manifest.

Chains ingest (external_mic_ingest) -> crosswalk resolve (external_crosswalk) ->
re-key native MIC keys to BioSamples -> per-drug `build_drug_labels` -> writes
`data/raw/oxford_extval_<drug>/selected_{strict,relaxed}.tsv` + `buckets_<drug>.json`
(BioSample-keyed) AND the single `cohort_manifest_external_<run_id>.json` — the exact
scored-cohort definition consumed by preflight (exact-set), scorer (drift guard), and
roll-up. ABORTS on any crosswalk HARD CONFLICT (no silent collapse).

NOTE: `cohort_manifest_external_<run_id>.json` is a per-run ARTIFACT, unrelated to the
FROZEN `dna_decode/eval/cohort_manifest.py` leakage-registry MODULE.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.external_cohort_genomes import is_biosample_key
from dna_decode.data.external_crosswalk import resolve_keys, write_crosswalk
from dna_decode.data.external_mic_ingest import ingest_mic_table
from dna_decode.data.external_mic_labels import (
    RELAXED_EXTRA,
    STRICT_TIERS,
    build_drug_labels,
    canonical_drug,
    tier_for_isolate,
    write_labels,
)


class CrosswalkConflictError(RuntimeError):
    """Raised when the crosswalk reports HARD CONFLICTS touching the cohort keys."""


def rekey_to_biosample(data: dict, resolved: dict[str, str]) -> tuple[dict, list[str]]:
    """Re-key {native_key: {drug: {...}}} to {biosample: {drug: {mics, calls}}}.

    Native keys absent from `resolved` are dropped (returned as `dropped`, not silent).
    Multiple native keys resolving to the SAME BioSample are merged (mics concatenated,
    calls unioned).
    """
    bs_data: dict[str, dict[str, dict]] = {}
    dropped: list[str] = []
    for native, per_drug in data.items():
        bs = resolved.get(native)
        if bs is None:
            dropped.append(native)
            continue
        dst = bs_data.setdefault(bs, {})
        for drug, slot in per_drug.items():
            d = dst.setdefault(drug, {"mics": [], "calls": set()})
            d["mics"].extend(slot["mics"])
            d["calls"] |= slot["calls"]
    return bs_data, sorted(dropped)


def manifest_rows_for_drug(bs_data: dict, drug: str) -> list[dict]:
    """Per-BioSample manifest rows for one drug (tier + label + censor_meta)."""
    rows = []
    for bs, per_drug in sorted(bs_data.items()):
        slot = per_drug.get(drug)
        if slot is None:
            continue
        tier = tier_for_isolate([mv.raw for mv in slot["mics"]], slot["calls"], drug)
        if tier in STRICT_TIERS:
            label, strict, relaxed = STRICT_TIERS[tier], True, True
        elif tier in RELAXED_EXTRA:
            label, strict, relaxed = RELAXED_EXTRA[tier], False, True
        else:
            label, strict, relaxed = "EXCLUDED", False, False
        rows.append({"biosample": bs, "drug": drug, "tier": tier, "label": label,
                     "strict": strict, "relaxed": relaxed,
                     "censor_meta": tier.startswith("CENSORED"), "conflict_status": "ok"})
    return rows


def _drug_inputs(bs_data: dict, drug: str) -> tuple[dict, dict]:
    """Build build_drug_labels inputs ({bs:[raw tokens]}, {bs:calls}) for one drug."""
    iso_mics = {bs: [mv.raw for mv in pd[drug]["mics"]] for bs, pd in bs_data.items() if drug in pd}
    iso_calls = {bs: pd[drug]["calls"] for bs, pd in bs_data.items() if drug in pd}
    return iso_mics, iso_calls


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", required=True)
    ap.add_argument("--mic-table", required=True)
    ap.add_argument("--key-col", required=True)
    ap.add_argument("--drug-col", action="append", default=[], metavar="COL=ALIAS")
    ap.add_argument("--call-col", action="append", default=[], metavar="COL=ALIAS")
    ap.add_argument("--run-id", default=f"run_{_date.today().isoformat()}")
    ap.add_argument("--out-root", default="data/raw")
    ap.add_argument("--wiki-dir", default="wiki")
    a = ap.parse_args()
    drug_cols = dict(kv.split("=", 1) for kv in a.drug_col)
    call_cols = dict(kv.split("=", 1) for kv in a.call_col)

    ingested = ingest_mic_table(a.mic_table, key_col=a.key_col, drug_cols=drug_cols, call_cols=call_cols)
    if ingested["unmapped_columns"]:
        print(f"NOTE unmapped columns: {ingested['unmapped_columns']}")

    from dna_decode.eval.biosample_resolver import BioSampleResolver
    records = BioSampleResolver().read_run_records_for_project(a.project)
    resolution = resolve_keys(list(ingested["data"].keys()), records)
    write_crosswalk(Path(a.wiki_dir) / f"oxford_crosswalk_{a.run_id}.json",
                    {"crosswalk": resolution["resolved"], "conflicts": resolution["conflicts"]})
    if resolution["conflicts"]:
        raise CrosswalkConflictError(
            f"{len(resolution['conflicts'])} crosswalk conflict(s); refusing to build labels. "
            f"See oxford_crosswalk_{a.run_id}.json.")
    if resolution["unresolved"]:
        print(f"NOTE {len(resolution['unresolved'])} unresolved MIC keys dropped (e.g. "
              f"{resolution['unresolved'][:5]})")

    bs_data, dropped = rekey_to_biosample(ingested["data"], resolution["resolved"])

    drugs = sorted({d for d in (canonical_drug(a) for a in drug_cols.values()) if d})
    manifest = []
    for drug in drugs:
        iso_mics, iso_calls = _drug_inputs(bs_data, drug)
        if not iso_mics:
            continue
        res = build_drug_labels(iso_mics, drug, iso_calls)
        assert all(is_biosample_key(k) for k in res["strict"])  # BioSample-keyed by construction
        write_labels(Path(a.out_root) / f"oxford_extval_{drug}", res)
        rows = manifest_rows_for_drug(bs_data, drug)
        manifest.extend(rows)
        print(f"  {drug}: strict {res['n_strict']} / relaxed {res['n_relaxed']} of {res['n_total']} BioSamples")

    man_path = Path(a.wiki_dir) / f"cohort_manifest_external_{a.run_id}.json"
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(
        {"_schema": "external-cohort-manifest-v1", "run_id": a.run_id, "project": a.project,
         "date": _date.today().isoformat(), "rows": manifest,
         "biosamples": sorted({r["biosample"] for r in manifest})}, indent=2), encoding="utf-8")
    print(f"manifest -> {man_path} ({len(manifest)} rows, {len({r['biosample'] for r in manifest})} BioSamples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
