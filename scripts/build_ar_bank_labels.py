"""Build BioSample-keyed R/S labels for the FROZEN external-cohort arm from the CDC AR Isolate Bank.

Enumerate the (public) AR Isolate Bank, filter to an organism, fetch each isolate's reference
broth-microdilution MIC, tier it through the FROZEN `build_drug_labels`, and emit — per drug —
`data/raw/<cohort>_extval_<drug>/selected_{strict,relaxed}.tsv` + `buckets_<drug>.json`, plus the
single `wiki/cohort_manifest_external_<run_id>.json` that `external_cohort_preflight` (exact-set),
`external_cohort_revalidate` (drift-guard), and the roll-up consume.

Unlike the Oxford ingester there is NO crosswalk — the AR Bank isolate page carries the NCBI
BioSample directly, so labels are BioSample-keyed by construction. The MIC is CDC reference BMD (a
real G1 phenotype) and the isolates are curated CDC outbreak/surveillance isolates -> provenance-
separable from the decoder's NCBI-PD tuning set (the BioSample-level preflight verifies disjointness).

Run (network; cache-first, resumable):
  uv run python -m scripts.build_ar_bank_labels --cohort ar_bank_ecoli \
      --organism Escherichia --drug ciprofloxacin --drug ceftriaxone --drug gentamicin
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import ar_isolate_bank as ab
from dna_decode.data.external_mic_labels import (
    RELAXED_EXTRA,
    STRICT_TIERS,
    build_drug_labels,
    canonical_drug,
    tier_for_isolate,
    write_labels,
)


def collect_details(cohort_rows, cache_dir: str, limit: int | None = None,
                    offline_ok: bool = False) -> list[ab.IsolateDetail]:
    """Fetch + parse the MIC detail page for each enumerated isolate (cache-first, resumable)."""
    details: list[ab.IsolateDetail] = []
    rows = cohort_rows[:limit] if limit else cohort_rows
    for i, row in enumerate(rows, 1):
        det = ab.fetch_isolate_detail(row, cache_dir, offline_ok=offline_ok)
        if det is not None and det.biosample:
            details.append(det)
        if i % 25 == 0:
            print(f"  ...{i}/{len(rows)} detail pages ({len(details)} parsed)")
    return details


def manifest_rows_for_drug(details, drug: str) -> list[dict]:
    """Per-BioSample manifest rows for one drug (tier + label + censor flag), BioSample-keyed."""
    canon = canonical_drug(drug)
    iso_mics, iso_calls = ab.to_label_inputs(details, drug)
    rows = []
    for bs, tokens in sorted(iso_mics.items()):
        tier = tier_for_isolate(tokens, iso_calls.get(bs, set()), canon)
        if tier in STRICT_TIERS:
            label, strict, relaxed = STRICT_TIERS[tier], True, True
        elif tier in RELAXED_EXTRA:
            label, strict, relaxed = RELAXED_EXTRA[tier], False, True
        else:
            label, strict, relaxed = "EXCLUDED", False, False
        rows.append({"biosample": bs, "drug": canon, "tier": tier, "label": label,
                     "strict": strict, "relaxed": relaxed,
                     "censor_meta": tier.startswith("CENSORED"), "conflict_status": "ok"})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", default="ar_bank_ecoli", help="cohort slug (drives output paths)")
    ap.add_argument("--organism", action="append", default=[], metavar="NEEDLE",
                    help="organism substring filter (repeatable), e.g. Escherichia")
    ap.add_argument("--drug", action="append", default=[], metavar="DRUG",
                    help="pilot drug (repeatable); default ciprofloxacin/ceftriaxone/gentamicin")
    ap.add_argument("--cache", default="data/raw/ar_isolate_bank/pages", help="HTML page cache dir")
    ap.add_argument("--run-id", default=f"arbank_{_date.today().isoformat()}")
    ap.add_argument("--out-root", default="data/raw")
    ap.add_argument("--wiki-dir", default="wiki")
    ap.add_argument("--limit", type=int, default=None, help="cap isolates (smoke)")
    ap.add_argument("--offline-ok", action="store_true", help="use cache only; skip network misses")
    ap.add_argument("--details-json", default=None,
                    help="optional: write/reuse parsed details here to avoid re-fetch")
    a = ap.parse_args()
    drugs = a.drug or ["ciprofloxacin", "ceftriaxone", "gentamicin"]
    needles = tuple(a.organism) or None

    print(f"Enumerating AR Bank (organism={needles})...")
    rows = ab.enumerate_all(a.cache, organism_needles=needles, offline_ok=a.offline_ok)
    print(f"  {len(rows)} unique-BioSample isolates")

    details = None
    dj = Path(a.details_json) if a.details_json else None
    if dj and dj.exists():
        raw = json.loads(dj.read_text(encoding="utf-8"))
        details = [ab.IsolateDetail(**d) for d in raw]
        print(f"  reused {len(details)} parsed details from {dj}")
    if details is None:
        details = collect_details(rows, a.cache, limit=a.limit, offline_ok=a.offline_ok)
        if dj:
            dj.parent.mkdir(parents=True, exist_ok=True)
            dj.write_text(json.dumps([{**d.__dict__} for d in details], indent=2), encoding="utf-8")
    print(f"  {len(details)} isolates with BioSample + MIC parsed")

    manifest: list[dict] = []
    for drug in drugs:
        canon = canonical_drug(drug)
        if canon is None:
            print(f"  SKIP non-pilot drug {drug!r}")
            continue
        iso_mics, iso_calls = ab.to_label_inputs(details, drug)
        if not iso_mics:
            print(f"  {canon}: no isolates carry this drug")
            continue
        res = build_drug_labels(iso_mics, canon, iso_calls)
        write_labels(Path(a.out_root) / f"{a.cohort}_extval_{canon}", res)
        manifest.extend(manifest_rows_for_drug(details, drug))
        nR = sum(1 for v in res["strict"].values() if v == "R")
        nS = sum(1 for v in res["strict"].values() if v == "S")
        print(f"  {canon}: strict {res['n_strict']} ({nR}R/{nS}S) / relaxed {res['n_relaxed']} "
              f"of {res['n_total']}  buckets={res['buckets']}")

    man_path = Path(a.wiki_dir) / f"cohort_manifest_external_{a.run_id}.json"
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(
        {"_schema": "external-cohort-manifest-v1", "run_id": a.run_id, "project": a.cohort,
         "source": "cdc_fda_ar_isolate_bank", "organism_filter": list(needles or []),
         "date": _date.today().isoformat(), "rows": manifest,
         "biosamples": sorted({r["biosample"] for r in manifest})}, indent=2), encoding="utf-8")
    print(f"manifest -> {man_path} ({len(manifest)} rows, "
          f"{len({r['biosample'] for r in manifest})} BioSamples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
