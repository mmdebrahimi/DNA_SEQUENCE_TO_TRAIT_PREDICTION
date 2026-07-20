"""Generalized AR Isolate Bank INT label builder — any organism in `ar_bank_registry`.

Enumerate the AR Bank for `--organism`'s needle, fetch each isolate's detail, and label each drug from the
CDC S/I/R INT column (R vs S; INTERMEDIATE excluded), emitting the selected_strict.tsv + buckets +
cohort_manifest the preflight + generalized scorer consume. Generalizes `build_ar_bank_gono_labels` via the
registry (drug -> AR-Bank INT column name). Frozen surface untouched.

  uv run python -m scripts.build_ar_bank_organism_labels --organism enterococcus_faecium
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data import ar_isolate_bank as ab
from dna_decode.organism_rules.ar_bank_registry import config_for, shown_name_map


def int_labels(details, shown: str) -> dict[str, str]:
    """{biosample: 'R'|'S'} from the CDC INT call for the AR-Bank column `shown`; I/SDD/NS excluded."""
    out: dict[str, str] = {}
    for d in details:
        if not d.biosample:
            continue
        c = (d.calls.get(shown) or "").strip().upper()
        if c in ("R", "S"):
            out[d.biosample] = c
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, help="registry key (e.g. enterococcus_faecium)")
    ap.add_argument("--cache", default="data/raw/ar_isolate_bank/pages")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--out-root", default="data/raw")
    ap.add_argument("--wiki-dir", default="wiki")
    ap.add_argument("--details-json", default=None)
    ap.add_argument("--exclude-biosamples", default=None)
    ap.add_argument("--offline-ok", action="store_true")
    a = ap.parse_args()
    cfg = config_for(a.organism)
    run_id = a.run_id or f"arbank_{a.organism}_{_date.today().isoformat()}"
    drug_shown = shown_name_map(a.organism)

    print(f"Enumerating AR Bank {cfg['needle']!r}...")
    rows = ab.enumerate_all(a.cache, organism_needles=(cfg["needle"],), offline_ok=a.offline_ok)
    print(f"  {len(rows)} isolates")

    details = None
    dj = Path(a.details_json) if a.details_json else None
    if dj and dj.exists():
        details = [ab.IsolateDetail(**d) for d in json.loads(dj.read_text(encoding="utf-8"))]
        print(f"  reused {len(details)} details from {dj}")
    if details is None:
        details = []
        for i, r in enumerate(rows, 1):
            d = ab.fetch_isolate_detail(r, a.cache, offline_ok=a.offline_ok)
            if d is not None and d.biosample:
                details.append(d)
            if i % 25 == 0:
                print(f"  ...{i}/{len(rows)} ({len(details)} parsed)")
        if dj:
            dj.parent.mkdir(parents=True, exist_ok=True)
            dj.write_text(json.dumps([{**d.__dict__} for d in details], indent=2), encoding="utf-8")
    print(f"  {len(details)} isolates with BioSample + MIC")

    excluded: set[str] = set()
    if a.exclude_biosamples:
        txt = Path(a.exclude_biosamples).read_text(encoding="utf-8")
        excluded = {t.strip() for t in re.split(r"[,\s]+", txt) if t.strip()}
        before = len(details)
        details = [d for d in details if d.biosample not in excluded]
        print(f"  leakage exclusion: dropped {before - len(details)} -> {len(details)} disjoint")

    manifest = []
    for drug, shown in drug_shown.items():
        labels = int_labels(details, shown)
        if not labels:
            print(f"  {drug} ({shown}): no INT R/S")
            continue
        outd = Path(a.out_root) / f"ar_bank_{a.organism}_extval_{drug}"
        outd.mkdir(parents=True, exist_ok=True)
        tsv = "".join(f"{bs}\t{rs}\n" for bs, rs in sorted(labels.items()))
        (outd / "selected_strict.tsv").write_text(tsv, encoding="utf-8")
        (outd / "selected_relaxed.tsv").write_text(tsv, encoding="utf-8")
        nR = sum(1 for v in labels.values() if v == "R")
        (outd / f"buckets_{drug}.json").write_text(json.dumps(
            {"drug": drug, "shown": shown, "label_source": "ar_bank_INT_cdc",
             "n_total": len(labels), "n_R": nR, "n_S": len(labels) - nR}, indent=2), encoding="utf-8")
        for bs, rs in labels.items():
            manifest.append({"biosample": bs, "drug": drug, "label": rs, "strict": True, "relaxed": True})
        print(f"  {drug} ({shown}): {len(labels)} ({nR}R/{len(labels) - nR}S)")

    mp = Path(a.wiki_dir) / f"cohort_manifest_external_{run_id}.json"
    mp.parent.mkdir(parents=True, exist_ok=True)
    mp.write_text(json.dumps(
        {"_schema": "external-cohort-manifest-v1", "run_id": run_id, "project": f"ar_bank_{a.organism}",
         "source": "cdc_fda_ar_isolate_bank", "organism_filter": [cfg["needle"]],
         "label_source": "ar_bank_INT_cdc", "leakage_excluded": sorted(excluded),
         "date": _date.today().isoformat(), "rows": manifest,
         "biosamples": sorted({r["biosample"] for r in manifest})}, indent=2), encoding="utf-8")
    print(f"manifest -> {mp} ({len(manifest)} rows, {len({r['biosample'] for r in manifest})} BioSamples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
