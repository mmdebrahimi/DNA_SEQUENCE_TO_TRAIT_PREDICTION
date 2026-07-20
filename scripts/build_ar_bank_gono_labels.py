"""Build BioSample-keyed R/S labels for the AR Bank N. gonorrhoeae cohort from CDC's OWN S/I/R calls.

Gonococcal CLSI breakpoints are SPECIES-SPECIFIC (NOT the Enterobacterales values in the frozen
`mic_tiers.py` -- e.g. gono ceftriaxone S<=0.25 vs Enterobacterales S<=1). Re-tiering gono MICs through
the frozen Enterobacterales breakpoints would mislabel; and `mic_tiers.py` is FROZEN so gono breakpoints
can't be added there. So this builder uses the **AR Bank isolate page's own INT (S/I/R) interpretation**
-- CDC applied the correct gonococcal breakpoints -- as the label (R vs S; INTERMEDIATE excluded). This is
categorical (no 4x strict-MIC margin); a non-frozen gono strict-breakpoint path is a later refinement.

Emits per drug `data/raw/ar_bank_gono_extval_<drug>/selected_{strict,relaxed}.tsv` + `buckets_<drug>.json`
(strict == relaxed here: both are the INT-derived R/S set) + the single
`wiki/cohort_manifest_external_<run_id>.json` consumed by the preflight (exact-set) + the gono scorer
(drift guard). Reuses the ingester (`dna_decode/data/ar_isolate_bank`) for enumerate/fetch/parse.

Run (network; cache-first, resumable):
  uv run python -m scripts.build_ar_bank_gono_labels --details-json data/raw/ar_isolate_bank/gono_details.json
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

# canonical gono drug -> the drug name as shown on the AR Bank isolate page.
GONO_DRUGS = {
    "azithromycin": "Azithromycin", "cefixime": "Cefixime", "ceftriaxone": "Ceftriaxone",
    "ciprofloxacin": "Ciprofloxacin", "penicillin": "Penicillin", "tetracycline": "Tetracycline",
}
# gentamicin is on the AR Bank panel but has NO validated determinant -> the cell ABSTAINS (not built here).


def int_labels_for_drug(details, drug: str) -> dict[str, str]:
    """{biosample: 'R'|'S'} from CDC's INT call for `drug`. INTERMEDIATE / SDD / NS are excluded."""
    shown = GONO_DRUGS[drug]
    out: dict[str, str] = {}
    for d in details:
        if not d.biosample:
            continue
        call = (d.calls.get(shown) or "").strip().upper()
        if call in ("R", "S"):
            out[d.biosample] = call
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", default="ar_bank_gono")
    ap.add_argument("--drug", action="append", default=[], help="canonical gono drug (repeatable); default all 6")
    ap.add_argument("--cache", default="data/raw/ar_isolate_bank/pages")
    ap.add_argument("--run-id", default=f"arbank_gono_{_date.today().isoformat()}")
    ap.add_argument("--out-root", default="data/raw")
    ap.add_argument("--wiki-dir", default="wiki")
    ap.add_argument("--details-json", default=None, help="write/reuse parsed details to avoid re-fetch")
    ap.add_argument("--exclude-biosamples", default=None, help="leakage-overlap BioSamples to drop")
    ap.add_argument("--offline-ok", action="store_true")
    a = ap.parse_args()
    drugs = a.drug or list(GONO_DRUGS)
    bad = [d for d in drugs if d not in GONO_DRUGS]
    if bad:
        print(f"REFUSE: non-gono-scorable drug(s) {bad} (scorable: {sorted(GONO_DRUGS)}; gentamicin ABSTAINS)")
        return 2

    print("Enumerating AR Bank gonococcus...")
    rows = ab.enumerate_all(a.cache, organism_needles=("gonorrhoeae",), offline_ok=a.offline_ok)
    print(f"  {len(rows)} unique-BioSample gono isolates")

    details = None
    dj = Path(a.details_json) if a.details_json else None
    if dj and dj.exists():
        details = [ab.IsolateDetail(**d) for d in json.loads(dj.read_text(encoding="utf-8"))]
        print(f"  reused {len(details)} parsed details from {dj}")
    if details is None:
        details = []
        for i, row in enumerate(rows, 1):
            det = ab.fetch_isolate_detail(row, a.cache, offline_ok=a.offline_ok)
            if det is not None and det.biosample:
                details.append(det)
            if i % 25 == 0:
                print(f"  ...{i}/{len(rows)} detail pages ({len(details)} parsed)")
        if dj:
            dj.parent.mkdir(parents=True, exist_ok=True)
            dj.write_text(json.dumps([{**d.__dict__} for d in details], indent=2), encoding="utf-8")
    print(f"  {len(details)} isolates with BioSample + MIC parsed")

    excluded: set[str] = set()
    if a.exclude_biosamples:
        txt = Path(a.exclude_biosamples).read_text(encoding="utf-8")
        excluded = {t.strip() for t in re.split(r"[,\s]+", txt) if t.strip()}
        before = len(details)
        details = [d for d in details if d.biosample not in excluded]
        print(f"  leakage exclusion: dropped {before - len(details)} of {len(excluded)} listed -> {len(details)} disjoint")

    manifest: list[dict] = []
    for drug in drugs:
        labels = int_labels_for_drug(details, drug)
        if not labels:
            print(f"  {drug}: no INT R/S calls")
            continue
        out_dir = Path(a.out_root) / f"{a.cohort}_extval_{drug}"
        out_dir.mkdir(parents=True, exist_ok=True)
        tsv = "".join(f"{bs}\t{rs}\n" for bs, rs in sorted(labels.items()))
        (out_dir / "selected_strict.tsv").write_text(tsv, encoding="utf-8")     # strict == relaxed (INT-derived)
        (out_dir / "selected_relaxed.tsv").write_text(tsv, encoding="utf-8")
        nR = sum(1 for v in labels.values() if v == "R"); nS = len(labels) - nR
        (out_dir / f"buckets_{drug}.json").write_text(json.dumps(
            {"drug": drug, "label_source": "ar_bank_INT_cdc_gonococcal_breakpoints",
             "n_total": len(labels), "n_R": nR, "n_S": nS}, indent=2), encoding="utf-8")
        for bs, rs in sorted(labels.items()):
            manifest.append({"biosample": bs, "drug": drug, "label": rs, "strict": True, "relaxed": True,
                             "label_source": "ar_bank_INT", "conflict_status": "ok"})
        print(f"  {drug}: {len(labels)} labels ({nR}R/{nS}S)")

    man_path = Path(a.wiki_dir) / f"cohort_manifest_external_{a.run_id}.json"
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(
        {"_schema": "external-cohort-manifest-v1", "run_id": a.run_id, "project": a.cohort,
         "source": "cdc_fda_ar_isolate_bank", "organism_filter": ["gonorrhoeae"],
         "label_source": "ar_bank_INT_cdc_gonococcal_breakpoints",
         "leakage_excluded": sorted(excluded), "n_leakage_excluded": len(excluded),
         "date": _date.today().isoformat(), "rows": manifest,
         "biosamples": sorted({r["biosample"] for r in manifest})}, indent=2), encoding="utf-8")
    print(f"manifest -> {man_path} ({len(manifest)} rows, {len({r['biosample'] for r in manifest})} BioSamples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
