"""Discovery: do the ready-made multimodal SOTA models (SaProt / VenusREM / ...) beat our ESM2+ProSST?

The forward cell's validated deployable lift is the 2-way `ESM2+ProSST` (+0.067 vs ESM2, 93%, N=56). The
directive: "discover ready-made multimodal SOTA (SaProt/VenusREM)". The DECISIVE discovery question is whether
adopting one of those off-the-shelf models would BEAT our hybrid — answerable PAIRED per-protein from
ProteinGym's PRECOMPUTED columns (SaProt/VenusREM/ESCOTT/S3F/ProtSSN are all on the leaderboard), with ZERO
install. The install matters only if the answer is "yes, adopt one".

Method: over every ProteinGym assay with all needed columns, abs-Spearman(col, DMS_score) with mid-rank ties.
Our hybrids are rank-averages of the PRECOMPUTED component columns (apples-to-apples vs the ready-made
columns). Each ready-made model is compared PAIRED (per-protein delta, sign test) vs our 2-way ESM2+ProSST —
NOT a difference of medians (the documented trap).

Run:  uv run python scripts/multimodal_sota_discovery.py
Exit: 0 = ran; 2 = substrate unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.forward_blosum_proteingym_sweep import _spearman           # noqa: E402
from dna_decode.forward.variant_effect import rank_average_hybrid        # noqa: E402

ZS_DIR = Path("D:/dna_decode_cache/proteingym/pg_zeroshot")

# our validated components + the ready-made SOTA columns to compare against
OUR_2WAY = ("ESM2_650M", "ProSST-2048")
OUR_3WAY = ("ESM2_650M", "ProSST-2048", "GEMME")
READY_MADE = ["VenusREM", "SaProt_650M_AF2", "ESCOTT", "S3F_MSA", "ProtSSN_ensemble", "MSA_Transformer_ensemble"]
NEEDED = set(OUR_3WAY) | set(READY_MADE)
MIN_N = 20


def _f(v):
    try:
        x = float(v); return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


def _hybrid_spearman(rows, cols, dms_by_mut):
    tables = []
    for c in cols:
        t = {r["mutant"]: _f(r.get(c)) for r in rows if _f(r.get(c)) is not None}
        tables.append(t)
    shared = set(tables[0])
    for t in tables[1:]:
        shared &= set(t)
    shared &= set(dms_by_mut)
    if len(shared) < MIN_N:
        return None
    combined = rank_average_hybrid([{k: t[k] for k in shared} for t in tables])
    ks = sorted(shared)
    return abs(_spearman([combined[k] for k in ks], [dms_by_mut[k] for k in ks]))


def _col_spearman(rows, col, dms_by_mut):
    xs, ys = [], []
    for r in rows:
        a = _f(r.get(col))
        if a is not None and r["mutant"] in dms_by_mut:
            xs.append(a); ys.append(dms_by_mut[r["mutant"]])
    return abs(_spearman(xs, ys)) if len(xs) >= MIN_N else None


def sign_p(wins, losses):
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    files = sorted(ZS_DIR.glob("*.csv"))
    if not files:
        print("no pg_zeroshot CSVs", file=sys.stderr); return 2

    per_assay = {}
    for fp in files:
        with fp.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        if not rows or not all(c in rows[0] for c in NEEDED):
            continue
        dms = {r["mutant"]: _f(r.get("DMS_score")) for r in rows if _f(r.get("DMS_score")) is not None}
        if len(dms) < MIN_N:
            continue
        rec = {"two_way": _hybrid_spearman(rows, OUR_2WAY, dms),
               "three_way": _hybrid_spearman(rows, OUR_3WAY, dms)}
        for m in READY_MADE:
            rec[m] = _col_spearman(rows, m, dms)
        if rec["two_way"] is not None:
            per_assay[fp.stem] = rec

    n = len(per_assay)
    if not n:
        print("no usable assays", file=sys.stderr); return 2

    def med(key):
        vals = [r[key] for r in per_assay.values() if r.get(key) is not None]
        return round(statistics.median(vals), 4) if vals else None

    # PAIRED: each ready-made model (and our 3-way) vs our 2-way ESM2+ProSST
    report = []
    for model in ["three_way", *READY_MADE]:
        deltas = [r[model] - r["two_way"] for r in per_assay.values()
                  if r.get(model) is not None and r.get("two_way") is not None]
        w = sum(1 for d in deltas if d > 1e-9); l = sum(1 for d in deltas if d < -1e-9)
        report.append({"model": model, "n": len(deltas), "median_self": med(model),
                       "median_delta_vs_2way": round(statistics.median(deltas), 4) if deltas else None,
                       "win_vs_2way": f"{w}/{len(deltas)}", "sign_p": round(sign_p(w, l), 4),
                       "beats_2way": bool(deltas and statistics.median(deltas) > 0 and w > l and sign_p(w, l) < 0.05)})
    report.sort(key=lambda r: -(r["median_self"] or -9))

    summary = {"generated": date.today().isoformat(), "n_assays": n,
               "our_2way_median": med("two_way"), "our_3way_median": med("three_way"),
               "baseline": "our ESM2+ProSST 2-way (rank-avg of precomputed columns)",
               "report_paired_vs_2way": report}
    out = Path(args.out) if args.out else REPO / "wiki" / f"multimodal_sota_discovery_{date.today().isoformat()}.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"N={n} assays | our 2-way ESM2+ProSST median |Spearman| = {med('two_way')} "
          f"| our 3-way = {med('three_way')}")
    print(f"{'model':26s} {'median':>7} {'Δ vs 2way':>10} {'win':>8} {'sign-p':>7}  beats-2way")
    for r in report:
        print(f"{r['model']:26s} {r['median_self']!s:>7} {r['median_delta_vs_2way']:>+10} "
              f"{r['win_vs_2way']:>8} {r['sign_p']!s:>7}  {'YES' if r['beats_2way'] else ''}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
