"""The modality-hybrid test: is there ANY move past the ESM2-650M molecular ceiling — and does a HYBRID beat
its parts?

The world model works in exactly ONE regime (molecular fitness / DMS). At scale its sequence baseline is
ESM2-650M (median |Spearman| ~0.49, N~194; scale is a DEAD END — 3B/15B regress). The measured headroom is
MODALITY, not parameters: +evolution (MSA) and +structure. This sweep answers three things from ProteinGym's
OWN per-variant precomputed scores (`pg_zeroshot/*.csv`, 99 model columns incl. GEMME / MSA-Transformer /
ProSST / SaProt / TranceptEVE), PAIRED per protein against the ESM2-650M baseline:

  1. Which single MODALITY beats ESM2-650M? (evolution vs structure vs retrieval-hybrid)
  2. Does a naive rank-average HYBRID (ESM2 (+) GEMME, ESM2 (+) ProSST) beat BOTH its components?
     -- the direct test of the hybrid-world-model thesis (orthogonal signals combine).
  3. Re-verify the scale dead-end (ESM2 8M..15B) on this project's own paired terms.

Pre-registered bar (a candidate "beats ESM2-650M"): PAIRED median delta > 0 AND win-rate >= 60% AND a
two-sided sign test p < 0.05 over the assays where both have a valid Spearman. Falsifiable: if nothing
clears it, the single-modality ceiling is real and the world model is at its measured max in its one regime.

Metric: abs-Spearman(model score, DMS_score) per assay, mid-rank ties (the documented tie-order trap),
>= MIN_N shared non-NaN variants. Hybrids orient each component by ProteinGym's fixed higher=fitter
convention (verified in-run), rank-average, then Spearman -- no label is fit, only the standard orientation.

Run:  uv run python scripts/forward_modality_hybrid_sweep.py
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

from scripts.forward_blosum_proteingym_sweep import _spearman  # noqa: E402  (mid-rank Spearman)

ZS_DIR = Path("D:/dna_decode_cache/proteingym/pg_zeroshot")
REF_CSV = Path("D:/dna_decode_cache/proteingym/pg_reference.csv")
MIN_N = 20                       # shared non-NaN variants for a valid per-assay Spearman
MIN_ASSAYS_PER_CATEGORY = 8      # below this a per-category paired median is too thin to report
BASELINE = "ESM2_650M"

# The candidates whose per-phenotype-category behaviour answers "which modality for which trait".
CATEGORY_CANDIDATES = ["ProSST-2048", "SaProt_650M_AF2", "GEMME", "MSA_Transformer_ensemble",
                       "HYB_ESM2+GEMME", "HYB_ESM2+ProSST", "HYB_ESM2+GEMME+ProSST"]

# Curated model set, grouped by modality. Every name is a real column in pg_zeroshot/*.csv.
MODELS = {
    "sequence":   [BASELINE],
    "evolution":  ["GEMME", "EVmutation", "Site_Independent", "MSA_Transformer_ensemble"],
    "structure":  ["ProSST-2048", "SaProt_650M_AF2", "ESM-IF1", "S3F"],
    "retrieval":  ["TranceptEVE_L", "VenusREM", "ESCOTT"],
    "esm2_scale": ["ESM2_8M", "ESM2_35M", "ESM2_150M", "ESM2_650M", "ESM2_3B", "ESM2_15B"],
}
# Naive rank-average hybrids: sequence (+) a second modality. The world-model thesis test.
HYBRIDS = {
    "HYB_ESM2+GEMME":        (BASELINE, "GEMME"),
    "HYB_ESM2+ProSST":       (BASELINE, "ProSST-2048"),
    "HYB_ESM2+MSA_T":        (BASELINE, "MSA_Transformer_ensemble"),
    "HYB_ESM2+SaProt":       (BASELINE, "SaProt_650M_AF2"),
    "HYB_GEMME+ProSST":      ("GEMME", "ProSST-2048"),
}
HYB3 = {"HYB_ESM2+GEMME+ProSST": (BASELINE, "GEMME", "ProSST-2048")}


def _f(v: str) -> float | None:
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


def _midrank(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        mid = (i + j) / 2.0
        for k in range(i, j + 1):
            r[order[k]] = mid
        i = j + 1
    return r


def _oriented_ranks(scores: list[float], dms: list[float]) -> list[float]:
    """Mid-ranks oriented so higher rank = higher DMS (ProteinGym's fixed higher=fitter convention).
    Orientation uses the SIGN of the in-assay signed Spearman only -- no magnitude is fit."""
    rk = _midrank(scores)
    sign = 1.0 if _spearman(scores, dms) >= 0 else -1.0
    n = len(rk)
    return rk if sign > 0 else [n - 1 - x for x in rk]


def load_assay(path: Path) -> dict:
    """Return {col: [values]} with DMS_score + all model cols, rows kept as strings->float|None."""
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return rows


def spearman_for(rows: list[dict], col: str) -> tuple[float, int] | None:
    xs, ys = [], []
    for r in rows:
        a, b = _f(r.get(col, "")), _f(r.get("DMS_score", ""))
        if a is not None and b is not None:
            xs.append(a); ys.append(b)
    if len(xs) < MIN_N:
        return None
    return abs(_spearman(xs, ys)), len(xs)


def signed_spearman_for(rows: list[dict], col: str) -> float | None:
    xs, ys = [], []
    for r in rows:
        a, b = _f(r.get(col, "")), _f(r.get("DMS_score", ""))
        if a is not None and b is not None:
            xs.append(a); ys.append(b)
    if len(xs) < MIN_N:
        return None
    return _spearman(xs, ys)


def hybrid_spearman(rows: list[dict], cols: tuple[str, ...]) -> tuple[float, int] | None:
    """Rank-average of >=2 oriented score vectors over the shared non-NaN variants, then Spearman vs DMS."""
    shared = []
    for r in rows:
        vals = [_f(r.get(c, "")) for c in cols]
        d = _f(r.get("DMS_score", ""))
        if d is not None and all(v is not None for v in vals):
            shared.append((vals, d))
    if len(shared) < MIN_N:
        return None
    dms = [d for _, d in shared]
    per_col_ranks = []
    for ci in range(len(cols)):
        col_scores = [vals[ci] for vals, _ in shared]
        per_col_ranks.append(_oriented_ranks(col_scores, dms))
    avg = [sum(per_col_ranks[ci][i] for ci in range(len(cols))) / len(cols) for i in range(len(shared))]
    return abs(_spearman(avg, dms)), len(shared)


def sign_test_p(wins: int, losses: int) -> float:
    """Two-sided binomial sign test under p=0.5 (ties dropped)."""
    n = wins + losses
    if n == 0:
        return 1.0
    k = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def load_categories(ref_csv: Path) -> dict[str, str]:
    """DMS_id -> ProteinGym coarse phenotype category (Activity/Binding/Expression/Stability/OrganismalFitness)."""
    if not ref_csv.exists():
        return {}
    with ref_csv.open(encoding="utf-8") as fh:
        return {row["DMS_id"]: row.get("coarse_selection_type", "") for row in csv.DictReader(fh)}


def paired_in_subset(per_assay: dict, model: str, dms_ids: set[str]) -> dict:
    """Paired deltas abs_spearman(model) - abs_spearman(BASELINE) over the given assay subset."""
    deltas = []
    for dms_id in dms_ids:
        m = per_assay.get(dms_id, {})
        if model in m and BASELINE in m:
            deltas.append(m[model] - m[BASELINE])
    wins = sum(1 for d in deltas if d > 1e-9)
    losses = sum(1 for d in deltas if d < -1e-9)
    return {
        "n_paired": len(deltas),
        "median_delta": round(statistics.median(deltas), 4) if deltas else None,
        "win_rate": round(wins / len(deltas), 3) if deltas else None,
        "sign_test_p": round(sign_test_p(wins, losses), 5) if deltas else None,
    }


def by_category(per_assay: dict, categories: dict[str, str]) -> dict:
    """Per-phenotype-category paired lift of each modality candidate vs ESM2-650M -- 'which modality for
    which trait'. Only categories with >= MIN_ASSAYS_PER_CATEGORY scored assays are reported."""
    cats: dict[str, set[str]] = {}
    for dms_id in per_assay:
        c = categories.get(dms_id, "")
        if c and BASELINE in per_assay[dms_id]:
            cats.setdefault(c, set()).add(dms_id)
    out = {}
    for cat, ids in sorted(cats.items(), key=lambda kv: -len(kv[1])):
        base_vals = [per_assay[i][BASELINE] for i in ids]
        row = {"n_assays": len(ids),
               "reportable": len(ids) >= MIN_ASSAYS_PER_CATEGORY,
               "baseline_median": round(statistics.median(base_vals), 4) if base_vals else None,
               "candidates": {}}
        for model in CATEGORY_CANDIDATES:
            row["candidates"][model] = paired_in_subset(per_assay, model, ids)
        out[cat] = row
    return out


def paired_vs_baseline(per_assay: dict, model: str) -> dict:
    """Paired deltas abs_spearman(model) - abs_spearman(BASELINE) over assays where both are valid."""
    deltas = []
    for dms_id, m in per_assay.items():
        if model in m and BASELINE in m:
            deltas.append(m[model] - m[BASELINE])
    wins = sum(1 for d in deltas if d > 1e-9)
    losses = sum(1 for d in deltas if d < -1e-9)
    return {
        "model": model,
        "n_paired": len(deltas),
        "median_delta": round(statistics.median(deltas), 4) if deltas else None,
        "mean_delta": round(statistics.mean(deltas), 4) if deltas else None,
        "win_rate": round(wins / len(deltas), 3) if deltas else None,
        "wins": wins, "losses": losses,
        "sign_test_p": round(sign_test_p(wins, losses), 5),
        "beats_baseline": bool(deltas and statistics.median(deltas) > 0
                               and (wins / len(deltas)) >= 0.60 and sign_test_p(wins, losses) < 0.05),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--zs-dir", default=str(ZS_DIR))
    ap.add_argument("--ref-csv", default=str(REF_CSV))
    ap.add_argument("--by-category", action="store_true",
                    help="also break the modality lift down by phenotype category (which modality for which trait)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    zs = Path(args.zs_dir)
    files = sorted(zs.glob("*.csv"))
    if not files:
        print(f"NO ZEROSHOT CSVs at {zs}", file=sys.stderr)
        return 2

    all_cols = sorted({c for grp in MODELS.values() for c in grp})
    per_assay: dict[str, dict[str, float]] = {}   # dms_id -> {model: abs_spearman}
    per_assay_n: dict[str, dict[str, int]] = {}
    orient_ok = {"positive": 0, "negative": 0}    # convention check on the baseline

    for fp in files:
        rows = load_assay(fp)
        dms_id = fp.stem
        m: dict[str, float] = {}
        mn: dict[str, int] = {}
        for col in all_cols:
            res = spearman_for(rows, col)
            if res:
                m[col], mn[col] = res
        # hybrids
        for hname, cols in {**HYBRIDS, **HYB3}.items():
            res = hybrid_spearman(rows, cols)
            if res:
                m[hname], mn[hname] = res
        # orientation audit on baseline
        s = signed_spearman_for(rows, BASELINE)
        if s is not None:
            orient_ok["positive" if s >= 0 else "negative"] += 1
        per_assay[dms_id] = m
        per_assay_n[dms_id] = mn

    # baseline headline
    base_vals = sorted(m[BASELINE] for m in per_assay.values() if BASELINE in m)
    base_median = statistics.median(base_vals) if base_vals else None

    # per-model median + paired-vs-baseline
    candidates = [c for c in all_cols if c != BASELINE] + list(HYBRIDS) + list(HYB3)
    report = []
    for model in [BASELINE] + candidates:
        vals = sorted(m[model] for m in per_assay.values() if model in m)
        row = {"model": model,
               "n_assays": len(vals),
               "median_abs_spearman": round(statistics.median(vals), 4) if vals else None,
               "mean_abs_spearman": round(statistics.mean(vals), 4) if vals else None}
        if model != BASELINE:
            row.update(paired_vs_baseline(per_assay, model))
        report.append(row)

    beats = [r for r in report if r.get("beats_baseline")]
    beats.sort(key=lambda r: (-(r["median_delta"] or -9), -(r["win_rate"] or 0)))

    summary = {
        "generated": date.today().isoformat(),
        "substrate": "ProteinGym pg_zeroshot per-variant scores",
        "n_assays_total": len(files),
        "min_variants_per_assay": MIN_N,
        "baseline": BASELINE,
        "baseline_median_abs_spearman": round(base_median, 4) if base_median else None,
        "orientation_audit_baseline": orient_ok,
        "prereg_bar": "median_delta>0 AND win_rate>=0.60 AND sign_test_p<0.05 (paired per-assay)",
        "n_candidates_beating_baseline": len(beats),
        "top_beaters": beats[:8],
        "report": sorted(report, key=lambda r: -(r["median_abs_spearman"] or -9)),
    }
    if args.by_category:
        cats = load_categories(Path(args.ref_csv))
        summary["by_category"] = by_category(per_assay, cats)

    out = Path(args.out) if args.out else REPO / "wiki" / f"forward_modality_hybrid_{date.today().isoformat()}.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"baseline {BASELINE} median abs-Spearman = {summary['baseline_median_abs_spearman']} "
          f"(N={report[0]['n_assays']})")
    print(f"orientation audit (baseline signed Spearman): {orient_ok}")
    print(f"candidates beating baseline (prereg bar): {len(beats)}")
    for b in beats[:8]:
        print(f"  {b['model']:26s} median={b['median_abs_spearman']}  "
              f"delta={b['median_delta']:+}  win={b['win_rate']}  p={b['sign_test_p']}  N={b['n_paired']}")
    if args.by_category and "by_category" in summary:
        print("\n== which modality for which trait (paired median delta vs ESM2-650M) ==")
        show = ["ProSST-2048", "GEMME", "HYB_ESM2+GEMME+ProSST"]
        print(f"{'category':20s} {'N':>3} {'baseline':>8}  " + "  ".join(f"{m[:14]:>14}" for m in show))
        for cat, row in summary["by_category"].items():
            if not row["reportable"]:
                continue
            cells = []
            for m in show:
                c = row["candidates"][m]
                cells.append(f"{c['median_delta']:+.4f}({c['win_rate']})" if c["median_delta"] is not None else "     -")
            print(f"{cat:20s} {row['n_assays']:>3} {row['baseline_median']:>8}  " + "  ".join(f"{x:>14}" for x in cells))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
