"""HIV epistasis + quantitative-effect model — Family A of the genome world-model plan (2026-07-11).

QUESTION: does a CONSTRAINED pairwise-interaction (epistasis) model of HIV drug-resistance log10
fold-change BEAT the additive mutation-presence baseline OUT-OF-SAMPLE, on the FREE, independent
Stanford HIVDB PhenoSense labels? The deployed decoder rule is additive (count/OR); reality has
higher-order structure (accessory mutations, TAMs, PI 82+54). This tests whether that structure is
learnable from the data we already hold. A NEGATIVE ("additive suffices for HIV DR") is a valid,
shippable world-model finding — the point is the honest pre-registered comparison, not that epistasis wins.

PRE-REGISTERED DESIGN (derived, not asserted — plan `plans/Genome_World_Model_Creative_Data_Reuse_Plan_2026-07-11.md`):
  * Data: data/raw/hiv/{NNRTI,NRTI,PI,INI}_DataSet.txt (gitignored; cite Rhee 2003 Nucleic Acids Res 31:298).
    Cols: SeqID, per-drug fold-change (auto-detected = non-SeqID non-P<digits>), P<pos> = amino acid at position.
  * Additive feature space: binary mutation-presence token `<pos><AA>`, kept iff present in >= MIN_MUTS
    isolates (=10; the DRMcv.R default, `hiv_nnrti_baseline.py`).
  * Interaction space: pairwise products among the TOP_K_FOR_PAIRS most-prevalent features, kept iff the
    pair CO-OCCURS in >= MIN_COOC isolates (=10). Elastic-net (L1+L2) zeros out unsupported pairs.
  * Model: ElasticNetCV (inner CV picks alpha + l1_ratio) inside a StandardScaler pipeline; nested 5-fold
    OUTER via cross_val_predict -> out-of-fold (OOF) predictions. SAME outer folds + SAME model class for
    BOTH additive and interaction => a fair paired comparison (never a null strawman).
  * Metric: Spearman rho + R2 of OOF vs log10 fold. "BEAT" = paired bootstrap (B=BOOT) on the OOF
    predictions: delta_rho (interaction - additive) 95%-CI lower bound > 0, per drug.
  * Powered drug-cell: n_isolates >= N_MIN (=30; the existing harness floor) AND non-degenerate fold spread.
  * PASS = CI-positive interaction gain on >= PASS_FRACTION (=0.5) of powered drug-cells across all classes.
  * verify-in-batch: the top interaction coefficients are reported per drug for a biology sanity check
    (are they real accessory/primary DRM pairs, or overfit noise?).

Honest scope: censored folds ('>'/'<') are kept at the numeric bound (rank-metric-safe; slightly biases R2 —
disclosed). This is a mechanism-FEATURE statistical model, NOT a sequence embedding (the closed 0-for-5 arm).
CPU-only, data in hand, no network / Docker / GPU / money. Frozen AMR surface is untouched (READ-only on data).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.hiv_nnrti_validate import _parse_fold, load_rows  # reuse the tested parsers

# --- frozen pre-registration constants (derived; see module docstring) ---
MIN_MUTS = 10            # feature kept iff present in >= this many isolates (DRMcv.R default)
MIN_COOC = 10            # interaction pair kept iff it co-occurs in >= this many isolates
TOP_K_FOR_PAIRS = 40     # form pairwise interactions among the top-K most-prevalent features
MAX_PAIRS = 400          # hard cap on interaction terms (dimensionality control)
N_MIN = 30               # powered drug-cell floor (matches hiv_nnrti_baseline's `< 30 -> too few`)
CV_OUTER = 5
BOOT = 1000
PASS_FRACTION = 0.5      # PASS iff >= this fraction of powered cells are CI-positive
SEED = 0
_L1_GRID = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
_CLASSES = {"NNRTI": "NNRTI_DataSet.txt", "NRTI": "NRTI_DataSet.txt",
            "PI": "PI_DataSet.txt", "INI": "INI_DataSet.txt"}


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho via rank-then-Pearson (numpy only; ties = average rank)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = _avg_rank(x); ry = _avg_rank(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def _avg_rank(v: np.ndarray) -> np.ndarray:
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v), float)
    sv = v[order]
    i = 0
    while i < len(v):
        j = i
        while j < len(v) and sv[j] == sv[i]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2.0 + 1.0
        i = j
    return ranks


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, float); y_pred = np.asarray(y_pred, float)
    ss_res = float(((y_true - y_pred) ** 2).sum())
    ss_tot = float(((y_true - y_true.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def _position_columns(header: list[str]) -> list[str]:
    return [h for h in header if len(h) > 1 and h[0] == "P" and h[1:].isdigit()]


def load_class(path: Path):
    """Return (rows, drug_cols, pos_cols). drug_cols auto-detected = header minus SeqID minus P<digits>."""
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
    rows = load_rows(path)
    pos_cols = _position_columns(header)
    drug_cols = [h for h in header if h != "SeqID" and h not in pos_cols and h.strip()]
    return rows, drug_cols, pos_cols


def build_presence(rows, pos_cols, min_muts=MIN_MUTS):
    """Binary mutation-presence design matrix. Feature `<pos><AA>` = 1 iff isolate carries AA at pos
    ('-'/'.' = consensus). Keep features present in >= min_muts isolates. Returns (X, feat_names, counts)."""
    counts: dict[str, int] = {}
    per_row: list[set[str]] = []
    for row in rows:
        feats = set()
        for c in pos_cols:
            cell = (row.get(c) or "").strip()
            if cell in ("", "-", ".", "NA"):
                continue
            pos = c[1:]
            for aa in cell:
                if aa.isalpha():
                    feats.add(f"{pos}{aa}")
        per_row.append(feats)
        for feat in feats:
            counts[feat] = counts.get(feat, 0) + 1
    feat_names = sorted(f for f, n in counts.items() if n >= min_muts)
    fidx = {f: j for j, f in enumerate(feat_names)}
    X = np.zeros((len(rows), len(feat_names)), float)
    for i, feats in enumerate(per_row):
        for feat in feats:
            j = fidx.get(feat)
            if j is not None:
                X[i, j] = 1.0
    return X, feat_names, {f: counts[f] for f in feat_names}


def build_pairs(X, feat_names, counts, top_k=TOP_K_FOR_PAIRS, min_cooc=MIN_COOC, max_pairs=MAX_PAIRS):
    """Constrained interaction pairs among the top-K most-prevalent features, co-occurring >= min_cooc.
    Returns (pairs[list of (i,j)], pair_names). Ranked by co-occurrence; capped at max_pairs."""
    top = sorted(range(len(feat_names)), key=lambda j: -counts[feat_names[j]])[:top_k]
    cand = []
    for a_i in range(len(top)):
        for b_i in range(a_i + 1, len(top)):
            i, j = top[a_i], top[b_i]
            cooc = int(np.sum((X[:, i] > 0) & (X[:, j] > 0)))
            if cooc >= min_cooc:
                cand.append((cooc, i, j))
    cand.sort(key=lambda t: -t[0])
    cand = cand[:max_pairs]
    pairs = [(i, j) for _, i, j in cand]
    pair_names = [f"{feat_names[i]}:{feat_names[j]}" for _, i, j in cand]
    return pairs, pair_names


def augment(X, pairs):
    if not pairs:
        return X
    extra = np.zeros((X.shape[0], len(pairs)), float)
    for k, (i, j) in enumerate(pairs):
        extra[:, k] = X[:, i] * X[:, j]
    return np.hstack([X, extra])


def _model():
    return make_pipeline(
        StandardScaler(),
        ElasticNetCV(l1_ratio=_L1_GRID, n_alphas=20, cv=3, max_iter=5000, tol=1e-3, random_state=SEED),
    )


def nested_oof(X, y, seed=SEED):
    """Nested-CV out-of-fold predictions: ElasticNetCV (inner CV picks lambda) inside a 5-fold outer."""
    kf = KFold(n_splits=CV_OUTER, shuffle=True, random_state=seed)
    return cross_val_predict(_model(), X, y, cv=kf)


def paired_bootstrap(y, add_oof, int_oof, boot=BOOT, seed=SEED):
    """Paired bootstrap on OOF predictions: delta_rho = rho(interaction) - rho(additive), 95% CI."""
    rng = np.random.default_rng(seed)
    n = len(y)
    add_rho = _spearman(add_oof, y)
    int_rho = _spearman(int_oof, y)
    deltas = np.empty(boot, float)
    for b in range(boot):
        idx = rng.integers(0, n, n)
        deltas[b] = _spearman(int_oof[idx], y[idx]) - _spearman(add_oof[idx], y[idx])
    deltas = deltas[~np.isnan(deltas)]
    lo, hi = (np.percentile(deltas, [2.5, 97.5]) if len(deltas) else (float("nan"), float("nan")))
    return {"add_rho": round(add_rho, 4), "int_rho": round(int_rho, 4),
            "delta_rho": round(int_rho - add_rho, 4),
            "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
            "ci_positive": bool(lo > 0)}


def top_interactions(X_int, y, pair_names, n_single, top=6):
    """Fit the interaction model on FULL data; report the largest-magnitude INTERACTION coefficients
    (verify-in-batch: are they real biology or overfit noise?)."""
    mdl = _model().fit(X_int, y)
    coef = mdl[-1].coef_
    inter = [(pair_names[k], float(coef[n_single + k])) for k in range(len(pair_names))]
    inter = [t for t in inter if abs(t[1]) > 1e-8]
    inter.sort(key=lambda t: -abs(t[1]))
    return [{"pair": p, "coef": round(c, 4)} for p, c in inter[:top]]


def run_drug(rows, drug_col, X_all, feat_names, pairs, pair_names, seed=SEED):
    keep = [i for i, r in enumerate(rows) if _parse_fold(r.get(drug_col, "")) not in (None,)
            and (_parse_fold(r.get(drug_col, "")) or 0) > 0]
    y = np.array([np.log10(_parse_fold(rows[i][drug_col])) for i in keep], float)
    n = len(keep)
    if n < N_MIN or float(np.var(y)) < 1e-6:
        return {"drug_col": drug_col, "n": n, "powered": False,
                "note": "too few isolates" if n < N_MIN else "degenerate fold spread"}
    idx = np.array(keep)
    Xa = X_all[idx]
    Xi = augment(Xa, pairs)
    add_oof = nested_oof(Xa, y, seed)
    int_oof = nested_oof(Xi, y, seed)
    boot = paired_bootstrap(y, add_oof, int_oof, seed=seed)
    return {"drug_col": drug_col, "n": n, "powered": True,
            "add_r2": round(_r2(y, add_oof), 4), "int_r2": round(_r2(y, int_oof), 4),
            "n_interaction_terms": len(pairs), **boot,
            "top_interactions": top_interactions(Xi, y, pair_names, Xa.shape[1])}


def run_class(path: Path, class_name: str, max_drugs=None, seed=SEED):
    rows, drug_cols, pos_cols = load_class(path)
    X_all, feat_names, counts = build_presence(rows, pos_cols)
    pairs, pair_names = build_pairs(X_all, feat_names, counts)
    if max_drugs:
        drug_cols = drug_cols[:max_drugs]
    per_drug = {}
    for dc in drug_cols:
        per_drug[dc] = run_drug(rows, dc, X_all, feat_names, pairs, pair_names, seed)
    return {"class": class_name, "n_isolates": len(rows), "n_single_features": len(feat_names),
            "n_interaction_candidates": len(pairs), "per_drug": per_drug}


def run_all(data_dir: Path, classes=None, max_drugs=None, seed=SEED):
    classes = classes or list(_CLASSES)
    per_class = {}
    for cname in classes:
        p = data_dir / _CLASSES[cname]
        if not p.exists():
            per_class[cname] = {"class": cname, "note": f"dataset absent at {p}"}
            continue
        per_class[cname] = run_class(p, cname, max_drugs=max_drugs, seed=seed)

    # verdict: across all powered drug-cells, fraction CI-positive
    powered, ci_pos = [], []
    for c in per_class.values():
        for dc, m in c.get("per_drug", {}).items():
            if m.get("powered"):
                powered.append(f"{c['class']}:{dc}")
                if m.get("ci_positive"):
                    ci_pos.append(f"{c['class']}:{dc}")
    n_pow, n_pos = len(powered), len(ci_pos)
    frac = (n_pos / n_pow) if n_pow else 0.0
    verdict = ("PASS_EPISTASIS_BEATS_ADDITIVE" if (n_pow > 0 and frac >= PASS_FRACTION)
               else ("FAIL_ADDITIVE_SUFFICES" if n_pow > 0 else "NO_POWERED_CELLS"))
    return {
        "artifact": "hiv_epistasis_vs_additive",
        "schema": "hiv-epistasis-v1",
        "question": "Does a constrained pairwise-interaction model beat the additive mutation-presence "
                    "baseline out-of-sample on Stanford HIVDB PhenoSense fold-change?",
        "prereg": {"MIN_MUTS": MIN_MUTS, "MIN_COOC": MIN_COOC, "TOP_K_FOR_PAIRS": TOP_K_FOR_PAIRS,
                   "MAX_PAIRS": MAX_PAIRS, "N_MIN": N_MIN, "CV_OUTER": CV_OUTER, "BOOT": BOOT,
                   "PASS_FRACTION": PASS_FRACTION, "seed": seed,
                   "beat_rule": "paired bootstrap delta_rho (interaction-additive) 95% CI lower bound > 0",
                   "model": "ElasticNetCV(l1_ratio grid, inner cv=3) in StandardScaler pipeline; nested 5-fold OOF"},
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent wet-lab; NOT Sierra interpretation)",
        "verdict": verdict,
        "n_powered_cells": n_pow, "n_ci_positive": n_pos, "fraction_ci_positive": round(frac, 3),
        "ci_positive_cells": ci_pos, "powered_cells": powered,
        "honest_caveats": [
            "A mechanism-FEATURE interaction model, NOT a sequence embedding (the closed 0-for-5 arm).",
            "Censored folds ('>'/'<') kept at the numeric bound: rank-metric-safe, slightly biases R2.",
            "Paired bootstrap on FIXED OOF predictions is the standard held-out paired test (mildly optimistic; disclosed).",
            "Interactions constrained to top-K prevalent co-occurring (>=MIN_COOC) pairs + L1 — not the full O(p^2) space.",
            "A FAIL ('additive suffices for HIV DR') is a valid world-model finding, not a bug.",
        ],
        "citation": "Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset per HIVDB Terms of Use",
        "per_class": per_class,
    }


def render_md(res: dict, generated: str) -> str:
    L = [f"# HIV epistasis vs additive — does interaction structure beat the additive rule? ({generated})", "",
         f"**Verdict: {res['verdict']}** — {res['n_ci_positive']}/{res['n_powered_cells']} powered drug-cells "
         f"show a bootstrap-CI-positive interaction gain (fraction {res['fraction_ci_positive']}; "
         f"PASS bar = {res['prereg']['PASS_FRACTION']}).", "",
         f"Label = {res['label_source']}. {res['question']}", "",
         "Both models = ElasticNetCV in a nested 5-fold OOF harness on the SAME folds (fair paired comparison). "
         "`beat` = paired bootstrap delta_rho (interaction-additive) 95%-CI lower bound > 0.", "",
         "| class:drug | n | add rho | int rho | **d_rho [95% CI]** | CI+ | add/int R2 |",
         "|---|---|---|---|---|---|---|"]
    for c in res["per_class"].values():
        for dc, m in c.get("per_drug", {}).items():
            if not m.get("powered"):
                L.append(f"| {c['class']}:{dc} | {m.get('n')} | — | — | {m.get('note','unpowered')} | — | — |")
                continue
            L.append(f"| {c['class']}:{dc} | {m['n']} | {m['add_rho']} | {m['int_rho']} | "
                     f"**{m['delta_rho']} [{m['ci_lo']}, {m['ci_hi']}]** | {'YES' if m['ci_positive'] else 'no'} | "
                     f"{m['add_r2']}/{m['int_r2']} |")
    L += ["", "## Verify-in-batch — top interaction terms (biology sanity check)"]
    for c in res["per_class"].values():
        for dc, m in c.get("per_drug", {}).items():
            if m.get("powered") and m.get("top_interactions"):
                terms = ", ".join(f"{t['pair']}({t['coef']:+})" for t in m["top_interactions"][:4])
                L.append(f"- **{c['class']}:{dc}** — {terms}")
    L += ["", "## Interpretation (how to read a FAIL)",
          "A FAIL means the ADDITIVE mutation-presence model is out-of-sample-competitive with the explicit "
          "interaction model for RANK-prediction of fold-change — NOT that epistasis is absent. Read the "
          "top-interaction terms above: if they are known synergy pairs recovered with LARGE coefficients "
          "(e.g. INI G140S+Q148H `140S:148H`, NRTI TAMs `41L:215Y`/`215Y:210W`, PI `82A:84V`/`46I:54V`), the "
          "epistasis is REAL and correctly localized — it just does not improve the population rank-metric, "
          "because the additive main effects already rank double-mutants high (Spearman is rank-based) and the "
          "genuine primary-weak+accessory-weak synergy cases are rare in the population. The world-model "
          "implication: the curated ADDITIVE catalog wins in the quantitative/continuous regime too, extending "
          "the deterministic-decoder thesis beyond binary R/S. Interactions can even HURT on small cells "
          "(variance/overfit) — see any negative delta_rho."]
    L += ["", "## Honest caveats"] + [f"- {x}" for x in res["honest_caveats"]]
    L += ["", f"Citation: {res['citation']}."]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=REPO / "data" / "raw" / "hiv")
    ap.add_argument("--classes", default=None, help="comma list subset of NNRTI,NRTI,PI,INI")
    ap.add_argument("--max-drugs", type=int, default=None, help="smoke: cap drugs/class")
    ap.add_argument("--boot", type=int, default=None, help="override BOOT (smoke)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if a.boot:
        globals()["BOOT"] = a.boot
    classes = [c.strip().upper() for c in a.classes.split(",")] if a.classes else None
    if not any((a.data_dir / f).exists() for f in _CLASSES.values()):
        print(f"ERROR: no HIV DataSets under {a.data_dir} (gitignored; see hiv_nnrti_validate for download)",
              file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run_all(a.data_dir, classes=classes, max_drugs=a.max_drugs)
    out = a.out or (REPO / "wiki" / f"hiv_epistasis_result_{today}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
