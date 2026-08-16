"""E-Flux bridge: does expression-constrained FBA improve conditional essentiality?

The bar is frozen in `wiki/fba_eflux_bridge_prereg_2026-08-16.md` BEFORE any solve. Both arms run on
the SAME conditions -- the shipped 25-condition 0.7368 is not a valid comparator for an 11-condition
panel, and comparing against it would repeat the regulatory-arm verdict bug.

Arm A: plain FBA (apply_carbon_condition only).
Arm B: identical + E-Flux bounds from PRECISE-1K expression (Colijn et al. 2009).

Usage:
    uv run python scripts/fba_eflux_bridge.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    conditionally_essential_genes,
    confusion_from_calls,
    mcc,
)
from dna_decode.fba.fitness_browser import (  # noqa: E402
    ESSENTIAL_FITNESS,
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

FRAC = 0.01
P1K_DIR = Path("D:/dna_decode_cache/precise1k")
EXPR_PCTL = 99.0


# ---------------------------------------------------------------- substrate name normalization
def _norm(s: str) -> str:
    s = str(s).lower().replace("_", " ")
    for j in ("monohydrate", "hydrate", "disodium salt", "sodium salt", "dibasic hexahydrate",
              "hexahydrate", "hydrochloride", "potassium ", "sodium ", "acid", "ic ", " salt"):
        s = s.replace(j, " ")
    return re.sub(r"[^a-z0-9]", "", s)


_SYN = {
    "dfructose": "fructose", "fructose": "fructose", "dgalactose": "galactose",
    "galactose": "galactose", "dglucon": "gluconate", "gluconate": "gluconate",
    "dglucosamine": "glucosamine", "dglucose": "glucose", "glucose": "glucose",
    "dglucose6phosphate": "g6p", "dmaltose": "maltose", "maltose": "maltose",
    "dmannitol": "mannitol", "dmannose": "mannose", "dribose": "ribose", "ribose": "ribose",
    "dserine": "dserine", "dsorbitol": "sorbitol", "sorbitol": "sorbitol",
    "dxylose": "xylose", "xylose": "xylose", "glycerol": "glycerol", "glycol": "glycolate",
    "lfucose": "fucose", "lmal": "lmalate", "nacetyldglucosamine": "nag",
    "nacetylglucosamine": "nag", "acetate": "acetate", "dlactate": "dlactate",
    "pyruvate": "pyruvate", "succinate": "succinate", "aketoglutar": "akg",
}


def canon(s: str) -> str:
    return _SYN.get(_norm(s), _norm(s))


# ---------------------------------------------------------------- GPR evaluation (AND=min, OR=sum)
_TOKEN = re.compile(r"\(|\)|\band\b|\bor\b|[A-Za-z0-9_.\-]+", re.I)


def eval_gpr(rule: str, expr: dict[str, float]) -> float | None:
    """Evaluate a GPR to an expression score. None if any referenced gene is unmeasured.

    Returning None (rather than 0.0) is load-bearing: an unmeasured gene must leave its reaction
    UNCONSTRAINED, never silently knocked out.
    """
    toks = [t for t in _TOKEN.findall(rule or "") if t.strip()]
    if not toks:
        return None
    pos = 0

    def peek() -> str | None:
        return toks[pos] if pos < len(toks) else None

    def parse_or() -> float | None:
        vals = [parse_and()]
        while peek() and peek().lower() == "or":
            nonlocal pos
            pos += 1
            vals.append(parse_and())
        if any(v is None for v in vals):
            return None
        return float(sum(vals))

    def parse_and() -> float | None:
        vals = [parse_atom()]
        while peek() and peek().lower() == "and":
            nonlocal pos
            pos += 1
            vals.append(parse_atom())
        if any(v is None for v in vals):
            return None
        return float(min(vals))

    def parse_atom() -> float | None:
        nonlocal pos
        t = peek()
        if t is None:
            return None
        if t == "(":
            pos += 1
            v = parse_or()
            if peek() == ")":
                pos += 1
            return v
        pos += 1
        return expr.get(t)

    try:
        return parse_or()
    except RecursionError:
        return None


def apply_eflux(model, expr: dict[str, float]) -> dict[str, int]:
    """Scale each GPR-carrying, non-exchange reaction's bounds by normalized expression. IN PLACE."""
    scores: dict[str, float] = {}
    exchanges = {r.id for r in model.exchanges}
    for rxn in model.reactions:
        if rxn.id in exchanges or not rxn.gene_reaction_rule.strip():
            continue
        s = eval_gpr(rxn.gene_reaction_rule, expr)
        if s is not None:
            scores[rxn.id] = s
    if not scores:
        return {"n_constrained": 0, "n_skipped_unmeasured": 0}
    ref = float(np.percentile(list(scores.values()), EXPR_PCTL)) or 1.0
    n = 0
    for rid, s in scores.items():
        rxn = model.reactions.get_by_id(rid)
        scale = min(1.0, s / ref)
        ub, lb = rxn.upper_bound, rxn.lower_bound
        if ub > 0:
            rxn.upper_bound = ub * scale
        if lb < 0:
            rxn.lower_bound = lb * scale
        n += 1
    n_gpr = sum(1 for r in model.reactions
                if r.id not in exchanges and r.gene_reaction_rule.strip())
    return {"n_constrained": n, "n_skipped_unmeasured": n_gpr - n}


# ---------------------------------------------------------------- expression per condition
def build_condition_expression(conds: dict[str, str]) -> tuple[dict[str, dict[str, float]], dict]:
    meta = pd.read_csv(P1K_DIR / "p1k_meta.csv", low_memory=False)
    col = [c for c in meta.columns if "Carbon" in c][0]
    tpm = pd.read_csv(P1K_DIR / "log_tpm_qc.csv", index_col=0)

    # The join key is whichever column actually overlaps the expression matrix columns -- NOT the
    # column named `sample_id`, which holds a descriptive label (`control__wt_glc__1`), not the id.
    cols = set(tpm.columns)
    sid = max(meta.columns, key=lambda c: len(set(meta[c].astype(str)) & cols))
    n_ov = len(set(meta[sid].astype(str)) & cols)
    if n_ov == 0:
        raise SystemExit("no metadata column joins to the expression matrix")
    print(f"join key: {sid!r} ({n_ov} samples matched)")

    meta = meta[meta[sid].isin(tpm.columns)]
    meta["_canon"] = meta[col].astype(str).str.replace(r"\(.*?\)", "", regex=True).map(canon)

    out, prov = {}, {}
    for name in conds:
        want = canon(name)
        samples = meta.loc[meta["_canon"] == want, sid].tolist()
        if not samples:
            continue
        sub = tpm[samples]
        # log2(TPM) -> linear TPM, then mean across replicates
        lin = (2.0 ** sub).mean(axis=1)
        out[name] = {g: float(v) for g, v in lin.items()}
        prov[name] = {"n_samples": len(samples)}
    return out, prov


# ---------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_eflux_bridge_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    conn = open_db()
    conds_all = carbon_conditions(conn, model)
    expr_by_cond, prov = build_condition_expression(conds_all)
    keys = sorted(expr_by_cond)
    print(f"conditions with PRECISE-1K expression: {len(keys)}/{len(conds_all)}")
    for k in keys:
        print(f"   {k:45} n_samples={prov[k]['n_samples']}")
    if len(keys) < 2:
        print("too few matched conditions", file=sys.stderr)
        return 2

    conds = {k: conds_all[k] for k in keys}
    all_ex = tuple(conds_all.values())
    model_genes = {g.id for g in model.genes}

    records = load_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    print(f"\nconditionally essential on this {len(keys)}-condition panel: {len(subset)}")
    if not subset:
        print("no two-sided genes on this panel", file=sys.stderr)
        return 2
    genes = [r.gene_id for r in subset]

    arms: dict[str, dict] = {}
    for arm in ("baseline", "eflux"):
        calls: dict[str, dict[str, bool]] = {}
        wt_by, eflux_stats = {}, {}
        print(f"\n=== ARM {arm}: {len(genes)} genes x {len(keys)} conditions ===", flush=True)
        for n, cond in enumerate(keys, 1):
            with model:
                apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
                if arm == "eflux":
                    eflux_stats[cond] = apply_eflux(model, expr_by_cond[cond])
                wt = wildtype_growth(model)
                wt_by[cond] = round(float(wt), 5)
                d: dict[str, bool] = {}
                if wt > 1e-9:
                    res = single_gene_deletion(
                        model, gene_list=[model.genes.get_by_id(g) for g in genes])
                    # Identical to the shipped validate script, NaN coding included: an infeasible
                    # solve is GENUINE essentiality (ATPM floor unmet), not a solver failure.
                    for _, row in res.iterrows():
                        gid = next(iter(row["ids"]))
                        g = row["growth"]
                        d[gid] = (g != g) or (g < FRAC * wt)
                else:
                    d = {g: True for g in genes}
                calls[cond] = d
            print(f"  [{n}/{len(keys)}] {cond:42} wt={wt_by[cond]:.4f}", flush=True)

        cells_right = cells_total = 0
        per_cond = {}
        for cond in keys:
            e = {r.gene_id: r.experimental[cond] for r in subset}
            p = {g: calls[cond].get(g, False) for g in genes}
            cm = confusion_from_calls(e, p)
            right = cm["tp"] + cm["tn"]
            cells_right += right
            cells_total += cm["n"]
            per_cond[cond] = {"agreement": round(right / cm["n"], 4) if cm["n"] else None,
                              "mcc": round(mcc(cm), 4), **cm}
        exact = sum(1 for r in subset
                    if {c for c in keys if calls[c].get(r.gene_id, False)}
                    == {c for c in keys if r.experimental[c]})
        arms[arm] = {
            "per_cell_agreement": round(cells_right / cells_total, 4) if cells_total else None,
            "n_cells": cells_total,
            "exact_set_match": exact,
            "exact_set_match_rate": round(exact / len(subset), 4),
            "wildtype_growth": wt_by,
            "per_condition": per_cond,
            **({"eflux": eflux_stats} if arm == "eflux" else {}),
        }
        print(f"  -> per-cell {arms[arm]['per_cell_agreement']} | exact-set {exact}/{len(subset)}")

    b, e = arms["baseline"]["per_cell_agreement"], arms["eflux"]["per_cell_agreement"]
    delta = round(e - b, 4)
    verdict = ("EFLUX_IMPROVES" if delta >= 0.02
               else "MACHADO_PRIOR_CONFIRMED_ON_ESSENTIALITY" if delta <= 0
               else "AMBIGUOUS_BELOW_PREREGISTERED_BAR")

    out = {
        "record": "fba-eflux-bridge-v1",
        "date": date.today().isoformat(),
        "prereg": "wiki/fba_eflux_bridge_prereg_2026-08-16.md",
        "model": "iML1515",
        "expression": "PRECISE-1K log_tpm_qc.csv (SBRG/precise1k)",
        "labels": f"Fitness Browser RB-TnSeq orgId=Keio, fit<{a.threshold}",
        "n_conditions": len(keys),
        "conditions": {k: prov[k] for k in keys},
        "n_conditionally_essential": len(subset),
        "arms": arms,
        "delta_per_cell": delta,
        "verdict": verdict,
        "prereg_bar": {"go": ">= +0.02", "machado_confirmed": "<= 0"},
        "caveats": [
            "PRECISE-1K is K-12 MG1655; Fitness Browser labels are orgId=Keio (BW25113 parent).",
            "Sample counts are extremely imbalanced (glucose 621, glycerol 111, others 2-8).",
            "This panel is 11 of the shipped 25 conditions; a null here does not close the bridge "
            "on the full panel.",
            "The shipped 0.7368 per-cell figure is a 25-condition number and is NOT the comparator; "
            "the baseline arm here is recomputed on the same conditions.",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nbaseline {b} | eflux {e} | delta {delta:+.4f}\nVERDICT: {verdict}")
    print(f"wrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
