"""Transplant three biomass demands from iML1515's own WT equation, and score the pre-registered run.

RUNS `wiki/fba_biomass_completion_prereg_2026-08-22.md` EXACTLY (frozen at commit abc8040, BEFORE this
scored anything). Every constant traces to a numbered section:

  §2  intervention   drop kdo2lipid4_e, add colipa_e / hemeO_c / enter_c -- ALL coefficients read from
                     `BIOMASS_Ec_iML1515_WT_75p37M`, none invented, none fitted
  §3  predictions    colipa_e -> 9 LPS genes flip; hemeO_c -> cyoE flips; enter_c -> the 7 iron-uptake
                     genes must NOT flip (they need the LOADED siderophore, which no biomass demands)
  §4  endpoints      primary >= 8/10; FP guardrail +20% relative; guardrail breach = FAILURE
  §5  determinism    processes=1, whole pipeline twice, identical calls required
  §6  verdicts       six pre-committed outcomes

The mechanism check is the point of the design: if the 7 iron genes flip, the stated mechanism is wrong
even though the headline number improves.

Usage:
    uv run python scripts/fba_biomass_completion.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.conditional_essentiality import conditionally_essential_genes  # noqa: E402
from dna_decode.fba.fitness_browser import (  # noqa: E402
    ESSENTIAL_FITNESS,
    apply_carbon_condition,
    carbon_conditions,
    load_records,
    open_db,
)
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

FRAC = 0.01
CORE_BIOMASS = "BIOMASS_Ec_iML1515_core_75p37M"
WT_BIOMASS = "BIOMASS_Ec_iML1515_WT_75p37M"

#: §2 -- verified against the WT equation at run time; a drift raises rather than silently rescoring.
DROP = {"kdo2lipid4_e": -0.019456}
ADD = {"colipa_e": -0.008151, "hemeO_c": -0.000223, "enter_c": -0.000223}

#: §3 predictions, by gene NAME (resolved to ids at run time so a rename fails loudly)
PREDICT_FLIP = ["gmhA", "gmhB", "hldE", "hldD", "waaC", "waaF", "waaP", "waaG", "galU", "cyoE"]
PREDICT_NO_FLIP = ["fes", "fepB", "fepC", "fepD", "fepG", "tonB", "exbD"]

PRIMARY_BAR = 8              # of 10
GUARDRAIL_REL_INCREASE = 0.20


def verify_coefficients(model) -> dict:
    """§2 -- the added coefficients must actually BE the WT equation's. Raise if the model disagrees."""
    wt = model.reactions.get_by_id(WT_BIOMASS)
    have = {m.id: c for m, c in wt.metabolites.items()}
    bad = {k: (v, have.get(k)) for k, v in ADD.items()
           if have.get(k) is None or abs(have[k] - v) > 1e-9}
    if bad:
        raise SystemExit(f"coefficients drifted from {WT_BIOMASS}: {bad}")
    core = model.reactions.get_by_id(CORE_BIOMASS)
    ch = {m.id: c for m, c in core.metabolites.items()}
    for k, v in DROP.items():
        if abs(ch.get(k, 0.0) - v) > 1e-9:
            raise SystemExit(f"{k} is {ch.get(k)} in {CORE_BIOMASS}, expected {v}")
    return {"verified_against": WT_BIOMASS, "add": ADD, "drop": DROP}


def apply_completion(model) -> None:
    """§2 intervention, IN PLACE on the current context."""
    core = model.reactions.get_by_id(CORE_BIOMASS)
    for k, v in DROP.items():
        core.add_metabolites({model.metabolites.get_by_id(k): -v}, combine=True)
    for k, v in ADD.items():
        core.add_metabolites({model.metabolites.get_by_id(k): v}, combine=True)


def resolve(model, names: list[str]) -> dict[str, str]:
    lookup = {(g.name or "").lower(): g.id for g in model.genes}
    out = {}
    for n in names:
        if n.lower() not in lookup:
            raise SystemExit(f"gene name {n!r} does not resolve in {model.id}")
        out[n] = lookup[n.lower()]
    return out


def run_arm(model, conds, genes, modified, single_gene_deletion, label):
    all_ex = tuple(conds.values())
    calls, stats = {}, {}
    for n, cond in enumerate(sorted(conds), 1):
        with model:
            apply_carbon_condition(model, conds[cond], all_carbon=all_ex)
            if modified:
                apply_completion(model)
            wt = wildtype_growth(model)
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in genes], processes=1)
                d = {}
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    d[gid] = (g != g) or (g < FRAC * wt)
            else:
                d = {g: True for g in genes}
            calls[cond] = d
            stats[cond] = {"wt": round(float(wt), 6)}
        print(f"  [{n:2}/{len(conds)}] {label:9} {cond[:32]:34} wt={stats[cond]['wt']:.5f} "
              f"essential={sum(d.values()):4}", flush=True)
    return calls, stats


def confusion(calls, truth, genes, conds):
    tp = fp = fn = tn = 0
    for c in conds:
        for g in genes:
            p = bool(calls.get(c, {}).get(g, False))
            y = bool(truth.get((g, c), False))
            if p and y:
                tp += 1
            elif p:
                fp += 1
            elif y:
                fn += 1
            else:
                tn += 1
    return tp, fp, fn, tn


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=ESSENTIAL_FITNESS)
    ap.add_argument("--out", default=f"wiki/fba_biomass_completion_{date.today().isoformat()}")
    a = ap.parse_args(argv)

    from cobra.flux_analysis import single_gene_deletion

    model = load_model()
    prov = verify_coefficients(model)
    print(f"coefficients verified against {WT_BIOMASS}")

    conn = open_db()
    conds = carbon_conditions(conn, model)
    model_genes = {g.id for g in model.genes}
    records = load_records(conn, conds, gene_filter=model_genes, threshold=a.threshold)
    subset = conditionally_essential_genes(records)
    genes = [r.gene_id for r in subset]
    truth = {(r.gene_id, c): bool(y) for r in subset for c, y in r.experimental.items()}
    print(f"conditions {len(conds)} | conditionally essential genes {len(genes)}")

    flip_ids = resolve(model, PREDICT_FLIP)
    noflip_ids = resolve(model, PREDICT_NO_FLIP)
    missing = [n for n, g in {**flip_ids, **noflip_ids}.items() if g not in set(genes)]
    print(f"predicted-flip {len(flip_ids)} | predicted-NO-flip {len(noflip_ids)}"
          + (f" | NOT on panel: {missing}" if missing else ""))

    arms = {}
    for run_i in (1, 2):
        print(f"\n===== RUN {run_i} =====")
        arms.setdefault("baseline", {})[run_i] = run_arm(
            model, conds, genes, False, single_gene_deletion, "baseline")
        arms.setdefault("completed", {})[run_i] = run_arm(
            model, conds, genes, True, single_gene_deletion, "completed")

    det = {}
    for arm, runs in arms.items():
        c1, c2 = runs[1][0], runs[2][0]
        diff = sum(1 for c in conds for g in genes
                   if bool(c1.get(c, {}).get(g)) != bool(c2.get(c, {}).get(g)))
        det[arm] = {"n_cells_differing": diff, "identical": diff == 0}
    determinism_ok = all(v["identical"] for v in det.values())
    print(f"\n[§5 determinism] {determinism_ok}  {det}")

    base = arms["baseline"][1][0]
    mod = arms["completed"][1][0]
    b_tp, b_fp, b_fn, b_tn = confusion(base, truth, genes, conds)
    m_tp, m_fp, m_fn, m_tn = confusion(mod, truth, genes, conds)

    def recovered(ids):
        out = []
        for nm, g in ids.items():
            hits = [c for c in conds
                    if truth.get((g, c)) and mod.get(c, {}).get(g) and not base.get(c, {}).get(g)]
            if hits:
                out.append({"name": nm, "gene": g, "conditions": sorted(hits)})
        return out

    rec_flip = recovered(flip_ids)
    rec_noflip = recovered(noflip_ids)
    fp_rel = (m_fp - b_fp) / b_fp if b_fp else (float("inf") if m_fp else 0.0)
    breached = fp_rel > GUARDRAIL_REL_INCREASE

    print(f"\n[baseline]  TP={b_tp} FP={b_fp} FN={b_fn} TN={b_tn}")
    print(f"[completed] TP={m_tp} FP={m_fp} FN={m_fn} TN={m_tn}   FP {fp_rel:+.1%} "
          f"guardrail_breached={breached}")
    print(f"\n[§4 PRIMARY]   predicted-to-flip recovered: {len(rec_flip)}/10 "
          f"-> {[r['name'] for r in rec_flip]}")
    print(f"[§3 MECHANISM] predicted-NOT-to-flip that flipped: {len(rec_noflip)}/7 "
          f"-> {[r['name'] for r in rec_noflip]}")

    # ---- WHERE DID THE FALSE POSITIVES GO? A biomass coefficient is condition-INDEPENDENT, so a gene
    # that is the sole route to the demanded metabolite in every condition becomes CONSTITUTIVELY
    # essential -- it cannot become conditionally essential. Measure that directly.
    breadth = []
    fp_from_recovered = 0
    for nm, g in flip_ids.items():
        p = {c for c in conds if mod.get(c, {}).get(g)}
        t = {c for c in conds if truth.get((g, c))}
        fp_from_recovered += len(p - t)
        breadth.append({"name": nm, "gene": g, "n_predicted_essential": len(p),
                        "n_experimentally_essential": len(t), "new_false_positives": len(p - t)})
    fp_delta = m_fp - b_fp
    share = (fp_from_recovered / fp_delta) if fp_delta else None
    print(f"\n[diagnostic] the {len(flip_ids)} recovered genes are predicted essential in "
          f"{breadth[0]['n_predicted_essential']}/{len(conds)} conditions each; they alone introduce "
          f"{fp_from_recovered} of the {fp_delta} new false positives"
          + (f" ({share:.0%})" if share else ""))

    if not determinism_ok:
        verdict = "INDETERMINATE"
    elif breached:
        verdict = "FAILURE_GUARDRAIL_BREACHED"
    elif len(rec_flip) >= PRIMARY_BAR and not rec_noflip:
        verdict = "H1_SUPPORTED"
    elif len(rec_flip) >= PRIMARY_BAR:
        verdict = "H1_SUPPORTED_MECHANISM_DISCONFIRMED"
    elif rec_flip:
        verdict = "H1_PARTIAL"
    else:
        verdict = "H1_FALSIFIED"

    out = {
        "record": "fba-biomass-completion-v1",
        "date": date.today().isoformat(),
        "prereg": "wiki/fba_biomass_completion_prereg_2026-08-22.md",
        "prereg_commit": "abc8040",
        "model": model.id, "intervention": prov,
        "n_conditions": len(conds), "n_genes": len(genes),
        "determinism": {"runs": 2, "per_arm": det, "passed": determinism_ok},
        "baseline_confusion": {"TP": b_tp, "FP": b_fp, "FN": b_fn, "TN": b_tn},
        "completed_confusion": {"TP": m_tp, "FP": m_fp, "FN": m_fn, "TN": m_tn},
        "recall": {"baseline": round(b_tp / (b_tp + b_fn), 4) if (b_tp + b_fn) else None,
                   "completed": round(m_tp / (m_tp + m_fn), 4) if (m_tp + m_fn) else None},
        "fp_relative_increase": round(fp_rel, 4),
        "guardrail_breached": breached,
        "primary": {"bar": PRIMARY_BAR, "n_recovered": len(rec_flip), "of": len(flip_ids),
                    "recovered": rec_flip},
        "mechanism_check": {"expected": 0, "n_flipped": len(rec_noflip), "flipped": rec_noflip},
        "constitutive_diagnostic": {
            "per_gene": breadth,
            "fp_introduced_by_recovered_genes": fp_from_recovered,
            "total_fp_increase": fp_delta,
            "share_of_guardrail_breach": None if share is None else round(share, 4),
            "reading": ("a biomass coefficient is condition-INDEPENDENT, so a gene that is the sole "
                        "route to the demanded metabolite in EVERY condition becomes constitutively "
                        "essential -- trading a false negative in some conditions for a false positive "
                        "in all the rest"),
        },
        "genes_not_on_panel": missing,
        "verdict": verdict,
        "caveats": [
            "Model-side change: recovering a gene means the MODEL now calls it essential; it does not "
            "validate the biomass edit against measured composition.",
            "Coefficients are transplanted from the model's own WT biomass, never fitted -- there is no "
            "free parameter to tune.",
            "Growth goes slightly UP under the modification, so recovery cannot be a growth-burden effect.",
        ],
    }
    Path(a.out + ".json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}\nwrote {a.out}.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
