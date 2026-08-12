"""Does LABEL-BLIND gap-filling move the conditional-essentiality metric?

    uv run python scripts/fba_gapfill_conditional_test.py

The conditional cell (`wiki/fba_conditional_essentiality_2026-08-12.md`) found FBA reproduces the
medium-dependent essentiality switch for ~5% of the genes where it matters. The obvious next lever is
gap-filling: the model is incomplete, so add the missing biochemistry.

PRE-REGISTERED PREDICTION, written before running: **it should make things WORSE.** Conditional
essentiality requires the ABSENCE of an alternative route in one medium. Gap-filling ADDS routes, so it
should push genes toward "dispensable everywhere" -- the constant pattern the model already over-produces.

Three label-blind arms, none of which consults the essentiality labels:

  RANDOM    k reactions drawn at random from a donor reconstruction, several doses x several seeds.
            The dose-response: if adding biochemistry mattered at all, more of it should matter more.
  TARGETED  only donor reactions touching one of iML1515's DEAD-END metabolites -- what a structural
            gap-filler would actually propose.
  MAXIMAL   every donor reaction absent from iML1515. The upper bound on "add more biochemistry".

Donor is the Salmonella pan-reactome (iYS1720): enterobacterial, so a biologically plausible source, and
large enough to be a real test.

The run reports BOTH the binary metric and the threshold-free one, because the two answer different
questions -- and it reports the number of BINARY CALL FLIPS, which is the mechanism.

Exit 0 always: this is an experiment, not a gate.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.conditional_essentiality import (  # noqa: E402
    CONDITIONS,
    apply_condition,
    conditionally_essential_genes,
    confusion_from_calls,
    continuous_readout,
    load_labels,
    mcc,
    switch_accuracy,
)
from dna_decode.fba.gapfill import model_dead_ends  # noqa: E402
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

ESSENTIAL_FRAC = 0.01


def knockout_ratios(model, gene_ids: list[str]) -> tuple[dict, dict]:
    """{condition: {gene: mutant/wild-type growth}} plus the wild-type growth per condition."""
    from cobra.flux_analysis import single_gene_deletion  # noqa: PLC0415

    ratios, wts = {}, {}
    for c in sorted(CONDITIONS):
        with model:
            apply_condition(model, c)
            wt = wildtype_growth(model)
            wts[c] = round(wt, 4)
            out = {}
            if wt > 1e-9:
                res = single_gene_deletion(
                    model, gene_list=[model.genes.get_by_id(g) for g in gene_ids])
                for _, row in res.iterrows():
                    gid = next(iter(row["ids"]))
                    g = row["growth"]
                    out[gid] = 0.0 if g != g else g / wt
            ratios[c] = out
    return ratios, wts


def score_arm(records, ratios: dict) -> dict:
    calls = {c: {g: r <= ESSENTIAL_FRAC for g, r in d.items()} for c, d in ratios.items()}
    sw = switch_accuracy(records, calls)
    mccs = []
    for c in sorted(CONDITIONS):
        cm = confusion_from_calls({r.gene_id: r.experimental[c] for r in records}, calls[c])
        mccs.append(mcc(cm))
    cont = continuous_readout(records, ratios)
    return {
        "exact_set_match": sw["exact_set_match"],
        "n_conditionally_essential": sw["n_conditionally_essential"],
        "per_condition_agreement": sw["per_condition_agreement"],
        "mean_per_condition_mcc": round(sum(mccs) / len(mccs), 4),
        "auroc_threshold_free": cont["auroc"],
        "deployed_mcc": cont["deployed_mcc"],
        "calls": calls,
    }


def count_flips(base_calls: dict, arm_calls: dict) -> int:
    """Binary calls that changed. THE mechanism number: ratios can move without any call moving."""
    n = 0
    for c in base_calls:
        for g, v in base_calls[c].items():
            if g in arm_calls.get(c, {}) and arm_calls[c][g] != v:
                n += 1
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="ecoli")
    ap.add_argument("--donor", default="salmonella")
    ap.add_argument("--doses", default="25,100,400")
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--date", default=str(date.today()))
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    base = load_model(organism=a.organism)
    present = {g.id for g in base.genes}
    records = conditionally_essential_genes([r for r in load_labels() if r.gene_id in present])
    genes = [r.gene_id for r in records]
    print(f"{base.id}: {len(base.reactions)} reactions | conditionally-essential genes scored: {len(genes)}")

    base_ratios, base_wts = knockout_ratios(base, genes)
    baseline = score_arm(records, base_ratios)
    base_calls = baseline.pop("calls")
    print(f"BASELINE: exact-set {baseline['exact_set_match']}/{baseline['n_conditionally_essential']} | "
          f"per-cell {baseline['per_condition_agreement']} | deployed MCC {baseline['deployed_mcc']} | "
          f"AUROC {baseline['auroc_threshold_free']}")

    donor = load_model(organism=a.donor)
    have = {r.id for r in base.reactions}
    pool = [r for r in donor.reactions if r.id not in have]
    dead = {d.metabolite for d in model_dead_ends(base)}
    targeted = [r for r in pool if any(m.id in dead for m in r.metabolites)]
    print(f"donor {donor.id}: {len(pool)} reactions absent from {base.id}; "
          f"{len(targeted)} touch a dead-end metabolite")

    def run_arm(label, rxns):
        m = base.copy()
        m.add_reactions([r.copy() for r in rxns])
        ratios, wts = knockout_ratios(m, genes)
        s = score_arm(records, ratios)
        calls = s.pop("calls")
        s["n_reactions_added"] = len(rxns)
        s["binary_call_flips_vs_baseline"] = count_flips(base_calls, calls)
        s["n_ratios_changed_vs_baseline"] = sum(
            1 for c in ratios for g in ratios[c]
            if abs(ratios[c][g] - base_ratios[c].get(g, ratios[c][g])) > 1e-6)
        s["wildtype_growth"] = wts
        print(f"   {label:28s} +{len(rxns):4d} rxns | exact-set {s['exact_set_match']}/"
              f"{s['n_conditionally_essential']} | per-cell {s['per_condition_agreement']} | "
              f"deployed MCC {s['deployed_mcc']} | CALL FLIPS {s['binary_call_flips_vs_baseline']} | "
              f"ratios moved {s['n_ratios_changed_vs_baseline']}")
        return s

    arms = {}
    print("\nRANDOM arm (label-blind dose-response):")
    for k in [int(x) for x in a.doses.split(",") if x.strip()]:
        for seed in [int(x) for x in a.seeds.split(",") if x.strip()]:
            pick = random.Random(seed).sample(pool, min(k, len(pool)))
            arms[f"random_k{k}_seed{seed}"] = run_arm(f"random k={k} seed={seed}", pick)
    print("\nTARGETED arm (dead-end-closing reactions only):")
    arms["targeted_dead_end"] = run_arm("targeted", targeted)
    print("\nMAXIMAL arm (every donor reaction):")
    arms["maximal"] = run_arm("maximal", pool)

    any_flip = any(v["binary_call_flips_vs_baseline"] for v in arms.values())
    verdict = ("GAP_FILLING_DOES_NOT_MOVE_THE_CONDITIONAL_METRIC" if not any_flip
               else "GAP_FILLING_CHANGES_SOME_CALLS")
    print(f"\nVERDICT: {verdict}")

    result = {
        "record": "fba-gapfill-conditional-test-v1",
        "date": a.date,
        "model": base.id,
        "donor": donor.id,
        "prediction_registered_before_running": (
            "gap-filling should make the conditional metric WORSE, because it adds alternative routes and "
            "conditional essentiality requires their ABSENCE in one medium"),
        "baseline": baseline,
        "baseline_wildtype_growth": base_wts,
        "n_donor_reactions_available": len(pool),
        "n_donor_reactions_touching_a_dead_end": len(targeted),
        "arms": arms,
        "verdict": verdict,
        "caveats": [
            "Adding a donor reaction also imports its GPR, so the augmented model gains genes; the SCORED "
            "gene set is held fixed to the conditionally-essential genes present in the base model.",
            "AUROC carries run-to-run and invocation-shape variation of a few points (degenerate LP "
            "optima); small AUROC deltas between arms are NOT evidence of anything.",
            "One donor only. A different donor could in principle supply a reaction that matters, though "
            "the maximal arm already adds every reaction this donor has.",
            "This tests whether ADDING REACTIONS helps. It does not test regulatory/uptake constraints, "
            "which is the other candidate explanation for the flat ratios.",
        ],
    }
    outdir = Path(a.out_dir) if a.out_dir else Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"fba_gapfill_conditional_test_{a.date}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
