"""Can the stress axis be tested with iML1515? A three-layer feasibility probe.

The carbon and nitrogen axes both work because the assay compound IS the medium component: opening its
exchange reconstructs the experimental condition. Stress is a different contract -- the compound is an
INHIBITOR added on top of a normal medium -- so mappability has to be tested, not assumed.

Three layers, each stricter than the last:

  L1  does an exchange exist at all?          (antibiotics/ionic liquids are not metabolites)
  L2  does adding it REDUCE growth?           (an exchange models supplementation, not toxicity)
  L3  is the molecular target in the model?   (the only remaining path: target-directed constraints)

Usage:
    uv run python scripts/fba_stress_feasibility_probe.py
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.fba.fitness_browser import ORG_ID, open_db  # noqa: E402
from dna_decode.fba.model import load_model, wildtype_growth  # noqa: E402

STRESS_EXP_GROUP = "stress"

#: Best-effort {stress compound -> candidate iML1515 exchange} for the panel members that are genuine
#: metabolites. Everything absent from this map is an antibiotic / ionic liquid / detergent / cytotoxic
#: with no metabolite representation whatsoever.
#:
#: NOTE: an earlier draft of this map paired "sodium fluoride" with EX_fe2_e. That is FERROUS IRON, not
#: fluoride -- a wrong mapping that would have manufactured a fake data point. It is dropped rather than
#: guessed at; iML1515 has no fluoride exchange.
CANDIDATE_EXCHANGES: dict[str, str] = {
    "Sodium acetate": "EX_ac_e",
    "Sodium Chloride": "EX_na1_e",
    "L-Lysine": "EX_lys__L_e",
    "Sodium nitrite": "EX_no2_e",
    "Dimethyl Sulfoxide": "EX_dmso_e",
    "Cobalt chloride hexahydrate": "EX_cobalt2_e",
    "Nickel (II) chloride hexahydrate": "EX_ni2_e",
    "copper (II) chloride dihydrate": "EX_cu2_e",
}

#: {compound -> (molecular target, candidate iML1515 gene/reaction ids)}. Empty id list = the target is
#: non-metabolic (ribosome / gyrase / membrane / cytoskeleton) and therefore outside a metabolic model
#: BY CONSTRUCTION, not by omission.
MOLECULAR_TARGETS: dict[str, tuple[str, list[str]]] = {
    "Phosphomycin disodium salt": ("MurA / UDP-GlcNAc enolpyruvyl transferase", ["murA", "UAGCVT"]),
    "D-Cycloserine": ("alanine racemase + D-Ala-D-Ala ligase",
                      ["alr", "ddlA", "ddlB", "ALARi", "ALAALAr"]),
    "Chloramphenicol": ("50S ribosome", []),
    "Tetracycline hydrochloride": ("30S ribosome", []),
    "Doxycycline hyclate": ("30S ribosome", []),
    "Nalidixic acid sodium salt": ("DNA gyrase", []),
    "Spectinomycin dihydrochloride pentahydrate": ("30S ribosome", []),
    "Fusidic acid sodium salt": ("EF-G translation elongation", []),
    "Bacitracin": ("undecaprenyl-PP carrier recycling (membrane)", []),
    "MreB Perturbing Compound A22": ("MreB cytoskeleton", []),
    "Cisplatin": ("DNA crosslinking", []),
    "Carbenicillin disodium salt": ("penicillin-binding proteins (transpeptidation)", []),
    "Cephalothin sodium salt": ("penicillin-binding proteins (transpeptidation)", []),
}

STRESS_TOLERANCE = 1e-6


def stress_conditions(conn) -> dict[str, int]:
    """{stress condition -> n experiments} for this organism."""
    return {c: n for c, n in conn.execute(
        "SELECT condition_1, COUNT(*) FROM Experiment WHERE orgId=? AND expGroup=? "
        "GROUP BY condition_1", (ORG_ID, STRESS_EXP_GROUP))}


def probe_growth_effect(model, exchange: str, baseline: float, uptake: float = 10.0) -> dict:
    """Does opening `exchange` REDUCE growth? A stress must; a nutrient does not.

    This is the load-bearing check. An exchange existing in the model says only that the compound is
    representable as a METABOLITE. Opening it models SUPPLEMENTATION -- the opposite of a stress.
    """
    with model:
        medium = dict(model.medium)
        medium[exchange] = uptake
        model.medium = medium
        g = float(wildtype_growth(model))
    delta = g - baseline
    return {"exchange": exchange, "growth": g, "delta": delta,
            "reduces_growth": delta < -STRESS_TOLERANCE}


def target_in_model(model, ids: list[str]) -> list[str]:
    """Which of `ids` resolve to a real reaction or gene in the model."""
    rxns = {r.id for r in model.reactions}
    gene_names = {(g.name or "").lower() for g in model.genes}
    gene_ids = {g.id for g in model.genes}
    found = []
    for i in ids:
        if i in rxns:
            found.append(f"rxn:{i}")
        elif i.lower() in gene_names or i in gene_ids:
            found.append(f"gene:{i}")
    return found


def main() -> int:
    model = load_model()
    conn = open_db()
    conds = stress_conditions(conn)
    baseline = float(wildtype_growth(model))

    # L1 -- exchange existence
    mappable = {k: v for k, v in CANDIDATE_EXCHANGES.items()
                if k in conds and v in {r.id for r in model.exchanges}}
    l1 = {"n_conditions": len(conds), "n_with_candidate_exchange": len(mappable),
          "fraction": round(len(mappable) / len(conds), 4) if conds else 0.0}

    # L2 -- does it actually behave as a stress?
    l2_rows = {k: probe_growth_effect(model, v, baseline) for k, v in mappable.items()}
    n_stress = sum(1 for r in l2_rows.values() if r["reduces_growth"])

    # L3 -- target-directed path
    l3_rows = {}
    for ab, (tgt, ids) in MOLECULAR_TARGETS.items():
        if ab not in conds:
            continue
        found = target_in_model(model, ids)
        l3_rows[ab] = {"target": tgt, "in_model": found, "metabolic_target": bool(found)}
    n_target = sum(1 for r in l3_rows.values() if r["metabolic_target"])

    out = {
        "record": "fba-stress-feasibility-v1",
        "date": date.today().isoformat(),
        "model": "iML1515",
        "labels": f"Fitness Browser RB-TnSeq orgId={ORG_ID}, expGroup='{STRESS_EXP_GROUP}'",
        "baseline_wildtype_growth": baseline,
        "L1_exchange_existence": l1,
        "L2_behaves_as_stress": {
            "n_probed": len(l2_rows), "n_reducing_growth": n_stress, "detail": l2_rows,
            "note": "an exchange models SUPPLEMENTATION; a stress must REDUCE growth",
        },
        "L3_target_directed": {
            "n_panel_antibiotics_checked": len(l3_rows), "n_with_metabolic_target_in_model": n_target,
            "detail": l3_rows,
            "note": "the only remaining path, and it needs a target-directed reaction constraint -- a "
                    "different experimental contract from the medium-swap used for carbon/nitrogen",
        },
        "verdict": ("STRESS_AXIS_NOT_REPRESENTABLE_IN_iML1515" if n_stress == 0
                    else "PARTIALLY_REPRESENTABLE"),
    }
    Path(f"wiki/fba_stress_feasibility_{date.today().isoformat()}.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")

    print(f"baseline wildtype growth: {baseline:.6f}")
    print(f"\nL1  stress conditions: {l1['n_conditions']} | with a candidate exchange: "
          f"{l1['n_with_candidate_exchange']} ({l1['fraction']:.0%})")
    print(f"L2  of those, REDUCING growth (i.e. actually a stress): {n_stress}/{len(l2_rows)}")
    for k, r in l2_rows.items():
        print(f"      {k[:33]:34} {r['exchange']:14} delta {r['delta']:+.6f}  "
              f"{'stress' if r['reduces_growth'] else 'NOT a stress'}")
    print(f"L3  panel antibiotics with a metabolic target in iML1515: {n_target}/{len(l3_rows)}")
    for k, r in l3_rows.items():
        if r["metabolic_target"]:
            print(f"      {k[:33]:34} {r['target'][:38]:40} {', '.join(r['in_model'])}")
    print(f"\nVERDICT: {out['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
