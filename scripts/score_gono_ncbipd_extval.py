"""External validation of the NON-FROZEN neisseria_amr cell on an INDEPENDENT NCBI Pathogen Detection
N. gonorrhoeae cohort — a large, provenance-disjoint (from the AR-Bank isolates) measured-AST cohort.

Why this matters: the AR-Bank Neisseria validation (2026-07-21) endorsed cipro + cefixime on N=20, but the
cefixime v0.1 rule was DERIVED + validated on that same 20-isolate cohort (an honest tension). This scores
the SAME cell on 170 independent NCBI-PD isolates (measured AST + NCBI's own pre-computed AMRFinderPlus
determinant calls) -> externally confirms (or refutes) the cell + the cefixime v0.1 narrowing at scale.

FEASIBILITY: no compute wall. NCBI-PD publishes its AMRFinderPlus calls in the `AMR_genotypes` field
(point mutations included, e.g. `gyrA_S91F=POINT`, `penA_I312M=POINT`), and this cohort's isolates have
downloadable assemblies -> NCBI ran AMRFinder for us. So this is a pure metadata JOIN + score, offline,
runnable in seconds from the two committed files (cohort.tsv + determinants.tsv).

HONEST CAVEATS:
- provenance-disjoint (different isolates than AR-Bank) but NOT methodology-independent (same AMRFinderPlus,
  same neisseria_amr cell) -- the standard caveat, same as the frozen cells.
- RAW sens/spec is CLONALITY-INFLATED (gonococci are clonally structured); lineage-collapse (Mash) is the
  follow-up (needs the assemblies). Read raw as a stress test, not a lineage-independent number.

  uv run python scripts/score_gono_ncbipd_extval.py
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules import neisseria_amr as NG  # noqa: E402
from scripts.amr_portal_score_independent import wilson_ci  # noqa: E402
from scripts.independent_cohort_validate import _conf  # noqa: E402

COHORT = "data/raw/gono_ncbipd_extval/cohort.tsv"
DETERMINANTS = "data/raw/gono_ncbipd_extval/determinants.tsv"
SPEC_FLOOR = 0.85
MIN_PER_CLASS = 5
# cell drugs -> call function (all scored from the SAME determinant set per isolate)
DRUGS = {
    "ciprofloxacin": NG.call_ng_ciprofloxacin,
    "cefixime": NG.call_ng_cefixime,
    "ceftriaxone": NG.call_ng_ceftriaxone,
    "azithromycin": NG.call_ng_azithromycin,
    "penicillin": NG.call_ng_penicillin,
    "tetracycline": NG.call_ng_tetracycline,
}


def main() -> int:
    labels = {r["biosample"]: r for r in csv.DictReader(open(COHORT), delimiter="\t")}
    dets = {r["biosample"]: [s for s in (r["determinants"] or "").split(";") if s]
            for r in csv.DictReader(open(DETERMINANTS), delimiter="\t")}

    results = {}
    print(f"Cohort: {len(labels)} NCBI-PD N. gonorrhoeae isolates (provenance-disjoint from AR-Bank)")
    for drug, fn in DRUGS.items():
        scored, records = [], []
        for bs, row in labels.items():
            rs = row.get(drug, "")
            if rs not in ("R", "S"):
                continue
            symbols = dets.get(bs, [])
            pred = fn(symbols)["prediction"]
            records.append({"biosample": bs, "pred": pred, "label": rs})
            if str(pred).upper() in ("R", "S"):
                scored.append((pred, 1 if rs == "R" else 0))
        conf = _conf(scored)
        n_R = conf["tp"] + conf["fn"]
        n_S = conf["tn"] + conf["fp"]
        powered = n_R >= MIN_PER_CLASS and n_S >= MIN_PER_CLASS
        # DEGENERATE guard: a cell that predicts all-one-class scores 1.0 on one metric while being
        # useless (azithromycin here: all-S -> sens 0.0, spec 1.0). Such a cell is NOT endorsed even
        # though spec passes the floor -- it does not discriminate. (Same class of bug as the
        # empty-assembly-as-S gate: a degenerate output must not read as a validated one.)
        degenerate = ((n_R >= MIN_PER_CLASS and conf["sens"] == 0.0)
                      or (n_S >= MIN_PER_CLASS and conf["spec"] == 0.0))
        endorsed = bool(powered and not degenerate
                        and conf["spec"] is not None and conf["spec"] >= SPEC_FLOOR
                        and (conf["sens"] is None or conf["sens"] >= 0.5))
        headline = ("DEGENERATE_NOT_ENDORSED" if degenerate
                    else "SCORED_ENDORSED" if endorsed
                    else "UNDERPOWERED" if not powered
                    else "SCORED_NOT_ENDORSED")
        results[drug] = {
            "drug": drug, "rule": fn(["_probe_"])["rule"], "binary": conf, "n_R": n_R, "n_S": n_S,
            "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
            "powered": powered, "headline": headline,
        }
        print(f"  {drug:14s}: n={conf['n_scored']:3d} ({n_R}R/{n_S}S) acc={conf['acc']} "
              f"sens={conf['sens']} spec={conf['spec']} -> {headline}")

    artifact = {
        "_schema": "gono-ncbipd-extval-v1", "date": _date.today().isoformat(),
        "organism": "Neisseria gonorrhoeae", "cell": "neisseria_amr (NON-FROZEN)",
        "cohort_source": "NCBI Pathogen Detection PDG000000032.673 (measured AST + NCBI AMRFinderPlus calls)",
        "n_isolates": len(labels),
        "independence": ("provenance-disjoint from the AR-Bank isolates (different isolates/BioSamples); "
                         "NOT methodology-independent (same AMRFinderPlus + same neisseria_amr cell)."),
        "clonality_caveat": ("RAW sens/spec is clonality-inflated (gonococci clonally structured); "
                             "lineage-collapse (Mash on assemblies) is the follow-up."),
        "spec_floor": SPEC_FLOOR, "min_per_class": MIN_PER_CLASS,
        "results": results, "frozen_surface_changed": False,
    }
    out = Path(f"wiki/gono_ncbipd_extval_{_date.today().isoformat()}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"artifact: {out}")
    return 0 if any(r["powered"] for r in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
