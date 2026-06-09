"""Build the calibrated per-organism AMR rule registry from cached AMRFinder cohorts.

Runs `calibrate_organism.calibrate` on each cached organism×drug cohort under data/raw and writes the
distilled configs to the COMMITTED registry `dna_decode/data/calibrated_amr_rules.json` (the cohorts
themselves are gitignored; the registry is the shippable output). `call_resistance(organism=...)` consults
this registry.

HONESTY: these calibrations are IN-SAMPLE (N~30 each, calibrated and evaluated on the same cohort). The
registry is therefore OPT-IN (DRUG_RULE stays the default when organism is unspecified) and stamped with
that provenance. The EXPRESSION_FLOOR abstain verdicts are safe regardless of in-sample-ness (abstaining is
conservative); the CALIBRATED configs need an independent cohort before becoming a default.

Usage: uv run python -m scripts.build_calibrated_registry   (or .venv/Scripts/python.exe scripts/...)
"""
from __future__ import annotations

import json
from pathlib import Path

from dna_decode.eval.calibrate_organism import calibrate, features_from_main_tsv
from scripts.organism_drug_validate import _run_dir

# (data/raw slug, registry organism key, drug)
COHORTS = [
    ("campylobacter_ciprofloxacin", "Campylobacter", "ciprofloxacin"),
    ("klebsiella_cipro", "Klebsiella", "ciprofloxacin"),
    ("salmonella_ciprofloxacin", "Salmonella", "ciprofloxacin"),
    ("acinetobacter_meropenem", "Acinetobacter", "meropenem"),
    ("pseudomonas_aeruginosa_meropenem", "Pseudomonas_aeruginosa", "meropenem"),
]
REGISTRY_PATH = Path("dna_decode/data/calibrated_amr_rules.json")


def _load_cohort(slug: str, drug: str):
    base = Path(f"data/raw/{slug}")
    sel = base / "selected.tsv"
    if not sel.exists():
        return None
    # Resolve each accession's AMRFinder run via the SAME reuse_glob the validator uses (own dir +
    # data/raw/<genus>_*/amrfinder_runs) — a strain's run is identical regardless of which drug cohort
    # downloaded it, and many runs live in a sibling cohort dir. Reading only base/amrfinder_runs
    # under-loads (the Pseudomonas degenerate-cohort bug: 7 of 30 found -> all-S -> bogus verdict).
    reuse_glob = f"data/raw/{slug.split('_')[0]}_*/amrfinder_runs"
    own_runs = base / "amrfinder_runs"
    strains, labels = [], []
    for ln in sel.read_text().splitlines():
        if "\t" not in ln:
            continue
        acc, lab = ln.split("\t")
        rd = _run_dir(acc, own_runs, reuse_glob)
        if rd is not None:
            strains.append(features_from_main_tsv(rd / "main.tsv", drug))
            labels.append(lab.strip())
    return (strains, labels) if strains else None


def build() -> dict:
    rules = {}
    for slug, organism, drug in COHORTS:
        loaded = _load_cohort(slug, drug)
        if loaded is None:
            print(f"  SKIP {organism}|{drug}: no cached cohort at data/raw/{slug}")
            continue
        strains, labels = loaded
        r = calibrate(strains, labels, drug)
        rules[f"{organism}|{drug}"] = {
            "verdict": r.verdict, "counter": r.counter, "threshold": r.threshold,
            "intrinsic_families_excluded": r.intrinsic_families_excluded,
            "loo_balanced_accuracy": r.loo_balanced_accuracy, "n": r.n, "n_R": r.n_R, "n_S": r.n_S,
        }
        print(f"  {organism}|{drug}: {r.verdict} {r.counter}@{r.threshold} "
              f"excl={r.intrinsic_families_excluded} LOO={r.loo_balanced_accuracy}")
    return {
        "_provenance": ("IN-SAMPLE calibration on cached NCBI/BV-BRC cohorts (N~30 each), 2026-06-08. "
                        "OPT-IN: call_resistance uses these ONLY when organism= is passed; DRUG_RULE "
                        "remains the default. CALIBRATED configs need an INDEPENDENT cohort before becoming "
                        "a default; EXPRESSION_FLOOR abstain verdicts are conservative and safe to ship. "
                        "See wiki/calibrate_organism_validation_2026-06-08.md."),
        "_schema": "calibrated-amr-rules-v1",
        "rules": rules,
    }


def main():
    reg = build()
    REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)
    REGISTRY_PATH.write_text(json.dumps(reg, indent=2))
    print(f"\nwrote {REGISTRY_PATH} ({len(reg['rules'])} rules)")


if __name__ == "__main__":
    main()
