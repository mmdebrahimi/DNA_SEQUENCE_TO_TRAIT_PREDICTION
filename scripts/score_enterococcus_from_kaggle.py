"""Score the AR-Bank E. faecium cell from the Kaggle AMRFinder factory output (no local Docker).

The Kaggle kernel `enterococcus-amrfinder` SRA-assembles + AMRFinder-scans the 17 AR-Bank E. faecium
isolates (the assembly path that wedges Docker locally) and emits one `<biosample>.amrfinder.tsv` per
isolate. This scorer pulls those, parses the Element symbols, runs the NON-FROZEN enterococcus_amr cell
(call_efm_tetracycline -> doxycycline; call_efm_ciprofloxacin -> levofloxacin), and scores vs the CDC
S/I/R labels -- the first Enterococcus AR-Bank validation (locally impossible: 0 downloadable assemblies).

  # 1. pull the Kaggle output (when the kernel finishes):
  uv run python scripts/kaggle_push_poll.py pull enterococcus-amrfinder data/raw/ar_bank_efm_kaggle
  # 2. score:
  uv run python scripts/score_enterococcus_from_kaggle.py --amrfinder-dir data/raw/ar_bank_efm_kaggle
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.organism_rules import enterococcus_amr as EFM  # noqa: E402
from scripts.amr_portal_score_independent import wilson_ci  # noqa: E402
from scripts.external_cohort_revalidate import _read_selected, parse_determinant_symbols  # noqa: E402
from scripts.independent_cohort_validate import _conf  # noqa: E402

SPEC_FLOOR = 0.85
# canonical drug -> (label dir, call function). Both scored from the SAME AMRFinder run per isolate.
DRUGS = {
    "doxycycline": ("data/raw/ar_bank_enterococcus_faecium_extval_doxycycline", EFM.call_efm_tetracycline),
    "levofloxacin": ("data/raw/ar_bank_enterococcus_faecium_extval_levofloxacin", EFM.call_efm_ciprofloxacin),
    # Same AMRFinder run per isolate -> two more powerable glycopeptide cells (van clusters):
    "vancomycin": ("data/raw/ar_bank_enterococcus_faecium_extval_vancomycin", EFM.call_efm_vancomycin),
    "teicoplanin": ("data/raw/ar_bank_enterococcus_faecium_extval_teicoplanin", EFM.call_efm_teicoplanin),
}


def _find_tsv(amr_dir: Path, biosample: str) -> Path | None:
    for name in (f"{biosample}.amrfinder.tsv", f"{biosample}.tsv"):
        p = amr_dir / name
        if p.exists():
            return p
    return None


# Assembly-quality floor: a genuine E. faecium assembly yields many determinants (13-24 in the good
# isolates here, incl. intrinsic aac(6')-Ii even in susceptible strains). An assembly that produces
# FEWER than this is a FAILED/empty assembly (single-end SKESA on the SAMN15040xxx block emitted 0-1
# determinants = header-only TSVs). Such an isolate MUST be INDETERMINATE (assembly_failed), NOT scored
# S-by-absence -- else an empty assembly masquerades as susceptible + inflates specificity. The floor is
# deliberately low (an isolate can be genuinely determinant-light) but 0-2 = the empty-assembly signature.
MIN_DETERMINANTS_FOR_VALID_ASSEMBLY = 3


def score_drug(drug: str, labels_dir: str, call_fn, amr_dir: Path, min_per_class: int) -> dict:
    labels = _read_selected(Path(labels_dir) / "selected_strict.tsv")
    records, n_missing, n_assembly_failed = [], 0, 0
    for bs, rs in sorted(labels.items()):
        tsv = _find_tsv(amr_dir, bs)
        if tsv is None:
            n_missing += 1
            continue
        symbols = parse_determinant_symbols(tsv)
        if len(symbols) < MIN_DETERMINANTS_FOR_VALID_ASSEMBLY:
            # empty/failed assembly -> INDETERMINATE, not S-by-absence
            n_assembly_failed += 1
            records.append({"biosample": bs, "prediction": "INDETERMINATE_ASSEMBLY_FAILED",
                            "label": rs, "y": 1 if rs == "R" else 0, "symbols": symbols,
                            "n_determinants": len(symbols)})
            continue
        pred = call_fn(symbols)["prediction"]
        records.append({"biosample": bs, "prediction": pred, "label": rs,
                        "y": 1 if rs == "R" else 0, "symbols": symbols,
                        "n_determinants": len(symbols)})
    scored = [(r["prediction"], r["y"]) for r in records if str(r["prediction"]).upper() in ("R", "S")]
    conf = _conf(scored)
    n_R = conf["tp"] + conf["fn"]
    n_S = conf["tn"] + conf["fp"]
    powered = n_R >= min_per_class and n_S >= min_per_class
    endorsed = bool(powered and conf["spec"] is not None and conf["spec"] >= SPEC_FLOOR)
    return {
        "drug": drug, "gene_rule": call_fn(["_probe_"])["rule"], "binary": conf,
        "n_R": n_R, "n_S": n_S, "n_missing_assembly": n_missing,
        "n_assembly_failed": n_assembly_failed,
        "sens_wilson95": wilson_ci(conf["tp"], n_R), "spec_wilson95": wilson_ci(conf["tn"], n_S),
        "powered": powered, "spec_floor": SPEC_FLOOR,
        "headline": "SCORED_ENDORSED" if endorsed else ("UNDERPOWERED" if not powered else "SCORED_NOT_ENDORSED"),
        "records": records,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--amrfinder-dir", default="data/raw/ar_bank_efm_kaggle",
                    help="dir of Kaggle-produced <biosample>.amrfinder.tsv files")
    ap.add_argument("--min-per-class", type=int, default=5)
    ap.add_argument("--run-id", default=f"efm_kaggle_{_date.today().isoformat()}")
    a = ap.parse_args()
    amr_dir = Path(a.amrfinder_dir)
    if not amr_dir.exists():
        print(f"REFUSE: AMRFinder dir {amr_dir} not found -- pull the Kaggle output first "
              f"(kaggle_push_poll.py pull enterococcus-amrfinder {amr_dir})")
        return 2
    n_tsv = len(list(amr_dir.glob("*.amrfinder.tsv")) + list(amr_dir.glob("*.tsv")))
    print(f"AMRFinder TSVs present: {n_tsv}")

    results = {}
    for drug, (labels_dir, fn) in DRUGS.items():
        r = score_drug(drug, labels_dir, fn, amr_dir, a.min_per_class)
        results[drug] = r
        c = r["binary"]
        print(f"  {drug}: n={c['n_scored']} ({r['n_R']}R/{r['n_S']}S) acc={c['acc']} "
              f"sens={c['sens']} spec={c['spec']} missing={r['n_missing_assembly']} "
              f"asm_failed={r['n_assembly_failed']} -> {r['headline']}")

    artifact = {
        "_schema": "ar-bank-efm-kaggle-validation-v1", "date": _date.today().isoformat(),
        "run_id": a.run_id, "organism": "Enterococcus faecium",
        "rule_status": "CURATED_NONFROZEN", "rule_scope": "scorer_local", "not_in_shipped_surface": True,
        "genotype_source": ("Kaggle-native bioconda SRA-assembly (skesa) + AMRFinderPlus "
                            "-O Enterococcus_faecium; local Docker-wedge path bypassed"),
        "label_source": "ar_bank_INT_cdc (S/I/R; I excluded)",
        "independence_tier": ("CDC AR Isolate Bank measured S/I/R; NON-FROZEN curated cell (hand rule, not "
                              "tuned on any cohort -> no tuning-leak); NOT methodology-independent."),
        "results": results, "frozen_surface_changed": False,
    }
    out = Path(f"wiki/ar_bank_efm_kaggle_validation_{a.run_id}.json")
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"artifact: {out}")
    powered_any = any(r["powered"] for r in results.values())
    return 0 if powered_any else 1


if __name__ == "__main__":
    raise SystemExit(main())
