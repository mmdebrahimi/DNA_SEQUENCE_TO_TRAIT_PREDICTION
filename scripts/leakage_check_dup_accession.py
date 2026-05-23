# leakage_check_dup_accession.py
# Run from: dna_decode_repo/
# Output: dna_decode_repo/reports/cipro_leakage_check_dup_accession_2026-05-22.{json,md}
# Runtime: < 5 seconds
# Gating: BLOCKS bounded falsifier interpretation if loso_leakage_present == True.

from pathlib import Path
import json
import pickle
import pandas as pd

REPO = Path(__file__).resolve().parent
COHORT = REPO / "data/processed/stage2_n150_cipro_cohort.parquet"
MODEL  = REPO / "data/processed/models/ciprofloxacin_nucleotide_transformer.pkl"
DUP_ACC = "GCA_025200635.1"
OUT_JSON = REPO / "reports/cipro_leakage_check_dup_accession_2026-05-22.json"
OUT_MD   = REPO / "reports/cipro_leakage_check_dup_accession_2026-05-22.md"

cohort = pd.read_parquet(COHORT)

# 1. Confirm the duplicate
dup_rows = cohort[cohort["assembly_accession"] == DUP_ACC]
strain_ids = sorted(dup_rows["strain_id"].astype(str).tolist())
labels = sorted(dup_rows["label"].astype(str).tolist())

# 2. Did LOSO put them in the same fold or different folds?
#    By LOSO construction (each strain_id IS its own fold), if two distinct strain_ids
#    share an assembly_accession, the OTHER duplicate is in the training set for each.
loso_leakage = len(strain_ids) >= 2

# 3. Inspect trained model provenance
training_cohort = None
trained_on = None
n_strains_trained = None
try:
    with open(MODEL, "rb") as f:
        model_obj = pickle.load(f)
    if isinstance(model_obj, dict):
        provenance = model_obj.get("provenance", {})
        training_cohort = provenance.get("training_cohort")
        trained_on = provenance.get("trained_on")
        n_strains_trained = provenance.get("n_strains")
except Exception as exc:
    training_cohort = f"<error loading model: {exc}>"

# 4. AUROC inflation upper bound (pairs over N)
n_total = len(cohort)
upper_bound_inflation = round(2.0 / n_total, 4) if n_total else None

# 5. Sanity: label consistency on the duplicate
labels_consistent = len(set(labels)) == 1

result = {
    "duplicate_accession": DUP_ACC,
    "strain_ids_sharing_accession": strain_ids,
    "labels_each": labels,
    "labels_consistent": labels_consistent,
    "loso_leakage_present": loso_leakage,
    "n_total_cohort": int(n_total),
    "auroc_inflation_upper_bound_pairs_over_N": upper_bound_inflation,
    "model_provenance_training_cohort": training_cohort,
    "model_provenance_trained_on": trained_on,
    "model_provenance_n_strains": n_strains_trained,
    "recommendation": (
        "DEDUP cohort (keep one of strain_ids) + retrain XGBoost + recompute LOSO AUROC "
        "BEFORE interpreting bounded falsifier results"
        if loso_leakage else
        "no same-accession leakage in cohort; proceed to bounded falsifier"
    ),
    "blocking_for_falsifier": loso_leakage,
}

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(result, indent=2))

md = []
md.append("# Cipro Leakage Check - Duplicate Accession")
md.append("")
md.append(f"- Duplicate accession: `{DUP_ACC}`")
md.append(f"- Sharing strain_ids: {strain_ids}")
md.append(f"- Labels: {labels} (consistent={labels_consistent})")
md.append(f"- LOSO leakage present (same genome in train + test fold): **{loso_leakage}**")
md.append(f"- AUROC inflation upper bound (pairs/N): {upper_bound_inflation}")
md.append(f"- Model trained_on: {trained_on}")
md.append(f"- Model training_cohort: {training_cohort}")
md.append(f"- Model n_strains_trained: {n_strains_trained}")
md.append(f"- Blocking for falsifier: **{result['blocking_for_falsifier']}**")
md.append(f"- Recommendation: {result['recommendation']}")
md.append("")
OUT_MD.write_text("\n".join(md))

print(json.dumps(result, indent=2))
