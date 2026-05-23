# Cipro Bounded Falsifier — Coordination & Execution Plan

**Date:** 2026-05-22
**Authors:** Claude (GTX 860M laptop) + Codex (Precision 7780, RTX 3500 Ada)
**Status:** Ready to execute. Codex owns falsifier runner; Claude owns subset, leakage check, doc template.
**Substrate:** N=147 cipro cohort, 67 audited R strains, NT v2 100M frozen + XGBoost classifier.
**Inputs read:** `current_cipro_interpretability_audit_2026-05-21.{json,md}`, Codex bounded-falsifier spec.

---

## 0. TL;DR — What changes vs Codex's spec

| Item | Codex's spec | This plan | Why |
|---|---|---|---|
| Subset choice | "12 strains" (unspecified) | 4 ERS + 4 ELX-family + 4 all-negative-delta, named below | Hits 3 falsifier modes: control / failure case / negative-delta fix |
| Method change | Positive-only Δ ranking | Same + max-abs-Δ side-channel for diagnostic | Cheap; disambiguates saturation from confound |
| Pre-check | None | Leakage check on `GCA_025200635.1` (parallel) | Same-genome appears as 562.109860 + 562.111036 in cohort — could be in train + test split |
| Decision rule | pass/fail → continue/ship scope-limit | Same, with explicit residual-uncertainty gate before declaring "pass" | One subset of 12 can't certify the method on N=147 |

Net: keep Codex's bounded falsifier as the *primary* gate. Run leakage check in parallel (cheap). Diagnostic exports gated on falsifier outcome. Mash-cluster deferred until after falsifier verdict.

---

## 1. Coordination protocol

### Roles
- **Codex (Precision 7780):** Owns model, cache, cohort parquet, ISM runner. Implements + runs falsifier. Implements diagnostic exports.
- **Claude (GTX 860M laptop):** Owns analysis memos, plan documents, subset selection logic, audit JSON interpretation. Cannot run model.

### File transfer cadence
- Codex → Claude: every artifact (JSON + MD sidecar) gets copied to `C:\Users\Farshad\Downloads\` for transfer (Gmail .tab workaround if needed for `.py`).
- Claude → Codex: every plan / template lands in `C:\Users\Farshad\Downloads\` with date stamp and is also written to `dna_decode/reports/` if applicable.

### Naming convention (mandatory — both sides)
```
cipro_<artifact>_<YYYY-MM-DD>.{md,json}     # narrative + machine-readable pair
cipro_<artifact>_runner_<YYYY-MM-DD>.py     # runner script if standalone
```

### Escalation triggers (stop + check with the other session)
1. Falsifier pass criterion ambiguous (e.g., 3 of 4 ELX strains recover) — don't auto-decide.
2. Leakage check returns "duplicate in same split" — model AUROC is suspect; halt falsifier interpretation.
3. Any single ISM run > 30 min on Precision 7780 — abort; likely a runtime regression vs the audit's ~95s/strain.
4. Any disagreement between this plan and Codex's spec — Codex's spec wins on *runner mechanics*; this plan wins on *subset + leakage + scope-limit doc*.

---

## 2. Bounded falsifier — deliberate 12-strain subset

### Selection logic
3 buckets × 4 strains each. Each bucket tests a distinct failure mode:

- **Bucket A (control / "method works"):** ERS-prefix strains where QRDR attribution already ranks ≤ top-50. If positive-only Δ still recovers these, the runner is sane.
- **Bucket B (test case / "method collapses on near-clonal block"):** ELX-family strains with QRDR rank > 500. Falsifier must improve these; if it doesn't, ranking change isn't the fix.
- **Bucket C (negative-delta fix / "all hits ≤ 0"):** Strains where every known-locus hit currently has Δ ≤ 0. Positive-only ranking is structurally undefined here — the runner must handle this case explicitly (and the outcome documents whether the all-negative pattern reflects saturation or genuine non-signal).

### The 12 strains

#### Bucket A — ERS control (4 strains)
QRDR rank ≤ 50; positive-Δ already present.

| strain_id | accession | locus_tag prefix | best_QRDR_rank | max_pos_Δ |
|---|---|---|---|---|
| 562.7693 | GCA_001284665.1 | ERS | 12 | 0.1520 |
| 562.7775 | GCA_001286305.1 | ERS | 31 | 0.1289 |
| 562.7660 | GCA_001284005.1 | ERS | 32 | 0.0077 |
| 562.7572 | GCA_001277535.1 | ERS | 41 | 0.2021 |

Pass criterion for Bucket A: ≥ 3 of 4 still rank a QRDR locus in top-10 under positive-only Δ.

#### Bucket B — ELX-family failure cases (4 strains)
QRDR rank > 500 (deep failure); represents the 30-strain ELX-family block.

| strain_id | accession | locus_tag prefix | best_QRDR_rank | max_pos_Δ |
|---|---|---|---|---|
| 562.45910 | GCA_004569165.1 | ELT | 1727 | 0.0001 |
| 562.45849 | GCA_004567165.1 | ELX | 3178 | 0.0017 |
| 562.50247 | GCA_004566895.1 | ELX | 2644 | 0.0009 |
| 562.50233 | GCA_004566595.1 | ELY | 1110 | 0.0010 |

Pass criterion for Bucket B: ≥ 2 of 4 move into top-10 (positive-only ranking) AND median rank shift is ≥ 100-fold improvement (e.g., 1727 → ≤ 17).

#### Bucket C — All-negative-Δ strains (4 strains)
Every known-locus hit currently has Δ ≤ 0; `best_known_locus_rank = None` for 4 of 5 (one has positive QRDR Δ but no positive `best_known_locus_rank`).

| strain_id | accession | locus_tag prefix | n_hits | original best_rank |
|---|---|---|---|---|
| 562.115106 | GCA_026420925.1 | OLZ | 5 | None |
| 562.17620 | GCA_002192275.1 | AM | 4 | None |
| 562.19632 | GCA_002180195.1 | AM | 4 | None |
| 562.34085 | GCA_003571665.1 | AM | 4 | None |

Pass criterion for Bucket C: positive-only ranking either (a) recovers at least 1 hit into top-50 OR (b) the runner explicitly emits `INDETERMINATE_ALL_NEGATIVE_DELTA` and the diagnostic exports show classifier baseline P(R) ≥ 0.95 (saturation) — in which case "all-negative" reflects saturation, not attribution failure, and the result is interpretable.

### Overall falsifier verdict

| Bucket A | Bucket B | Bucket C | Verdict |
|---|---|---|---|
| pass | pass | pass-or-saturated | **PASS** → continue method refinement on full N=67 |
| pass | fail | any | **FAIL — ranking is not the bottleneck** → ship v0 with scope-limit doc |
| fail | any | any | **RUNNER REGRESSION** → halt; debug runner against audit baseline before interpreting |
| any | any | "all positive Δ become indeterminate" | **METHOD CHANGE BREAKS WORKING STRAINS** → revert; ship v0 with scope-limit |

---

## 3. Parallel train/test leakage check (Claude-owned snippet, Codex-executed)

### Why
`GCA_025200635.1` appears twice in the cohort as `562.109860` AND `562.111036` (byte-identical audit results confirm same genome). If LOSO CV put these two strain_ids in different folds, the model was trained on the test sample → AUROC is inflated. This invalidates *any* downstream interpretability conclusion.

### Snippet (Codex runs on Precision 7780)

```python
# leakage_check_dup_accession.py
# Run from: dna_decode_repo/
# Output: dna_decode_repo/reports/cipro_leakage_check_dup_accession_2026-05-22.{json,md}

from pathlib import Path
import json
import pandas as pd
import pickle

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
#    (LOSO = leave-one-strain-out, so each strain_id IS its own fold; the
#    question is whether the OTHER duplicate was in the training set for each.)
loso_leakage = len(strain_ids) >= 2  # by LOSO construction: yes, each is in the other's train fold
# Treat this as the formal "is there same-genome leakage between train and held-out?" check.

# 3. Reload the trained model + inspect fold metadata if persisted
with open(MODEL, "rb") as f:
    model_obj = pickle.load(f)
provenance = model_obj.get("provenance", {}) if isinstance(model_obj, dict) else {}
training_cohort = provenance.get("training_cohort")
trained_on = provenance.get("trained_on")
n_strains = provenance.get("n_strains")

# 4. Estimate AUROC inflation upper bound:
#    if duplicates always co-occur (train + test), the LOSO AUROC overcounts
#    by 2/N pairs at most.
n_total = len(cohort)
upper_bound_inflation = round(2.0 / n_total, 4) if n_total else None

result = {
    "duplicate_accession": DUP_ACC,
    "strain_ids_sharing_accession": strain_ids,
    "labels_each": labels,
    "labels_consistent": len(set(labels)) == 1,
    "loso_leakage_present": loso_leakage,
    "n_total_cohort": int(n_total),
    "auroc_inflation_upper_bound_pairs_over_N": upper_bound_inflation,
    "model_provenance_training_cohort": training_cohort,
    "model_provenance_trained_on": trained_on,
    "model_provenance_n_strains": n_strains,
    "recommendation": (
        "DEDUP cohort + retrain before any interpretability conclusion is drawn"
        if loso_leakage else
        "no same-accession leakage in cohort; proceed"
    ),
}

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(result, indent=2))
OUT_MD.write_text(
    "# Cipro Leakage Check — Duplicate Accession\n\n"
    f"- Duplicate accession: `{DUP_ACC}`\n"
    f"- Sharing strain_ids: {strain_ids}\n"
    f"- Labels: {labels} (consistent={result['labels_consistent']})\n"
    f"- LOSO leakage present (same genome in train + test fold): **{loso_leakage}**\n"
    f"- AUROC inflation upper bound (pairs/N): {upper_bound_inflation}\n"
    f"- Model trained_on: {trained_on}\n"
    f"- Recommendation: {result['recommendation']}\n"
)
print(json.dumps(result, indent=2))
```

### Interpretation matrix

| labels_consistent | loso_leakage_present | Decision |
|---|---|---|
| True | True | DEDUP cohort + retrain. Falsifier results stand only if model rebuilt on dedup'd cohort. |
| True | False | proceed (duplicate doesn't span train/test) |
| False | * | hard error: same genome has conflicting R/S labels → audit AST pipeline |
| * | * | if upper bound > 1.5 pp on a 75/75 cohort → flag in scope-limit doc even if dedup retrain unfeasible |

Run this **before or alongside** the falsifier — it's < 5 seconds and gates whether the falsifier result is interpretable at all.

---

## 4. Diagnostic exports (Codex-implemented, gated on falsifier outcome)

Trigger: falsifier returns PASS or "method change breaks working strains" (anything except clean FAIL on Bucket B).

For each of the 12 strains in the falsifier subset, export:

| Field | Definition | Use |
|---|---|---|
| `baseline_proba_R` | classifier output on full mean-pool, no ISM | Detect saturation (≥0.95 → all-negative Δ reflects saturation, not failure) |
| `baseline_logit` | pre-sigmoid logit on full mean-pool | Larger dynamic range than probability for Δ interpretation |
| `per_known_locus_proba_drop` | proba(full) − proba(drop locus) | Currently `prediction_delta`; keep |
| `per_known_locus_logit_drop` | logit(full) − logit(drop locus) | New; disambiguates saturated from true-zero contribution |
| `max_abs_delta_all_genes` | max over all ~5000 cached genes of \|Δproba\| | If max-abs-Δ on a strain is < 0.01 globally → classifier is saturated; no gene matters; "ELX failure" is calibration not lineage |
| `top_k_positive_delta_aliases` | aliases of top-10 genes by positive Δ | The falsifier's ranking output |

JSON schema (one entry per strain):
```json
{
  "strain_id": "562.50226",
  "accession": "GCA_004566415.1",
  "label": "R",
  "baseline_proba_R": 0.9871,
  "baseline_logit": 4.331,
  "max_abs_delta_all_genes": 0.0046,
  "saturation_flag": true,
  "per_known_locus": [
    {"alias": "gyrA", "proba_drop": 0.0046, "logit_drop": 0.412, "rank_pos_delta": 6}
  ],
  "top_10_positive_delta": [...]
}
```

Saturation flag rule: `baseline_proba_R >= 0.95 AND max_abs_delta_all_genes < 0.01` → saturation_flag=true. Use this to discriminate the 4 hypotheses I named in the original analysis (lineage confound vs annotation quality vs genome size vs classifier saturation). Cipher: if every Bucket B + C strain has saturation_flag=true, the headline finding is **classifier saturation, not lineage confound** — and the architectural fix is calibration, not Mash-clade stratification.

---

## 5. Mash-cluster work — deferred, gated on falsifier verdict

| Falsifier verdict | Mash-cluster work | Why |
|---|---|---|
| PASS | run full N=147 (R+S both) Mash sketch + dist; re-stratify LOSO; report per-clade AUROC | Method works; the question is whether the *training set* is over-represented on one clade |
| FAIL (ship scope-limit) | skip Mash-cluster entirely | Method bottleneck is mechanism class (distributed mobile-element), not lineage |
| RUNNER REGRESSION | skip | Debug runner first |
| Saturation diagnosed | run Mash-cluster only if calibration fix lands and AUROC restored | Saturation supersedes lineage confound |

Mash command (when triggered, on Precision 7780):
```bash
# requires Docker Desktop running + the pinned image
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "$PWD/data/processed/fastas":/data \
  quay.io/biocontainers/mash:2.3--hb105d93_10 \
  bash -c "mash sketch -o /data/sketch /data/*.fna && mash dist /data/sketch.msh /data/sketch.msh > /data/dist.tsv"
```

Use existing `dna_decode/eval/phylogeny.py::compute_mash_distances(..., use_docker=True)` rather than calling docker directly.

---

## 6. v0 ship-or-document-scope-limit fork

### If falsifier PASSES
1. Apply method change (positive-only Δ ranking) to full N=67 R strains; persist new audit JSON.
2. Run Mash-cluster diagnostic (Section 5).
3. Update `wiki/decoder_v0_ux_and_success_criterion.md` interpretability success criterion if QRDR top-10 recovery rate ≥ 50%.
4. Re-fire `scripts/pipeline.py predict` on 3 sample strains; check the new attribution surfaces in v0 JSON output.
5. Commit + tag v0.

### If falsifier FAILS
1. Skip method-refinement track entirely (don't burn cycles on ranking variants).
2. Ship v0 anyway, BUT with a scope-limit doc explicit about the interpretability failure mode.
3. The decoder still emits a top-K attribution table — but with a calibrated honesty caveat string in JSON + markdown sidecar.
4. v0 still passes Functional HARD + Predictive + Documentation criteria; **Interpretable criterion ships at "partial" tier with named scope-limit**.

The point: **shipping v0 with a documented scope-limit is a success**, not a failure. The north star is "AI DNA decoder tool, not papers" — a v0 that's honest about where attribution works and where it doesn't is exactly the artifact the project should produce.

---

## 7. v0 scope-limit doc template (fill-in-the-blank, used only if falsifier FAILS)

File: `wiki/cipro_v0_attribution_scope_limit_<YYYY-MM-DD>.md`

```markdown
# Cipro v0 — Attribution Scope Limit

**Status:** v0 SHIPS with documented attribution scope limit.
**Date:** <YYYY-MM-DD>
**Trigger:** Bounded falsifier (positive-only Δ ranking) failed on Bucket B (ELX-family) — see `cipro_bounded_falsifier_results_<DATE>.{md,json}`.

## What works
- **Predictive performance:** LOSO AUROC = <FILL> on N=147; meets v0 Predictive criterion (≥ 0.70).
- **Attribution on diverse-lineage cipro-R strains:** ERS-prefix subset (~12% of cohort) — QRDR loci rank in top-10 for X of 4 strains.
- **Audit framework:** mechanism × MIC × opacity merge + SUSPEND gate propagates to every prediction's `audit_verdict` field.

## What does NOT work
- **Attribution on near-clonal cipro-R blocks:** ELX/ELY/ELV/ELU/ELT-family strains (~45% of cohort, GCA_004566xxx–004570xxx batch) — QRDR loci rank > 500 in median; positive-only Δ ranking does NOT recover this. Median rank under failed-falsifier ranking: <FILL>.
- **Strains with all-negative Δ on known loci:** N=5 strains (562.115106 / 562.17620 / 562.19632 / 562.34085 / 562.50269) — `best_known_locus_rank=None`. May reflect classifier saturation rather than absent mechanism signal — diagnostic exports `baseline_proba_R` and `max_abs_delta_all_genes` reported in v0 JSON.
- **Distributed mobile-element resistance** (tet phenotype): mechanism class is OUT-OF-SCOPE for v0; cef + tet substrates deferred per cross-drug architectural finding (`wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`).

## v0 decoder behavior on scope-limited cases
- **predict CLI output:** every prediction includes `top_k_attribution` AND a new field `attribution_scope_confidence` ∈ {HIGH, PARTIAL, INDETERMINATE}.
  - HIGH: strain's `best_known_locus_rank` ≤ 50 in the bounded-falsifier subset *or* the closest-by-Mash-distance falsifier-PASS strain.
  - PARTIAL: strain is in a clade with at least one PASS exemplar but not directly verified.
  - INDETERMINATE: strain matches the all-negative-Δ pattern (`saturation_flag=true` OR `n_positive_delta_known_loci=0`).
- Markdown sidecar shows the same flag prominently.

## What unblocks "HIGH" attribution for all strains
1. Calibration fix (if saturation is the dominant cause) — temperature scaling on validation logits.
2. Per-gene NT window approach (replacing mean-pool with sliding-window concatenation + attention) — Phase 2 redesign work.
3. Larger + more diverse cohort with the BV-BRC `assembly_accession` bottleneck resolved (e.g., source-database swap to PATRIC / NARMS / EuSCAPE).

## Honest-criterion compliance
Per v0 spec (`wiki/decoder_v0_ux_and_success_criterion.md`), the **Honest HARD criterion** requires that the decoder never claim higher attribution confidence than the evidence supports. This scope-limit doc + the `attribution_scope_confidence` field satisfy that requirement.

## Reproducing the scope limit
```bash
uv run python scripts/cipro_bounded_falsifier.py \
  --cohort data/processed/stage2_n150_cipro_cohort.parquet \
  --model data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --cache <cache-path> --subset reports/cipro_bounded_falsifier_subset_2026-05-22.json
```
```

---

## 8. Decision tree summary (one-pass execution)

```
┌─ Step 1: Run leakage_check_dup_accession.py (Codex, < 5 s)
│      │
│      ├─ loso_leakage_present == True AND labels_consistent == True
│      │      → halt; dedup cohort + retrain BEFORE running falsifier
│      │      → fork: light dedup → keep ONE of (562.109860, 562.111036), retrain XGBoost,
│      │              recompute LOSO AUROC, then proceed
│      │
│      └─ otherwise → proceed
│
├─ Step 2: Run bounded falsifier on the 12-strain subset above (Codex, ~10-20 min)
│      │
│      ├─ Bucket A pass + Bucket B pass + Bucket C handled → PASS
│      │      → Step 3a (diagnostic exports + Mash-cluster + ship v0 standard)
│      │
│      ├─ Bucket A pass + Bucket B fail → FAIL
│      │      → Step 3b (scope-limit doc + ship v0 partial-attribution)
│      │
│      ├─ Bucket A fail → RUNNER REGRESSION
│      │      → halt; debug runner against audit baseline
│      │
│      └─ Method change breaks Bucket A → REVERT
│             → Step 3b (scope-limit doc + ship v0 with current attribution)
│
├─ Step 3a (PASS path):
│      │   Diagnostic exports → Mash-cluster N=147 → per-clade LOSO AUROC →
│      │   if any clade < 0.70 AUROC, scope-limit that clade explicitly →
│      │   update decoder_v0 spec, tag v0
│      │
│      └─ DONE
│
└─ Step 3b (FAIL/REVERT path):
       │   Fill scope-limit doc template (Section 7) →
       │   Add attribution_scope_confidence field to predict CLI →
       │   Tag v0 with scope-limit referenced in README + LESSONS_LEARNED
       │
       └─ DONE
```

---

## 9. Open questions (carried into next session)

1. Does the trained XGBoost classifier's `provenance` dict actually populate `training_cohort` + `trained_on`? (Claude's predict rewrite assumed yes; the leakage check will reveal it.) If empty → minor regression in `cmd_train` — file a bug.
2. Is the falsifier runner's positive-only Δ ranking equivalent to ranking by `max(0, prediction_delta)` or by `prediction_delta > 0` filter then re-rank? Codex's spec should pin this. The two are not identical when ties exist.
3. For Bucket C (all-negative-Δ), if diagnostic exports show `baseline_proba_R < 0.95` + `max_abs_delta_all_genes > 0.01`, what does the all-negative pattern mean? Possibly: the classifier learned a *negative-evidence* feature (presence of a gene tags the strain as S). Worth a Phase 2 brainstorm.
4. Is there a cheaper subset-of-12 — e.g., 3 buckets × 4 vs 4 × 3 — that would catch more failure modes? (Current 4×3 prioritizes per-bucket statistical confidence; 3×4 would add an "intermediate-rank" bucket but split power.)
5. If we PASS the falsifier but Mash-cluster reveals one clade < 0.70 AUROC, do we scope-limit that clade in v0 or block v0 ship? Default: scope-limit, ship — per north star.

---

## 10. Concrete deliverables Claude produces from this laptop (no model access)

- ✅ This plan document (you're reading it).
- ✅ The 12-strain subset list (Section 2) — committable as `reports/cipro_bounded_falsifier_subset_2026-05-22.json`:
  ```json
  {
    "bucket_A_control_ERS": ["562.7693","562.7775","562.7660","562.7572"],
    "bucket_B_ELX_failure": ["562.45910","562.45849","562.50247","562.50233"],
    "bucket_C_all_negative_delta": ["562.115106","562.17620","562.19632","562.34085"]
  }
  ```
- ✅ The leakage-check snippet (Section 3) — committable as `scripts/leakage_check_dup_accession.py`.
- ✅ The v0 scope-limit doc template (Section 7) — committable as `wiki/_template_cipro_v0_attribution_scope_limit.md`.
- ✅ The diagnostic exports spec (Section 4) — Codex implements the export emitter inside the falsifier runner.

## 11. Concrete deliverables Codex produces from Precision 7780

- The falsifier runner (`scripts/cipro_bounded_falsifier.py`) per Codex's existing spec at `reports/cipro_bounded_falsifier_experiment_spec_2026-05-22.md`.
- The leakage-check result (Section 3 snippet, runtime < 5 s).
- The diagnostic exports JSON (Section 4, runtime ≈ 12 × 95 s ≈ 20 min if reusing audit per-gene Δ; instant if exporting from cached audit JSON).
- The falsifier results: `reports/cipro_bounded_falsifier_results_2026-05-22.{md,json}`.

---

## 12. Pre-commit safety check

Before either side commits any artifact:
- [ ] Leakage check result included in commit message (one line: `leakage_present=<bool>`).
- [ ] Falsifier verdict included in commit message (one line: `verdict=PASS|FAIL|REGRESSION|REVERT`).
- [ ] If FAIL or REVERT path: scope-limit doc filled out + linked in README.
- [ ] If PASS path: diagnostic exports JSON + Mash-cluster TSV both included.
- [ ] LESSONS_LEARNED.md gets a one-liner regardless of outcome.

---

**End of plan.** Codex: ack this in your next turn before running the falsifier. Claude: monitor `Downloads/` for falsifier output.
