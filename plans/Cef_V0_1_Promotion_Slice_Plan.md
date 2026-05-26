# Cef v0.1 Promotion Slice — Plan

> Promote cef cached-strain from "internal viability proven" (Codex 2026-05-25 handoff: 49 strains, AUROC 0.895, 2 real prediction examples) to a **formal v0.1 product packet** matching the cipro v0/v0.1 release-packet pattern. Per Codex's explicit ask: "recommend the smallest credible next slice to promote cef from internal viability to a formal v0.1 product packet, including whether duplicate-accession audit, audit-sidecar design, or genome-input cef should come next."

**Status:** DRAFT 2026-05-26.
**Anchors on:** `Downloads/dna_decoder_v0_1_cef_cached_handoff_2026-05-25.md` (Codex's actual cef cached-strain ship), `plans/Cef_Cached_Strain_V0_1_Slice_Plan.md` (the slice plan that informed Codex's work but was partially superseded by Codex going bigger).

**Supersedes:** `plans/Cef_Cached_Strain_V0_1_Slice_Plan.md` partially — Codex went bigger (dedicated 67-strain cohort + cache, NOT the N=12 mini reuse I scoped). The mini-cohort-based slice is no longer the right framing; promotion = closing the gaps on the actual cef substrate Codex built.

---

## What Codex shipped overnight (the actual current state on Precision 7780)

| Asset | Value |
|---|---|
| Cef cohort | `data/processed/gate_b_cohort.parquet` — 50 cef-labeled strains; 26R / 24S |
| Cef NT cache | `nt_gate_b_cohort_67.h5` — 67 strains; 312,585 embeddings |
| Cef trained classifier | `data/processed/models/ceftriaxone_nucleotide_transformer.pkl` — AUROC 0.895; 49 usable strains; 25R/24S; **CV strategy: `loso` with `strain_id` grouping** |
| Cef predict examples | `reports/dna_decoder_v0_1_cef_cached_example_{R,S}_2026-05-25.{md,json}` |

NONE of this is on origin per `git fetch`. Same drift pattern continues.

## Critical finding (Codex flagged this themselves)

**CV used `strain_id` grouping, NOT `leave_one_accession_out`.** AUROC 0.895 is suspect until duplicate-accession audit runs. The cipro N=147 cohort had a known duplicate (`GCA_025200635.1` = `562.109860` + `562.111036`); the cef cohort hasn't been audited the same way. If cef has analogous duplicates, the 0.895 is inflated by same-genome train/test leakage — same shape as the cipro 2026-05-22 LESSON.

This is the single most-important gap blocking promotion to a "formal v0.1 product packet."

## Design Decisions

### D1: Duplicate-accession audit is the FIRST promotion step

**Decision:** Before writing any cef release packet, run `scripts/leakage_check_dup_accession.py` (already exists in repo, drug-agnostic by accident — checks the COHORT, not just cipro) against the cef cohort. If duplicates found, RETRAIN with `leave_one_accession_out` CV before any promotion artifact is written.

**Rationale:** Per Codex's own recommendation #1; matches the 2026-05-22 LESSON on duplicate-accession leakage. Promoting a model with potentially-inflated AUROC to a release packet = repeat of the 2026-05-22 sin.

**Trade-off:** Adds 1-2 hr if dupes found + retrain needed. Acceptable — alternative is shipping a number that could be 0.05-0.10 pp inflated.

### D2: Audit-sidecar discipline is REQUIRED for canonical reporting (not optional)

**Decision:** Cef v0.1 promotion REQUIRES a cef merge-gate audit JSON (`--audit-merge-json` source). Until that exists, predict outputs are debug-mode + cannot be in a formal release packet. Per the RELOCKED v0 spec: "Canonical v0 reporting requires `--audit-merge-json`."

**Rationale:** v0.1-cipro shipped WITH audit propagation; v0.1-cef should match for product-contract consistency. Codex's overnight examples used `--allow-missing-audit --no-attribution`; that's fine for internal viability proof but NOT for the formal product packet.

**Trade-off:** Building the cef audit sidecar = run `scripts/drug_mechanism_audit.py --drug ceftriaxone` on the 50 cef strains (~50 × 95s = ~80 min on Precision 7780 Docker) + write a cef mechanism × phenotype merge JSON. ~2 hr total.

### D3: Cef genome-input is DEFERRED (don't change two axes at once)

**Decision:** Cef genome-input is v0.2+ work. Promotion stays at cached-strain only.

**Rationale:** Same "one axis at a time" discipline that worked for v0.1 cipro genome-input. Promote cef cached-strain first; do cef genome-input only after the cached promotion is clean.

**Trade-off:** Slower to "full feature parity" with cipro. Acceptable — the marginal value of cef genome-input over cef cached-strain is small until we have multiple cef-using customers.

### D4: AUROC re-report after leakage check, with explicit `_post_dedup` field if retrain happens

**Decision:** The cef release packet reports BOTH the original strain-id-CV AUROC AND the leakage-safe `leave_one_accession_out` AUROC (if retrain fires). Annotate the change explicitly. Match the cipro v0 pattern (cipro v0 RELOCKED spec ships `cv_strategy: leave_one_accession_out` + `cv_auroc: 0.8697`).

**Rationale:** Honest reporting + matches cipro provenance schema. Don't quietly replace one number with another; show the delta.

**Trade-off:** Slightly longer provenance block. Acceptable.

### D5: Release packet matches the cipro v0 packet structure verbatim

**Decision:** `reports/dna_decoder_cef_v0_1_release_candidate_2026-05-26.md` uses the same section structure as `reports/dna_decoder_v0_release_candidate_2026-05-24.md`. Differences are content, not format.

**Rationale:** Cross-drug consistency = future-me reading either packet finds the same layout. Galaxy-style "reproducibility as feature" discipline.

**Trade-off:** None.

---

## Implementation Steps (Codex executes on Precision 7780 after origin sync close)

### Step 0 (PRECONDITION): Sync the v0.1 cipro slice 1 + EP-0 close artifacts to origin

Codex's prior 2026-05-24 EP-0 close work + 2026-05-25 v0.1 cipro slice 1 + 2026-05-25 cef cached-strain work are ALL on Precision 7780 only. `scripts/cross_machine_sync_check.py` on this side reports 2/5 missing markers per the 2026-05-26 run.

```bash
git push origin main   # Codex on Precision 7780
```

Then this side: `git pull --ff-only origin main` → `scripts/cross_machine_sync_check.py` → confirm IN SYNC.

Until this happens, nothing else can ship.

### Step 1: Duplicate-accession audit on the cef cohort

```bash
# Adapt scripts/leakage_check_dup_accession.py to take --cohort arg (or fork it)
uv run python scripts/leakage_check_dup_accession.py \
  --cohort data/processed/gate_b_cohort.parquet \
  --model data/processed/models/ceftriaxone_nucleotide_transformer.pkl \
  --output reports/cef_leakage_check_2026-05-26.json
```

Decision tree on output:
- `loso_leakage_present: false` → proceed to Step 2.
- `loso_leakage_present: true` → retrain with `leave_one_accession_out` CV (Step 1.5).

### Step 1.5 (CONDITIONAL): Retrain cef classifier with leakage-safe CV

If Step 1 finds duplicates:

```bash
uv run python -m scripts.pipeline train \
  --drug ceftriaxone \
  --model nucleotide_transformer \
  --cohort data/processed/gate_b_cohort.parquet \
  --cache <nt_gate_b_cohort_67.h5> \
  --cv-strategy leave_one_accession_out \
  --output data/processed/models/ceftriaxone_nucleotide_transformer_leakage_safe.pkl
```

Record BOTH AUROCs (original strain-id-CV 0.895 AND new leakage-safe N). Decide which is canonical: if new AUROC ≥ 0.70, ship the leakage-safe model as the canonical. If it drops below 0.70, ship the cef v0.1 with documented scope-limit ("cef predictive performance dropped to <0.70 after leakage correction; same architectural class as tet failure mode; mechanism may be more distributed than initially measured").

### Step 2: Cef merge-gate audit sidecar

```bash
# Run AMRFinder mechanism audit on the 50 cef strains
uv run python -m scripts.drug_mechanism_audit \
  --drug ceftriaxone \
  --cohort data/processed/gate_b_cohort.parquet \
  --refseq-cache <refseq-cache-path> \
  --output wiki/ceftriaxone_mechanism_audit_2026-05-26.md
# ~80 min on Precision 7780 Docker.

# Then build the merge-gate JSON (cef mechanism × MIC × opacity).
# Per cef slice plan R2: write scripts/drug_mechanism_phenotype_merge.py
# parallel to the existing cipro version (do NOT refactor cipro per CLAUDE.md gotcha).
uv run python scripts/drug_mechanism_phenotype_merge.py \
  --drug ceftriaxone \
  --mechanism-audit wiki/ceftriaxone_mechanism_audit_2026-05-26.json \
  --mic-audit <if-exists> \
  --output wiki/cef_mechanism_phenotype_merge_2026-05-26.json
```

If MIC audit doesn't exist for cef yet, scope the merge gate to mechanism-only (no SUSPEND condition; just opacity flag + per-strain primary mechanism).

### Step 3: Re-run cef predict examples with audit propagation

```bash
# Same 2 strains Codex used overnight (562.12960 R + 562.7572 S):
uv run python -m scripts.pipeline predict \
  --drug ceftriaxone \
  --strain-id 562.12960 \
  --model-path data/processed/models/ceftriaxone_nucleotide_transformer.pkl \
  --cache <nt_gate_b_cohort_67.h5> \
  --annotations <562.12960.gff3> \
  --audit-merge-json wiki/cef_mechanism_phenotype_merge_2026-05-26.json \
  --output reports/dna_decoder_cef_v0_1_canonical_example_R_2026-05-26.json
# Same for 562.7572 S
```

Difference from Codex's overnight run: NO `--allow-missing-audit`, NO `--no-attribution`. Full canonical_audit_aware mode.

### Step 4: Cef release packet

Write `reports/dna_decoder_cef_v0_1_release_candidate_2026-05-26.md` matching the cipro v0 packet structure. Include:
- Sample command (canonical mode, audit-aware)
- 2 example outputs (R + S from Step 3)
- Provenance: `cv_strategy: leave_one_accession_out` (if Step 1.5 fired) OR `cv_strategy: loso (caveat: duplicate-accession audit returned no leakage)` (if Step 1 found nothing)
- AUROC: leakage-safe number (or the post-dedup number if retrain fired)
- Honest scope: N=49 usable strains; smaller than cipro N=146; cef strict-MIC infeasible at scale per 2026-05-18 census
- v0.2 horizon: cef genome-input OR cef-S label backfill OR external benchmark

### Step 5: Tag `v0.1-cef`

Local-only tag at the post-Step-4 commit. Don't push the tag until user confirms.

### Step 6: Push to origin + this side runs post-push checklist

```bash
git add data/processed/models/ceftriaxone_*.pkl reports/cef_*.{md,json} wiki/cef_*.json
git commit -m "ship(cef-v0.1): cef cached-strain release candidate"
git push origin main
```

This side: `git pull --ff-only origin main` → `scripts/cross_machine_sync_check.py` (expect KNOWN_DIVERGENCE_TARGETS to add cef markers in a follow-up commit) → `scripts/preflight_runnable.py` (cef_slice should flip to "RUNNABLE: yes" if Codex's pickle is leakage-safe).

---

## Success gate (the actual ship criteria)

A cef v0.1 promotion-slice ships when ALL of:

1. **Leakage check ran** + recorded in `reports/cef_leakage_check_2026-05-26.json`. If duplicates found, retrain with `leave_one_accession_out` CV happened.
2. **Canonical-mode predict outputs exist** (no `--allow-missing-audit`; full audit_verdict propagation).
3. **Release candidate doc** matches the cipro v0 packet structure + names honest scope.
4. **CV AUROC ≥ 0.70** on the canonical (leakage-safe) CV strategy. Per v0 spec.
5. **Pushed to origin** + this side's `scripts/preflight_runnable.py` reports `cef_slice: runnable=YES`.

---

## Codex's specific axis-choice question

> "decide whether the next axis is: cef genome-input / cef audit-aware packet / or a second-drug validation panel"

**Pick: cef audit-aware packet** (this plan's Steps 2-5). Reasons:
- It closes the v0 spec's HARD gate (audit-aware required for canonical reporting).
- It catches the strain_id-CV leakage risk THIS WEEK, not after a v0.2 ships on top of it.
- Cef genome-input is v0.2+ per D3 (one axis at a time).
- A second-drug validation panel (e.g., tet smoke) is interesting but is EP-1.5 architecture-decision territory, not cef promotion.

---

## Risks

- **R1 (HIGH):** If cef has duplicate accessions + retrain drops AUROC < 0.70, the cef v0.1 release ships with documented scope-limit. Same FAIL-branch ship pattern that cipro v0 used 2026-05-24. **Mitigation:** acknowledge in advance that this is acceptable + match the v0 north-star (honest failure-tolerant iteration).
- **R2 (MEDIUM):** Cef MIC audit doesn't exist; merge gate would be mechanism-only without MIC signal. **Mitigation:** scope the gate to mechanism + opacity only; cef MIC audit can be a v0.2 follow-up.
- **R3 (MEDIUM):** `scripts/drug_mechanism_phenotype_merge.py` doesn't exist yet (called out in `Cef_Cached_Strain_V0_1_Slice_Plan.md` R2). Needs to be written. ~1 hr.
- **R4 (LOW):** Cross-machine sync drift recurs again. `scripts/cross_machine_sync_check.py` + the 5 KNOWN_DIVERGENCE_TARGETS markers should catch it.

---

## What this plan deliberately does NOT cover

- Cef genome-input (v0.2+).
- Tet / gent / 4th-mechanism-class drugs (EP-1.5 architecture decision territory).
- Cef-S label backfill from PATRIC / NARMS / EuSCAPE (v0.2+).
- Per-gene NT windows architecture rework (EP-1.5).
- External benchmark for cef vs AMRFinderPlus / RGI (EP-1B-equivalent for cef; v0.2+).

---

## Estimated effort (Codex on Precision 7780)

- Step 0 (push): ~5 min
- Step 1 (leakage check): ~5 min
- Step 1.5 (retrain if needed): ~30 min
- Step 2 (mechanism audit + merge): ~2-3 hr
- Step 3 (canonical examples): ~10 min
- Step 4 (release packet): ~30 min
- Step 5 (tag): ~1 min
- Step 6 (push): ~5 min

**Total: ~3-4 hr if Step 1.5 fires; ~2-3 hr if no duplicates found.**

---

## Bottom line

The cef cached-strain substrate is REAL on Precision 7780 + ready for promotion. Three concrete gaps must close before formal v0.1: **leakage check**, **audit sidecar**, **canonical-mode examples + release packet**. Genome-input cef + second-drug validation panel are explicitly deferred.

Single most important step: **Codex on Precision 7780 must push to origin first** (same drift pattern as 2026-05-24 + 2026-05-25). Until that lands, this side cannot validate, can't tag, can't run the post-push checklist.
