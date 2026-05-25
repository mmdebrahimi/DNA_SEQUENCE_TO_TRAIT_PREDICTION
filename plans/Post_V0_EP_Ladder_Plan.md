# Post-v0 Evidence Packet Ladder — Plan

> Path from cipro v0 (cached-strain predictor, shipped 2026-05-24) toward the long-term goal: a DNA-input → phenotype + gene-level trait identification tool. Uses Evidence Packet labels per the 2026-05-15 repo framing reset; "Phase" labels stay retrospective-only.

**Status:** LOCKED 2026-05-25 after multiple brainstorm + probe + review rounds in this session.
**Supersedes:** the misframed "Commercial Discovery Prerequisites" technical plan drafted earlier (NOT saved; framing corrected by user).
**Long-term goal (re-confirmed by user 2026-05-25):** DNA input → phenotype + trait identification at gene level, per the original 2026-05-11 `/project-init` goals + CLAUDE.md L9 long-term vision. "Sell to a large company" is a possible downstream outcome IF the tool works — NOT the framing for current development.

---

## Problem Statement

After v0's FAIL-branch ship 2026-05-24 (cached-strain cipro predictor with documented scope-limit), the project needs an articulated multi-EP path from "current narrow v0" toward "DNA-input → phenotype tool at maturity." Without that ladder:

1. Per-EP decisions feel ad-hoc (e.g., "do cef next or Klebsiella?")
2. Architectural debt accumulates silently (the 2026-05-17 cross-drug finding showed mean-pooling fails on distributed mobile-element mechanisms; this is unaddressed)
3. "Decode any DNA" stays unbounded research-program scope (the 2026-05-11 `/project-init` verdict that triggered the bounded refinement) — phases must each have a **terminal claim** to bound the ladder

Three brainstorm rounds + 1 /probe + 1 /review converged on this ladder. Saving it before further drift.

---

## Design Decisions

### D1: Use Evidence Packet labels, not "Phase" labels

**Decision:** EP-0, EP-1A, EP-1B, EP-1C, EP-1.5, EP-2, EP-2.5, EP-3, EP-4+. NOT "Phase 1/2/3."

**Rationale:** Per CLAUDE.md "Project framing: Phase 1 / 2 / 3 labels are now retrospective-only" (2026-05-15 framing reset). Reintroducing "Phase" labels for forward work conflicts with the established convention.

**Trade-off:** EP labels are less linearly clear than numbered Phases. Acceptable — consistency with the existing ledger matters more.

### D2: Each EP must have a terminal claim

**Decision:** Every EP in the ladder states an explicit, testable claim that bounds its scope. E.g., EP-1A = "E. coli cipro genome-input predictor with same-strain parity ≤ ε."

**Rationale:** The original 2026-05-11 `/project-init` flagged "decode any DNA" as RESEARCH-PROGRAM (unbounded). Terminal claims per EP prevent the unbounded-scope failure mode from recurring per-EP.

**Trade-off:** Defining terminal claims upfront adds planning friction. Mitigation: terminal claims can be revised post-completion (not pre-committed forever); they exist to bound execution, not to lock down forever.

### D3: EP-2 (multi-drug) is forked by mechanism class, NOT drug list

**Decision:** Cef (concentrated β-lactamase signal) reuses the v0 architecture. Tet (distributed mobile-element mechanisms) does NOT — must use the EP-1.5 architectural decision.

**Rationale:** The 2026-05-17 cross-drug architectural finding (`wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`) showed NT whole-genome mean-pooling PASSES on concentrated-signal mechanisms (cipro QRDR + cef plasmid β-lactamases) and FAILS on distributed mobile-element mechanisms (tet efflux + ribosomal protection: AUROC 0.400 anti-predictive). Treating multi-drug as a drug-list expansion would inherit the architectural failure.

**Trade-off:** Adds an architecture decision packet (EP-1.5) before EP-2 can ship for distributed-mechanism drugs. Cef can still ship as EP-2a using current architecture; tet/gent waits for EP-1.5.

### D4: Honest-output framework generalization is its own EP (EP-1C), not background work

**Decision:** Lift the audit verdict + SUSPEND gate + `attribution_scope_confidence` + INDETERMINATE propagation from cipro-shaped to drug-agnostic + organism-agnostic, in its own EP.

**Rationale:** This is the project's actual differentiator vs existing AMR tools (Galaxy provenance-as-feature analog). Doing it as an explicit EP forces consistency across all future EPs; doing it as background work risks per-EP variation.

**Trade-off:** Adds an EP that has no new prediction capability. Net win because every downstream EP inherits the framework for free.

### D5: External benchmark is its own EP (EP-1B), not optional

**Decision:** EP-1B runs v0 against a small public AST panel (e.g., NARMS / PATRIC subset) side-by-side with AMRFinderPlus + RGI, before broader expansion.

**Rationale:** v0's moat is honest output, but no external comparison has been done. Without EP-1B, all moat claims are internal. Cheap (1 small public panel) + high information.

**Trade-off:** Delays multi-drug expansion by 1-2 weeks. Acceptable — protects against building EP-2 + EP-3 on an unvalidated v0.

### D6: EP-4+ (non-AMR / eukaryotic / multimodal) require separate /idea-anchor + /project-init cycles

**Decision:** Do not commit non-AMR or eukaryotic phenotypes as roadmap items at this ladder. Each requires its own `/idea-anchor + /project-init` invocation when a paired dataset materializes.

**Rationale:** They are research-program-shaped at the moment (no concrete dataset, no terminal claim). Including them in the committed ladder would re-introduce the unbounded-scope failure.

**Trade-off:** The "DNA → phenotype tool at maturity" framing implies these will happen eventually. Mitigation: track as horizon options in `plans/v1_Horizon_Framing_Plan.md` (already drafted), not as commitments here.

### D7: "Commercialize" is removed from the ladder entirely

**Decision:** Commercialize is a downstream OUTCOME (if/when the tool actually works on a useful slice), not a development phase. Not numbered, not committed.

**Rationale:** A development phase changes technical capability; commercialization changes distribution / support / validation / legal posture. The two are different shapes of work. The user clarified 2026-05-25 that "sell to a large company" is a possible downstream outcome, NOT the framing for current dev.

**Trade-off:** None — eliminates the over-weighting that the misframed brainstorm chain introduced.

---

## Implementation Plan

### EP-0 close — Reproducibility freeze (BLOCKER — happens FIRST)

**Owner:** Codex on Precision 7780 + Claude on GTX 860M (sync verification).

**Terminal claim:** v0.0-cipro is a fully reproducible reference point. Pulling origin HEAD + checking out tag `v0.0-cipro` reproduces the exact model + DBs + Docker images + example output that was shipped 2026-05-24.

**Sub-actions:**
1. Codex on Precision 7780: `git push origin main` with the 5 outstanding artifacts (retrained `ciprofloxacin_nucleotide_transformer.pkl` with `leave_one_accession_out` CV / runtime `pipeline.py` changes / `reports/cipro_v0_scope_limit_decision_2026-05-23.md` / `reports/dna_decoder_v0_release_candidate_2026-05-24.{md,json}` / `wiki/cipro_bounded_falsifier_results_2026-05-23.{md,json}`).
2. Claude on GTX 860M: `git pull --ff-only origin main`; run `scripts/cross_machine_sync_check.py`; confirm 0 spec-divergence + 0 Downloads/-recent artifacts (or all resolved).
3. Tag `v0.0-cipro` on the post-pull commit. Tag locally only; don't push the tag until user confirms.
4. Lock the AMRFinder DB version (already pinned in CLAUDE.md to `ncbi/amr:4.2.7-2026-03-24.1`); document Bakta DB version in `wiki/v0_reproducibility_freeze_<DATE>.md`.

**Success criteria:**
- `scripts/cross_machine_sync_check.py` reports IN SYNC.
- `git tag --list 'v0.0-cipro'` returns the tag.
- A new `wiki/v0_reproducibility_freeze_<DATE>.md` document lists exact model file + DBs + Docker images + Python version + key environment.

**Out-of-scope:** anything beyond making v0 reproducible. EP-0 is hygiene; the project's actual capability does not change.

---

### EP-1A — Genome-input contract

**Owner:** Codex on Precision 7780 (compute) + Claude on GTX 860M (planning + non-GPU code).

**Terminal claim:** `scripts/pipeline.py predict --genome-fasta <novel-E-coli-FASTA> --drug ciprofloxacin ...` runs end-to-end on a public held-out E. coli genome (not in the N=147 training cohort) and emits the v0 JSON+MD schema. Same-strain parity test: when run with the FASTA of a cohort strain, the output matches the cached-strain prediction within tolerance ε (TBD, candidate ε = 0.01 on calibrated probability).

**Sub-actions:**
1. Choose the public held-out E. coli genome — candidate sources: NCBI Pathogen Detection (recent BV-BRC-listed strain with AST); GenomeTrakr; NARMS isolate.
2. Define ε for the same-strain parity test.
3. Implement `dna_decode/data/genome_ingest.py` — FASTA loader + CDS extraction + transient cache.
4. Extend `scripts/pipeline.py predict` to accept `--genome-fasta` + `--annotations` (parallel to `--strain-id`).
5. Run the same-strain parity test against 3 cohort strains; record results.
6. Run the public-held-out genome through the new path; record results.
7. Document the EP in `wiki/ep1a_genome_input_contract_<DATE>.md`.

**Success criteria:**
- Same-strain parity ≤ ε on all 3 cohort strains.
- 1 public held-out E. coli genome runs to completion with all v0 JSON fields populated.
- `attribution_scope_confidence` field renders correctly (PARTIAL or INDETERMINATE for novel input; HIGH only for ERS-prefix strains per locus-tag-prefix proxy).

**Out-of-scope:** multi-drug; Bakta-on-the-fly annotation; persistence of novel-genome embeddings to the canonical cache.

---

### EP-1B — External benchmark smoke

**Owner:** Codex on Precision 7780 (runs benchmarks) + Claude on GTX 860M (table rendering).

**Terminal claim:** A side-by-side comparison of v0 + AMRFinderPlus + RGI on a small public AST-labeled panel (suggested N = 10 E. coli genomes from NARMS or PATRIC with known cipro phenotype) exists in `reports/external_benchmark_cipro_<DATE>.{md,json}`. The comparison surfaces whether v0's calibrated probability + audit verdict carry meaningful additional information vs AMRFinderPlus's gene call.

**Sub-actions:**
1. Select 10 public E. coli genomes with known cipro AST labels. Diversity goal: at least 3 with known QRDR mutations + 3 plasmid-carrying + 3 wild-type.
2. Run v0 (genome-input via EP-1A) + AMRFinderPlus + RGI on all 10.
3. Render side-by-side: cohort label / v0 prediction + probability / AMRFinderPlus mechanism calls / RGI call.
4. Compute: agreement matrix; v0-additional-info score (does v0's calibrated probability or audit verdict provide signal AMRFinderPlus doesn't?).
5. Document in `reports/external_benchmark_cipro_<DATE>.{md,json}`.

**Success criteria:**
- 10/10 genomes complete without errors.
- v0 prediction agrees with AST label on at least 7/10 (≈ 0.70 accuracy floor matching v0 spec).
- The side-by-side renders a "what does v0 add vs AMRFinderPlus" answer one way or the other.

**Out-of-scope:** clinical-grade validation; multi-drug benchmark.

---

### EP-1C — Honest-output framework generalization

**Owner:** Claude on GTX 860M (mostly refactor; can run locally).

**Terminal claim:** The audit verdict + `attribution_scope_confidence` + SUSPEND propagation + INDETERMINATE handling all operate via drug-agnostic + organism-agnostic interfaces. Adding a new drug or organism does NOT require touching `scripts/pipeline.py` predict logic; it requires adding to `dna_decode/data/mic_tiers.py` (already drug-agnostic).

**Sub-actions:**
1. Audit current `scripts/pipeline.py predict` for cipro-shaped assumptions. Specifically: locus-tag-prefix proxy constants (`_ATTRIBUTION_ERS_CONTROL_PREFIXES`, `_ATTRIBUTION_ELX_FAMILY_PREFIXES`) — currently hardcoded from cipro audit; need a generic prefix-family registry.
2. Refactor `_classify_attribution_scope` to accept a drug parameter + look up prefix families from `mic_tiers.py` (extending the catalog as needed).
3. Update tests to exercise cef + tet + gent on the helper (even if no real model exists yet).
4. Document the generalization in `wiki/ep1c_honest_output_framework_<DATE>.md`.

**Success criteria:**
- `_classify_attribution_scope` no longer references cipro-specific constants directly; all per-drug data lives in `mic_tiers.py`.
- New tests for cef + tet + gent pass.
- v0 cipro behavior unchanged (regression-safe).

**Out-of-scope:** adding any new drug's actual mechanism catalog (that's EP-2 territory).

---

### EP-1.5 — Architecture decision packet

**Owner:** Claude on GTX 860M (decision-only) + Codex on Precision 7780 (any small-N proof-of-concept).

**Terminal claim:** A written decision exists at `plans/EP1_5_Architecture_Decision_Plan.md` choosing the architectural fix for distributed-mechanism resistance (tet, gent). Options: per-gene NT windows + attention pooling; k-mer + AMRFinder-feature fusion; lightweight transformer head on top of pooled NT. The decision is grounded in a small-N proof-of-concept (e.g., re-fire tet 12-strain smoke with each candidate; pick the one with highest AUROC).

**Sub-actions:**
1. Read `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md` to re-anchor.
2. Run 3 candidate architectures on the tet 12-strain smoke (already populated cache, low cost).
3. Pick the candidate with highest tet AUROC AND that doesn't regress cipro/cef AUROC.
4. Document the decision packet.

**Success criteria:**
- 1 candidate architecture beats current mean-pool on tet AUROC by ≥ 10 pp (current 0.400 → ≥ 0.500).
- Same architecture matches or beats current on cipro (≥ 0.75) + cef (≥ 0.80).
- Decision documented + linked from this plan.

**Out-of-scope:** full rework. EP-1.5 is decision-only; implementation lives in EP-2/3.

---

### EP-2 — Multi-drug AMR for E. coli (BY MECHANISM CLASS)

**Owner:** Codex on Precision 7780 (cohort + model) + Claude on GTX 860M (planning + tests).

**Terminal claim:** Per-drug models exist for cef (concentrated mechanism, reuses v0 arch) AND for tet (distributed mechanism, uses EP-1.5 architecture). Both predict via `scripts/pipeline.py predict --drug <cef|tet>`. CV AUROC ≥ 0.70 per drug.

**Sub-actions:**
1. **Cef sub-track:** BV-BRC categorical-MIC feasibility check; build cef cohort; populate cef NT cache on Precision 7780; train cef classifier; LOSO CV.
2. **Tet sub-track:** same BV-BRC feasibility check; if cohort feasible, build + populate + train using EP-1.5's chosen architecture; LOSO CV.
3. Validate honest-output framework (EP-1C) works for both drugs without code change.

**Success criteria:**
- Cef CV AUROC ≥ 0.70.
- Tet CV AUROC ≥ 0.70 OR documented INDETERMINATE with explicit reason (cohort infeasibility / architecture still insufficient).
- Per-drug provenance + audit verdict propagate correctly.

**Out-of-scope:** multi-organism; non-AMR phenotypes.

---

### EP-2.5 — Cohort acquisition automation

**Owner:** Claude on GTX 860M.

**Terminal claim:** `scripts/cohort_discovery.py` (new) automates: BV-BRC AST + assembly availability join; per-drug per-class label tiering (strict-MIC / categorical-MIC / decisive); per-class balancing; output cohort parquet. No more manual per-drug feasibility checks.

**Sub-actions:**
1. Generalize existing `scripts/build_stage2_n150_cohort.py` (cipro-hardcoded) → `scripts/cohort_discovery.py` (drug-agnostic).
2. Reuse `dna_decode/data/mic_tiers.py` per-drug catalogs.
3. CLI: `--drug X --target-total N --label-quality {strict|categorical|decisive}`.
4. Run against cef + tet + gent + colistin (4th-mechanism-class falsifier) for parallel cohort discovery.

**Success criteria:**
- 1 CLI invocation builds a working cohort parquet for ≥ 3 drugs.
- Per-drug feasibility report at the same time.

**Out-of-scope:** anything beyond E. coli.

---

### EP-3 — Multi-organism AMR

**Owner:** TBD.

**Terminal claim:** Per-organism models exist for at least Klebsiella pneumoniae (next priority after E. coli per the long-term vision in CLAUDE.md L9). Same predict CLI + JSON schema.

**Out-of-scope until E. coli ladder is stable:** Pseudomonas, Acinetobacter, mycobacteria, etc.

---

### EP-4+ — Non-AMR + eukaryotic + multimodal

**These are NOT roadmap commitments.** Each requires its own `/idea-anchor + /project-init` cycle when a paired dataset exists. Tracked in `plans/v1_Horizon_Framing_Plan.md`.

---

## Execution Preview

```
Wave 0 (BLOCKER):     EP-0 close (Codex push → sync verify → tag v0.0-cipro)
Wave 1 (1 EP):        EP-1A (genome-input contract)
Wave 2 (3 parallel):  EP-1B (external benchmark), EP-1C (honest-output generalization), EP-1.5 (architecture decision)
Wave 3 (1 EP):        EP-2 (multi-drug by mechanism class) — cef + tet sub-tracks
Wave 4 (1 EP):        EP-2.5 (cohort acquisition automation)
Wave 5 (1 EP):        EP-3 (multi-organism AMR)
```

**Critical path:** EP-0 → EP-1A → EP-1.5 → EP-2 (tet sub-track) → EP-3.

**Parallelism:** Max 3 EPs at Wave 2 (EP-1B, EP-1C, EP-1.5 are independent).

**Wall-clock estimate (rough):**
- EP-0 close: 1 day (gated on Codex push)
- EP-1A: 3-5 days post-EP-0
- EP-1B + EP-1C + EP-1.5 in parallel: 1-2 weeks elapsed
- EP-2: 2-4 weeks (per-drug cohort populate is the bottleneck)
- EP-2.5: 3-5 days
- EP-3: 4-8 weeks per new organism

---

## Risk Flags

- **R1 (HIGH):** EP-0 close is blocked on Codex push from Precision 7780. If Codex doesn't push for > 1 week, the entire ladder stalls. Mitigation: this side keeps running planning + non-GPU work in parallel; sync check runs daily until push lands.
- **R2 (MEDIUM):** EP-1.5 architecture decision could surface that NO simple fix works for distributed mechanisms within solo capacity. Mitigation: if all 3 candidate architectures fail, EP-2 tet sub-track is scoped to "documented INDETERMINATE" rather than "shipped with ≥ 0.70 AUROC."
- **R3 (MEDIUM):** Per-drug cohort feasibility (the 2026-05-18 BV-BRC census problem) may block EP-2 even after EP-2.5 automation. Mitigation: EP-2.5 surfaces feasibility per-drug; EP-2 explicitly skips infeasible drugs.
- **R4 (LOW):** Cross-machine sync drift recurs. Mitigation: `scripts/cross_machine_sync_check.py` runs at every EP boundary.
- **R5 (LOW):** EP-1A same-strain parity test fails (cached-strain prediction ≠ FASTA-path prediction for the same strain). Likely root causes: CDS extraction differences; mean-pool ordering; transient cache float drift. Mitigation: ε tolerance set generously (0.01 calibrated probability); if even that fails, investigate ingest path before declaring EP-1A done.
- **R6 (Behavioral):** Analysis-fatigue is a real risk after 3 brainstorms + 1 probe + 1 review in one session. This plan is the convergent recommendation across all rounds; don't re-brainstorm before executing EP-0 close.

---

## Verification

After this plan is committed + pushed:

1. `git log --oneline origin/main -1` shows this plan's commit.
2. `wiki/plans-index.md` has an entry for `plans/Post_V0_EP_Ladder_Plan.md`.
3. `project_state/dna-decode-2026-05-11.md` Bellman frame Current State references this plan.
4. `scripts/cross_machine_sync_check.py` continues to flag the 5 outstanding Codex artifacts as drift (no change here on this side; EP-0 close requires Codex push to resolve).
5. Future sessions can read this plan + know where to pick up without re-running the brainstorm chain.

---

## What this plan deliberately does NOT cover

- **Commercialization.** Removed per D7. Track as a downstream outcome, NOT a development phase.
- **/idea-anchor + /project-init invocation on the long-term goal.** Both belong AFTER EP-1B's external benchmark produces evidence; running them earlier would be premature commitment + would likely fail the same gates the 2026-05-11 init failed.
- **Customer-discovery calls.** Was over-weighted in the misframed brainstorm chain; not in scope for a tool-building plan.
- **EP-4+ horizon items.** Tracked in `plans/v1_Horizon_Framing_Plan.md` as horizon options.
- **Concrete dates.** Solo developer with day job — wall-clock estimates are rough; user picks pace.

---

## Bottom line

Single most-important next action: **Codex on Precision 7780 pushes the 5 outstanding v0-closeout artifacts to origin.** Everything else is gated on that.

After Codex pushes: EP-0 close → tag → EP-1A. The plan covers the path from there onward.
