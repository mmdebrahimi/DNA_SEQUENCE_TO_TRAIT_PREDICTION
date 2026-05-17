# EP1 + EP2 Cross-Drug Synthesis Plan

> Single-deliverable plan to write the cross-drug architectural-finding synthesis from EP1 (cipro) + EP2 (cef + tet). Closes Phase 1 evidence collection internally. External publication remains deferred per EP1 closeout discipline.

---

## Problem Statement

The 3-experiment evidence chain across cipro + cef + tet is complete. EP1 cipro was internally closed 2026-05-17 as an adversarial audit/infrastructure packet. EP2 cef + tet smoke fired 2026-05-17 and produced the cross-drug architectural finding (cef PASS at 0.833 = k-mer 0.833; tet FAIL at 0.400 anti-predictive vs k-mer 0.722; H17 falsified).

Per /brainstorm round-3 framing critique on `plans/Cef_Mechanism_Audit_Plan.md`: **the cross-drug architectural claim is already on paper** in `wiki/EP2_cef_tet_verdict_2026-05-17.md`. Further mechanism-audit packets are corroborating, not decision-bearing. The project's near-term value is closing a coherent research packet, not maximizing additional diagnostics.

This plan writes the synthesis. Single deliverable. Listed residual uncertainty for any future EP. Defines the stopping rule for Phase 1 evidence collection.

## Design Decisions

### D1: One deliverable; no new code

**Decision:** Single artifact = `wiki/ep1_ep2_cross_drug_architectural_finding_<date>.md`. No new scripts, no new tests, no new modules. The synthesis cites + integrates existing artifacts.

**Rationale:** /brainstorm round 3 caught that "more evidence packets ≠ better" once the theory of the case is stable. Synthesis writing is the right next move; framework-expansion isn't.

### D2: Internal scope only

**Decision:** This synthesis is INTERNAL. Not arXiv. Not blog. PC1 (Phase 1 framing) remains `internal_closeout` per the EP1 closeout's pre-conditions artifact.

**Rationale:** EP1 closeout already locked `internal_closeout`. External publication requires either (a) BV-BRC-wide cohort rebuild for at least one drug or (b) more drugs' worth of evidence. Neither is in scope; both are heavy work for marginal scientific gain at the current evidence level.

### D3: Stopping rule defined

**Decision:** This synthesis closes Phase 1 evidence collection. No further audits + no additional mechanism work + no per-gene NT diagnostics + no Stage 1 N=40 escalations on cipro/cef/tet UNLESS the synthesis surfaces a specific contradiction or open hypothesis worth investigating.

**Rationale:** /brainstorm round 3 surfaced "without a stopping rule, every corroborating audit can spawn another audit." The synthesis is the stopping point. Phase 2 work (other drugs, larger cohorts, multimodal, etc.) is a separate strategic decision that requires its own /idea-anchor + /project-init cycle.

### D4: Residual uncertainty explicitly listed

**Decision:** The synthesis has a dedicated "Residual uncertainty" section listing what was NOT resolved by EP1 + EP2:
- Cef mechanism question (N=12 mini-cohort; was the 0.833 mechanism-driven or lineage-tracking? — UNRESOLVED)
- Tet failure mode disambiguation (architectural mismatch vs N=12 calibration vs label noise — UNRESOLVED)
- Cipro Stage 1 N=40 architectural verdict (cohort label noise was the dominant confound — partial diagnosis, not full architectural test)
- BV-BRC-wide cef cohort feasibility (NOT YET CHECKED)
- Per-gene NT windows on CLEAN_R cipro strains (DEFERRED architectural diagnostic)
- Bakta cohort-wide annotation (DEFERRED)
- Mean+max preflight v3 closeout falsifier (DEFERRED conditional)

These items become Phase 2 candidates if Phase 2 fires.

### D5: Cef cross-tab artifact (if produced via `plans/Cef_Mechanism_Audit_Plan.md` Step 3) cited as corroborating evidence only

**Decision:** If the user fires the cef NT-vs-mechanism cross-tab from the reduced cef plan, the synthesis cites it in a "Corroborating Diagnostic" subsection — NOT as part of the load-bearing claim. The synthesis stands whether or not the cef cross-tab fires.

**Rationale:** Round-3 framing critique: "if the synthesis is defensible without cef cross-tab, the cef audit collapses to 'deferred residual uncertainty' entirely." Honor that — synthesis must be writable + defensible WITHOUT the cef cross-tab.

## Implementation Plan

### Step 1: Write the synthesis
Files: wiki/ep1_ep2_cross_drug_architectural_finding_<date>.md (new)
Depends on: none

**What changes:**
- NEW file `wiki/ep1_ep2_cross_drug_architectural_finding_<date>.md`.
- Sections:
  1. **Executive summary (1 paragraph).** The cross-drug architectural finding in one sentence: "Frozen-NT-whole-genome-pooling localizes concentrated-signal AMR mechanisms (QRDR point mutations + plasmid acquired-gene β-lactamases) at smoke fidelity but fails on distributed mobile-element mechanisms (tet efflux + ribosomal protection). The architecture's failure mode is data-shape-dependent, not drug-or-cohort-dependent in isolation."
  2. **Evidence chain table.** 3-experiment summary:
     | Drug | Cohort | NT-XGBoost AUROC | Mechanism class | Verdict |
     |---|---|---:|---|---|
     | Cipro (smoke N=12) | Cipro mini | 0.750 | QRDR point mutations | PASS (2026-05-14) |
     | Cipro (N=38 Stage 1) | Cipro N=38 | 0.568-0.615 | QRDR + label noise | FAIL → EP1 audit closeout |
     | Cef (smoke N=12) | Cef mini | **0.833** | plasmid β-lactamases | **PASS** |
     | Tet (smoke N=12) | Tet mini | **0.400** | tet-family efflux + ribosomal protection | **FAIL anti-predictive** |
  3. **Architectural pattern.** The "concentrated-signal vs distributed-mobile-element" mechanism class division.
  4. **EP1 cipro audit-infrastructure findings.** Summarize the 6-experiment chain (Stage 1 + 1b + preflight v2 + Experiment 1 mechanism + Experiment 2 MIC + merge SUSPEND_CONDITION_4) + the 4-tier audit infrastructure shipped (mechanism × MIC merge + opacity flag + structurally-enforced gate).
  5. **EP2 H17 falsification.** The cef-PASS / tet-FAIL split.
  6. **Residual uncertainty (per D4).** Enumerated list.
  7. **Architectural implications.** What this finding means for production AMR-prediction pipelines: per-drug architecture selection or targeted per-gene windows for mobile-element-dominant resistance.
  8. **Phase 1 stopping point.** Declare evidence collection complete.
  9. **Corroborating diagnostic (OPTIONAL section).** If the cef cross-tab from `plans/Cef_Mechanism_Audit_Plan.md` Step 3 has fired, cite its findings as corroborating evidence — does NT-XGBoost's per-strain prediction agree with AMRFinder β-lactamase detection? Does NOT change the synthesis verdict.
  10. **Lessons (cross-cutting).** Reusable lessons for future EPs:
      - Pre-conditions discipline (PC1/PC2 lock) catches statistical bugs (binomial-threshold-below-uniform-error-null in cipro).
      - Mechanism × phenotype merge with opacity flag prevents "labels are wrong" / "tool is wrong" conflation.
      - Smoke-tier infrastructure (N=12, calibrate=False, k-mer baseline) reusable across drugs with cosmetic patches.
      - Anti-predictive AUROC at N=12 with non-calibrated classifier is data-shape divergence, NOT a plumbing bug.
      - Plan-language drift (e.g., "smaller alternative" preserving prior reduction authority while authorizing similar expansion) is a framing trap; require LIT in-place edits BEFORE /execute-plan.
- Cross-references:
  - `wiki/cipro_ep1_closeout_2026-05-17.md`
  - `wiki/cipro_decision_bundle_pre_conditions_2026-05-17.md`
  - `wiki/EP2_cef_tet_verdict_2026-05-17.md`
  - `wiki/smoke_gate_12strain_cipro_2026-05-14.md`
  - `wiki/smoke_gate_12strain_ceftriaxone_2026-05-17.md`
  - `wiki/smoke_gate_12strain_tetracycline_2026-05-17.md`
  - `wiki/stage1_n40_cipro_2026-05-15.md`
  - `wiki/stage1_n40_cipro_mean-plus-max_2026-05-16.md`
  - `wiki/cipro_attribution_preflight_2026-05-16.md`
  - `wiki/cipro_mechanism_audit_2026-05-17.md`
  - `wiki/cipro_mic_audit_2026-05-17.md`
  - `wiki/cipro_mechanism_phenotype_audit_2026-05-17.md`

**Verification:** The synthesis stands on its own as an internal closeout artifact. Reader unfamiliar with the project's intermediate evidence can navigate from the synthesis to the source packets.

### Step 2: Update project ledger
Files: project_state/dna-decode-2026-05-11.md
Depends on: Step 1

**What changes:**
- Append Action Log row recording EP1+EP2 synthesis written + Phase 1 evidence collection closed.
- Refresh Bellman current-state to "Phase 1 evidence collection CLOSED; cross-drug architectural finding written internally; external publication deferred; Phase 2 strategic decision queued (not yet a planned action)."
- Update H17 entry: status remains `falsified` (no change); the EP2 verdict packet + this synthesis are now the audit trail.
- Mark Phase 1 EP1 + EP2 as `CLOSED_INTERNAL` in the Mid-term EP table.

**Verification:** `grep -E "EP1.*closed\|Phase 1 evidence.*closed" project_state/dna-decode-2026-05-11.md` returns the closure entries.

### Step 3 (optional): Decide whether to fire the cef cross-tab
Files: (decision only)
Depends on: Step 1

**What changes:**
- AFTER Step 1's synthesis is written, re-read it.
- If the synthesis is defensible without the cef cross-tab corroboration → mark `plans/Cef_Mechanism_Audit_Plan.md` as DEFERRED (don't fire); update its plans-index entry.
- If the synthesis explicitly relies on the cef cross-tab to ground a claim → fire `plans/Cef_Mechanism_Audit_Plan.md` as a follow-up.

**Verification:** decision recorded in project ledger as an Action Log row.

## Execution Preview

```
Wave 0 (1 sequential): Step 1 — write synthesis
Wave 1 (1 sequential): Step 2 — update project ledger
Wave 2 (1 sequential, optional): Step 3 — decide on cef cross-tab firing

Critical path: Step 1 → Step 2 → (optional Step 3) (2-3 waves)
Max parallelism: 1 agent (sequential writing work)
```

This is intentionally a write-mode plan, NOT a build-mode plan. /execute-plan in parallel mode is overkill; sequential is the right shape.

## Verification (end-to-end)

After Step 1 + Step 2:

- `wiki/ep1_ep2_cross_drug_architectural_finding_<date>.md` exists with all 10 required sections.
- Synthesis cross-references 10+ source artifacts.
- Project ledger Action Log has the closeout entry + Bellman current-state refreshed.
- Phase 1 evidence collection has a defined stopping point.

After Step 3 (if fired):
- Either `plans/Cef_Mechanism_Audit_Plan.md` is marked DEFERRED in the plans-index, OR the cef cross-tab Step 3 has been executed and its findings cited in the synthesis's optional "Corroborating Diagnostic" subsection.

## Out of scope

- External publication (arXiv / blog) — deferred per EP1 closeout PC1=`internal_closeout`.
- Stage 1 N=40 cipro re-runs — Stage 1 + 1b results stand; no more cipro experiments.
- Stage 1 N=40 cef rebuild — would be Phase 2.
- Stage 1 N=40 tet rebuild — would be Phase 2.
- Per-gene NT windows on CLEAN_R cipro strains — Phase 2 architectural diagnostic.
- Bakta cohort-wide annotation — Phase 2.
- BV-BRC-wide cef cohort feasibility census — Phase 2.
- 2nd-organism portability, MIC continuous head, pan-genome graph, multimodal — parked horizon tracks per framing reset.
