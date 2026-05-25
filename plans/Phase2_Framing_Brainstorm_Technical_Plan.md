# Phase 2 Framing Brainstorm — Technical Plan

> Apply 5 /review edits to `plans/Phase2_Framing_Brainstorm_Plan.md` to convert the document from "framing brainstorm with implicit preference" → "framing brainstorm with explicit 3-candidate slate + decision gates." User-confirmed scope: add classifier-tier Candidate 3; Candidate 1 merge sentence permits Months-including-Databricks budget.

---

## Problem Statement

The /review of `plans/Phase2_Framing_Brainstorm_Plan.md` surfaced 5 concrete edits that align the document with its actual function:

1. **`(preferred)` label pre-decides what `/idea-anchor` should decide.** Strip it.
2. **Candidate 1 oversells current state.** Audit scripts are cipro-hardcoded at module-constant level (CIPRO_LOCI_BY_MECHANISM, CLSI cipro breakpoints, frozen STAGE1B_*_AUROC literals); "audit-first AMR evidence-packet system" implies a multi-drug factory that requires a generalization refactor. Merge with Candidate 2's verb + permit full-budget scope.
3. **Slate is incomplete.** No classifier-tier candidate; user picked classifier-tier (post-feasibility-census A1) as the third candidate.
4. **Open Questions Q1-Q3 are footnote, not gate.** They determine which candidate gets selected; ordering must be explicit.
5. **"Honest Read on User Intent" is encoded as premise, not hypothesis.** Reframe.

Deliverable = edited `plans/Phase2_Framing_Brainstorm_Plan.md` + refreshed entry in `wiki/plans-index.md`. No code changes.

## Design Decisions

### D1: Documentation-only plan

**Decision:** No code changes. Edits target two markdown files only.

**Rationale:** The framing brainstorm's deliverable is candidate anchor sentences; review feedback is applied to the deliverable, not to source code. Downstream code work happens after `/idea-anchor` + `/project-init` fire in a separate session.

### D2: Third candidate is classifier-tier, not exit-tier

**Decision:** Candidate 3 = "After the BV-BRC strict-MIC feasibility census proves drug-X has ≥150 clean labels, ship a single deployable classifier for drug-X with full LOMO-clade-out validation, using whichever architecture the audit framework selects."

**Rationale:** User-confirmed via AskUserQuestion. Surfaces the deferred ship-a-deployable-classifier path so /idea-anchor sees the real tradeoff space (audit-tier ×2 + classifier-tier ×1).

**Trade-off:** Exit-tier ("wind down dna-decode") was offered + rejected by user; classifier-tier is the chosen hedge.

### D3: Candidate 1 merge sentence permits full-budget scope

**Decision:** Candidate 1 body explicitly names "generalizing cipro-hardcoded mechanism catalogs + breakpoints + Stage 1b reference baselines" + "audit-gated cross-drug experiments (cipro + cef + tet + a 4th-mechanism-class smoke) on BV-BRC-wide cohorts" + "before committing each drug to Databricks burst + N=150 NT cache populate."

**Rationale:** User-confirmed budget = "Months including Databricks (full)" via AskUserQuestion. Candidate 1 honors that without committing every drug to Databricks simultaneously (audit gates fire per-drug first).

**Trade-off:** 1-2 weeks minimum (Tradeoff B lower bound) was offered + rejected. Months-extended without Databricks was offered + rejected.

### D4: Strip `(preferred)` label; framing brainstorm exposes choice rather than pre-ranks

**Decision:** Remove the `(preferred)` annotation on Candidate 1 + delete the "Why Candidate 1 preferred" rationale paragraph. Replace with neutral "How the 3 candidates differ" + explicit "No pre-ranking" note.

**Rationale:** /idea-anchor is where selection happens; encoding a preference in the brainstorm pre-commits the project to one candidate before Open Questions are resolved.

### D5: Open Questions Q1-Q3 promoted to gate

**Decision:** Section heading changes to "## Open Questions — gate /idea-anchor on these"; intro paragraph names them as required-resolved-before-anchor + annotates each with which candidate dimension it affects.

**Rationale:** Q1-Q3 (external artifact type, portability requirement, replication threshold) determine which candidate gets selected. As footnote, they were skippable; as gate, they discipline the anchor selection.

### D6: Honest Read on User Intent reframed as hypothesis

**Decision:** Prepend a sentence noting the read is "hypothesis worth confirming at /idea-anchor input, not a precommitment."

**Rationale:** Codifying the read as candidate-preference shaping = self-fulfilling loop. Surfacing as a check lets the user confirm, refine, or reject.

## Codebase Context

**Files affected:**
- `plans/Phase2_Framing_Brainstorm_Plan.md` (exists; written earlier this session)
- `wiki/plans-index.md` (exists; current entry inserted at top)

**No source code changes.**

**Verified facts from /review (preserved as constraints for edit content):**
- 4 audit scripts (`scripts/cipro_mechanism_audit.py`, `cipro_mic_audit.py`, `cipro_mechanism_phenotype_merge.py`, `cipro_curated_baseline.py`) cipro-hardcoded at module-constant level (CLAUDE.md gotchas + eng review)
- Smoke runner (`scripts/smoke_gate_12strain_cipro.py`) drug-parameterized via `--drug`
- N=147 cipro Stage 2 cohort exists at `data/processed/stage2_n150_cipro_cohort.parquet`; NT cache NOT populated (only `nt_n40_cipro.h5` at 456 MB)
- AMRFinder / Bakta / Mash Docker toolchain installed + smoke-validated 2026-05-15

## Implementation Plan

### Step 1: Apply 5 review edits to Phase2_Framing_Brainstorm_Plan.md

**Files:** `plans/Phase2_Framing_Brainstorm_Plan.md`
**Depends on:** none

**What changes:**

- **Edit 1 — strip `(preferred)` label.** Heading `**Candidate 1 (preferred):**` → `**Candidate 1:**`. Delete the prose paragraph beginning `**Why Candidate 1 preferred:** absorbs A2 (audit-as-product)...`.

- **Edit 2 — replace Candidate 1 body with full-budget merge sentence.** Replace the existing Candidate 1 blockquote with:
  > Validate and operationalize the Phase 1 mechanism-class architectural finding by extending the audit infrastructure across drug families — generalizing the cipro-hardcoded mechanism catalogs + breakpoints + Stage 1b reference baselines into a per-drug configuration layer — and running audit-gated cross-drug experiments (cipro + cef + tet + a 4th-mechanism-class smoke) on BV-BRC-wide cohorts, using strict-MIC, mechanism, opacity, and smoke-model gates to determine whether frozen NT whole-genome pooling is appropriate per drug × mechanism class before committing each drug to Databricks burst + N=150 NT cache populate.

- **Edit 3 — add Candidate 3 (classifier-tier, post-census).** Append below Candidate 2:
  > **Candidate 3 (classifier-tier, post-feasibility-census):**
  > > After the BV-BRC strict-MIC feasibility census proves drug-X has ≥150 clean labels, ship a single deployable classifier for drug-X with full LOMO-clade-out validation, using whichever architecture (frozen NT whole-genome, per-locus NT windows, or mechanism-only curated baseline) the audit framework selects per the EP1+EP2 cross-drug architectural finding.

- **Edit 4 — replace `Why Candidate 1 preferred:` rationale with slate-framing paragraph.** Add after Candidate 3:
  > **How the 3 candidates differ:**
  > - **C1 (audit-tier, forward-looking):** anchors on building the multi-drug audit factory; classifier-shipping is a downstream artifact of audit gates passing.
  > - **C2 (audit-tier, backward-looking):** anchors on validating Phase 1's specific finding via replication; lower scope, narrower claim.
  > - **C3 (classifier-tier, post-census):** anchors on shipping a deployable model; conditional on feasibility census proving labels exist; commits to Databricks burst for whichever drug clears the census.
  >
  > **No pre-ranking.** /idea-anchor is where the selection happens; this slate exposes the choice rather than encoding a preference.

- **Edit 5 — promote Open Questions to gate.** Replace section heading + intro with:
  > ## Open Questions — gate /idea-anchor on these
  >
  > These 3 questions must be answered before `/idea-anchor` selects a candidate. The anchor sentence will encode implicit answers to all three regardless; making them explicit prevents the project from inheriting unstated constraints.
  >
  > 1. Is the desired external artifact a paper/blog, GitHub tool, or portfolio case study? *(Affects D1/D2 publication threshold + C1 toolkit packaging.)*
  > 2. Strictly E. coli, or is portability a required success criterion? *(Affects whether audit framework refactor goes per-drug-only or per-drug × per-organism.)*
  > 3. Minimum replication for calling the audit framework real: cipro + one more drug, or cipro + cef + tet + 4th class? *(Affects Candidate 1 scope width.)*

- **Edit 6 — reframe Honest Read on User Intent as hypothesis.** Prepend after the section heading:
  > **Hypothesis worth confirming at /idea-anchor input, not a precommitment.** Codifying this read as candidate-preference shaping would be a self-fulfilling loop; surfacing it as a check for the user lets it be confirmed, refined, or rejected.

  Keep existing "Most likely" + "Second-order" bullets unchanged.

**Key details:**
- All edits are local markdown surgery. No section reorders, no body deletions outside what's explicitly named.
- Preserve existing tier structure (Product/Research/Horizon/Meta), axis critique, missing anchors, cross-domain analogs, sub-problem decomposition, tradeoffs, assumptions sections.
- Date stamps + "Next Step" footer stay; "Next Step" still says fresh /idea-anchor + /project-init.

### Step 2: Update plans-index.md entry to reflect 3-candidate slate

**Files:** `wiki/plans-index.md`
**Depends on:** Step 1

**What changes:**

Locate `## [plan_file: Phase2_Framing_Brainstorm_Plan.md] 2026-05-17` entry. Replace its `**Summary:**` + `**Key decisions:**` body with:

```
**Summary:** Pre-/idea-anchor ideation pass capturing the Phase 2 anchor search space, axis critique (A1-D3), missing anchors, and 3 candidate /idea-anchor sentences (audit-tier forward, audit-tier backward, classifier-tier post-census) after Phase 1 closeout 2026-05-17. Revised post-/review 2026-05-17 to strip preference pre-ranking + add classifier-tier candidate + gate /idea-anchor on Open Questions Q1-Q3.
**Key decisions:**
- D1: Axis menu structured by tier (Product / Research / Horizon / Meta), not by topic
- D2: Generative-ideation v2.1 payload pattern — critique + generation in one round
- D3: Stopping discipline preserved — no Phase 2 experiment fires; fresh /idea-anchor + /project-init required
- D4 (post-/review): Slate is 3 candidates with no pre-ranking; /idea-anchor is where selection happens
- D5 (post-/review): Open Questions Q1-Q3 are gate, not footnote
- D6 (post-/review): Honest Read on User Intent is hypothesis to confirm, not premise shaping candidate preference
- Candidate 1 (audit-tier, forward): multi-drug audit factory + generalization refactor + audit-gated cross-drug experiments before Databricks burst per drug
- Candidate 2 (audit-tier, backward): validate Phase 1 finding via low-cost replication
- Candidate 3 (classifier-tier, post-census): ship deployable classifier conditional on feasibility census
```

**Key details:**
- Preserve `---` separator + position at top of index.
- No restructuring of other index entries.

## Execution Preview

```
Wave 0 (1 step):     Step 1 — Apply 5 review edits to plan
Wave 1 (1 step):     Step 2 — Update plans-index.md entry

Critical path: Step 1 → Step 2 (2 waves)
Max parallelism: 1 agent
```

## Risk Flags

- **Documentation-only plan; no code + no tests fire.** Verification is visual + structural.
- **Step 1's Edit 2 (Candidate 1 merge sentence) is opinionated.** Full-budget framing was user-confirmed; if budget assumption changes before /idea-anchor, this sentence will need re-edit.
- **Cheaper-compute alternative under discussion (Google Colab / personal GPU rig).** If user adopts a non-Databricks compute strategy, the Candidate 1 sentence's `before committing each drug to Databricks burst + N=150 NT cache populate` phrasing becomes inaccurate; flagged for revision before applying.
- **No file overlap between steps.** Step 1 → plans/, Step 2 → wiki/.

## Verification

1. Open `plans/Phase2_Framing_Brainstorm_Plan.md`; scroll to "## Candidate /idea-anchor Sentences". Confirm:
   - 3 candidate blockquotes labeled Candidate 1 / 2 / 3
   - No `(preferred)` label on any candidate
   - "How the 3 candidates differ" paragraph present between Candidate 3 and next section
2. Scroll to "## Open Questions — gate /idea-anchor on these"; confirm heading change + ordering-note intro present.
3. Scroll to "## Honest Read on User Intent"; confirm hypothesis-framing sentence present as opening line.
4. Open `wiki/plans-index.md`; confirm the `[plan_file: Phase2_Framing_Brainstorm_Plan.md]` entry mentions 3 candidates + has D4-D6 decision bullets.
5. No other plan or index entry modified.
