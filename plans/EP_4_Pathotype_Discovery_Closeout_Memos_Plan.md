# EP-4 Pathotype Discovery Closeout Memos Plan

> Two-step sequential execution of the remaining discovery-machine actions for the EP-4 pathotype project: NCBI Pathogen Detection `host_disease` facet audit + EP-1 SUSPEND-gate reuse feasibility memo. User has overridden the prior "deferred until Gate A signal" recommendation.

---

## Problem Statement

The EP-4 pathotype project ledger at `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md` has 5 open Short-term actions on the discovery machine. Three are already completed in this session (substrate-survey, T1+T2 resolution, architecture-fork final lock). Two remain queued behind the workhorse's Gate A signal:

1. **Action 4-new — NCBI Pathogen Detection `host_disease` facet query.** Required to close the COMMENSAL Tier-3 substrate-density gap. The Horesh 2021 independent-label subset (audit shipped 2026-05-27) has only N=2 commensal records; Whittam DECA per-pathotype counts unknown until contact returns. NCBI Pathogen Detection is the only remaining candidate Tier-3 commensal source.

2. **Action 5 — EP-1 SUSPEND-gate audit-gate reuse feasibility verdict.** Required by H5 (`VirulenceFinder/ETECFinder → resolver concordance is stable under pinned versions`). Produces a feasibility-grade memo on whether the existing AMR-mechanism SUSPEND-gate pattern (cipro-specific at `scripts/cipro_mechanism_phenotype_merge.py`; drug-agnostic generalization at `scripts/drug_mechanism_audit.py` + `drug_mechanism_phenotype_merge.py` + `dna_decode/data/mic_tiers.py`) generalizes to pathotype-call opacity.

Both deliverables are research/audit memos, NOT code changes. The v0 implementation contract belongs on the workhorse per the existing handoff doc at `research_outputs/phase4_pathotype_discovery_handoff_to_workhorse_2026-05-27.md`.

## Design Decisions

### D1: Override the "deferred until Gate A" recommendation

**Decision:** execute both items NOW on the discovery machine, rather than wait for Gate A signal from the workhorse.

**Rationale:** the items retain partial value even if Gate A fails. Item 1 (NCBI facet) hardens substrate strategy regardless of architecture — the COMMENSAL gap is a property of the dataset landscape. Item 2 (SUSPEND-gate feasibility) informs the pivot direction if Gate A fails (e.g., direct-reuse strengthens the case for the deterministic-resolver architecture even on a pivoted scope like an ECTyper PR).

**Trade-off:** lower information-gain-per-unit-cost than running these post-Gate-A. If Gate A surfaces install / DB-schema / decision-table friction, the memos become reference documents rather than direct inputs to the next decision. Acceptable risk per user.

### D2: Discovery memos are feasibility-grade, NOT implementation specs

**Decision:** Item 2 produces a 1-page feasibility verdict + analogy mapping table + risks/caveats. NOT a full design spec for `pathotype_tiers.py` or implementation pseudo-code for `pathotype_mechanism_phenotype_merge.py`.

**Rationale:** the workhorse handoff doc explicitly warned "produce feasibility-grade memos, not broad roadmap sprawl." The v0 implementation contract belongs on the workhorse. The discovery machine's job is to surface whether the reuse pattern is viable + at what adaptation cost — not to pre-specify the workhorse's code.

**Trade-off:** the workhorse will need to do the catalog enumeration + classify_noise() rewrite + test rewrite when v0 implementation fires (estimated 4-6 hrs engineering per the Phase 1 exploration). Discovery machine could pre-write some of that, but doing so over-extends into workhorse territory.

### D3: Single commit covering both memos

**Decision:** commit + push both deliverables in one commit with a clear "discovery(ep-4): NCBI Tier-3 facet audit + EP-1 SUSPEND-gate reuse feasibility" message.

**Rationale:** narrow commit boundary; both are closeout-of-discovery artifacts; logically one work unit.

**Trade-off:** if either memo turns out to be wrong / needs revision, the rollback affects both. Acceptable given the small surface (2 markdown files, ~3 pages total).

### D4: Conditional project-ledger update

**Decision:** update `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md` ONLY if either memo produces a clean go/no-go verdict that retires a Pending Decision or updates a Hypothesis status. Skip if verdicts are partial / honest-gap-flagged.

**Rationale:** the ledger is currently in a stable v3 state. Stale updates that don't actually retire pending items add noise. The conditional pattern matches discipline established earlier in the session (only update ledger when decisions actually fire).

## Implementation Plan

### Step 1 — NCBI Pathogen Detection `host_disease` facet query (~45 min)

**Phase 1a — Discover query mechanism (~10 min)**

- WebSearch for: NCBI Pathogen Detection isolate browser SOLR / facet / EUtils API documentation.
- WebFetch the isolate browser landing page + any documented API endpoint to identify the query URL pattern.
- Confirm: how to facet on `host_disease` for organism=Escherichia coli.

**Phase 1b — Run the facet query (~15 min)**

- Execute the discovered query via WebFetch (one-shot, no integration code).
- Capture the `host_disease` value distribution + counts where extractable.
- Identify commensal-indicating values (e.g., "asymptomatic", "healthy carriage", "screening", "non-clinical", "stool surveillance", absent host_disease).
- Cross-check with `isolation_source` (free-text per substrate-survey Row 9) for soft-overlap.

**Phase 1c — Write audit memo (~20 min)**

- Output file: `research_outputs/ncbi-pathogen-detection-host-disease-facet-2026-05-27.md` (~1-2 pages).
- Tier-3 substrate-feasibility verdict: does the host_disease + isolation_source combination close the H2 commensal floor (N≥75)?
- Per-value isolate counts where extractable; honest gaps where blocked.
- Recommendation: confirm Tier-3 as commensal source OR fall back to alternative (e.g., curated commensal panel from K-12 reference + healthy-isolate sub-studies).

**Fallback if API blocked:** document what failed + the manual query path the user can run via the public NCBI isolate browser UI. Treat as legitimate honest-gap.

### Step 2 — EP-1 SUSPEND-gate pathotype reuse feasibility memo (~45 min)

**Phase 2a — Distill verdict structure (~10 min)**

- The Phase 1 exploration agent already produced a detailed feasibility analysis (~80% direct reuse, ~20% adaptation cost; full notes in conversation thread). Distill into a 1-page feasibility memo.
- Verdict shape: `DIRECT_REUSE` / `ADAPTED_REUSE` / `NEW_DESIGN_NEEDED`. Expected verdict: `ADAPTED_REUSE`.

**Phase 2b — Map the analogy (~15 min)**

- Per-strain `noise_class` mapping table (AMR → pathotype). Examples:
  - `CLEAN_R_primary_mechanism` → `CLEAN_pathotype_call_with_primary_cluster_present`
  - `OPAQUE_R_no_mechanism` → `OPAQUE_pathotype_label_no_cluster_marker_detected`
  - `SUSPECT_S_silent_primary_mechanism` → `SUSPECT_commensal_call_with_primary_cluster_present`
  - `NOISY_R_borderline` → `NOISY_pathotype_call_with_partial_cluster_match`
- Catalog mapping table: `DRUG_LOCI_BY_MECHANISM` → `PATHOTYPE_CLUSTER_SIGNATURES`; `DRUG_PRIMARY_MECHANISMS` → `PATHOTYPE_PRIMARY_CLUSTERS`; `CO_RESISTANCE_MECHANISMS` → `CO_PATHOTYPE_MODIFIERS` (efflux/regulatory analogs).
- Verdict thresholds: confirm 0.70/0.40 reusability as starting points + flag retuning expected on real pathotype cohort.
- Input-audit redesign: `mechanism_audit` → `cluster_audit` (VirulenceFinder/Bakta gene-presence matrix + call confidence); `mic_audit` → `phenotype_audit` (serotype/MLST/independent-label).

**Phase 2c — Write feasibility memo (~20 min)**

- Output file: `research_outputs/ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md` (~1 page).
- Sections: feasibility verdict + analogy mapping table + risks/caveats + recommendation to the workhorse.
- Recommendation: re-use `drug_mechanism_phenotype_merge.py` as template; design `pathotype_tiers.py` analog to `mic_tiers.py`; retune verdict thresholds on real pathotype cohort.

### Step 3 — Commit + push both memos (~5 min)

Single commit on `main`, message:

```
discovery(ep-4): NCBI Tier-3 facet audit + EP-1 SUSPEND-gate reuse feasibility

Closes the two remaining Short-term actions on the discovery machine. Both
were deferred behind Gate A signal; user override to proceed now.

Step 1: NCBI Pathogen Detection host_disease facet audit → Tier-3 substrate
feasibility verdict for COMMENSAL class.

Step 2: EP-1 SUSPEND-gate (drug-agnostic version at scripts/drug_mechanism_
phenotype_merge.py) reuse feasibility verdict for pathotype-call opacity.
Verdict: ADAPTED_REUSE; ~80% code reuse + ~20% adaptation per exploration.
```

Push to `origin/main`.

### Step 4 — Conditional project-ledger update (~10 min, conditional)

If Step 1's verdict is FAVORABLE (Tier-3 closes commensal gap with documented per-value counts ≥75):
- Add Decision Made to `project_state/ecoli-pathotype-prediction-cli-2026-05-26.md`.
- Retire the corresponding Pending Decision row.

If Step 2's verdict is FAVORABLE (`ADAPTED_REUSE` confirmed with named adaptation deltas):
- Add Decision Made.
- Update H5's substrate-requirement note.

If either memo lands on `honest-gap` or `partial` verdict, leave ledger as-is.

## Files to be created

- `C:\Users\Farshad\PythonProjects\dna_decode\research_outputs\ncbi-pathogen-detection-host-disease-facet-2026-05-27.md` (~1-2 pages)
- `C:\Users\Farshad\PythonProjects\dna_decode\research_outputs\ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md` (~1 page)

## Files NOT to be modified

- `scripts/cipro_mechanism_phenotype_merge.py` (EP-1 closed)
- `scripts/drug_mechanism_audit.py` / `scripts/drug_mechanism_phenotype_merge.py` (drug-agnostic generalization closed)
- `dna_decode/data/mic_tiers.py` (per-drug catalog frozen)
- `dna_decode/data/bvbrc_genome.py` / `cohort.py` / `refseq.py` (data-layer adapters)
- Any test file (workhorse-owned)
- `research_outputs/phase4_pathotype_discovery_handoff_to_workhorse_2026-05-27.md` (already gate-order-corrected in commit `140eb74`)

## Existing patterns to reuse (NOT to re-implement)

- **For Step 1's WebFetch pattern:** same WebFetch + parse + audit-memo pattern already used in the `/research` substrate survey on 2026-05-27 (Mission Control L1 run `2026-05-27-0140-research-ecoli-pathotype-substrate`). No new integration code.
- **For Step 2's analogy mapping:** the existing drug-agnostic generalization pattern in `scripts/drug_mechanism_audit.py` + `dna_decode/data/mic_tiers.py` (per-drug catalogs swapped at call site). Document the analogous swap for pathotype; do NOT propose a new architecture variant beyond Candidate 5.

## Verification

After Step 3 push:

- `git log --oneline -n 3` on `origin/main` shows the discovery memos commit.
- `git status --short` shows no uncommitted changes from this session.
- The two deliverable memos exist at the named paths and are pushed.
- Workhorse can `git pull` and see both memos alongside the existing handoff.

End-to-end signal: the project ledger's Short-term actions table has 5 of the original 6 actions marked completed. Remaining open items at end of plan:

- Whittam STEC Center direct-contact follow-up (user-owned; email draft ready at `research_outputs/whittam-stec-center-contact-draft-2026-05-27.md`)
- Gate A on workhorse (workhorse-owned)
- Gate B cold-email (user-owned; packet ready at `research_outputs/gate_b_cold_email_packet_2026-05-27.md`)

## Scope guardrail

Per the handoff doc's "What NOT to do on the workhorse" mirroring, corresponding "what NOT to do on this discovery machine" applies:

- Do not write the v0 cluster resolver code (workhorse-owned).
- Do not write a full `pathotype_tiers.py` catalog (workhorse-owned; only enumerate the analogy shape in the feasibility memo).
- Do not run `/research-verify` on the new memos (they're discovery-tier, not audit-tier).
- Do not propose a new architecture variant beyond what's already locked (Candidate 5 from `## Refinement Candidates` in the project ledger).
- Do not modify the existing project ledger v0 Output Contract section.
