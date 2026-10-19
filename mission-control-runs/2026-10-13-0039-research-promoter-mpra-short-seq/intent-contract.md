<!-- run-format: 1.0 -->
<!-- intent-contract: 1.0 -->

# Intent Contract — promoter-mpra-short-seq — 2026-10-13T00:39Z

## Verbatim Input

```
Predicting promoter/regulatory activity from short DNA sequences (100-300 bp), especially bacterial promoter MPRA: what model families beat composition baselines, whether multi-scale/dilated convolution helps on SHORT sequences or is long-context-specific, bacterial promoter MPRA state of the art vs replicate ceiling, named techniques with measured effect sizes for escaping a sequence-only plateau near a composition baseline, and whether fixed-background mutational series are easier to predict than diverse-sequence panels
```

### Caller-supplied scoping context (NOT to be re-derived or verified)

- Substrate: GSE144621 bacterial promoter tiles; n=44,106 tiles; every tile exactly 150 bp; leave-peak-out CV.
- Replicate-derived ceiling: Spearman 0.941 (two lab replicates of the same measurement).
- Best local model: Spearman 0.3234 (3-width CNN, 13,217 params). GC-content-only baseline: 0.3112.
  4-scale pyramid (widths 3/12/48/150, 1.70 decades span): 0.3167. All three CNN arms have
  paired CIs including zero vs GC.
- The caller's priority: sub-question 2 (is span/multi-scale gain short-sequence-absent in the
  literature too?) and sub-question 4 (what the field does at exactly this plateau).

## Restated Decomposition

Standard `/research` v0.3 workflow: web search + synthesis → intake validation → followup queue
update. Research is organised around the caller's five numbered sub-questions, with sub-questions
2 and 4 prioritised for depth.

Sub-tasks:
1. Web research (Steps 1-3 produce `<slug>.raw.md`)
2. Intake validation (Step 4 produces `<slug>.md` + `<slug>_unsupported.md`)
3. Followup queue update (Step 5 updates `_followup_queue.md`)

## Method Used

Implicit decomposition by `/research` orchestrator pipeline (v0.3), with the caller's five
sub-questions as the search-angle enumeration.

## Verification Criteria

| Sub-task | Verification criterion | Evidence shape |
|---|---|---|
| Web research | ≥5 audit-grade rows OR honest-gap declared | `<slug>.raw.md` exists with V1 13-column table; row count + honest-gap status visible |
| Intake validation | Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan + source-identity advisory | `<slug>.md` (≥1 row OR explicit empty-memo marker) + `<slug>_unsupported.md` exist |
| Followup queue update | `_followup_queue.md` exists + touched within this run | mtime of `_followup_queue.md` ≥ run-start timestamp |

## Out-of-Scope Boundary

This run explicitly will NOT:

- Write outside `<cwd>/research_outputs/` or `<cwd>/mission-control-runs/2026-10-13-0039-research-promoter-mpra-short-seq/`
- Write anything into `wiki/`, `project_state/`, `features/`, or any project ledger
  (explicit caller prohibition, reinforcing the standing path gate)
- Modify rules YAML, ADSP units, or any code under `dna_decode/` or `scripts/`
- Fabricate sources to hit row count targets
- Auto-promote candidates into rules/wiki (Promotion Gate remains manual)
- Re-derive or attempt to verify the caller-supplied substrate figures (n=44,106; 150 bp;
  ceiling 0.941; 0.3234 / 0.3112 / 0.3167) — these are given, not research targets
- Run a model, train anything, or touch the GSE144621 data
- Spend more than 15% of daily token budget
- Run longer than 30 minutes wall-clock
- Make more than 100 tool calls

## Escalation Conditions

- Zero useful web search results after 3 attempts → halt + ask user to refine topic
- Intake rejects ≥80% of rows → halt + surface likely topic-shape failure
- Slug validation fails → halt + ask user for manual slug
- Followup queue file write fails → halt; supported memo still preserved
- Token budget ≥15% daily exhausted
- Wall-clock ≥30 minutes elapsed
- Tool-call count ≥100
- Unresolved-uncertainty count ≥5

## Run Metadata

- **Run ID:** 2026-10-13-0039-research-promoter-mpra-short-seq
- **Started:** 2026-10-13T00:39:00Z
- **Autonomy Level:** L1
- **Department(s):** Research
- **Skill invoked:** `/research`
- **Caller:** user (relayed via parent agent)
- **Budget caps:**
  - Token: 15% daily
  - Wall-clock: 30 minutes
  - Tool-calls: 100
  - Unresolved-uncertainty: 5

## Contract version

v1.0
