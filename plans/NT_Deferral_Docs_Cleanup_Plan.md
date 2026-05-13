# NT Deferral Docs Cleanup Plan

> Docs-only follow-up to commit `d4a4652` — apply the three issues surfaced by `/brainstorm` against the just-shipped NT AutoModel refactor deferral: stale plans-index, conflated TODOS scope, missing reproducibility metadata.

---

## Problem Statement

Commit `d4a4652` deferred Commit 2 of `plans/Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` after the equivalence test failed at NT v2 100M model-load time (state_dict shape mismatch under `trust_remote_code`). The cleanup committed alongside that deferral captured the right technical facts but missed three things, identified by `/brainstorm`:

1. **`wiki/plans-index.md` is stale on a load-bearing disposition.** Both index entries for the audit/NT plans still describe Commit 2 / the `AutoModel` switch as in-flight key decisions. A future reader would think Commit 2 is pending, not blocked.
2. **The TODOS entry conflates two separable scopes.** The current entry marks the NT `output_hidden_states=True` drop as `BLOCKED` — but the failed gate only proves `AutoModel.from_pretrained` is unsafe. The loaded `AutoModelForMaskedLM` may still expose a base-encoder accessor that avoids `output_hidden_states=True`. The broader hidden-state hardening category should stay open as a separate diagnostic, not be closed by the AutoModel-specific failure.
3. **Failure record lacks reproducibility metadata.** For a `trust_remote_code` failure, future retries against changed remote code can silently behave differently. TODOS currently records only the error text — not `transformers` version, model revision, torch version, or loader class attempted.

## Design Decisions

### D1: Re-run `/save-plan` before manual index edit

**Decision:** Re-invoke `/save-plan` first to regenerate the plans-index entries for the audit/NT plans. If the regenerated summaries still omit the deferral disposition, then manually append a `**Status:**` line under each affected entry.

**Rationale:** `wiki/plans-index.md` header says "Auto-maintained by /save-plan. Do not edit manually." The auto-maintenance convention should be preserved when possible. Manual edits are justified only when the auto-tool can't capture load-bearing project state.

**Trade-off:** Considered going straight to manual edit since the disposition note is high-signal — rejected because preserving the index's auto-maintenance pattern reduces drift risk for other entries.

### D2: Split the TODOS entry into specific + general

**Decision:** Replace the single "NT model: drop `output_hidden_states=True` ... **BLOCKED**" entry with two:
- **Specific (blocked):** `AutoModel.from_pretrained` swap blocked indefinitely with state_dict evidence.
- **General (open):** NT hidden-state hardening category — diagnostic needed to inspect whether the loaded `AutoModelForMaskedLM` exposes `.base_model` / `.encoder` / `.nucleotide_transformer` or equivalent accessor that returns final hidden states without `output_hidden_states=True`.

**Rationale:** Conflating the two closes the broader optimization investigation on the basis of one specific failure. The general entry preserves the optimization opportunity for a future bounded diagnostic spike.

**Trade-off:** Considered keeping a single entry with nested sub-bullets — rejected because TODOS structure is flat and sub-bullets get missed in scans.

### D3: Lean reproducibility metadata, one line

**Decision:** Append one line to the specific TODOS entry:
`Environment at failure: transformers==5.8.0, torch==2.11.0, NT revision=<hash or "not recorded">, loader=AutoModel.from_pretrained(trust_remote_code=True)`.

**Rationale:** Minimum needed to distinguish "we fixed our code" from "remote code moved" on future retry. Don't dump full environment — that's noise. Don't include local cache paths — irrelevant to reproducibility.

**Trade-off:** Considered adding it to the plan status header instead — rejected because the plan header is human-prose, and the metadata is structured-lookup material; TODOS is the right surface.

### D4: Diagnostic spike deferred, NOT killed

**Decision:** Park the MLM-wrapper base-encoder accessor diagnostic as a follow-up under the new general TODOS entry. Do not run it now. Trigger condition: NT becomes critical-path during Phase 2 smoke (i.e., `cache.populate` runs NT over many genes).

**Rationale:** Optimization on a speculative hot path is scope creep. The general TODOS entry preserves the opportunity without burning current focus.

## Implementation Plan

### Step 1: Re-run `/save-plan` for stale index entries

**Files:** `wiki/plans-index.md`

Re-invoke `/save-plan` against `plans/Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` (the deferred plan) so the auto-maintained summary picks up the new `**Status:**` header on that plan. Verify the regenerated index entry mentions the deferral.

If the regenerated entry still omits "Commit 2 deferred 2026-05-13 — gate failed," then manually append a single `**Status:**` line under both affected entries:
- `[plan_file: Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md]` entry (around line 4-11)
- `[plan_file: Audit_Calibration_NT_AutoModel_Plan.md]` entry (around line 15-22, since D4 there advertises the AutoModel switch)

Note in the eventual commit message: "index correction after a failed gate."

### Step 2: Split the TODOS entry

**Files:** `TODOS.md`

Under `## Phase 2.5 perf hardening (deferred from Gate B prep)`, replace the current single bullet:

> `**NT model: drop output_hidden_states=True** (foundation.py:259) — low-risk speedup; current code requests all hidden states then takes only the last. Use outputs.last_hidden_state if available. **BLOCKED 2026-05-13:** AutoModel.from_pretrained cannot load the NT v2 100M checkpoint — trust_remote_code modeling code defines a Linear(512, 2048) layer where the checkpoint has [4096, 512]. ...`

With two bullets:

> `**[BLOCKED] NT AutoModel swap** (`foundation.py:239`) — `AutoModel.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-100m-multi-species", trust_remote_code=True)` fails at load time with state_dict shape mismatch (`Linear[4096, 512]` checkpoint vs `Linear[512, 2048]` AutoModel architecture). InstaDeep's trust_remote_code defines architecturally distinct AutoModel vs AutoModelForMaskedLM variants, not just head-differential. Equivalence test attempted per `plans/Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` Commit 2; deferred indefinitely per D3 gating rule. Revisit only if InstaDeep ships an AutoModel-compatible checkpoint or we drop trust_remote_code in favor of a manual reimplementation. Environment at failure: transformers==5.8.0, torch==2.11.0, NT revision=<not recorded>, loader=AutoModel.from_pretrained(trust_remote_code=True).`

> `**[OPEN] NT hidden-state hardening** (`foundation.py:259`) — current code calls `model(**inputs, output_hidden_states=True)` then takes only `hidden_states[-1]`. Drop the materialization of intermediate layers for a low-risk speedup. AutoModel swap is blocked (see above), but the loaded `AutoModelForMaskedLM` may expose a base-encoder accessor (`.base_model`, `.encoder`, `.nucleotide_transformer`, or similar) that returns the final encoder hidden state directly. Diagnostic spike deferred until NT becomes critical-path during Phase 2 smoke (i.e., `cache.populate` runs NT over many genes). Premature now.`

### Step 3: Commit the doc deltas

**Files:** `wiki/plans-index.md`, `TODOS.md`, `plans/NT_Deferral_Docs_Cleanup_Plan.md` (this file)

Single commit, message: `docs: clean up NT deferral aftermath — index entry status + TODOS split + reproducibility metadata`.

### Step 4: Push

`git push origin main` — currently 1 commit ahead (`d4a4652`); after Step 3 will be 2 ahead.

## Verification

- `wiki/plans-index.md` entry for the ship-path plan mentions Commit 2 deferral with date.
- `TODOS.md` has TWO bullets under Phase 2.5 perf hardening covering the AutoModel swap (BLOCKED) and the hidden-state hardening (OPEN) as separable items.
- Specific bullet includes the four-field environment line.
- `git status` clean after Step 3.
- `git log --oneline -2` shows the cleanup commit on top of `d4a4652`.
- `origin/main == main` after Step 4.

## Scope boundaries

- **Out of scope:** running the MLM-wrapper diagnostic spike (D4 — deferred until NT becomes critical-path).
- **Out of scope:** any change to `dna_decode/models/foundation.py` (the production NT code stays on `AutoModelForMaskedLM`).
- **Out of scope:** Phase 2 entry real-data smoke — separate decision after this docs cleanup ships.
