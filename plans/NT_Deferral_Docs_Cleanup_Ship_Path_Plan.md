# NT Deferral Docs Cleanup — Ship Path Plan

> Scope-tightened delta from `NT_Deferral_Docs_Cleanup_Plan.md` after `/review` synthesis (2026-05-13). Drops D1's `/save-plan` hedge (already disproven on disk), trims the over-prose `[BLOCKED]` bullet, sharpens the deferral annotation, attempts cheap NT revision retrieval before recording the gap, and adds an explicit untracked-file staging note.

---

## Problem Statement

The technical plan `NT_Deferral_Docs_Cleanup_Plan.md` (3 steps + 2 metadata steps, Wave 0 parallel) was correctly sized and structurally sound, but `/review` (CEO + Eng lenses, 2026-05-13) surfaced selective tightening opportunities:

1. **D1 hedging is moot.** The technical plan's Step 1 still hedges with "Re-run `/save-plan` first; if regenerated entry still omits 'Commit 2 deferred', then manually annotate." Eng verified via the on-disk unstaged diff that `/save-plan` was already run (the `NT_Deferral_Docs_Cleanup_Plan.md` entry was prepended) and did NOT update the two stale `Audit_Calibration_NT_*` entries. The "if /save-plan didn't capture it" conditional is already resolved — go straight to manual edit.

2. **`[BLOCKED]` bullet text exceeds surrounding TODOS density.** Proposed bullet is ~150 words in one prose blob; existing Phase 2.5 perf hardening bullets are 2-3 sentences each. Scannability cost is real when TODOS is the entry point for "what's blocked vs open."

3. **Repro metadata is half-done.** `NT revision=<not recorded>` is honest gap-marking, but the gap remains permanent if nobody retrieves the actual revision when it's cheap (the model was downloaded into `~/.cache/huggingface/hub/` during the equivalence test run; the snapshot hash is sitting on disk).

4. **"Gate failed" is opaque six months out.** Three-word annotation forces a reader to open the saved plan to learn what gate. One-word expansion gives skim-readability.

5. **Step 3 staging gap.** `plans/NT_Deferral_Docs_Cleanup_Plan.md` is currently untracked. The original Step 3 lists it in `Files:` but doesn't explicitly call out `git add` for the untracked path — easy to miss during execution.

The review's clarifying questions (Q1 NT-revision retrievable / Q2 failure-trace pointer / Q3 wave framing keep-or-drop) are resolved inline below.

---

## Design Decisions

### D1: Drop the `/save-plan` hedge — direct manual edit only

**Decision:** Step 1 (annotate stale plans-index entries) is unconditional manual edit. No `/save-plan` first-prong attempt.

**Rationale:** Eng lens verified that `/save-plan` already ran this conversation (the `NT_Deferral_Docs_Cleanup_Plan.md` entry exists in the unstaged `wiki/plans-index.md` diff) and did NOT update the two stale entries. The conditional ("if /save-plan didn't capture it") is moot — it didn't. Per the `/save-plan` skill spec Step 6, existing entries are explicitly skipped. Keeping the hedge in the ship-path plan invites ambiguity for the executor.

**Trade-off:** Considered keeping the hedge as defensive scaffolding — rejected because confirmed evidence resolves the conditional, and ambiguity-as-defense costs executor clarity.

### D2: Trim `[BLOCKED]` bullet to ≤3 sentences + separate Environment line

**Decision:** Replace the long single-prose `[BLOCKED]` bullet with: (a) ≤3 sentences capturing failure mode + plan reference + deferral rationale, (b) separate `Environment at failure:` line, (c) separate `Revisit when:` line.

**Rationale:** Existing Phase 2.5 perf hardening bullets are 2-3 sentences each (`load_bvbrc_ast` slow, `cache.populate` per-sequence, `output_hidden_states=True`). A 150-word prose blob breaks scannability when TODOS is the "what's blocked vs open" entry point.

**Trade-off:** Considered keeping the long form for completeness — rejected because the saved plan file (`plans/NT_Deferral_Docs_Cleanup_Plan.md`) is the canonical full record; TODOS is for at-a-glance status.

### D3: Attempt NT revision retrieval before recording the gap (with env-var precedence + multi-snapshot handling)

**Decision:** Step 2 begins with a retrieval attempt that respects HF cache env-var precedence:
1. `HF_HUB_CACHE` (if set, that IS the cache root)
2. else `HF_HOME` + `/hub`
3. else `XDG_CACHE_HOME` + `/huggingface/hub`
4. else `TRANSFORMERS_CACHE` (legacy)
5. else `%USERPROFILE%\.cache\huggingface\hub` (Windows default; equivalent to `~/.cache/huggingface/hub`)

Within the resolved cache root, look for `models--InstaDeepAI--nucleotide-transformer-v2-100m-multi-species/snapshots/`. If exactly one snapshot dir exists, record its hash. If multiple exist, record `NT revision=<ambiguous: hash1, hash2, ...>` rather than guessing. If the cache directory doesn't exist or `snapshots/` is empty, keep `NT revision=<not recorded>` AND add a third `[OPEN]` TODO bullet to capture revision next time NT loads.

**Rationale:** CEO lens: "writing 'not recorded' into a permanent record adds noise without signal." Eng lens: "don't make up a hash." Default-path-only listing would falsely record `<not recorded>` if any HF env var redirected the cache. Multi-snapshot ambiguity must be surfaced, not silently resolved by picking one. Retrieval is still a single directory listing per candidate root — sub-second cost.

**Trade-off:** Considered dropping the env line entirely if revision is unrecoverable — rejected because `transformers==5.8.0` + `torch==2.11.0` + `python==3.11.5` are still useful even without the model hash.

### D4: Sharpen "gate failed" → "equivalence test failed at model load"

**Decision:** Replace the `**Status:**` line phrasing in Step 1's `wiki/plans-index.md` annotation. New form:

> `**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.`

**Rationale:** "Gate failed" forces a reader six months from now to open the saved plan to learn which gate. One-phrase expansion ("equivalence test failed at model load") gives skim-readability without bloating the entry.

**Trade-off:** Considered including the full state_dict mismatch line — rejected because that detail belongs in the saved plan + TODOS entry, not the index summary.

### D5: Step 3 explicit `git add` reminder for untracked plan file

**Decision:** Step 3 instructions explicitly call out `git add plans/NT_Deferral_Docs_Cleanup_Plan.md` (currently untracked) as a separate sub-bullet from the modified-file staging. Add a pre-commit verification: `git status --short` shows ALL THREE files in the staging area (one new + two modified).

**Rationale:** `plans/NT_Deferral_Docs_Cleanup_Plan.md` is in the untracked file set (`?? plans/NT_Deferral_Docs_Cleanup_Plan.md` per `git status`). The original Step 3's `Files:` list includes it, but `git commit -a` doesn't stage untracked files. An executor running `git add wiki/plans-index.md TODOS.md && git commit` would silently miss it.

**Trade-off:** Considered using `git add -A plans/ wiki/ TODOS.md` — rejected because `-A` semantics are too broad (could accidentally stage other untracked files in those directories).

### D6: Keep wave framing (Q3 resolved)

**Decision:** Preserve the 3-wave / max-parallelism-2 structure from the technical plan.

**Rationale:** CEO lens flagged wave framing as over-spec for the size; Eng lens confirmed it's structurally correct. The framing has near-zero cost (a few lines of metadata) and preserves `/execute-plan` compatibility. Removing it would require executor judgment to re-derive dependencies.

**Trade-off:** Considered dropping waves in favor of a single sequential commit narrative — rejected because keeping the wave block costs nothing and preserves option value for `/execute-plan`.

### D7: Drop failure-trace pointer addition (Q2 resolved)

**Decision:** Do NOT add a `Failure trace:` line to the deferral annotation.

**Rationale:** The equivalence-test failure output ran inline via the Bash tool — it was not persisted to a separate tool-result file. The closest concrete trace is the LESSONS_LEARNED.md 2026-05-13 entry (in `d4a4652`) which captures the salient error. A `Failure trace: see LESSONS_LEARNED.md 2026-05-13` pointer is one degree of indirection away from the BLOCKED bullet itself, which already cites the state_dict mismatch verbatim. Adds words, not signal.

**Trade-off:** Considered persisting the bash output to a `wiki/failures/` log retroactively — rejected because the error message is already captured in TODOS BLOCKED bullet text + LESSONS_LEARNED entry; reconstructing a "trace file" after the fact is artifact-creation noise.

### D8: Drop verification greps (CEO suggestion, light scope reduction)

**Decision:** Replace the grep-based verification of Steps 1-2 with eyeball check + `git diff --stat`. Keep only the post-Step-4 `git status --short --branch` verification.

**Rationale:** Greps that test for the exact strings just typed are tautological — they'd only fail on typos the executor would notice immediately during the edit. `git diff --stat` confirms file shape; visual inspection confirms semantic correctness. Removing the greps cuts noise without losing signal.

**Trade-off:** Considered keeping the greps as defensive — rejected because the greps don't catch realistic failure modes for a manual edit of this size.

---

## Implementation Plan

### Step 1: Annotate stale plans-index entries with sharpened Status
Files: wiki/plans-index.md
Depends on: none

Under both `[plan_file: Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md]` and `[plan_file: Audit_Calibration_NT_AutoModel_Plan.md]` entries, append a `**Status:**` line **between the `**Key decisions:**` block and the trailing `---` separator**:

```
**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.
```

Do NOT touch the prepended `NT_Deferral_Docs_Cleanup_Plan.md` entry (already correct).

### Step 2: NT revision retrieval + split TODOS bullet
Files: TODOS.md
Depends on: none

**Sub-step 2a — Retrieve NT revision (cheap, env-var-aware):**
- Resolve HF cache root in precedence order: `HF_HUB_CACHE` → `HF_HOME/hub` → `XDG_CACHE_HOME/huggingface/hub` → `TRANSFORMERS_CACHE` → `%USERPROFILE%\.cache\huggingface\hub` (Windows default).
- Within the resolved root, list `models--InstaDeepAI--nucleotide-transformer-v2-100m-multi-species/snapshots/`.
- Outcomes:
  - Exactly one snapshot dir → record its hash (40-char hex).
  - Multiple snapshot dirs → record `NT revision=<ambiguous: hash1, hash2, ...>` (do NOT guess).
  - Cache root or model dir missing / `snapshots/` empty → record `<not recorded>` (triggers sub-step 2c).

**Sub-step 2b — Replace TODOS line 40 single bullet with:**

```
- [ ] **[BLOCKED] NT AutoModel swap** (`foundation.py:239`) — `AutoModel.from_pretrained` fails at load on the NT v2 100M `trust_remote_code` checkpoint with state_dict shape mismatch (`Linear[4096, 512]` vs `Linear[512, 2048]`). InstaDeep's `trust_remote_code` defines architecturally distinct AutoModel vs AutoModelForMaskedLM variants. Deferred indefinitely per `plans/Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` Commit 2 gating rule.
  - Environment at failure: `transformers==5.8.0`, `torch==2.11.0`, `python==3.11.5`, os=`Windows 10`, cuda=`available`, NT revision=`<HASH | "ambiguous: ..." | "not recorded">`, loader=`AutoModel.from_pretrained(trust_remote_code=True)`.
  - Revisit when: InstaDeep ships an AutoModel-compatible checkpoint, OR project drops `trust_remote_code` in favor of a manual reimplementation.
- [ ] **[OPEN] NT hidden-state hardening** (`foundation.py:259`) — current code calls `model(**inputs, output_hidden_states=True)` then takes only `hidden_states[-1]`. Drop intermediate-layer materialization for a low-risk speedup. AutoModel swap is blocked (see above); diagnostic spike needed to check whether the loaded `AutoModelForMaskedLM` exposes a base-encoder accessor (`.base_model` / `.encoder` / `.nucleotide_transformer`) returning final hidden state directly. Defer until NT becomes critical-path during Phase 2 smoke.
```

**Sub-step 2c — IF NT revision was `<not recorded>` in 2a, append a THIRD bullet:**

```
- [ ] **[OPEN] Capture NT revision next NT load** — current failure record at `[BLOCKED] NT AutoModel swap` has `NT revision=<not recorded>`. Next time NT v2 100M loads (Phase 2 smoke or any future test), capture the snapshot hash from `~/.cache/huggingface/hub/models--InstaDeepAI--nucleotide-transformer-v2-100m-multi-species/snapshots/` and back-fill the BLOCKED bullet's Environment line.
```

### Step 3: Commit doc deltas (with explicit untracked staging)
Files: wiki/plans-index.md, TODOS.md, plans/NT_Deferral_Docs_Cleanup_Plan.md, plans/NT_Deferral_Docs_Cleanup_Ship_Path_Plan.md
Depends on: Step 1, Step 2

**Sub-step 3a — Stage explicitly (one command, all four files):**
```
git add wiki/plans-index.md TODOS.md plans/NT_Deferral_Docs_Cleanup_Plan.md plans/NT_Deferral_Docs_Cleanup_Ship_Path_Plan.md
```

Note: the two `plans/` files are currently untracked (`??` status). Without explicit `git add`, a `git commit -a` would silently skip them.

**Sub-step 3b — Pre-commit verification:**
- `git status --short` shows exactly four entries in the staging area (one `M` for `wiki/plans-index.md`, one `M` for `TODOS.md`, two `A` for the new plan files).

**Sub-step 3c — Commit:**
```
docs: clean up NT deferral aftermath — index entries + TODOS split + reproducibility

Per /brainstorm + /review against d4a4652:
- wiki/plans-index.md: append Status: lines to the two stale Audit_Calibration_NT_*
  entries recording Commit 2 deferral (manual edit; /save-plan cannot update
  existing entries per its skip-on-duplicate rule)
- TODOS.md: split the conflated "NT output_hidden_states drop" bullet into
  [BLOCKED] AutoModel swap (with reproducibility metadata) + [OPEN] hidden-state
  hardening (broader optimization category) [+ [OPEN] capture NT revision, if
  retrieval failed]
- plans/NT_Deferral_Docs_Cleanup_Plan.md: new technical plan saved earlier this
  conversation
- plans/NT_Deferral_Docs_Cleanup_Ship_Path_Plan.md: review-delta plan (this file)
```

### Step 4: Push to origin
Files: (none)
Depends on: Step 3

`git push origin main`. Currently 1 commit ahead (`d4a4652`); after Step 3 will be 2 ahead.

---

## Verification

- `git diff --stat HEAD~1 HEAD` lists exactly four files changed (eyeball check that they're the expected files).
- `git log --oneline -2` shows `<new-sha> docs: clean up NT deferral aftermath ...` on top of `d4a4652`.
- `git status --short --branch` shows `## main...origin/main` with no `[ahead N]` indicator after Step 4.
- Eyeball: `wiki/plans-index.md` shows the `**Status:**` line under both audit/NT entries, positioned between `**Key decisions:**` and `---`.
- Eyeball: `TODOS.md` Phase 2.5 perf hardening section has TWO bullets (BLOCKED + OPEN) OR THREE bullets (BLOCKED + OPEN + capture-revision) depending on Step 2a retrieval outcome.

## Scope Boundaries

- **In scope:** doc edits to `wiki/plans-index.md`, `TODOS.md`; commit + push of both untracked plan files.
- **Out of scope:** the MLM-wrapper base-encoder diagnostic spike (deferred per saved plan D4 until NT becomes Phase 2 critical path).
- **Out of scope:** any change to `dna_decode/models/foundation.py` (production NT code stays on `AutoModelForMaskedLM`).
- **Out of scope:** Phase 2 entry real-data smoke — separate decision after this cleanup ships.
- **Out of scope:** retroactive "failure trace" log creation (D7 — already captured in LESSONS_LEARNED).
