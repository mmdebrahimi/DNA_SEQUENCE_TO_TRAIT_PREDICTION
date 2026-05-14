# Sidework Sequence — Ship Path Plan

> **Status:** ✅ EXECUTED 2026-05-14 — all 6 steps complete via `/execute-plan` (sequential). Commits: A `66dfde2`, C `2ca1799`, E `9538ef8`, test-epilogue `8f3dbe0`. N=40 cohort built at `data/processed/gate_b_n40_cipro_cohort.parquet`. Auto-memory updated with GPU reality. Final test count: 368 passed / 1 skipped / 0 regressions.

> Delta from `Sidework_Sequence_Plan.md` after `/review` (CEO + Eng lenses, 2026-05-13). Resolves the load-bearing B-scope problem the Eng lens surfaced (narrow fix is degenerate on N=12 unique MLST), expands Step C scope (missed ARCHITECTURE.md + GATE_A_REPORT.md), tightens code-quality choices (helper relocation, test count, drop over-engineered dict counter), corrects TODOS.md hunk reality (line 17 needs new author edit not staging), and adds CEO-flagged process discipline (numerical before/after snapshot in B, time-box C, post-populate slow tests).

---

## Problem Statement

`/review` of `Sidework_Sequence_Plan.md`'s technical-plan expansion (CEO + Eng agents, 2026-05-13) surfaced one critical issue and several material refinements:

**Critical (Eng lens, grounded):** Step B's "narrow fix" doesn't deliver its stated goal on the actual 12-strain cohort. Verified by reading `data/processed/gate_b_mini_cohort.parquet`: all 12 strains have 12 *unique* MLST strings. Replacing `hash(s.mlst) % 10` with a real MLST parser maps each strain to a singleton clade ID. Under LOSO, every held-out strain is in a clade unseen during training → `predict_clade_only` falls back to `global_positive_rate` → clade-only AUROC ≈ 0.5 regardless of foundation model quality. The OLD hash placeholder at least produced groups of size >1 (artificially); the "fix" is numerically *worse* for the smoke gate's purposes.

The deferred `per_clade_baseline` strain-keying semantics fix (originally Step B's deeper-scope option) is the actually load-bearing one. Three resolution options identified.

**Material (Eng lens, grounded):**
- Step C scope misses `docs/ARCHITECTURE.md` (lines 79, 118 contain `4-bit` references) and possibly `wiki/GATE_A_REPORT.md` (line 60: "Evo 4-bit needs Linux/CUDA via WSL2"). Plan claimed 11 files but listed 12 in the table — count drift.
- Step B test count over-scoped: 5 unit tests for a single helper with one call site is unnecessary. 1 parametrized test (~6 cases) + 1 integration assertion in `tests/test_pipeline_cli.py` is sufficient.
- `mlst_to_clade_id` belongs in `dna_decode/data/cohort.py` (data-layer concern; CandidateStrain.mlst is defined there), NOT `eval/clade_baseline.py`. Avoids back-reference from eval to data semantics.
- `fallback_counter: dict | None` parameter is over-engineered. Use stable per-input hash (e.g., `-(hash(s) & 0xffff)`) for singleton IDs; rely on `predict_clade_only`'s existing `global_positive_rate` fallback for unseen clades.
- Comma-joined MLST format `MLST.A.X,MLST.B.Y` (BV-BRC emits multi-scheme records) is not handled. Trailing-integer-only parser produces silently-wrong results.
- Scheme collision: `MLST.ecoli_achtman_4.131` and `MLST.Escherichia_coli_1.131` are different lineages but trailing-integer parser collapses them.
- TODOS.md hunk discipline broken: Step C planned to take "the line 17 GPU-sanity-check hunk via `git add -p`", but line 17 in the working tree is **not modified** — Step C must author a new edit, not stage an existing one. Plan misframed C as staging-only.
- Step D safety margin is zero: source cohort has exactly 25 R + 42 S for cipro; `--per-class 25` consumes the entire R class. Any future cohort regen that drops a single R strain breaks D silently.
- `pytest -m "not slow"` emits `PytestUnknownMarkWarning` because the marker isn't registered. Works but noisy.

**Material (CEO lens, inferred):**
- M (auto-memory edit) doesn't belong in the wave graph as a "step" — it's a different repo with no git op + no dependency. Demote to "side errand, do whenever."
- E sandwiched between code changes (Wave 2, between C and B) is mildly awkward — planning docs that describe work-not-yet-done get committed mid-implementation. Minor optics issue for solo project; not load-bearing.
- B is the only commit that can shift the smoke-gate acceptance bar. Add a commit-message line item with numerical before/after snapshot of clade-only baseline on a fixed input.
- Time-box C to 30 min; if it overruns, commit partial + defer rest. Don't let doc-correction balloon past populate completion.
- After populate finishes, run the GPU/slow tests A's pre-flight skipped, BEFORE running the smoke gate. Don't proceed to smoke gate on a HEAD whose tests were partially validated.

---

## Design Decisions

### D1: B-resolution path is an open decision pending user choice

**Decision (LOCKED 2026-05-14 via post-save /brainstorm): B-B (Skip clade-only from smoke gate).** Smoke gate runs **4 variants** — AMRFinder + k-mer + gene-presence + NT-XGBoost. Clade-only column is dropped because 12 unique MLST strings × LOSO = guaranteed singleton-clade fallback regardless of helper correctness. Deferred clade-only baseline picks up at Stage 1 N=50 where MLST replication exists.

Other options considered + rejected:
- **B-A (Widen)** — would have added the strain-keying fix in `pipeline.py:294-300`: rename `per_clade_*` → `per_strain_*` OR replace LOSO with `leave_one_mlst_out_cv` for the clade-baseline path. ~1-2 hours added. Rejected — still doesn't answer much at N=12; fixes a baseline whose information content is structurally bounded by cohort size.
- **B-C (Defer)** — would have kept B as a placeholder swap (helper + tests, no behavioral change on N=12), documenting the degenerate outcome explicitly. Rejected — keeps a known-degenerate result in the smoke path; mask-without-fix.

**B-B implementation scope (cheap, ~15 min):**
- Edit `scripts/pipeline.py` to drop the `--include-clade-baseline` flag from the smoke gate invocation (or simply don't pass it).
- Update `plans/Phase2_Decision_Gate_Plan.md` D5 smoke acceptance bar: "NT-XGBoost not obviously worse than k-mer + gene-presence + AMRFinder" (drop clade-only-fixed reference; 4 variants).
- No helper extraction needed for smoke. The `mlst_to_clade_id` helper still gets built at Stage 1 N=50 time when clade-only becomes meaningful.

**Rationale for lock:** Post-save /brainstorm verified the 12-unique-MLST finding empirically. B-B is cheapest (15 min vs 30-120 min) + most honest (no degenerate baseline in the result packet) + best preserves the deferred fix for the right N.

**Trade-off:** B-A and B-C are both technically viable; B-B is the right choice for N=12 specifically. At Stage 1 (N=50), if MLST replication exists, the deferred work activates with full context.

### D2: Step C scope = 11 + `docs/ARCHITECTURE.md` with per-line judgment; GATE_A_REPORT.md NOT in scope

**Decision (LOCKED 2026-05-14):** Add `docs/ARCHITECTURE.md` to Step C, but **inspect each `4-bit` reference individually**; rewrite only lines that imply current-feasibility claims about 4-bit quantization on the actual hardware. Leave future-state quantization-workflow descriptions intact (they may remain accurate as future/optional check descriptions).

**`wiki/GATE_A_REPORT.md` confirmed OUT of scope:** Eng lens verified line 60 reads "Quantization — Evo 4-bit needs Linux/CUDA via WSL2; not exercised here" — this is a historical Gate A caveat (correctly stating Evo quantization was NOT exercised), NOT a stale RTX 4090 claim. Different document type than the rest of C's scope.

**Final Step C scope: 12 files** — `CLAUDE.md`, `TODOS.md`, `FUTURE_FEATURES.md`, `README.md`, `wiki/phase1_ship_report.md`, `wiki/PHASE2_PREFLIGHT.md`, `pyproject.toml`, `config/datasources.yaml`, `scripts/quantize_fidelity_check.py`, `dna_decode/data/pilot.py`, `dna_decode/models/foundation.py`, `LESSONS_LEARNED.md` (append GPU-spec lesson) + `docs/ARCHITECTURE.md` (per-line judgment). 13 files total. Note: today's earlier /documentation invocation already added 2 NEW lessons to `LESSONS_LEARNED.md`; Step C will add a third (the GPU-spec one).

**Rationale:** Per-line judgment for ARCHITECTURE.md avoids mechanically replacing valid future-state descriptions. Excluding GATE_A_REPORT.md keeps Step C scope on docs that mislead future-self; gate reports are historical artifacts.

**Trade-off:** Considered mechanically rewriting all `4-bit` mentions in ARCHITECTURE.md (rejected — would alter potentially-correct future-state language). Considered including GATE_A_REPORT.md for thoroughness (rejected — line 60 is a correct historical caveat, not a current-feasibility claim).

### D3: Step C is edit-then-stage, not stage-only

**Decision:** Reframe Step C as: author new edits across 12-14 files (including TODOS.md line 17 rewrite), THEN stage + commit. Plan previously framed it as "take the line 17 hunk via `git add -p`" which incorrectly assumed the line was already modified.

**Rationale:** `git diff TODOS.md` confirms line 17 is unchanged in working tree. The line says "confirm RTX 4090 has ≥24GB VRAM" — needs rewriting, not just staging.

**Trade-off:** None; the original framing was simply wrong.

### D4: `mlst_to_clade_id` helper lives in `dna_decode/data/cohort.py`

**Decision:** Move the helper from the originally-proposed `dna_decode/eval/clade_baseline.py` to `dna_decode/data/cohort.py`.

**Rationale:** `CandidateStrain.mlst` is defined in `cohort.py`. `_mlst_balanced_selection` (existing helper) also lives there. A future `leave_one_mlst_out_cv` consumer (referenced in `Ecoli_G2P_Platform_Technical_Plan.md:236`) would import from `cohort.py` not `eval/`. Putting it in `eval/` creates an inverted dependency.

**Trade-off:** Considered keeping it in `eval/clade_baseline.py` to minimize file count (rejected — wrong architectural home).

### D5: Test scope reduced from 5 unit tests to 1 parametrized + 1 integration

**Decision:** Replace the 5 individual unit tests with:
- 1 parametrized test in `tests/test_data_cohort.py` (extend — file confirmed existing 2026-05-14 via /brainstorm) — covers all 6 MLST format cases (achtman_4, Escherichia_coli_1, comma-joined multi-scheme, empty string, None, garbage)
- 1 integration assertion in `tests/test_pipeline_cli.py` (extend — file confirmed existing 2026-05-14) — `pipeline.py train --include-clade-baseline` runs without raising on the 12-strain mini cohort

**Note:** B-B was selected (clade-only dropped from smoke); the integration test still has value as a regression guard for the `--include-clade-baseline` flag path, which remains in `pipeline.py` for future N=50 use. If B-B's lighter scope is taken (just dropping the smoke-time flag), the helper + tests are deferred to Stage 1 prep — see D1.

**Rationale:** Eng lens flagged 5-unit-test/single-helper as over-budget. Parametrized cases cover the same input domain with less code. Integration assertion catches signature drift in the call site.

**Trade-off:** Considered the original 5-test approach (rejected — over-engineered for one helper). Considered no integration test (rejected — call site is the actual risk surface).

### D6: Drop `fallback_counter: dict` parameter from `mlst_to_clade_id`; use deterministic hashing

**Decision:** Helper signature is `mlst_to_clade_id(mlst_str: str | None) -> int`. No mutable counter. For None / unparseable inputs, use a **deterministic** stable per-input hash via `hashlib.blake2b` or `zlib.crc32` (NOT Python's built-in `hash()` — see D7 for why). Let downstream `predict_clade_only` handle unseen-clade fallback via existing `global_positive_rate` logic.

**Rationale:** Eng lens flagged the counter dict as over-engineered side-effect coupling. The downstream code already handles unseen-clade fallback. Singleton IDs only need to be distinct per-input AND stable across process runs. Python's `hash()` is process-salted (PYTHONHASHSEED randomized per process) — using it would break reproducibility across runs.

**Trade-off:** Considered keeping the counter for "stable distinct IDs per call site" (rejected — adds coupling for zero current benefit). Considered Python's built-in `hash()` (rejected — non-deterministic across processes, breaks ML reproducibility).

### D7: Multi-scheme MLST + scheme-collision semantics — explicit decision

**Decision:** For multi-scheme records (`MLST.A.X,MLST.B.Y`), take the first scheme's ST + scheme name as a tuple, hash to int via **deterministic** hashing. For scheme-collision (same trailing ST in different schemes), use `(scheme_name, ST)` as the key so they collapse only when scheme matches. Pseudocode (uses `zlib.crc32`, NOT Python's process-salted `hash()`):

```python
import zlib

def _stable_hash(s: str) -> int:
    """Deterministic 32-bit hash. Stable across processes/platforms."""
    return zlib.crc32(s.encode("utf-8"))

def mlst_to_clade_id(mlst_str: str | None) -> int:
    if not mlst_str:
        return -1  # singleton via downstream fallback
    first_scheme = mlst_str.split(",")[0]  # take first if multi
    parts = first_scheme.split(".")
    if len(parts) < 3:
        return -(_stable_hash(mlst_str) & 0xffff)  # negative singleton
    scheme, st = parts[1], parts[-1]  # e.g., "ecoli_achtman_4", "410"
    try:
        st_int = int(st)
    except ValueError:
        return -(_stable_hash(mlst_str) & 0xffff)
    return _stable_hash(f"{scheme}.{st_int}")  # positive 32-bit int, scheme-aware, deterministic
```

**Rationale:** Eng lens flagged "trailing integer of whole string" as silently wrong on (a) comma-joined multi-scheme records and (b) cross-scheme ST collisions. Scheme-aware hashing fixes both. Post-save /brainstorm 2026-05-14 caught that the earlier pseudocode used Python's built-in `hash()` which is **process-salted** (PYTHONHASHSEED is randomized per process) → clade IDs would change across runs → ML reproducibility broken silently. ChatGPT cross-engine review 2026-05-14 preferred `zlib.crc32` over `hashlib.blake2b` for this use case (MLST cardinality is ~thousand-scale; CRC32's 32-bit space is more than sufficient + smaller + faster). Both are deterministic stdlib; CRC32 wins on simplicity.

**Trade-off:** Considered `hashlib.blake2b` (more bits, marginally heavier — chose CRC32 for MLST grouping where collision robustness isn't load-bearing at this cardinality). Considered an explicit `{mlst_string: int}` persisted mapping in cohort parquet metadata (deferred to `FUTURE_FEATURES.md` Phase 2.5 entry — "cleanest scientific approach" per ChatGPT but adds schema work). Considered (a) take-first-scheme-only, (b) explode-into-membership-in-both-groups (rejected — overcomplicated for the smoke gate's purpose).

### D8: Step M demoted out of the wave graph

**Decision:** Auto-memory update on `~/.claude/projects/C--Users-Farshad/memory/user_environment.md` is no longer a "Step." Labelled as a "side errand — do whenever, no dependency, no commit." Removed from Wave 0.

**Rationale:** Wave graph is about commit dependencies through shared files. M has neither (different repo, no git op). Listing it as a step inflated the apparent work.

**Trade-off:** None.

### D9: Time-box C to 30 min; commit partial if it overruns

**Decision:** Step C has a 30-min hard time-box. At minute 30: if all 12-14 files done, commit. If not, commit what's done + leave the rest with a `# TODO(gpu-docs)` marker + file a `[OPEN]` TODO for the unfinished files.

**Rationale:** CEO lens flagged C as the "doc-correction balloon" risk. Time-boxing forces explicit deferral over implicit scope creep.

**Trade-off:** Considered no time box (rejected — open-ended docs work has historically expanded). Considered 60 min (rejected — 30 is tight enough to force focus, loose enough to handle the realistic file count).

### D10: B commit message includes numerical before/after snapshot

**Decision:** Step B's commit body includes a one-line numerical before/after of the clade-only baseline AUROC on the 12-strain mini cohort:

```
Clade-only AUROC on gate_b_mini_cohort:
  before (hash placeholder): <X>
  after (real MLST parse):   <Y>
```

**Rationale:** CEO lens noted B is the only commit that can shift the smoke-gate acceptance bar. Capturing the before/after in the commit message makes the bar-shift explicit and auditable.

**Trade-off:** Considered logging to a separate file (rejected — commit message is the right durability surface for a one-line numerical fact).

### D11: Post-populate slow tests gate the smoke gate

**Decision:** When NT populate completes, run `uv run pytest tests/test_models_foundation.py -m slow -v` BEFORE running the smoke gate. Don't run smoke on a HEAD whose tests were only partially validated.

**Rationale:** CEO lens flagged that A's pre-flight used `-m "not slow"` to avoid GPU contention with the running populate. Once populate finishes, the slow tests are runnable; should be cleared before the next gate.

**Trade-off:** Considered skipping (rejected — minor cost, large confidence gain).

### D12: Step D `--per-class 20` (locked — leaves 5 R margin)

**Decision (LOCKED 2026-05-14):** Step D uses `--per-class 20` → N=40 cohort (20R / 20S). Leaves 5-strain margin in the R class for future cohort regenerations that might drop a strain.

**Rationale:** Eng lens caught that the source cohort has exactly 25 R for cipro; `--per-class 25` has zero margin and `build_mini_cohort.py:58` hard-fails when R count < requested. N=40 is still a substantial improvement over N=12 for the Stage 1 engineering screen (AUROC 95% CI width ≈ ±0.10 at N=40 vs ±0.19 at N=12). The 2-percentage-point noise-floor difference between N=40 and N=50 is negligible compared to the brittleness cost of zero-margin builds.

**Trade-off:** Considered `--per-class 25` for exactly-N=50 (rejected — brittle to regen). Considered `--per-class 22` (rejected — splits the margin difference; clean N=40 is easier to communicate).

### D13: Register `slow` pytest marker

**Decision:** Add to `pyproject.toml` under `[tool.pytest.ini_options]`:

```toml
markers = [
    "slow: tests that require GPU + large model downloads (excluded by default in CI)",
]
```

**Rationale:** Eliminates the `PytestUnknownMarkWarning` that A's pre-flight emits. Cosmetic but warns get easier to spot.

**Trade-off:** None; small mechanical fix.

---

## Revised Implementation Plan

### Step 1: Validate batching changes (A pre-flight)

```bash
cd C:/Users/Farshad/PythonProjects/dna_decode
uv run pytest tests/test_models_foundation.py -m "not slow" -v
```

If the `slow` marker isn't registered yet (per D13), this emits a warning but still works.

### Step 2: Register `slow` marker (D13 small fix, lands with A)

Edit `pyproject.toml` `[tool.pytest.ini_options]` to add the markers list.

### Step 3: Commit A — batching refactor + tests + marker registration

```bash
git add dna_decode/models/foundation.py dna_decode/models/cache.py tests/test_models_foundation.py pyproject.toml
git add -p TODOS.md  # ONLY the [x] RESOLVED hunk for batching
git commit -m "perf(cache): batched embedding populate; 25× max speedup on Ada+ GPUs"
```

### Step 4: Commit C — GPU-docs correction (D2 expanded scope; D3 edit-then-stage; D9 time-boxed)

Author new edits across 12-14 files (replace RTX 4090 / 24 GiB / 4-bit Evo claims with GTX 860M reality). Time-box: 30 min. If overrun at minute 30: commit what's done + add `# TODO(gpu-docs)` markers + file `[OPEN]` TODO for unfinished files.

Files (12 confirmed + 2 conditional):
- `CLAUDE.md`, `TODOS.md` (line 17 rewrite), `FUTURE_FEATURES.md`, `README.md`, `wiki/phase1_ship_report.md`, `wiki/PHASE2_PREFLIGHT.md`, `pyproject.toml` (dep comments), `config/datasources.yaml`, `scripts/quantize_fidelity_check.py:3`, `dna_decode/data/pilot.py`, `dna_decode/models/foundation.py` (EvoModel docstring), `LESSONS_LEARNED.md` (NEW lesson entry), `docs/ARCHITECTURE.md` (D2 added)
- Conditional: `wiki/GATE_A_REPORT.md` (D2: only if grep confirms stale GPU claims)

```bash
git add <12-14 files>
git commit -m "docs: correct GPU spec across docs — was RTX 4090, actually GTX 860M"
```

### Step 5: Commit E — planning docs + research outputs

```bash
git add plans/Phase2_Decision_Gate_Plan.md plans/Sidework_Sequence_Plan.md plans/Sidework_Sequence_Ship_Path_Plan.md wiki/plans-index.md
git add -p TODOS.md  # ONLY the Phase 2 prerequisites section hunk
git add research_outputs/
git commit -m "docs: Phase 2 decision-gate plan + research synthesis"
```

### Step 6: Commit B — clade-only fix (D1: USER CHOOSES B-A vs B-B vs B-C)

**Path B-A (Widen — strain-keying fix included):**
- Move helper `mlst_to_clade_id` to `dna_decode/data/cohort.py` (D4) with scheme-aware tuple-hash (D7)
- Fix `pipeline.py:294-300` strain-keying — either rename `per_clade_*` → `per_strain_*` OR use `leave_one_mlst_out_cv` for the clade-baseline path
- 1 parametrized test in `tests/test_data_cohort.py` + 1 integration in `tests/test_pipeline_cli.py` (D5)
- TODOS.md: no new `[OPEN]` needed (this option subsumes the deferred TODO)
- Commit body includes numerical before/after (D10)
- Estimated: 1-2 hours

**Path B-B (Skip — drop clade-only from smoke):**
- Edit `scripts/pipeline.py` to skip the `--include-clade-baseline` path for the smoke gate (or just don't pass the flag)
- No helper extraction needed; defer to Stage 1
- Update `plans/Phase2_Decision_Gate_Plan.md` to reflect 4 smoke variants instead of 5; remove the "clade-only-fixed" reference from D5's smoke-gate acceptance bar
- Estimated: 15 min

**Path B-C (Defer — narrow placeholder swap + acknowledged degeneracy):**
- Move helper to `dna_decode/data/cohort.py` (D4) with scheme-aware tuple-hash (D7)
- Swap `pipeline.py:285` call site
- 1 parametrized test + 1 integration test (D5)
- TODOS.md: append `[OPEN] per_clade_baseline strain-keying fix (load-bearing for clade-only at N=12)` per source plan D4
- Add a comment at `pipeline.py:285` documenting the known degenerate outcome on N=12 unique-MLST cohorts
- Estimated: 30-45 min

### Step 7 (out-of-wave): Auto-memory update (D8 demoted)

Edit `~/.claude/projects/C--Users-Farshad/memory/user_environment.md` — replace RTX 4090 claim with GTX 860M reality. No commit. Do whenever.

### Step 8 (out-of-wave): Optional N=50 cohort prep

```bash
uv run python scripts/build_mini_cohort.py \
  --source data/processed/gate_b_cohort.parquet \
  --output data/processed/gate_b_n50_cipro_cohort.parquet \
  --drug ciprofloxacin --per-class 25  # or 20 for safety margin per D12
```

### Step 9: Post-populate slow tests (D11)

After NT populate completes:

```bash
uv run pytest tests/test_models_foundation.py -m slow -v
```

Should pass on the post-A HEAD. Run before the smoke gate.

---

## Verification

After Steps 1-6 complete (and 7-9 as applicable):

- `git log --oneline -4` shows four commits (A, C, E, B-X) in order.
- `git status --short` clean except optionally `data/processed/gate_b_n50_cipro_cohort.parquet` (D output).
- `grep -rE "RTX 4090|24 GiB" CLAUDE.md TODOS.md FUTURE_FEATURES.md README.md wiki/PHASE2_PREFLIGHT.md wiki/phase1_ship_report.md pyproject.toml config/datasources.yaml dna_decode/ scripts/quantize_fidelity_check.py LESSONS_LEARNED.md docs/ARCHITECTURE.md` returns no hits.
- `uv run pytest tests/test_data_cohort.py tests/test_pipeline_cli.py -m "not slow"` passes (D5 new tests).
- `uv run pytest tests/test_models_foundation.py -m slow -v` passes (D11 post-populate).
- `grep "hash(s.mlst)" scripts/pipeline.py` returns no matches (B-A and B-C only; B-B leaves it).
- B commit message includes the numerical before/after of clade-only AUROC (B-A and B-C only).
- Auto-memory file at `~/.claude/...user_environment.md` has GTX 860M reality recorded.

---

## Resolved Decisions (2026-05-14 via post-save /brainstorm + user lock)

All three "Open Decisions for User" from the initial save are now locked:

1. **B-resolution path** — **B-B (Skip clade-only from smoke; 4 variants).** See D1.
2. **`wiki/GATE_A_REPORT.md` in C scope** — **NO.** Line 60 is a correct historical caveat, not a stale RTX claim. See D2.
3. **Step D `--per-class` value** — **20** (N=40 with 5-strain R margin). See D12.

## Open Question (non-blocking, follow-up)

- If B-B (clade-only dropped from smoke), should `scripts/pipeline.py --include-clade-baseline` flag stay enabled at code level (just unused for smoke) or be disabled with a warning on singleton-MLST cohorts? **Current call:** keep enabled, document as "skipped for N=12 smoke." Reconsider at Stage 1 if a defensive warning would prevent foot-shooting. Not blocking for /execute-plan.
