# Audit Calibration + NT AutoModel — Ship Path Plan

> **Status:** Commit 1 shipped 2026-05-13 at `473b8eb`. Commit 2 (NT AutoModel refactor) **deferred indefinitely** 2026-05-13 — equivalence gate failed: `AutoModel.from_pretrained` cannot load the NT v2 100M checkpoint (state_dict shape mismatch on `Linear[4096, 512]` vs `Linear[512, 2048]` in InstaDeep's trust_remote_code modeling code). Per D3 gating rule, Commit 2 ships only after the equivalence test passes locally. See `TODOS.md` "Phase 2.5 perf hardening" for the deferral context.

> Scope-reduced delta from `Audit_Calibration_NT_AutoModel_Plan.md` after `/review` synthesis (2026-05-13). Reduces 5 steps / 4 waves → 2 commits, 2 waves. Drops dual-verdict columns (institutionalization risk) and the wiki update (no text to replace). Splits the bundled NT refactor into a separate, gated commit.

---

## Problem Statement

The original technical plan (`Audit_Calibration_NT_AutoModel_Plan.md`) bundled two unrelated fixes into one commit and over-engineered the audit-fix presentation:

1. **Dual-verdict columns** (Step C in original) risked institutionalizing the relaxed Gate-B-infra threshold profile as a "supported configuration" rather than a one-time concession. Naming `GATE_B_INFRA_RULES` as a module constant elevates "temporary infra-bringup concession" to "documented profile." The mixed CLI semantics (flags affect Phase 1 column but not Gate B) is "the worst of both" — users misread the Gate B column as "what they asked for."
2. **Wiki update** (Step E in original) targets `wiki/GATE_B_REPORT.md` which contains zero `"GO"` / `"WARN"` / `"verdict"` references — its Results section is all TBD placeholders. The step was busywork.
3. **NT AutoModel refactor** (Step A in original) is cosmetic cleanup with no shipping urgency; bundling violates single-purpose-commit discipline. The cache-compat claim in D5 (numerical equivalence between `AutoModelForMaskedLM.hidden_states[-1]` and `AutoModel.last_hidden_state`) is unverified for `trust_remote_code=True` models — author-controlled modeling code can diverge. Without an empirical equivalence test, shipping the refactor creates hidden hybrid-cache risk.

The /review synthesis (CEO + Eng lenses, 2026-05-13) converged on scope reduction with two specific changes:
- Replace dual-column rendering with an asymmetric warning banner on relaxed flags.
- Split or defer the NT refactor; require an equivalence test before it ships.

---

## Design Decisions

### D1: Asymmetric warning banner replaces dual-verdict columns

**Decision:** Drop the `PHASE1_PRODUCTION_RULES` / `GATE_B_INFRA_RULES` constants. Single canonical verdict against passed `rules`. If `rules.target_per_drug < 150` OR `rules.min_minority_class < 30`, prepend a top-of-report banner: `⚠️ WARNING: NON-DEFAULT THRESHOLDS — VERDICT NOT COMPARABLE TO PHASE 1 GATE`.

**Rationale:** A symmetric two-column presentation legitimizes the relaxed profile as a "supported configuration" — exactly the failure mode the audit was built to prevent. An asymmetric warning makes the deviation unmissable without elevating the relaxed profile to first-class status. The threshold-disclosure block from Step B already makes any verdict self-interpreting.

**Trade-off:** Considered keeping dual columns for "Gate B reviewers who want a quick reference" — rejected because the named constant + test + documentation surface area institutionalizes a one-time concession.

### D2: Drop wiki/GATE_B_REPORT.md update entirely

**Decision:** Remove Step E. Add a single Reliability Notes bullet retracting the `b646fc9` GO verdict instead of trying to "replace single-verdict references" that don't exist in the file.

**Rationale:** Grep verified the file has no GO/WARN/verdict text to replace. The original step would have been pure overhead. The bad-verdict-already-shipped concern is addressed by the commit body + a single appended line.

**Trade-off:** Considered nothing — the original step was based on a wrong assumption about the file contents.

### D3: Split NT refactor into a separate, gated commit

**Decision:** Commit 1 covers the audit calibration only. NT refactor (original Step A) becomes Commit 2 and ships ONLY IF an `np.allclose` equivalence test passes between `AutoModelForMaskedLM.hidden_states[-1]` and `AutoModel.last_hidden_state` on a fixed 200-bp input. If GPU + storage aren't available to run the test, Commit 2 is deferred indefinitely.

**Rationale:** Bundling unrelated changes hurts `git bisect`. The NT refactor's cache-compat claim is unverified for `trust_remote_code` models. Without empirical equivalence verification, the refactor trades a cosmetic cleanup for hidden hybrid-cache risk on future NT runs. Equivalence test is the trigger for Commit 2 being safe to ship at all.

**Trade-off:** Considered shipping NT refactor unverified for the "all sibling models already use AutoModel" reason — rejected because sibling pattern alignment isn't urgent and the verification cost is low (~10 min once storage returns).

### D4: Keep Step D's default-semantics test as the regression lock

**Decision:** Add `test_default_thresholds_warn_on_undersized_cohort` to `tests/test_audit_cohort.py`. Also add `test_audit_report_includes_thresholds_in_header` and `test_warning_banner_on_relaxed_flags`. Drop the dual-column test (no dual columns).

**Rationale:** The original miscalibration happened because no test exercised default semantics. Locking the documented Phase 1 thresholds prevents silent drift.

**Trade-off:** None — universally accepted.

### D5: thresholds_block(rules) helper, not inline string list

**Decision:** Extract the threshold-formatting code into a `thresholds_block(rules: VerdictRules) -> list[str]` helper consistent with the existing `*_section()` helpers in `audit_cohort.py`.

**Rationale:** Matches existing module conventions. Easier to test in isolation.

**Trade-off:** Adds one function. Marginal — but the inline-list pattern was already inconsistent with the rest of the file.

### D6: stdout echo must also carry threshold context (post-/brainstorm patch C1)

**Decision:** `audit_cohort.py:main()` line ~500 currently emits `[audit_cohort] verdict: {verdict.verdict}` — bare verdict with zero threshold context. The original ship-path plan added a banner only to the markdown report; the stdout escape hatch remained. Patch: emit a SECOND line before the verdict echo when any threshold deviates from Phase 1 defaults: `[audit_cohort] WARNING: relaxed thresholds (target_per_drug=50, min_minority_class=10) — verdict not comparable to Phase 1 gate`.

**Rationale:** Stdout consumers (CI / pipeline / terminal viewers) propagate the bare verdict; report-only banners don't reach them. Stdout coverage must match report coverage for the calibration fix to be complete.

**Trade-off:** Slightly noisier CLI output when non-default thresholds are used. Acceptable — the noise is the point.

### D7: Check all 6 VerdictRules fields, not just 2 (post-/brainstorm patch C2)

**Decision:** Original plan checked only `target_per_drug < 150` OR `min_minority_class < 30`. `VerdictRules` has 6 fields. Patch: introduce `PHASE1_DEFAULT_RULES = VerdictRules()` constant. Helper `non_default_threshold_fields(rules, default=PHASE1_DEFAULT_RULES) -> list[str]` returns the field names where `rules` is more permissive than default (e.g., higher `max_pct_missing_metadata`, lower `min_pct_broth_microdilution`). Banner triggers if non-empty AND lists the deviated fields by name with canonical defaults.

**Rationale:** CLI exposes `--max-pct-missing-metadata` and `--min-pct-broth-microdilution`; relaxing those would silently bypass the original banner. Bug-class completeness requires covering all six.

**Trade-off:** Marginally more complex banner logic. The complexity is appropriate — the banner exists to prevent silent threshold relaxation, so it must cover all relaxable fields.

### D8: Retraction target is the superseded plan + commit body, NOT GATE_B_REPORT.md (post-/brainstorm patch C3)

**Decision:** `wiki/GATE_B_REPORT.md` is explicitly the 12-strain mini-cohort plan (per file header line 5). The bad `b646fc9` GO verdict was on the 67-strain audit cohort. Placing 67-strain retraction in a 12-strain plan creates audit-trail confusion. Patch: append a "Superseded by" closing addendum to `plans/Audit_Calibration_NT_AutoModel_Plan.md` (the original plan, now superseded by this ship-path plan) noting the b646fc9 misstep + the fix-forward commit. Drop the wiki/GATE_B_REPORT.md update entirely.

**Rationale:** Annotation lives next to the artifact that produced the bug (the original technical plan). Audit trail flows: original tech plan → superseded addendum → ship-path plan → execute-plan → commit. Future readers can trace the misstep back to its origin without confusion.

**Trade-off:** None — original wiki target was the wrong file.

### D9: NT equivalence test must run on CPU when no GPU available (post-/brainstorm patch M1+M2)

**Decision (Commit 2 only):** Original plan used `@pytest.mark.skipif(not torch.cuda.is_available())`, which silently skips on CPU CI / no-GPU machines — the gate becomes invisible. Patch: REMOVE skipif. Test runs on whichever device is available (`torch.cuda.is_available() ? "cuda" : "cpu"`). Mark `@pytest.mark.slow` only.

**Test scope expansion:** Multi-input matrix (5 inputs) covering edge cases: short (`"ACGT"*50`), long (`"A"*200`), ambiguous bases (`"N"*50 + "ACGT"*38`), GC-skewed (`"GCGCGC"*34`), and one input exceeding `max_context` to exercise sliding-window. Assertion starts strict: `torch.allclose(a, b, rtol=1e-6, atol=1e-7)`. Loosen ONLY on observed divergence with a documented reason. Test also records max-abs-diff + cosine similarity in stdout for audit visibility.

**Rationale:** A silently-skippable gate doesn't prevent the cache-hybridization bug it exists to catch. CPU mode is slow (~30 s for NT 100M load + 200-bp inference) but feasible. Multi-input matrix catches tokenization edge cases the original single-input test missed.

**Trade-off:** Slower test (CPU + 5 inputs). Acceptable — the test runs only when Commit 2 is being shipped, which is gated on this test passing anyway.

---

## Implementation Plan

### Commit 1 — `feat(audit): emit thresholds + warn-on-relaxed-flags + lock default semantics`

#### Step 1: `thresholds_block` helper + report header
Files: `scripts/audit_cohort.py`

- Add `thresholds_block(rules: VerdictRules) -> list[str]` that returns 8 lines: "**Thresholds applied:**" header + 6 bullet rows (target_per_drug, min_minority_class, max_pct_missing_metadata, min_pct_broth_microdilution, n50_min, contig_count_max) + blank trailer.
- `build_report()`: insert `lines += thresholds_block(rules)` after the timestamp lines, before the cohort overview.
- No change to `argparse`.

#### Step 2: Warn-on-relaxed-flags banner
Files: `scripts/audit_cohort.py`

- Add `relaxed_flags_warning(rules: VerdictRules) -> list[str]` that returns the warning banner if `rules.target_per_drug < 150` OR `rules.min_minority_class < 30`, else empty list.
- `build_report()`: insert `lines += relaxed_flags_warning(rules)` immediately after the title and before the timestamp lines.
- Warning text: `⚠️ **WARNING: NON-DEFAULT THRESHOLDS APPLIED — VERDICT NOT COMPARABLE TO PHASE 1 GATE.**` plus a line noting which threshold(s) are relaxed and what the canonical defaults are.

#### Step 3: Lock default semantics + new helper tests
Files: `tests/test_audit_cohort.py`

- `test_default_thresholds_warn_on_undersized_cohort` — 50R/50S cohort + `VerdictRules()` defaults → verdict == "WARN" with at least one rule referencing `target 150`.
- `test_audit_report_includes_thresholds_in_header` — assert `"target_per_drug: 150"` AND `"min_minority_class: 30"` present in default-rules report.
- `test_warning_banner_on_relaxed_flags` — pass `VerdictRules(target_per_drug=50)`; assert banner text present in report. Also pass `VerdictRules()` defaults; assert banner absent.
- `test_warning_banner_lists_relaxed_thresholds` — banner text mentions which threshold(s) deviated and the canonical default value.

#### Step 4: Retract the bad `b646fc9` GO verdict
Files: `wiki/GATE_B_REPORT.md`

- Add a single line to the "Reliability notes" section: `2026-05-13: audit_cohort.py commit b646fc9 produced a misleading "GO" verdict on this cohort by silently relaxing thresholds via CLI flags. Fix-forward applied in <new-commit-sha> — see commit body. The 67-strain cohort verdict against Phase 1 defaults is WARN (50 strains/drug < 150 target; minority 20-24 < 30 min).`

### Commit 2 — `refactor(foundation): NT uses AutoModel instead of AutoModelForMaskedLM` (gated)

#### Step 5: NT AutoModel refactor
Files: `dna_decode/models/foundation.py`

- Line 239: `from transformers import AutoModel, AutoTokenizer` (was: `AutoModelForMaskedLM`).
- Line 246-248: `AutoModel.from_pretrained(self.metadata.huggingface_id, trust_remote_code=True).to(self._device).eval()`.
- Line 252-261 (`_embed_window`): drop `output_hidden_states=True`; use `outputs.last_hidden_state.squeeze(0).mean(dim=0)`.

#### Step 6: NT numerical-equivalence test (gates Commit 2)
Files: `tests/test_models_foundation.py`

- New test `test_nt_automodel_vs_maskedlm_numerical_equivalence` marked `@pytest.mark.slow` and `@pytest.mark.skipif(not torch.cuda.is_available(), reason="GPU + storage required")`.
- Load NT under both APIs on a fixed 200-bp input ("ACGT" * 50); assert `np.allclose(autoModel_embed, autoModelForMaskedLM_embed, rtol=1e-4)`.
- If this test fails (numerical divergence), DO NOT SHIP Commit 2 — update D5's pooling-strategy tag instead and treat it as a new pooling strategy.

**Gating rule:** Commit 2 ships only after Step 6 passes locally. If GPU + storage aren't available, Commit 2 is deferred indefinitely.

---

## Verification

### Commit 1 verification

1. `uv run pytest tests/ -q` → 341 + 4 new tests = 345 passed, 1 skipped.
2. Regenerate the Gate B audit report:
   ```bash
   uv run python scripts/audit_cohort.py \
     --cohort data/processed/gate_b_cohort.parquet \
     --output reports/gate_b_audit.md
   ```
   With default flags: report contains "Thresholds applied" block, NO warning banner, single verdict line = "WARN", stdout echo `[audit_cohort] verdict: WARN`.
3. Re-run with relaxed flags:
   ```bash
   uv run python scripts/audit_cohort.py \
     --cohort data/processed/gate_b_cohort.parquet \
     --output reports/gate_b_audit_relaxed.md \
     --target-per-drug 50 --min-minority-class 10
   ```
   Output: warning banner present at top of report, threshold block shows the relaxed values, verdict line = "GO".
4. `git log --oneline -n 1` shows Commit 1 on top of `b646fc9`.

### Commit 2 verification (only if Commit 2 ships)

5. `uv run pytest tests/test_models_foundation.py -m slow -v` → equivalence test PASSES.
6. NT shape sanity:
   ```bash
   uv run python -c "
   from dna_decode.models.foundation import model_factory
   m = model_factory('nucleotide_transformer', device='cuda')
   m._ensure_loaded()
   emb = m.embed('ACGT' * 50)
   assert emb.shape == (1, 512)
   print('NT AutoModel embed OK:', emb.shape)
   "
   ```
7. `git log --oneline -n 2` shows both commits on top of `b646fc9`.

### Scope boundaries

- This plan does NOT add the dual-verdict columns (deliberately dropped per D1).
- This plan does NOT touch `load_bvbrc_ast` perf (Phase 2.5 TODO #1; separate brainstorm needed first).
- This plan does NOT depend on F: drive reconnect or WD Passport reformat for Commit 1. Commit 2 is gated on storage + GPU returning.
