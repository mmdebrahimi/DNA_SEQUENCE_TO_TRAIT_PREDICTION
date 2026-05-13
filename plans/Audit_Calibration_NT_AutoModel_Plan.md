# Audit Verdict Calibration + NT AutoModel Refactor

> **Status:** superseded 2026-05-13 by `Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md` after `/review` synthesis (scope reduction: dual-verdict columns dropped, wiki step dropped, NT refactor split into a separately-gated commit) and `/brainstorm` patches (D6/D7: stdout coverage + all-six-fields check; D8: retraction relocated; D9: equivalence test runs on CPU when no GPU). See the ship-path plan for the actual implementation.

> **Retraction of `audit_cohort.py` commit `b646fc9`:** that commit's stdout-emitted "GO" verdict on the 67-strain Gate B cohort was misleading — it was produced under `--target-per-drug 50 --min-minority-class 10`, well below the Phase 1 canonical thresholds (150 / 30). Under defaults, the correct verdict is **WARN** on 6 rules. The bug was caught by the post-ship adversarial review the same day. Fix-forward landed in the commit applying this ship-path plan.

> Fix a credibility bug in the just-shipped audit cohort generator (`scripts/audit_cohort.py`) — the "GO" verdict was emitted under silently-relaxed thresholds; defaults produce "WARN" — AND simultaneously replace `NucleotideTransformerModel`'s `AutoModelForMaskedLM` with `AutoModel` to eliminate the `output_hidden_states=True` workaround.

---

## Problem Statement

Two issues surfaced during the post-Gate-B-prep `/brainstorm` adversarial review (2026-05-13):

1. **Audit verdict miscalibration (C1, grounded).** `scripts/audit_cohort.py` (commit `b646fc9`) emitted "GO" for the 67-strain Gate B cohort. Verified: this was achieved by silently relaxing CLI flags to `--target-per-drug 50 --min-minority-class 10`. Under the documented Phase 1 defaults (`target_per_drug=150`, `min_minority_class=30`) the actual verdict is **WARN** on 6 rules (3 strain-count + 3 minority-class) because each drug pool has only 50 strains and the minority class is 20-24. The "GO" headline was misleading without threshold context. Tests don't lock the default semantics.

2. **NT uses `AutoModelForMaskedLM` when `AutoModel` suffices (M1/D1, grounded).** `dna_decode/models/foundation.py:239` constructs the model with `AutoModelForMaskedLM.from_pretrained(...)`. The MaskedLM head returns logits — wasted compute and VRAM for our use case (we discard logits, only use the encoder's hidden states). To access hidden states we currently set `output_hidden_states=True` (`foundation.py:259`) and index `outputs.hidden_states[-1]`. Switching to `AutoModel` returns the base encoder, which exposes `last_hidden_state` directly without the workaround.

Both issues are Phase 2.5 hardening surfaced before Gate B re-runs; neither requires F:-drive availability so they can ship while user decides about the WD Passport reformat.

---

## Design Decisions

### D1: Audit report header MUST surface threshold values

**Decision:** Add a "Thresholds applied" block in the report header next to the verdict line. Capture `target_per_drug`, `min_minority_class`, `max_pct_missing_metadata`, `min_pct_broth_microdilution`, `n50_min`, `contig_count_max`.

**Rationale:** A GO/WARN/NO-GO verdict without threshold context is uninterpretable; readers (including future-self) easily over-trust "GO" if they don't see the bar that was cleared.

**Trade-off:** Considered tagging the verdict label itself (e.g., `GO (lenient)`) — rejected because qualifiers degrade fast; explicit values are unambiguous.

### D2: Two verdict columns — Phase 1 production + Gate B infra-only

**Decision:** Compute and display TWO verdict rows side-by-side:
- **Phase 1 production thresholds:** `target_per_drug=150, min_minority_class=30` (the canonical Phase 1 ship thresholds documented in `Ecoli_G2P_Platform_Technical_Plan.md`).
- **Gate B infra-only thresholds:** `target_per_drug=12, min_minority_class=6` (matches the mini-cohort + Gate B's explicit "infrastructure dry-run, not model-quality" framing in `wiki/GATE_B_REPORT.md`).

CLI flags override both column thresholds simultaneously when explicitly set.

**Rationale:** Phase 1 has 150-strain target as canonical; Gate B mini-cohort target is much smaller. A single verdict label can't honestly cover both. Two columns make the dual-context explicit.

**Trade-off:** Considered a single verdict + a footnote — rejected because the footnote is easily missed when only the verdict label appears in summary tools or eye-tracking.

### D3: Pin default semantics in tests

**Decision:** Add `test_default_thresholds_warn_on_50_strain_cohort` (and equivalents for other thresholds) that locks the documented Phase 1 semantics. If anyone bumps `target_per_drug` default without intent, the test fires.

**Rationale:** The mis-headline happened because no test exercised default semantics. Tests today use custom-threshold cases. Regression discipline.

**Trade-off:** Considered "leave defaults flexible per project" — rejected; documented thresholds are documented for a reason, drift detection has value.

### D4: NT switches to `AutoModel`, not `AutoModelForMaskedLM`

**Decision:** Change `NucleotideTransformerModel._load_weights` (`foundation.py:237-250`) to use `AutoModel.from_pretrained(self.metadata.huggingface_id, trust_remote_code=True)`. Update `_embed_window` (`foundation.py:252-261`) to use `outputs.last_hidden_state.squeeze(0).mean(dim=0)` and drop `output_hidden_states=True`.

**Rationale:** `AutoModel` exposes the base encoder directly. Lighter VRAM (no MaskedLM head), cleaner code, no dance with `output_hidden_states` workaround. Embedding extraction has no use for the masking logits.

**Trade-off:** Considered keeping `AutoModelForMaskedLM` + dropping the flag — but the grounded review found `MaskedLM.outputs[0]` is logits, not hidden states; the flag was load-bearing for the existing path. `AutoModel` is the right shape for embedding extraction.

### D5: Pooling-strategy tag stays "single_seq_mean"

**Decision:** No change to `pooling_strategy` field added in commit `b384c96`. Both `AutoModelForMaskedLM` + `output_hidden_states[-1]` mean AND `AutoModel` + `last_hidden_state.mean` produce the same conceptual operation (mean over all token positions of the final encoder layer). The numerical output IS expected to be identical or near-identical.

**Rationale:** Avoid invalidating future caches if the refactor preserves the mean-pool semantics. The cache strategy tag captures the conceptual pooling, not the model class.

**Trade-off:** Considered tagging this as a new pooling strategy (e.g., "encoder_last_hidden_mean") — rejected because conceptually identical operations shouldn't fragment cache history. If the numerical output diverges in practice, revisit with a real equivalence test.

---

## Implementation Plan

### Step 1: Switch NT to `AutoModel`

**Files:** `dna_decode/models/foundation.py`

**Changes:**
- Replace `AutoModelForMaskedLM` import + usage (line 239 + 246) with `AutoModel`.
- Update `_embed_window` (line 252-261) to use `outputs.last_hidden_state.squeeze(0).mean(dim=0)` and drop `output_hidden_states=True`.
- No change to NT class metadata (huggingface_id, embedding_dim, max_context).

### Step 2: Verify NT refactor

**Files:** `tests/test_models_foundation.py` (extend), or new ad-hoc test invocation.

**Changes:**
- Existing tests (`test_factory_nucleotide_transformer_from_config`, `test_mock_does_not_eagerly_load_weights`) still pass — they cover construction + lazy-load, not real inference.
- Optional: add a "shape sanity" test that loads NT in CPU mode + runs `embed("ACGT" * 50)` → asserts shape `(512,)`. Slow test; mark with `pytest.mark.slow` if added.

### Step 3: Audit report header — emit thresholds

**Files:** `scripts/audit_cohort.py`

**Changes:**
- `build_report()` adds a "Thresholds applied" section after the timestamp:
  ```
  **Thresholds applied:**
  - target_per_drug: 150
  - min_minority_class: 30
  - max_pct_missing_metadata: 20.0
  - min_pct_broth_microdilution: 80.0
  - n50_min: 50000
  - contig_count_max: 500
  ```
  Values come from the `VerdictRules` instance passed in (so CLI flags reflect through automatically).
- No change to argparse; CLI flags already override these via `VerdictRules` construction in `main()`.

### Step 4: Dual-verdict columns

**Files:** `scripts/audit_cohort.py`

**Changes:**
- Define a new `PHASE1_PRODUCTION_RULES = VerdictRules()` (canonical 150/30 defaults).
- Define a new `GATE_B_INFRA_RULES = VerdictRules(target_per_drug=12, min_minority_class=6, ...)`.
- `evaluate_verdict()` already accepts arbitrary rules — reused for both contexts.
- `verdict_section()` becomes a side-by-side table comparing both:
  ```
  | Rule | Phase 1 prod | Gate B infra |
  |---|---|---|
  | strain count (cipro) | WARN (50 < 150) | PASS (50 ≥ 12) |
  ```
- CLI flags override BOTH columns when explicitly set (so users running the audit at custom thresholds still get one column under their custom rules + one under the OTHER documented context — keep both columns always visible).

**Constraint:** the existing single-verdict CLI behavior + `print(f"[audit_cohort] verdict: {verdict.verdict}")` stays — the CLI echo defaults to the Gate B verdict (since this audit tool's primary user is Gate B prep). Phase 1 verdict is in the report body only.

### Step 5: Lock default semantics in tests

**Files:** `tests/test_audit_cohort.py`

**Changes:**
- New test: `test_default_thresholds_warn_on_50_strain_cohort` — build a 50R/50S cohort and verify `evaluate_verdict(... VerdictRules())` returns WARN on strain-count + minority-class rules.
- New test: `test_audit_report_includes_threshold_values_in_header` — assert the report body contains the threshold lines.
- New test: `test_dual_verdict_columns_render_correctly` — assert both Phase 1 and Gate B columns appear in the rendered table.

### Step 6: Re-run audit + update GATE_B_REPORT.md

**Files:** `reports/gate_b_cohort_audit.md` (regenerated; gitignored), `wiki/GATE_B_REPORT.md` (update results section).

**Changes:**
- Re-run `audit_cohort.py` against the 67-strain cohort with default CLI flags (no threshold relaxation).
- Update `wiki/GATE_B_REPORT.md`'s pre-existing references to "GO" → updated dual-column verdict.

### Step 7: Commit + push

Single commit covering all of the above:
- `feat(phase2): audit dual-verdict columns + threshold disclosure + NT AutoModel refactor`
- Mention the C1 credibility fix from /brainstorm in the commit body.

---

## Verification

1. `uv run pytest tests/ -q` → 341 + 3 new tests (Step 5) = 344, 0 failed, 1 skipped.
2. `uv run python scripts/audit_cohort.py --cohort data/processed/gate_b_cohort.parquet --output reports/gate_b_audit.md` (no `--target-per-drug` etc.) → report contains "Thresholds applied" block AND dual-column verdict table AND Gate B verdict = "GO", Phase 1 verdict = "WARN".
3. Manual: open `reports/gate_b_audit.md`, confirm Phase 1 column WARN ≠ Gate B column GO.
4. NT smoke (requires storage + GPU available — defer until F: / new drive):
   ```bash
   uv run python -c "from dna_decode.models.foundation import model_factory; m = model_factory('nucleotide_transformer', device='cuda'); m._ensure_loaded(); import numpy as np; emb = m.embed('ACGT'*50); assert emb.shape == (1, 512); print('NT AutoModel embed OK')"
   ```
5. `git log --oneline -n 3` → single commit ahead of `b646fc9`.

**Scope boundaries:**
- This plan does NOT fix `load_bvbrc_ast` performance (Phase 2.5 TODO #1; needs benchmarking first per /brainstorm).
- This plan does NOT touch batched inference (`cache.populate` + `embed_batch`); that's the bigger /brainstorm Option B refactor.
- This plan does NOT depend on F: drive reconnect or WD Passport reformat — fully storage-independent.
