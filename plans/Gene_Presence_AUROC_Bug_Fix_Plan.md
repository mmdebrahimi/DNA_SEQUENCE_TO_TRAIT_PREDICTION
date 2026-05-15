# Gene-Presence AUROC=0.000 Bug Fix Plan

> Strengthen the diagnostic, confirm the strain-unique-identifier-domination hypothesis on real data, then add a `gene_symbol` column to `AnnotationTable` so the gene-presence smoke variant returns a non-degenerate AUROC at N=12.

---

## Problem Statement

The 12-strain cipro smoke gate (`scripts/smoke_gate_12strain_cipro.py`) reports gene-presence + XGBoost AUROC=0.000 / AUPRC=0.394 under LOSO at N=12 (6R/6S). NT-XGBoost (0.750) and k-mer (0.694) on the same cohort behave normally. `calibrate=False` is already set for all three variants, so the earlier `CalibratedClassifierCV` overcorrection bug is not the cause.

AUPRC near the 0.5 base rate combined with AUROC=0.000 indicates the predictions are **systematically rank-inverted**, not random. A synthetic check (12×5000 random Bernoulli, no class signal) returned AUROC=0.611, disproving the LOSO base-rate-inversion structural hypothesis. The bug is data-content-specific.

Top remaining hypothesis (verified-likely after reading `parse_gff3`): **strain-unique GFF3 `ID=` attribute domination of the gene-presence vocabulary**.

Evidence chain:
- `_extract_gene_ids` falls back `gene_id or locus_tag`.
- `parse_gff3` extracts `gene_id` from `attr_map.get("ID", attr_map.get("Name", ""))` — the GFF3 `ID=` attribute, NOT the `gene=` (gene-symbol) attribute.
- RefSeq GFF3 conventionally emits `ID=gene-b0001;Name=thrL;gene=thrL;locus_tag=b0001` → `gene_id` becomes `gene-b0001` (strain-unique ID).
- Under LOSO, training vocab is dominated by strain-unique IDs of the 11 training strains; held-out strain's `gene-*` IDs don't appear in training vocab → test rows are near-empty → XGBoost predictions collapse onto leaves systematically positioned against the held-out label.

Non-goal: making gene-presence "work well" at N=12. Goal: get a meaningful (≠ 0.000) signal so the smoke comparison surfaces a real story, and surface any irreducible-at-N=12 limitation honestly. Smoke verdict (PASS, +5.6 pp NT lift over k-mer) is unaffected by this fix.

## Design Decisions

### D1: Add `gene_symbol` column to AnnotationTable; do NOT rewrite `gene_id`

**Decision:** Extend `ANNOTATION_COLUMNS` with a new `gene_symbol` field. Populate it from GFF3 `gene=` attribute in `parse_gff3` and from `quals.get("gene", [""])[0]` in `parse_genbank` (which currently puts gene symbol into `gene_id` — fix that asymmetry while we're here). Have the smoke runner's `_extract_gene_ids` use `gene_symbol or gene_id or locus_tag`.

**Rationale:** `gene_id` is used as the per-CDS embedding cache key in `extract_cds_sequences` at `dna_decode/data/annotations.py:168`. Rewriting `gene_id` to emit gene symbol (e.g., `thrL`) would invalidate the existing HDF5 embedding cache AND create cache-key collisions across paralogs / multi-copy genes that share a symbol. Adding a separate column keeps cache semantics intact and consolidates the existing GFF3-vs-GenBank parser asymmetry.

**Trade-off:** Considered a smoke-only inline reparse helper (parse `gene=` directly from raw GFF3 lines inside the smoke runner). Rejected — adds parser duplication and the `parse_genbank` asymmetry remains a latent bug. Schema change is the cleaner consolidation.

### D2: Strengthen the diagnostic before mounting F: drive

**Decision:** Update `scripts/diagnose_gene_presence_auroc.py` to print:
- Absolute counts: `len(test_gene_set)`, `te_in_vocab`, `X_test.sum()`, `n_all_zero_test_rows`
- Per-prefix namespace breakdown of identifier patterns: counts matching `gene-`, `cds-`, `WP_`, `b\d+`, symbol-shaped `[a-z]{3,4}[A-Z]`
- First 10 example IDs per strain
- **Side-by-side overlap**: run the same fold loop twice — once with current `_extract_gene_ids`, once with a local raw-GFF3 reparse extracting `gene=` directly. Report overlap fraction + AUROC for both extractors in the same script.

**Rationale:** "Low train-vocab overlap" alone cannot uniquely confirm RefSeq ID/locus-tag domination — it could also come from parse failure, mismatched CDS feature types, plasmid/accessory skew, or paralog renaming. The side-by-side AUROC under the two extractors is the disambiguating signal: if `gene=` extraction restores AUROC and current extraction doesn't, the hypothesis is confirmed.

**Trade-off:** Adds ~50 LOC to the diagnostic. Rejected the "run minimal diagnostic, then iterate" approach because F: drive is removable and mounting it for multiple iterations costs more than a single thorough run.

### D3: Strengthen the synthetic falsifier too

**Decision:** Extend `scripts/diagnose_gene_presence_synthetic.py` (or add a second function) to model the suspected failure shape: strain-unique feature blocks per strain + a small shared core, then LOSO. Also test the explicit all-zero held-out row case against varied train matrices.

**Rationale:** The current synthetic (dense-ish random Bernoulli with stable train/test vocabulary) does NOT model the near-OOV-on-LOSO failure mode. The AUROC=0.611 result from the current synthetic rules out one structural artifact but cannot rule out XGBoost-leaf/default-direction artifacts under near-empty test rows.

**Trade-off:** Adds work before the real-data diagnostic. Accepted because the strain-unique-blocks synthetic provides a clean reproducer that doesn't depend on F: drive availability.

### D4: Add `INDETERMINATE_IDENTIFIER_OOV` smoke verdict as a guardrail

**Decision:** In `scripts/smoke_gate_12strain_cipro.py`, if median fold test-vocab overlap is below a threshold (e.g., < 0.20), tag the gene-presence variant verdict as `INDETERMINATE_IDENTIFIER_OOV` in the result packet instead of reporting AUROC as a biological/model result.

**Rationale:** Defense-in-depth. Even after the fix, future cohorts could re-hit this trap (e.g., Bakta-annotated genomes mixed with RefSeq-annotated genomes). Loud failure beats silent zero.

**Trade-off:** None significant. Threshold needs to be picked once and may need tuning later.

## Implementation Plan

1. **Strengthen the diagnostic script.**
   - File: `scripts/diagnose_gene_presence_auroc.py` (modify in place)
   - Add: absolute counts (`len(test_gene_set)`, `te_in_vocab`, `X_test.sum()`, `n_all_zero_test_rows`)
   - Add: per-prefix namespace breakdown (`gene-`, `cds-`, `WP_`, `b\d+`, symbol-shaped `[a-z]{3,4}[A-Z]`)
   - Add: first 10 example IDs per strain printed once at startup
   - Add: side-by-side dual-extractor loop — current `_extract_gene_ids` vs raw-GFF3-`gene=` reparse — reporting overlap fraction + LOSO AUROC for each extractor

2. **Extend the synthetic falsifier.**
   - File: `scripts/diagnose_gene_presence_synthetic.py` (modify in place; add a second `main()`-style function or `--mode` flag)
   - Add: strain-unique-feature-blocks synthetic — each strain gets ~500 unique identifiers + ~50 shared core — then LOSO with the same `train_xgboost_classifier(..., calibrate=False)` path
   - Add: explicit all-zero held-out row test case
   - Report AUROC + inverted AUROC + prediction value distribution for each synthetic variant

3. **Mount F: drive; run both diagnostics.**
   - `uv run python scripts/diagnose_gene_presence_synthetic.py` (no F: needed — confirms / refines structural hypothesis)
   - `uv run python scripts/diagnose_gene_presence_auroc.py` (F: needed — confirms data-content hypothesis)
   - Decision branch based on side-by-side AUROC under the two extractors:
     - If `gene=` extractor restores AUROC > 0.5 → D1 fix path; proceed to Step 4
     - If neither extractor restores AUROC → investigate next-most-likely causes (XGBoost overfitting; lineage-confounded distributions); revise plan
     - If overlap is structurally too low under both extractors → escalate to D4 `INDETERMINATE_IDENTIFIER_OOV` verdict only; no schema change

4. **(Conditional on Step 3) Schema change: add `gene_symbol` column.**
   - File: `dna_decode/data/annotations.py`
     - Extend `ANNOTATION_COLUMNS` tuple with `"gene_symbol"` (after `gene_id`, before `locus_tag` or at the end — pick a position consistent with downstream readers)
     - In `parse_gff3` (around line 92): populate `gene_symbol` from `attr_map.get("gene", "")`
     - In `parse_genbank` (around line 124): populate `gene_symbol` from `quals.get("gene", [""])[0]`; leave `gene_id` semantics intact (re-derive from a stable per-CDS identifier — `quals.get("protein_id", [""])[0]` or feature index, per existing convention)
   - File: `scripts/smoke_gate_12strain_cipro.py`
     - Modify `_extract_gene_ids` to prefer `gene_symbol` then fall back: `gid = row.get("gene_symbol") or row.get("gene_id") or row.get("locus_tag")`

5. **(Conditional on Step 3) Add `INDETERMINATE_IDENTIFIER_OOV` guardrail.**
   - File: `scripts/smoke_gate_12strain_cipro.py`
     - In `run_gene_presence_xgboost`: compute median per-fold `te_in_vocab_frac`; if below 0.20, return a result dict with verdict `INDETERMINATE_IDENTIFIER_OOV` instead of AUROC
     - In `write_packet`: detect this verdict and render `INDETERMINATE_IDENTIFIER_OOV` row in the per-variant results table; do not include it in gap analysis

6. **Re-run the smoke gate; verify gene-presence AUROC moves off 0.000.**
   - `uv run python scripts/smoke_gate_12strain_cipro.py`
   - Expected: gene-presence AUROC in roughly the 0.4-0.7 range (smoke noise floor is ±0.19 at N=12) — anywhere except perfectly anti-predictive
   - Expected: NT-XGBoost AUROC unchanged at 0.750 (no upstream changes)
   - Expected: PASS verdict preserved

7. **Update tests + docs.**
   - File: `tests/` — add unit test for `parse_gff3` populating `gene_symbol` from `gene=` attribute; add unit test for `parse_genbank` populating `gene_symbol` from `gene` qualifier; add regression test for the smoke runner's `_extract_gene_ids` preferring `gene_symbol`
   - File: `TODOS.md` — close the `[OPEN] Gene-presence + XGBoost returns AUROC=0.000 at N=12 even without calibration` entry with link to this plan + final AUROC
   - File: `LESSONS_LEARNED.md` — add an entry under 2026-05-14: "GFF3 `ID=` attribute is strain-unique by construction; do not use as gene-family vocabulary for cross-strain presence/absence features. Use `gene=` (symbol) with `locus_tag` fallback."
   - File: `CLAUDE.md` — update the GFF3 annotation-source gotcha to mention the new `gene_symbol` column

## Verification

- `scripts/diagnose_gene_presence_auroc.py` side-by-side run: current extractor shows AUROC ≈ 0.000 AND median `te_in_vocab_frac` is low (< 0.10); `gene=` extractor shows AUROC > 0.5 AND median `te_in_vocab_frac` rises substantially (e.g., > 0.50).
- `scripts/diagnose_gene_presence_synthetic.py` strain-unique-blocks variant reproduces AUROC ≈ 0.000 or AUROC ≈ 0.5 under the OOV failure shape, depending on XGBoost leaf-default behavior — confirms the structural mechanism.
- Smoke gate re-run: NT 0.750 preserved; gene-presence AUROC ∈ (0.0, 1.0) non-degenerate; PASS verdict preserved.
- All existing tests pass after schema change (`uv run pytest tests/ -v`).
