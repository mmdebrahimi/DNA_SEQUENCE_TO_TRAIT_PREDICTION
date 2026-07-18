## Lens status
Grounded via /research-department (modality-not-scale, confirmed vs current ProteinGym SOTA) + direct codebase read. No /brainstorm or /review run on this plan yet.

## Problem Statement
The forward variant-effect hybrid's biggest measured, unrealized lever is STRUCTURE: our own N=95 sweep found `ESM2 + ProSST` +0.05 paired vs ESM2-650M (win 87%) — ~4x the evolution/MSA path built overnight (+0.013), and the current ProteinGym leaderboard + the field agree structure tokens are the top modality (ProSST zero-shot 0.504). But that +0.05 exists only as ProteinGym's PRECOMPUTED `ProSST-2048` column; our own pipeline cannot produce a ProSST score. This plan builds a ProSST structure scorer (fed by free AlphaFold structures, reusing the existing `fetch_alphafold_pdb`) that emits a `{mutation: score}` table into the already-built-and-tested `rank_average_hybrid` slot — closing the structure path the same way the evolution path was closed, and validating whether our own ProSST reproduces the column + realizes the lift.

## Codebase Context
The forward package (`dna_decode/forward/`) already has the exact pattern to mirror and the AlphaFold-fetch half done:
- `structure_scorer.py` — `fetch_alphafold_pdb(uniprot, out_dir)` + `alphafold_pdb_url()` (free AlphaFold DB download, cached) ALREADY EXIST and work; `StructureMethodUnavailable` + `_load_esm_if()` lazy-import-or-raise + `esm_if_variant_table()` + `esm_if_tier()` are the template for a table-fed structure method.
- `variant_effect.py:130 rank_average_hybrid(tables)` accepts ANY list of `{mutation: score}` tables (higher=preserved) — the hybrid slot is DONE and tested; `predict_effect(...)` dispatches methods at lines 208-243 (esm2 / alphamissense / esm_if / hybrid), each table-fed.
- `msa_evolution.py:164 evolution_table_from_scores(scores)` is the pass-through adapter (any modality table -> hybrid).
- `msa_transformer.py` + `scripts/msa_transformer_lift.py` — the just-built evolution analog (lazy model load, per-variant table, restartable JSONL validation harness, reproduce-the-column-then-test-lift) — the structure scorer + its harness mirror these exactly.
- Data: `D:/dna_decode_cache/proteingym/pg_reference.csv` has `UniProt_ID` (AlphaFold key) + `coarse_selection_type` (phenotype category); each `pg_zeroshot/*.csv` has the `ProSST-2048` validation column (confirmed) + `DMS_score` + `ESM2_650M`. ESM2-650M masked-marginal tables for 201 proteins are cached at `D:/dna_decode_cache/esm/`.

External surface (verified via ai4protein/ProSST repo + HF, 2026-07-18): structure quantization = `from prosst.structure.quantizer import PdbQuantizer; PdbQuantizer(structure_vocab_size=2048)(pdb) -> [int tokens]` (GVP encoder — needs `torch_geometric`); model = `AutoModelForMaskedLM.from_pretrained("AI4Protein/ProSST-2048", trust_remote_code=True)` + `AutoTokenizer`; scoring = masked-marginal log-ratio `score(mut)=sum logP(mut|seq,struct) - logP(wt|seq,struct)`, higher=preserved. ProteinGym ships PRE-QUANTIZED structure tokens (transformer-only path avoids `torch_geometric`). `torch_geometric` is absent on this Windows/CPU host (the same blocker that keeps ESM-IF seam-only) — so the PdbQuantizer path runs on Kaggle T4 (where ESM-IF/ESM2 ran); the pre-quantized-token path is transformer-only and CPU-feasible.

### Reusable-Code Survey
- `dna_decode/forward/structure_scorer.py::fetch_alphafold_pdb` / `alphafold_pdb_url` / `StructureMethodUnavailable` — REUSE directly (AlphaFold fetch + the raise-if-absent pattern).
- `dna_decode/forward/variant_effect.py::rank_average_hybrid` + `predict_effect` dispatch — REUSE (hybrid slot + method seam; add a `prosst` branch mirroring `esm_if`).
- `dna_decode/forward/msa_transformer.py` + `scripts/msa_transformer_lift.py` — REUSE as the structural template (module shape + validation harness).
- `dna_decode/forward/msa_evolution.py::evolution_table_from_scores` — REUSE the adapter.
- Searched: `dna_decode/forward/*.py`, `scripts/msa_transformer_lift.py`, `scripts/forward_modality_hybrid_sweep.py`, `graphify-out/GRAPH_REPORT.md` (absent).

## Pre-Change Baseline
The forward hybrid has methods {blosum62, esm2, alphamissense, esm_if, hybrid}. `esm_if` (structure via GVP-GNN) is a complete but locally-unavailable seam (torch_geometric absent) AND underperforms (PTEN 0.479 < ESM2 0.518). NO ProSST scorer exists. `fetch_alphafold_pdb` exists + works. The measured structure lift ESM2+ProSST +0.05 (N=95 sweep, win 87%) comes ONLY from ProteinGym's precomputed `ProSST-2048` column — our own pipeline cannot produce a ProSST score today. `rank_average_hybrid` accepts any table (slot ready). Frozen decoder surface: `verify_lock OK` (must stay byte-unchanged). Current forward test count for the evolution work: 33.

## Verification Signal
- **Correctness (primary):** `prosst_variant_table` produces a `{mutation: score}` table for >=1 ProteinGym protein that REPRODUCES ProteinGym's own `ProSST-2048` column at Spearman >= 0.70 (mirrors the MSA-T 0.84 correctness bar; a correct implementation of the same model+scoring should agree strongly).
- **Lift (secondary, honest):** `ESM2 + our-ProSST` hybrid beats ESM2 alone paired on the Stability/Expression subset (the categories where the sweep showed structure helps: ProSST +0.10 on Expression / +0.07 on Stability). Reported per-category with the same honest-limit framing as the MSA-T run (report final n; the powered evidence is the N=95 precomputed-column sweep).
- **Seam:** `predict_effect(method='prosst', prosst_table=...)` returns a tiered ForwardPrediction; missing deps raise `StructureMethodUnavailable`; `rank_average_hybrid([esm2, prosst])` runs; offline mock tests green.
- **Frozen surface:** `uv run python -m scripts.prospective_lock_validate` prints `verify_lock OK` after all changes.

## Implementation Steps

### Step 1: ProSST structure scorer module
Files: dna_decode/forward/prosst_scorer.py
Depends on: none

**What changes:**
- New module mirroring `structure_scorer.py`. Reuse `fetch_alphafold_pdb` + `StructureMethodUnavailable` (import from `structure_scorer`).
- `_load_prosst(vocab=2048)` — lazy-import `transformers` (`AutoModelForMaskedLM`/`AutoTokenizer` from `AI4Protein/ProSST-{vocab}`, `trust_remote_code=True`, pinned `revision`); raise `StructureMethodUnavailable` with an install hint if absent.
- `quantize_structure(pdb_path, vocab=2048)` — lazy-import `prosst.structure.quantizer.PdbQuantizer`; return the structure-token list; raise `StructureMethodUnavailable` (torch_geometric) if absent.
- `prosst_variant_table(wt_seq, mutants, structure_tokens=None, pdb_path=None, vocab=2048, model_bundle=None) -> {mutation: score}` — accept EITHER pre-quantized `structure_tokens` (transformer-only) OR a `pdb_path` (calls `quantize_structure`). Feed sequence + structure tokens to the masked LM; per mutated position score = logP(alt)-logP(wt) (masked-marginal), higher=preserved. WT-vs-sequence coordinate check (fail loudly, mirroring esm_if). Skip multi-mutants / non-standard.
- `prosst_tier(delta)` — preserved/damaging/uncertain thresholds (mirror `esm_if_tier`; on the log-ratio scale).

**Test strategy:**
- Unit: `prosst_tier` thresholds; `prosst_variant_table` raises `StructureMethodUnavailable` when `transformers`/`prosst` absent (monkeypatch import to fail). Real forward is exercised by Step 6 (Kaggle), not CI.

### Step 2: Wire method='prosst' into predict_effect
Files: dna_decode/forward/variant_effect.py
Depends on: Step 1

**What changes:**
- Add a `prosst_table: dict | None = None` param to `predict_effect`.
- Add `elif method == "prosst":` mirroring the `esm_if` branch — read `prosst_table[mutation]` (raise a clear ValueError if the table or mutation is absent), tier via `prosst_scorer.prosst_tier`, set `raw_score` + regime `B_molecular`.
- Extend the trailing `NotImplementedError` method list to include `'prosst'`.

**Test strategy:**
- Unit: `predict_effect(method='prosst', prosst_table={...})` returns the right tier on a stub table; missing table/variant raises; no other method behavior changes.

### Step 3: Package exports
Files: dna_decode/forward/__init__.py
Depends on: Step 1

**What changes:**
- Export `prosst_variant_table`, `quantize_structure`, `prosst_tier` from `prosst_scorer` (alongside the existing `structure_scorer` exports).

**Test strategy:**
- Import-smoke: `from dna_decode.forward import prosst_variant_table, prosst_tier` resolves without importing torch/transformers at package import (lazy).

### Step 4: ProSST validation harness
Files: scripts/prosst_lift.py
Depends on: Step 1

**What changes:**
- Mirror `scripts/msa_transformer_lift.py`. For ProteinGym proteins with a `UniProt_ID` + an ESM2 table + a zeroshot CSV: obtain ProSST scores via our `prosst_variant_table` (pre-quantized structure tokens when supplied via `--structure-dir`, else `fetch_alphafold_pdb` + `quantize_structure`), compute abs-Spearman vs `DMS_score` + reproduction vs the `ProSST-2048` column, and the `ESM2 + our-ProSST` hybrid via `rank_average_hybrid`. Restartable JSONL checkpoint; per-category (coarse_selection_type) aggregation; `--limit`, `--vocab`, `--only`. Focus the lift read on Stability/Expression.
- Reuse `_spearman` from `scripts.forward_blosum_proteingym_sweep`.

**Test strategy:**
- Unit on the pure helpers (per-category bucketing / paired-delta / reproduction join) with synthetic records, mirroring `tests/test_forward_msa_transformer.py`. The model forward is a Step-6 real-surface run.

### Step 5: Offline tests
Files: tests/test_forward_prosst.py
Depends on: Step 1, Step 2, Step 3

**What changes:**
- `StructureMethodUnavailable` raised when `transformers`/`prosst`/`torch_geometric` absent (monkeypatch). `predict_effect(method='prosst')` tiering on a stub table + refusal on missing table/variant. `prosst_tier` thresholds. `alphafold_pdb_url` shape (reused). `rank_average_hybrid([esm2_stub, prosst_stub])` composes. `NotImplementedError` names 'prosst'. All offline — heavy deps mocked, no network/GPU.

**Test strategy:**
- These ARE the tests; run `uv run pytest tests/test_forward_prosst.py -q` green.

### Step 6: Real-surface validation run
Files: wiki/prosst_lift_RESULTDATE.md
Depends on: Step 1, Step 4

**What changes:**
- ATTENDED run (structure stack): install `prosst` + `transformers` (+ `torch_geometric` only for the PdbQuantizer path) on Kaggle T4 — the host where ESM-IF/ESM2 ran — OR use ProteinGym pre-quantized structure tokens + the ProSST transformer on CPU. Run `scripts/prosst_lift.py` over a protein set; write `wiki/prosst_lift_<date>.{md,json}`: the reproduction Spearman (vs ProSST-2048 column) + the ESM2+ProSST hybrid per-category lift, with the same honest-limit framing as the MSA-T artifact (report final n; powered evidence = the N=95 sweep).
- This step produces the empirical verdict; it is not run by /execute-plan's local loop (Kaggle/structure-stack), same posture as the ESM-IF Kaggle validation.

**Test strategy:**
- The run IS the integration test (real model, real structures, real artifact). Success = reproduction >= 0.70 on >=1 protein + a written verdict; do not report a passing build without inspecting the numbers.

### Step 7: Docs + closeout
Files: dna_decode/forward/README.md, project_state/dna-decode-2026-05-11.md
Depends on: Step 1, Step 2, Step 3, Step 4, Step 5, Step 6

**What changes:**
- README structure section: ProSST is the strong structure model (>ESM-IF); the +0.05 structure lever (biggest single modality); the phenotype routing (structure -> stability/expression). Note the deps (prosst + transformers + torch_geometric for novel proteins; pre-quantized tokens for ProteinGym).
- Ledger row recording the build + the Step-6 verdict; note frozen surface `verify_lock OK`; honest scope (gain small, one cell, labels are the binding constraint per the negative-results map).

**Test strategy:**
- Doc-only; re-run `verify_lock` + the full forward suite green as the closeout gate.

## Execution Preview
- Wave 0: Step 1 (foundation module).
- Wave 1: Step 2, Step 3, Step 4 (all depend only on Step 1; distinct files variant_effect.py / __init__.py / scripts/prosst_lift.py — no overlap).
- Wave 2: Step 5 (deps 1-3), Step 6 (deps 1,4; attended Kaggle/CPU).
- Wave 3: Step 7 (docs, deps all).
- Total waves: 4. Max parallelism: 3 (Wave 1). Critical path: Step 1 -> Step 2 -> Step 5 -> Step 7 (also Step 1 -> Step 4 -> Step 6 -> Step 7), length 4.

## Risk Flags
- **`trust_remote_code=True`** — ProSST's HF model executes remote modeling code. PIN `revision=<commit sha>` and review the modeling file once; a supply-chain consideration, not a blocker.
- **`torch_geometric` blocked on this host** [verified] — the PdbQuantizer (novel-protein path) will not run locally (same wall as ESM-IF); mitigated by the pre-quantized-token path (transformer-only, CPU) for ProteinGym validation + Kaggle for novel proteins. Step 6 carries the run.
- **MSA-T lesson (attenuation)** — a correct scorer may reproduce the column (rank) yet the hybrid lift can attenuate with a fast approximation. Use masked-marginal + ProSST's exact protocol; the PRIMARY bar is reproduction, the lift is secondary + honestly framed.
- **Pre-quantized ProteinGym structures obtainability** [unverified] — the repo points to a Google-Drive/HF download for quantized ProteinGym structures; if unobtainable, fall back to `fetch_alphafold_pdb` + PdbQuantizer on Kaggle. Step 6 resolves which path.
- **Strategic (north-star):** the gain is small (+0.05 on ONE cell); the project's binding constraint is LABELS not models (`wiki/negative_results_map_2026-06-13.md`). This is a nice-to-have, not a north-star move — flagged for the go/no-go.

## Open Questions
- Validation substrate: use ProteinGym PRE-QUANTIZED structures (lighter, transformer-only, possibly CPU) vs `fetch_alphafold_pdb`+PdbQuantizer (heavy, Kaggle)? Recommend pre-quantized for the reproduction bar.
- Inference mode: masked-marginal (L forwards, matches ProteinGym, reproduces +0.05) vs wt-marginal (1 forward, fast, risks the MSA-T attenuation)? Recommend masked-marginal.
- Scope: build the ESM2+ProSST two-modality hybrid only, or go straight to a ready-made multimodal SOTA (SaProt / VenusREM, which already fuse structure+MSA+retrieval)? This plan does ProSST; VenusREM is a follow-on.

## Verification
- `uv run pytest tests/test_forward_prosst.py -q` green (offline seam + tiering + refusal).
- `uv run pytest tests/ -q` no regressions in the forward suite.
- `uv run python -m scripts.prospective_lock_validate` prints `verify_lock OK` (frozen surface byte-unchanged).
- Step 6 artifact `wiki/prosst_lift_<date>.md` records reproduction Spearman (>= 0.70 on >=1 protein) + the honest per-category hybrid-lift verdict.
- `dna-decode forward` / existing CLI commands unchanged (additive only).

## Save-time amendments

Captured at: 2026-07-18
Source: `/save-plan` arguments
(Audit-notes only — /execute-plan reads ONLY `## Implementation Steps`; amendments are provenance, not executable instructions.)

- ProSST structure scorer -> rank_average_hybrid
- reuse fetch_alphafold_pdb
- validate reproduce ProSST-2048 column then per-category lift
- pre-quantized+masked-marginal recommended
<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,test-strategy-leak -->
