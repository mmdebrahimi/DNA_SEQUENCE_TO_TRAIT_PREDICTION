# HANDOFF — EP-4 representation bake-off — 2026-06-02

> Continues `HANDOFF_ep4_pathotype_model_2026-05-31.md`. Today resolved the stalled k-mer arm, ran a full-spectrum diagnostic, and launched the windowed-NT arm. Read the 2026-05-31 handoff for substrate/environment facts (ENA bare-URL fetch, cached genomes on C:, venv path). D: is BACK ONLINE today (`HF_HOME=D:\hf_cache`, NT weights cached there).

## Results landed today (ExPEC vs EPEC, N=24 LOSO, study-confounded Salipante-vs-Hazen)

| Representation | LOSO AUROC | note |
|---|---|---|
| k-mer k=8 **top-10k counts** (canonical bake-off) | **0.514** | chance; non-degenerate (24 distinct scores) |
| k-mer k=8 **full-spectrum counts** | 0.604 | +0.090 |
| k-mer k=8 **full-spectrum presence/absence** | **0.729** | +0.215 — BEST CLASSICAL = the bar to beat |

- Artifacts: `research_outputs/pathotype_bakeoff_kmer_expec_epec_2026-05-31.json` (0.514) + `research_outputs/pathotype_kmer_fullspectrum_diag_2026-06-02.json` (diagnostic).
- **Why the k-mer arm "stalled twice":** NOT a crash — `run_kmer_xgboost_loso` does within-fold vocab rebuild → ~5.5e9 pure-Python iters, longer than a session. Fixed with a cached-per-genome LOSO in `scripts/pathotype_kmer_bakeoff.py` (~45× faster, feature-level equivalence asserted vs the canonical runner). Use that path going forward.
- **Diagnostic verdict:** the 0.514 was UNDER-POWERED by top-N-most-frequent vocab (selects conserved core genome; discards accessory-genome discriminators). Full-spectrum recovers signal. **presence/absence ≫ counts (0.729 vs 0.604) ⇒ the signal is gene PRESENCE, not abundance.**

## Live decisions / re-pinned contract
- **The FM bar is 0.729 (best classical), NOT 0.514.** Per CLAUDE.md ">=3pp over classical", NT must reach **≥~0.76** LOSO to justify the foundation-model premise. The 2026-05-31 handoff's "beat k-mer by ≥3pp" is re-pinned to 0.729.
- **Dilution prediction:** counts(0.604) ≪ presence(0.729). NT *mean*-pool is a soft-average (analogue of counts) → predict ~0.60; NT *max*-pool is presence-aligned → predict higher. Running BOTH to test this.
- **Confound NOT resolved:** 0.514 only ruled out trivial bulk-composition batch. Accessory-genome presence/absence still carries lineage/study batch. Any NT win is necessary-not-sufficient — needs lineage/ST-aware splits before claiming biology.

## ✅ RESULT: windowed-NT arm FAILED (task `bljrvirar` completed 2026-06-02)

| NT pooling | LOSO AUROC | vs classical 0.729 |
|---|---|---|
| mean-pool | **0.382** | −0.347 |
| max-pool | **0.375** | −0.354 |

- **VERDICT: FM_FAILS.** Pooled whole-genome NT is ≤ chance and ~0.35 BELOW the classical presence/absence bar (0.729). Both arms fail; **max-pool did NOT recover the signal** (the "max tracks presence" sub-hypothesis is refuted — max ≈ mean).
- **The 0.38 (<0.5) is the no-signal LOSO artifact, not genuine anti-prediction:** with balanced 12/12 LOSO, holding out a class leaves the training set with the opposite-class majority → predictions biased away from the held-out label → AUROC dips just under 0.5 when there's no real signal. So read it as "no usable signal," not "inversely predictive." The k-mer arms escaped this because they HAD signal (0.60–0.73) overwhelming the prior-flip.
- **tet precedent CONFIRMED:** pooling (mean OR max) dilutes the localized gene-presence signal that the classical presence/absence baseline captures.
- **Scope of the claim:** this condemns POOLED whole-genome NT only. It does NOT test **per-CDS NT** (gene-level embeddings, aligned with the presence signal) — but that needs Bakta (~40 hr/250 genomes) and the classical presence/absence already gives 0.729 cheaply, so per-CDS NT is low-ROI unless there's a specific reason to expect gene-level FM embeddings to beat simple gene-presence.
- Artifact: `research_outputs/pathotype_nt_windowed_arm_2026-06-02.json`. Per-genome embeddings cached at `data/nt_windows/*.npy` (pooling variants are now free — no GPU re-pay).

### Decision (pre-committed branch 3): pooled NT < 0.729 → PIVOT to gene-targeted virulence detection
The pragmatic tool baseline is now **classical full-spectrum k-mer presence/absence (0.729)**. NT pooled adds nothing. Before trusting 0.729 (or anything) as *biology*, resolve the Salipante-vs-Hazen study confound with lineage/ST-aware splits.

<details><summary>(historical) in-flight notes for the run that produced the above</summary>

## windowed-NT arm config (task `bljrvirar`, 2026-06-02)
- Script: `scripts/pathotype_nt_windowed_arm.py`. NT v2 100M, whole-genome windowed (NO Bakta). **Run command:** `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True .venv/Scripts/python.exe scripts/pathotype_nt_windowed_arm.py` (from dna_decode/).
- **Final TDR-safe config: WINDOW=6144, STRIDE=6144, GPU_BATCH=2, MAX_N_FRAC=0.01** (skip >1%-N windows). ~748 clean windows/genome.
- **Two GTX 860M constraints learned the hard way (both will recur for any FM-embedding work here):**
  1. **4 GB VRAM OOM** — attention is O(B·H·T²). 12288 bp (~2048 tok) @ batch 8 OOMs; even 6144 bp @ batch 8 OOMs on N-inflated windows. Fix: skip >1%-N windows (the 100-N contig separators break 6-mer tokenization into long, memory-spiky sequences) so survivors are clean intra-contig ACGT (~1024 tok, uniform).
  2. **Windows TDR** — the 860M is ALSO the display GPU; any CUDA kernel >2 s triggers Windows Timeout Detection & Recovery (driver reset → `CUDA error: the launch timed out`). Batch 8 kernels exceed 2 s. **Only batch ≤2 is safe.** `expandable_segments:True` is a no-op on Windows (warns + ignores).
- **~1 window/s → ~12 min/genome → ~5-8 hr for 24.** Display GPU saturated throughout (user accepted 2026-06-02).
- **Per-genome embeddings cache to `data/nt_windows/<acc>_w6144_s6144.npy`** (C:, gitignored) as each genome finishes → RESTARTABLE. If it dies (session boundary / TDR), just re-run the same command; cached genomes skip the GPU. **First action in a fresh window: check `ls data/nt_windows/` for how many of 24 are done, then re-run to resume.**
- Pools mean + max across windows → LOSO XGBoost (same `train_xgboost_classifier`, calibrate=False, as the k-mer arm → fair).

</details>

## Next move when the NT number lands
1. **max-pool ≥0.76** → FM premise holds on this contrast → next: resolve the study-confound with lineage/ST-aware splits before trusting it as biology; then consider per-CDS NT.
2. **0.729 ≤ best < 0.76** → FM marginal, not worth the complexity → presence/absence k-mer is the pragmatic tool baseline.
3. **best < 0.729** → pooled NT dilutes localized signal (tet precedent confirmed) → pivot to gene-targeted virulence detection (the v0 deterministic-rules direction).
4. If max ≫ mean (as predicted), record that mean-pool is the wrong aggregation for pathotype — informs all downstream FM work.

## Uncommitted at this handoff (ahead of origin/main; commit when user approves)
- `scripts/pathotype_kmer_bakeoff.py` (cached LOSO + load_strains refactor), `scripts/pathotype_kmer_fullspectrum_diag.py` (new), `scripts/pathotype_nt_windowed_arm.py` (new), the two `research_outputs/*.json` results, this handoff.
- Stray empty-ish `2026-05-30-1200-ep4-pathotype/` Soraya-run dir at repo root still present (non-empty: approvals/audit/intent/recommendation/result .md stubs).
