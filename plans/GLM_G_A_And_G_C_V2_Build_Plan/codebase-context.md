### Reusable-Code Survey

- **`dna_decode/glm/genomewide.py`** — `load_peak_tiles`, `leave_peak_out(tiles, *, frac, seed)`,
  `position_blocked_split(frags, *, n_blocks, held_blocks, seed)`, `noise_ceiling`, `condition_effect(a, b,
  *, coordinate_key=True)`, `kmer_features`, `gc_features`, `positional_onehot`. **Both splits already take a
  `seed`**, so multi-seed is a loop rather than new machinery. `noise_ceiling` is duck-typed on `.rep1/.rep2`
  and already serves both `Fragment` and `Tile`.
- **`dna_decode/glm/expression.py`** — `load_pairs`, `leave_element_out`, `fit_predict_ridge(..., feature,
  alpha, log_target, return_pred)`. The grid arm reuses this unchanged.
- **`dna_decode/glm/falsifier.py`** — `heldout_auroc`, `self_control_ceiling`, `build_features`. Reused for
  G-A's optional active/inactive diagnostic and already reused by the reference-geometry triage.
- **`scripts/glm_genomewide_oracle.py`** — `_fit` is the arm-return shape every new arm mirrors.
- **`dna_decode/glm/generate.py`** — the **lazy-import convention to copy**: torch is imported only inside
  functions (lines 48 / 62 / 135 / 174 / 219), which is why the existing suite runs on a default install.
- **`tests/conftest.py`** — one `fixtures_dir` fixture; no model fixtures to reuse.
- **DUPLICATION FOUND:** there are **two independent k-mer implementations** —
  `falsifier.kmer_features` and `genomewide.kmer_features` — plus `expression.kmer_feature` which delegates to
  the falsifier one. Verified they currently return byte-identical output (width 256 on a test sequence), but
  that agreement is **coincidental, not enforced**. G-A's gate and the geometry triage would silently use
  different features if either drifts. Step 2 adds the guard.
- **None — searched:** no existing conv/CNN/encoder module anywhere under `dna_decode/` (searched for
  `Conv1d`, `torch.nn`, `Sequential`, `encoder`); no `graphify-out/GRAPH_REPORT.md`; no `src/utils`-style
  shared-helper directories. The encoder is genuinely new code.

### External surfaces — verified live, not inherited

- **torch 2.7.1+cu118**, sklearn 1.8.0 (`uv run python -c "import torch"`).
- **`torch.nn.Conv1d` expects `(N, C, L)`** — verified: input `(2, 4, 150)` with `kernel_size=6` returns
  `(2, 8, 145)`, i.e. L shrinks by `k-1` with no padding. Any multi-width branch therefore produces
  *different* output lengths and must pool before concatenation.
- **`torch.manual_seed` and `torch.use_deterministic_algorithms` both exist** on this version, so determinism
  is pinnable.
- **Local GPU is usable and this was measured, not assumed:** `torch.cuda.is_available()` is True on a GTX
  860M, compute capability **(5, 0)**, and a `Conv1d(4,32,6)` forward on a `(256,4,150)` batch succeeds
  (9.4 s first call, which is CUDA context init). No tensor cores, no bitsandbytes — irrelevant for fp32
  convolutions this small.
- **Kaggle escalation path exists:** `~/.kaggle/access_token` is present and `scripts/kaggle/` already holds
  working kernels (`pear_prosst_kernel.py`, `stage2_enformer_kernel.py`, …). The `kaggle` CLI is **not** on
  PATH, so the push path must follow the existing kernel pattern rather than assume a CLI.
- **torch is NOT in default deps** — present only in the `forward` and `ml` extras (verified by parsing
  `pyproject.toml`), which is why Step 6 lazy-imports.
- pytest markers `slow` and `integration` already exist for the heavier arms.
