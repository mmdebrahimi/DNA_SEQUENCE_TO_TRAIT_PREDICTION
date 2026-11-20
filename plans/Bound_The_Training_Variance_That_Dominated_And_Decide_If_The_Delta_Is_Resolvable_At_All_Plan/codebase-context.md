- `scripts/kaggle/build_glm_finetune_kernel.py` — the builder. Since the predecessor plan it exposes
  **`build_body(*, diagnostic=False, paired=False, payload="x")`** as the SINGLE body constructor,
  used by `main()` and by all three body-building tests, plus a guard that **REFUSES any body still
  carrying an unsubstituted `__PLACEHOLDER__`** (Python dunders excluded). A new placeholder must be
  added to `build_body` or the guard fires — that guard exists because adding `__PAIRED__` and
  `__SEEDS__` broke an exec-based test while two parse-only tests silently kept passing.
- `SEED_SET = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)` is defined once at builder level and injected via
  `__SEEDS__`. `run_arm(cli, label, model_path, device, outfile, *, seed=None)` appends `--seed`
  only inside `if seed is not None:` (0 unconditional occurrences).
- `step4_lora_shapecheck` at line 352 does `from peft import LoraConfig, get_peft_model` — **the
  insertion point for training-seed labelling**, because the adapter's random `lora_A` init is where
  the run-to-run divergence originates.
- `step6_finetune(peft_model, tok, seqs)` — hand-rolled PyTorch loop, saves a per-epoch adapter to
  `ADAPTER + "/ep%d"`. `ADAPTER` is under `/kaggle/working` and IS collected; `MERGED` is `/tmp` and
  is not.
- `paired_sweep(cli, device, dev, seqs, tok, model)` and `score_model_at_every_seed(cli, label,
  path, device, grid, epoch)` — the latter writes the grid after EVERY cell, which is the only
  reason the predecessor's spoiled verdict cost nothing.
- Line 447: `score_model_at_every_seed(cli, "zeroshot", MODEL_ID, device, grid, None)` scores the
  **pristine** `MODEL_ID`, so the zero-shot row is **training-independent** and the predecessor's 10
  committed zero-shot cells are reusable.
- `paired_reading` / `anchor_check` — the variant's own scale-free reading; the frozen `verdict()`
  is never called in that path.
- `scripts/glm_paired_sweep_verdict.py` — the offline recompute, with the corrected null gate:
  reproduction is anchored at **seed 7 only**, because the falsifier's single `--seed` also drives
  `build_natural` (which at line 57 does `random.Random(seed).shuffle(wins)` and so selects the
  fit/eval split) and `MarkovGenerator(seed=...).fit`, making the null refit per seed — it spanned
  **0.5674 to 0.7024** in the committed grid.
- Measured costs from two real runs: **138 s per scoring, 578 s per 3-epoch training**, whole
  predecessor run 40 scorings + 1 training in ~105 min.
- `scripts/kaggle_push_poll.py` — `push(code_py, slug, gpu, internet, dataset=None)`; imports the
  Kaggle SDK as a LIBRARY, so it must be importable in the project venv. The SDK is declared in
  `[dev]`, so `uv sync --extra forward` deletes it; restore with `uv pip install kaggle`.

### Reusable-Code Survey
- `build_body` + its placeholder guard — reused for the new variant; no second construction path.
- `run_arm(..., seed=...)` — reused verbatim; no new scoring code.
- `score_model_at_every_seed` — reused for the per-training-seed rows, keeping the per-cell grid write.
- `step6_finetune` + the per-epoch adapter saves — reused; this plan adds labelling around them, not
  a new training loop.
- The predecessor's 10 committed zero-shot cells in `wiki/glm_generator_paired_sweep_2026-11-20.json`
  — reused as data, since the row is training-independent.
- `scripts/glm_paired_sweep_verdict.py`'s `load_grid` / `cell` / `null_for` helpers — reused by the
  decomposition script rather than reimplemented.
- Searched: `scripts/kaggle/`, `scripts/glm_paired_sweep_verdict.py`, `dna_decode/glm/`, `tests/`;
  no `graphify-out/GRAPH_REPORT.md` on this host.
