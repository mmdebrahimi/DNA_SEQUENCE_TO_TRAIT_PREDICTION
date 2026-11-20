- `scripts/kaggle/build_glm_finetune_kernel.py` — the BUILDER. Emits a kernel by substituting
  `__DIAG__` / `__PAYLOAD__` into a `KERNEL` template string, then **AST-parses the emitted body and
  refuses to write one that does not parse**. Its own comment records why: on 2026-10-31 the builder
  parsed, the kernel did not, and Kaggle reported ERROR with an EMPTY failureMessage and zero output
  files — "indistinguishable from an infrastructure failure". Current args: only `--diagnostic`.
  **Its file-writing entry point's default output path IS `scripts/kaggle/glm_finetune_kernel.py`,**
  so any build-then-compare test would overwrite the artifact it is checking.
- `run_arm(cli, label, model_path, device, outfile)` at line 362 — scores ONE arm through the
  falsifier CLI; "only `--model` differs between arms". Passes `--modes kmer,both`,
  `--markov-order 3`, `--context`, `--n`, `--length`, `--limit-natural`, `--bins`, `--cache-dir`,
  `--out`. It does NOT pass `--seed`.
- `step6_finetune(peft_model, tok, seqs)` — a HAND-ROLLED PyTorch loop (`for ep in range(EPOCHS)`),
  NOT a HuggingFace Trainer, so `save_strategy="epoch"` does not apply and per-epoch adapter saving
  is an explicit call inside the existing loop.
- `verdict(base, ft, null_ok)` — labelled "The FROZEN reading. Ordered; first match wins." MUST NOT
  be edited; this diagnostic needs its own separate reading function.
- Frozen method constants: `EPOCHS = 3`, `LR = 1e-4`, `BATCH = 8`, `DTYPE = "float32"`,
  `LORA = {r:16, lora_alpha:32, lora_dropout:0.05, target_modules:[q,k,v,o]_proj}`, under
  `# ---- the frozen LoRA method ----` and mirrored in the pre-registration's `method_frozen`.
- `ADAPTER = OUT + "/lora_adapter"` where `OUT = /kaggle/working`, written at line 473 — so the
  adapter IS collected as kernel output. `MERGED = /tmp/glm_finetuned_merged` at line 475, so the
  4.8 GB merge is deliberately NOT collected.
- `dna_decode/glm/falsifier.py` — `heldout_auroc(a, b, *, k=4, seed=0, mode, bins)` (seeded, so the
  discriminator half is deterministic given the sequences) and
  `self_control_ceiling(natural, seeds=(0,1,2,3,4))` returning `(max_observed, all_values)`; the
  MAXIMUM is the ceiling "so a lucky split cannot manufacture a finding".
- `scripts/kaggle_push_poll.py` — `push(code_py, slug, gpu, internet, dataset=None)`; imports
  `from kaggle.api.kaggle_api_extended import KaggleApi`, so the SDK must be importable INSIDE the
  project venv. Verified live: SDK installed, `api.authenticate()` returns OK.
- Documented sampling trap, `generate.py:181-188`: passing `temperature=` to `generate()` is a NO-OP
  on this path (measured T=0.7 vs T=1.5, same seed, 6/6 byte-identical). Explicit
  `TemperatureLogitsWarper` / `TopKLogitsWarper` are the working route. Both committed arms ran at
  recorded defaults with `knob_wired: false`, so this plan changes NO sampling knob.

### Reusable-Code Survey
- `run_arm` — reused verbatim for every scoring; the only change is threading an explicit `--seed`.
  Nothing new is written to score a model.
- `falsifier.self_control_ceiling` — reused as the already-derived distinguishability floor. This
  plan ADDS a paired generation-level delta beside it rather than replacing it, because the two
  answer different questions.
- `step6_finetune`'s existing epoch loop — extended with one `save_pretrained` call rather than
  replaced with a Trainer.
- `scripts/kaggle_push_poll.py` — reused for push/status/pull; no new Kaggle plumbing.
- `tests/test_glm_finetune_kernel_builder.py` — existing builder tests, extended rather than
  duplicated (siblings: `test_finetune_kernel_builder.py`,
  `test_glm_sweep_kernel_attests_the_knob.py`).
- Searched: `scripts/kaggle/`, `dna_decode/glm/`, `tests/`, `scripts/kaggle_push_poll.py`; no
  `graphify-out/GRAPH_REPORT.md` on this host.
