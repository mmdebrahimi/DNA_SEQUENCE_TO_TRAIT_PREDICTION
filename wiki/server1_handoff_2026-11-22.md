# Server1 GPU handoff — activity-conditioned human generator, size curve + bar-mandated controls

**Date:** 2026-11-22 · **Target:** Server1, GTX 1070 (Pascal, sm_61, 8 GB), shared gaming box
**Bar:** `wiki/human_design_bar_2026-11-22.json` (FROZEN) · **Reservation:** `wiki/human_final_shard_reservation_2026-11-22.json`

---

## DECISIONS NEEDED

**None — and one blocker is dissolved rather than answered.** The request was held up confirming
whether the lentiMPRA deposit sits at `E:\dna_decode_cache\mpra_human`. That no longer matters:
Server1 receives a pre-exported training shard and never reads the deposit. Do not go looking for
it on that machine; if you find it, do not point this job at it.

One item is a user-authority call and is **not** being taken here: whether the `xl` point
(17.8 M params, ~5–8 h) is worth a slot at all. It is listed as Tier 3 and left unrun by default.

---

## CORRECTIONS to the plan as proposed

**1. Train on 200 bp, not 230.** The request specified 230 bp. The frozen bar's primary geometry is
`element_200nt_adaptors_stripped`. The 15 nt flanks are **measured constant** — one distinct value
each, share 1.0000, identical in both cell types, and an order-3 Markov chain reproduced them 0 of
120 times (`wiki/human_adaptor_probe_2026-11-21.json`). Training on 230 hands the model 30 constant
positions, which is the exact thing that inflated a distinguishability margin by up to **0.4296**.
The shard is already trimmed; the script hard-refuses anything that is not 200 bp.

**2. Do not copy the deposit to Server1 — and the reason is leakage, not convenience.** The deposit
carries folds 1–10. The generator trains on 1–6; 7–8 are the oracle's, 9 is working-held-out, and
**fold 10 is reserved FINAL_ONLY**. Copying the workbook puts the reserved final shard on a third
machine running a live session that will be edited later. The exported shard contains folds 1–6 and
nothing else, so the reserved shard **physically cannot be touched there**. Three further reasons
are in `scripts/glm_export_server1_shard.py`'s docstring: data vintage silently changes the
conditioning label, re-implementing a 200-line parser with two measured traps has no runner to gate
against, and Table S3 is 56 MB per openpyxl pass.

**3. Server1 must not score, and the script refuses to.** The judge is an oracle trained on folds
7–8 and the decision quantities are frozen. A box scoring its own output with its own re-derived
oracle would make every run a different generator seen through a different judge, and the two
machines' numbers would not be comparable. Server1 emits **samples + a manifest**; all scoring
happens on the repo host with the one frozen oracle. Do not add a metric to this script.

**4. A real bug in the shipped sampler, which affects the Kaggle kernel too.**
`dna_decode/glm/diffusion.py::sample` builds `torch.Generator(device="cpu")` and then calls
`torch.rand(x.shape, generator=g)` with no device, so on CUDA it raises a device mismatch. **The
shipped sampler is CPU-only by construction.** The standalone carries a device-correct fix. The
Kaggle kernel needs the same fix or it will fail at the sampling step after training completes —
the worst possible place to discover it. Consequence accepted: CUDA and CPU RNG streams differ, so
a seed reproduces on the **same device** only, which is all seed-spread needs.

**5. The curve is the right frame, but it was aimed at the wrong payload.** More curve points are
not what the second GPU is scarce for. Two quantities are **bar-mandated and Kaggle structurally
cannot supply them alone**:

- The bar deliberately freezes **no numeric tolerance**: *"The tolerance is the generated set's own
  seed spread, derived in the run."* A single Kaggle run produces one seed and therefore **no
  tolerance at all** — the secondary decision quantity (top-target lift vs the stratified null)
  is not evaluable without it.
- The bar's `NO_TARGET_CONTROL` reading says to *"check that the conditioning signal MOVES the
  output before concluding anything about capacity — an invariance from an unverified knob is
  plumbing, not a null result."* That needs a shuffled-label arm **at matched capacity**, which a
  size curve does not contain.

So Tier 1 below is seed replicates and the matched no-information arm at the headline capacity, and
the size curve is Tier 2. A size curve measured on an unverified knob is a curve of zeros.

---

## Install (free, no cost step)

Pascal needs a **cu126 or older** wheel; CUDA 12.8+ builds dropped sm_61. torch only — no numpy, no
project dependency.

```
python -m pip install --index-url https://download.pytorch.org/whl/cu126 "torch==2.7.1+cu126"
```

If that resolves badly, `torch==2.6.0+cu126` from the same index is the fallback. **Do not verify by
reading the version string** — the script calls `torch.cuda.get_arch_list()` and **refuses at
startup** if the wheel carries no kernel for sm_61, naming the cause. That refusal is the check;
discovering it mid-run is what it exists to prevent.

fp32 throughout, no AMP: sm_61 has no usable tensor cores and fp16 is not faster there.

## The package

Copy this directory to Server1 (12 MB):

`C:\Users\Laptop Portable\dna_decode_server1_handoff_2026-11-22\`

| file | what |
|---|---|
| `glm_human_diffusion_server1.py` | the job. stdlib + torch. imports nothing from any checkout |
| `train_K562_folds1-6.tsv` | 58,998 rows, `bin<TAB>sequence`, 200 bp ACGT. md5 `94f2eef17a11900cd0d6fda901381860` |
| `shard_manifest_K562.json` | provenance: source md5s, per-bin counts, activity ranges, the binning rule |

**Do not point the job at the stale `dna_decode` checkout on that box** (128 commits behind, with
another session's uncommitted edit). Nothing in the package reads it.

## Getting the package onto Server1 - the shard is NOT in git

**Verified, not assumed:** commit `31d68d2f` contains exactly three files -- the export script, the
trainer and this memo. `git ls-files` matches `train_K562_folds1-6.tsv` nowhere in the repo.

**It is deliberately not committed, and the repo's own convention is the reason.** No multi-MB file
is tracked under `wiki/` anywhere in this repository; every generated corpus
(`glm_training_corpus_*.txt` and siblings) is left untracked. Committing a 12 MB derived table would
break that convention rather than follow it. So a fresh clone of `31d68d2f` gives Server1 the
scripts and this memo and **no data** -- the clone is necessary and not sufficient.

**It also cannot be solved by letting Server1 build the shard itself.** Running the export there
needs the full deposit on that box, which puts the reserved FINAL_ONLY fold back on a third machine
-- the one control this whole shard design exists to enforce. Do not take that shortcut because the
copy is awkward.

**Push from the repo host is not available.** Measured 2026-11-22 against 192.168.2.197:
port 22 closed, 445 closed, 139 closed, **11434 open** (Ollama only), ICMP no reply. No scp, no SMB.

**So Server1 pulls.** Both boxes are on the same `/24` (repo host is `192.168.2.25`). Serve the
package from the repo host:

```
python -m http.server 8765 --bind 192.168.2.25 --directory "C:\Users\Laptop Portable\dna_decode_server1_handoff_2026-11-22"
```

Then on Server1:

```
mkdir C:\s1
cd C:\s1
curl -o train_K562_folds1-6.tsv       http://192.168.2.25:8765/train_K562_folds1-6.tsv
curl -o shard_manifest_K562.json      http://192.168.2.25:8765/shard_manifest_K562.json
curl -o glm_human_diffusion_server1.py http://192.168.2.25:8765/glm_human_diffusion_server1.py
certutil -hashfile train_K562_folds1-6.tsv MD5
```

`certutil` must print `94f2eef17a11900cd0d6fda901381860`. Passing `--shard-md5` to the job re-checks
it anyway and refuses on mismatch, so a truncated transfer cannot reach training.

**Limitation, stated plainly:** that HTTP server lives only as long as the repo-host session serving
it. It is a transfer window, not infrastructure. If it is down, the fallback is an ordinary manual
copy of the same 12 MB directory by any means the user prefers -- nothing about the job depends on
how the bytes arrive, only on the md5 matching.

## Run command

One line, with the data path explicit — there is no default data path, because a defaulted one
reports absent for data that is present.

```
python glm_human_diffusion_server1.py ^
  --shard      C:\s1\train_K562_folds1-6.tsv ^
  --shard-md5  94f2eef17a11900cd0d6fda901381860 ^
  --out-dir    C:\s1\runs ^
  --size l --seed 1 --steps 20000 --batch 64 --device cuda
```

Substitute the real paths. `--shard-md5` is optional but **pass it**: it is what catches a wrong
shard vintage, which would silently change the conditioning label and make these runs incomparable
to the headline.

### Yielding to a game — the whole contract

- `touch <out-dir>\PAUSE` → the job checkpoints at the next step boundary, **exits 0, and releases
  the GPU completely**. Delete the file and relaunch the identical command; `--resume auto` picks up
  from the checkpoint. Full release is why `exit` is the default rather than an in-process wait.
- Ctrl-C / SIGTERM → same: checkpoint, exit 0. Closing the window loses at most `--ckpt-every`
  (default 250) steps.
- `--on-pause wait` offloads to CPU and polls every 15 s, if an unattended hold is preferred.
- `--gpu-sleep-ms 20` throttles to leave the GPU responsive for light use.
- Checkpoints are written `.tmp` then `os.replace`d — atomic, because a torn checkpoint is worse
  than no checkpoint.
- Resuming refuses if the checkpoint was written under a different config or a different shard, so
  one experiment cannot be resumed into another.

## The curve, and who runs what

Real parameter counts, measured not estimated. Every point runs the **same `--steps`** so the curve
isolates capacity — varying capacity and token budget together cannot attribute a difference to
either.

| size | d_model / layers / heads / d_ff | params | tok/param | est. 20k steps on a 1070 |
|---|---|---|---|---|
| `s` | 128 / 4 / 4 / 512 | 820,740 | 14.38 | ~15–25 min |
| `m` | 192 / 6 / 6 / 768 | 2,710,660 | 4.35 | ~45–70 min |
| `l` | 256 / 8 / 8 / 1024 | **6,373,380** | **1.85** | ~1.7–2.7 h |
| `xl` | 384 / 10 / 8 / 1536 | 17,827,588 | 0.66 | ~5–8 h |

`l` sits at 1.85 tok/param against the borrowed 1.74 ratio, so **`l` is the sized-to-data point**
and the curve brackets it 7.8× below and 2.8× above. That is the shape that lets the sweep falsify
the borrowed ratio instead of assuming it.

**Runtimes are FLOP-based extrapolations, not measured on a 1070** — I have no 1070 to measure on.
The script prints measured `step/s`, peak VRAM and an ETA every 100 steps, so the real figure
arrives within the first minute and supersedes this table. VRAM at batch 64 should stay near
~0.3 GB (`s`) to ~2 GB (`xl`), well inside 8 GB; the printed peak is the authority.

### Kaggle keeps

`l` seed 0, conditioned — the headline run, already queued on fresh slug `glm-human-design-v1`.
It gets the large uninterruptible point because an interruption there costs the most.

### Server1 runs, in tier order

Tiers are ordered so that if the box is reclaimed for gaming, the bar-critical work is already done.
**Run Tier 1 first even though Tier 2 is cheaper.**

**Tier 1 — bar-mandated, ~5–8 h.** Without these the headline is not readable.
```
--size l --seed 1                    # with Kaggle's seed 0, starts the 3-seed spread
--size l --seed 2                    # completes it -> this IS the bar's tolerance
--size l --seed 0 --shuffle-labels   # matched-capacity no-information arm -> NO_TARGET_CONTROL
```

**Tier 2 — the size curve, ~2–3 h.** Two seeds each, so a curve point has an error bar.
```
--size s --seed 0      --size s --seed 1
--size m --seed 0      --size m --seed 1
```

**Tier 3 — only if the box stays idle.** `--size l --seed 1 --shuffle-labels` (a second control
point), then `--size xl --seed 0` if the user wants the top of the curve.

`--shuffle-labels` is **not a debug flag**. It permutes the bin column across rows, preserving the
label marginal while destroying the sequence–label association, so a model trained on it must show
~zero target monotonicity. Measuring against it is *differs from random*, not *differs from zero* —
the distinction that caught a near-harmless-looking prompt leak on the E. coli arm.

## What comes back

Per run, from `--out-dir`:

- `samples_K562_<size>_seed<k>[_shuffled].tsv` — 1,200 sequences, `target_bin<TAB>sequence`, 240 per
  bin × 5 bins. ~250 KB.
- `run_manifest_K562_<size>_seed<k>[_shuffled].json` — config, real param count, tokens seen, loss
  history with measured step/s, shard md5, environment, distinct-sample fraction, honest limits.

**Send back the `.tsv` and `.json` only.** Checkpoints stay on Server1 — they exist for resume and
are not needed for scoring.

## Verified before handoff

Each guard was made to fire, not merely written:

| guard | verified behaviour |
|---|---|
| wrong `--shard-md5` | REFUSE, exit 1, names the conditioning-label consequence |
| `--resume auto` | resumed step 20 → 30 and reported it |
| checkpoint config mismatch | REFUSE rather than resume one experiment into another |
| `PAUSE` file | checkpointed at step 0, exit 0, GPU released, named the file to delete |
| `--shuffle-labels` | ran, announced the arm, tagged its output `_shuffled` |
| end-to-end | trained, sampled 15/15 distinct, wrote samples + manifest |

Also cross-checked: the shard is 58,998 × 200 = **11,799,600 bp = 11.80 Mbp**, matching the Kaggle
sizing figure independently; bins are equal-count at 11799 / 11800 / 11799 / 11800 / 11800; zero
rows dropped for a non-ACGT base; `s` measures 820,740 params against the 813 K figure the Kaggle
plan quoted.

## Honest limits

- Runtimes are extrapolations. The box's own printed `step/s` supersedes them.
- `l` was chosen as the headline because it is the sized-to-data point under the **borrowed** 1.74
  ratio. The curve exists precisely because that ratio is unverified here; if `m` matches `l` on the
  decision quantities, the headline should move down and the Kaggle run was oversized.
- Eval loss on Server1 is a held-out slice of folds 1–6, not fold 9 — an in-corpus loss, **not** a
  design number. Neither decision quantity is computed on that box.
- Three seeds is the minimum that yields a spread at all. It is a thin tolerance and should be
  reported as such.
- The conditioning bins travel as data, so Server1 cannot disagree with the host about a label. It
  can still disagree about everything downstream of sampling, which is why it samples and does not
  score.
- `distinct_fraction` in the manifest is a cheap degeneracy tripwire, not the copying control. The
  copying check is longest-shared-substring against the training folds and runs on the host.
