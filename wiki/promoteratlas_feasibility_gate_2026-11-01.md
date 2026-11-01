# PromoterAtlas base-model feasibility gate — PASSES, and it corrected four of my own claims

<!-- doc-status: LIVE as of 2026-11-01. This is the FEASIBILITY GATE for the ceiling arm
     `promoteratlas-base-embeddings-frozen-150bp-padded-to-200`, run BEFORE any acceptance bar was
     frozen. It is NOT the arm and carries no arm result. Companion to
     wiki/promoteratlas_identity_resolved_2026-10-31.md (identity + circularity). -->

## Why this ran before the bar

This project has repeated one lesson four times in one week: a threshold specified before asking
what the arithmetic and the substrate permit is a threshold that gets amended mid-run. The
2026-10-31 self-distillation probe named four such thresholds in a single session — a control that
could not collapse, a drift tolerance below the float32 floor, a beat margin a control could clear,
and a memory layout that got the job killed. So the gate ran first, and it earned its place: it
**killed or corrected four claims**, three of them mine.

## CORRECTIONS

**1. A path in the committed identity memo is wrong.** I wrote the weights live at
`trained_model_weights/` with flat `.pt` files, taking the README's table at face value. The actual
layout is `trained_weights/` with three subdirectories:

```
trained_weights/base_model/promoteratlas-base.pt                        22,313,814 B
trained_weights/expression_prediction/promoteratlas-trspred-urtecho2018.pt  11,310,034 B
trained_weights/expression_prediction/promoteratlas-trspred-lafleur2022.pt  11,310,034 B
trained_weights/expression_prediction/promoteratlas-trspred-hossain2020.pt  11,310,034 B
trained_weights/expression_prediction/promoteratlas-trspred-yu2021.pt       11,309,354 B
trained_weights/expression_prediction/promoteratlas-tslpred-kosuri2013.pt   11,309,690 B
trained_weights/segmentation/promoteratlas-annotation.pt                 8,813,208 B
```

The fetch had explicitly warned me it *"can't confirm the exact folder name or contents"* and that
the live listing showed `trained_weights`. I wrote the README's version as fact anyway. **No LFS
pointers** — every file is real and usable.

**2. The withdrawal of the "worse lead" is now grounded in the real directory, and it holds.** The
head set is exactly the five the README named — `urtecho2018`, `lafleur2022`, `hossain2020`,
`yu2021`, `kosuri2013` (translation) — plus `annotation`. **No genome-wide / tiling Urtecho head
exists.** Previously this rested on a README table; it now rests on the file list.

**3. A new inference of mine died on measurement, before it reached a memo.** `NC_000913.1`
(*E. coli* K12) ships as one of **13 benchmark genomes**, each with `_PA_precision` / `_PA_recall`
annotated variants, and I took that as strengthening the "K-12 was held out of pretraining" reading.
**It does not.** Resolving the 12 benchmark BioSamples to assembly accessions and testing membership
in the pretraining corpus gives **4 PRESENT / 8 absent** — `P. putida` KT2440, `A. hydrophila`,
`Y. pestis` A1122 and `A. baylyi` ADP1 are all in the corpus *and* benchmarked. So benchmark
membership says nothing about held-out status, and the inference is **withdrawn**. (K-12's own
assembly did not resolve: `NC_000913.1` is an old RefSeq record carrying `DBLINK Project: 225` and
no BioSample.)

**4. My own probe used the wrong mask count, and then my verdict function was too crude.** I masked
20 positions per sequence, the argparse default; the checkpoint's recovered config says pretraining
used **60** (`wm0_pm60`). Worse, the length ranking *inverted* between the two settings — see below.

## What the gate measured

### A. The architecture inferred from source loads the published weights EXACTLY

`DNATransformer(dim=128, n_blocks=8, n_heads=8)`, strict `load_state_dict`: **0 missing, 0
unexpected**, **1,847,428 parameters** — matching the paper's "1.8 million parameter transformer"
independently. Two routes agreeing on the same number is the validation; a strict load alone could
pass on a mis-specified architecture only by coincidence.

### B. The base model is length-agnostic — by construction AND by measurement

`self.embedding = nn.Conv1d(4, dim, kernel_size=1)`, rotary positional embeddings computed at
forward time, no learned positional table, and no flattening anywhere in the base. Latent shape is
exactly `(batch, 128, L)` at L = 86 / 120 / 150 / 200 / 250 / 300 / 400.

**The `seq_length=86` default I first flagged as a blocker constrains only `ExpressionPredictor`,**
because that class flattens (`fc_input_dim = base_dim * seq_length`). The arm uses the base, so it
is not constrained. My initial reading of 86 as "the model's window" was wrong.

### C. In-distribution length — the confound, and what survives it

First attempt held mask **count** fixed, which lets density vary 70% → 20% across lengths. The
ranking then inverted with the count:

| mask count | 86 | 150 | 200 | 300 | argmax |
|---|---|---|---|---|---|
| 20 (argparse default) | 0.3271 | **0.4208** | 0.4188 | 0.4042 | 150 |
| 60 (true pretraining) | 0.2785 | 0.3514 | 0.3882 | **0.4299** | 300 |

**A ranking that inverts when a nuisance parameter changes is measuring the nuisance parameter.** So
density was pinned instead (random baseline 0.25):

| density | 86 | 120 | 150 | 200 | 250 | 300 | 400 | argmax |
|---|---|---|---|---|---|---|---|---|
| 10% | 0.3576 | 0.3776 | 0.4104 | 0.4391 | **0.4562** | 0.4146 | 0.4117 | 250 |
| 20% | 0.3290 | 0.3477 | 0.3708 | 0.3906 | 0.4031 | **0.4302** | 0.4074 | 300 |
| 30% | 0.2933 | 0.3472 | 0.3597 | 0.3708 | 0.3871 | **0.4076** | 0.3841 | 300 |

**The argmax is unstable (250 / 300 / 300); the SHAPE is stable at every density** — monotone rise
from 86 to a shallow plateau at ~200–300, then decline at 400. **86 bp is the worst length at all
three densities**, by 0.03–0.12.

So the arm's declared `150bp padded to 200` sits **inside the favourable plateau**, and the 86 the
expression heads use would have been the worst available choice. The automated verdict string in
`pa_density_probe.py` reads `UNSTABLE ... treat length as UNRESOLVED`; that is right about the
argmax and **too pessimistic about the implication**, because it tested argmax stability and spread
when the stable signal was the curve's form. The 0.10 "spread" it flagged is almost entirely the
86-vs-peak gap.

**Honest limit, and it is not small:** reconstruction accuracy measures the *pretraining objective*.
A length at which masked-nucleotide prediction does well is not guaranteed to be the length at which
frozen embeddings are most useful to a ridge. This validates the geometry as not-pathological; it
does not establish it as optimal for the arm's actual task.

### D. Non-vacuity control — passes everywhere

Real sequence beats shuffled at every length and every density (real − shuffled = +0.031 to +0.202),
with shuffled pinned near 0.25–0.27 as it should be. Without this, a flat accuracy curve would be
indistinguishable from a model doing nothing. One mild outlier: shuffled 0.3125 at 200bp/10%, on 640
scored positions — noise, not signal.

### E. Deterministic

Max `|latent delta|` over a repeated forward pass: **0.000e+00**, bit-identical. Worth stating
explicitly because the conv encoder in this project has a recorded unreproduced jitter anomaly; this
model shows none.

## The pretraining config, recovered from the checkpoint

```
epoch 238 (of n_epochs 2000, patience 25 -> early-stopped)    val_loss 0.29824
n_window_masks 0    window_size 6    n_point_masks 60    batch_size 1024    lr 1e-3
data_path data/processed/sequence_dataset.h5    architecture DNATransformer
```

**This independently re-confirms the circularity finding by a second route.** The objective is
masked-nucleotide prediction over genomic sequence; there is no expression term, no label, no MPRA
input. The identity memo established that from the corpus file list and the dataset card's own
provenance sentence; the checkpoint's own training config now says the same thing from inside the
artifact.

**No sequence length is recorded in the config** — it comes from the h5 (`batch, channels,
seq_length = sequence.size()`), which is 2.3 GB and was not downloaded. So the true pretraining
window remains unmeasured, which is exactly why section C matters and why its limit is stated.

## A safety decision, made on evidence rather than on trust

`weights_only=True` refuses this checkpoint — it embeds a `pathlib.PosixPath`, which is both off the
allowlist and un-instantiable on Windows. Aliasing `PosixPath` does not help, because `safe_globals`
keys on the *name* `pathlib.PosixPath`, not the resolved object.

Rather than disable the guard on a downloaded file or abandon the probe, the pickle was
**disassembled without executing it** (`pickletools.genops` over `best_model/data.pkl` inside the
zip). It references exactly four globals:

```
collections.OrderedDict   pathlib.PosixPath   torch.FloatStorage   torch._utils._rebuild_tensor_v2
```

No callable that could execute arbitrary logic is named, so `weights_only=False` is safe **for this
file, on this evidence**. The audit script refuses and exits 3 if any global falls outside a
torch/pathlib/stdlib allowlist. Nothing was installed at any point: `dna_transformer.py` imports
only torch, so it was loaded as a standalone module file and the downloaded package's `__init__`
never ran.

**The isolation mattered for a measured reason:** PromoterAtlas pins `numpy<2` and this project runs
**numpy 2.4.4**, so `pip install -e .` would have downgraded numpy under the entire 6,700-test suite
— the same hazard that forced Kaptive into an isolated install.

## Gate verdict

| check | result |
|---|---|
| A. architecture loads exactly | **PASS** (0/0, 1,847,428 params) |
| B. length-agnostic | **PASS** (by construction and measurement) |
| C. 150→200 in a favourable length regime | **PASS, with the proxy limit stated** |
| D. non-vacuity control | **PASS** at every length × density |
| E. determinism | **PASS** (0.000e+00) |

**The arm is feasible and its declared geometry is validated.** What remains before freezing a bar
is the pooling decision (the latent is `(B, 128, L)` and a ridge needs a fixed vector, so masked
mean-pool vs flatten must be declared pre-run) and the sequence-exposure question — whether any
*E. coli* assembly sits in the 3,371-genome pretraining corpus, which is being resolved separately.

**Reproduce:** `pickle_audit.py` → `pa_arch_probe.py [n_masks]` → `pa_density_probe.py`.
Artifacts `pa_arch_probe.json`, `pa_density_probe.json`.

## Honest limits

- Every number here is on **E. coli K-12 genomic windows** parsed from the repo's own
  `NC_000913.1.gb`. Not the tile substrate, not promoter-enriched — generic genomic sequence.
- n = 24–32 windows per cell. Fine for the gross shape of the length curve, not for fine ranking
  inside the plateau; the 250-vs-300 argmax difference is within that noise.
- The length curve is a **proxy** for in-distribution-ness, and the true pretraining window is still
  unmeasured.
- A feasibility gate says the arm *can* run and its geometry is sane. It says nothing about whether
  the arm's number will be any good — the standing inferred merit-negative still argues it will not.
