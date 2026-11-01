# The PromoterAtlas arm closes ATLAS_BEATS_THERMODYNAMIC — best on the rung, and it refutes a standing claim

<!-- doc-status: LIVE as of 2026-11-01. Graded against wiki/promoteratlas_arm_acceptance_bar.json
     (schema v2, frozen BEFORE the run, all seven outcomes and the check order fixed in advance).
     Artifact wiki/promoteratlas_tile_arm_2026-11-01.json. Orthogonality
     wiki/promoteratlas_orthogonality_2026-11-01.json. Feasibility gate
     wiki/promoteratlas_feasibility_gate_2026-11-01.md. Identity + circularity
     wiki/promoteratlas_identity_resolved_2026-10-31.md. Ungated first run preserved at
     wiki/promoteratlas_tile_arm_run1_ungated_2026-11-01.json. -->

## CORRECTIONS FIRST

**1. A standing project claim is REFUTED, not confirmed.** The recorded inference against this arm
was that PromoterAtlas would be a third *sequence-learned* channel beside the conv already in the
stack, and therefore *"the next channel should carry a different PRIOR, not a different neural net."*
Measured on the incumbents' own held-out tiles: **rank-average(max_tx, atlas) = 0.5402 against the
committed max_tx+conv hybrid at 0.4254, +0.1148.** A different neural net carries something the conv
does not. The inference was explicitly labelled INFERRED and is now measured — against itself.

**2. My forecast was badly wrong, and this is the sixth score.** The frozen bar put **p=0.55 on
ATLAS_LOSES_ALONE** and **p=0.08 on ATLAS_BEATS_THERMODYNAMIC**. The run landed on the 0.08 branch.
The orthogonality sub-prediction failed differently: I predicted atlas-vs-conv HIGH (both
sequence-learned) and atlas-vs-max_tx LOWER (different prior), and they came back
**indistinguishable** (0.4089 vs 0.4076). **My forecast record in this project is now 2 of 6 — do
not defer to these priors.**

**3. My verdict code let the strongest outcome through while a gate was failing.** Run 1 emitted this
same verdict while the shuffled-label control sat at 0.0203 against a 0.02 bar, because the clause
read `abs(med(null) - med(gc)) < MARGIN` — i.e. |0.0203 − 0.3112| = 0.2909, never below 0.02. So
`INVALID_NULL_SCORED` was **unreachable** and the most consequential line in a validity-first order
was dead code. Fixed, proven non-vacuous, and run 1 is preserved with its verdict unedited under a
loud `verdict_ungated` block. **The bar's own wording was also wrong** — *"|median gain over GC| <
0.02"* would require the null to be as *good* as GC — and the corrected kernel applies the reading
its derivation plainly means (null near zero), which is the **stricter** one. The bar's text was not
re-worded; a frozen bar is not edited after a result.

## The result

| arm | median Spearman | note |
|---|---|---|
| **promoteratlas base, frozen, native 150** | **0.5076** | **new best on this rung** |
| rank-avg(max_tx, atlas) | 0.5402 | vs committed 2-channel hybrid 0.4254 |
| promoteratlas, padded 200 (declared geometry) | 0.4983 | padding costs 0.0093 |
| orthogonal-modality stack (previous best) | 0.4404 | 3 channels |
| promcalc + conv hybrid | 0.4254 | 2 channels |
| promoter-calculator `max_tx` (thermodynamic) | 0.3968 | previous best single |
| promcalc full decomposition | 0.3925 | |
| conv, 13,217 params | 0.3234 | the other sequence-learned channel |
| pyramid 4-scale 150bp | 0.3167 | |
| GC, single feature | 0.3112 | |

**Paired against conv on byte-identical test sets: +0.1956, `median_ci95` [0.1704, 0.2182] at
achieved coverage 0.9785, 10 of 10 seeds.**

**Verdict `ATLAS_BEATS_THERMODYNAMIC`** — 0.5076 ≥ 0.3973 + 0.02, the bar's strongest outcome and
the one it flagged in advance as refuting a standing claim.

## Every gate, and the one that failed first

| gate | result |
|---|---|
| split equality, 10 seeds | **10/10** byte-identical to the committed `tile_indices` (sha256 per seed); `n_test` exact |
| GC pipeline anchor | **0.3000 at seed 0** against 0.3000 ± 0.01, median 0.3112 against the ledger's 0.3112 |
| embedding non-degeneracy | 0 constant dims of 128; mean inter-tile cosine 0.6473 |
| determinism | **0.00e+00 on CUDA** |
| shuffled-label null | **−0.0014** over **270 pooled seeds**, ci95 [−0.0114, +0.0098], achieved 0.9556 → **PASS** |

**The null failed first, at 0.0203 on ten seeds, and the run was not graded on it.** The bar's
pre-commitment was *"if the null median lands outside the bar, the honest move is to POOL MORE
SEEDS, not to re-specify"*, and the arithmetic backed that rather than preference: per-seed nulls
spanned −0.0909 to +0.1293 with sd 0.0600, so the median's own standard error (~0.024) was **larger
than the bar it was judged against**. Pooling to 270 moved the estimate to −0.0014 with an interval
comfortably inside 0.02. The ten-seed figure was noise, exactly as the sibling fragment arm's was.

**The GC anchor is the load-bearing gate, not the null.** Reproducing a known number to four decimal
places validates the loader, the `category=='tile' and len==150` filter, `leave_peak_out`, the
`log1p` target, the ridge protocol and the scorer *together*. A pipeline that cannot reproduce a
known number cannot establish a new one.

**Prediction kept, stated before the corrected run:** the arm's numbers would reproduce because the
arm is deterministic and the splits fixed. They did, identically — 0.5076 / 0.4983 / +0.1841 across
two independent Kaggle runs.

## Orthogonality — and the arm is not "another conv"

| atlas vs | median Spearman |
|---|---|
| decomposition | 0.4227 |
| conv | 0.4089 |
| max_tx | 0.4076 |
| gc | 0.3862 |

All four inside **0.386–0.423**. The arm sits roughly equidistant from every existing channel rather
than clustering with the other sequence-learned one, which is neither branch of what the bar
predicted. That is why the rank-average with the thermodynamic channel gains so much.

## Honest limits — the first one is the one that matters

**1. CORPUS-SELECTION EXPOSURE IS THE LEADING ALTERNATIVE EXPLANATION AND IS NOT EXCLUDED.** The
pretraining corpus holds **66 Escherichia/Shigella assemblies** (BL21(DE3), B str. REL606, UTI89,
IAI39, O157:H7, four *Shigella* species), and its own dataset card describes the corpus as regulatory
sequences *"extracted from GenBank files of gammaproteobacterial genomes"* — plausibly
upstream-of-gene regions. If so, *"resembles what was selected into my pretraining corpus"* is itself
a proxy for *"is a real promoter"*, which correlates with measured expression.

This is **not label leakage**: pretraining is masked-nucleotide prediction with no expression term,
the encoder is frozen, and the split is verified identical to the incumbents'. But it is a concrete
route by which a selection signal could inflate this number, and separating it needs the corpus
sequences (a 2.3 GB h5, or the 3,371 FASTAs), which were not downloaded. **Until that is done, read
0.5076 as an upper bound on what this encoder contributes.**

**2. Sequence exposure is real and measured.** The substrate is K-12 MG1655; that exact assembly
(`GCF_000005845`) is absent from the corpus, but BL21 and REL606 share most of K-12's core genome.
**The arm may not be described as evaluating on unseen sequence.**

**3. A 128-dim mean-pooled embedding is the weakest reasonable use of this encoder** — mean-pooling
discards the −35/−10 spacing promoter biology turns on. A positional or fine-tuned variant is a NEW
ARM under R6, never a re-read of this one. Stated in the bar *before* the run, so it cannot be used
as a post-hoc escape in either direction.

**4. 0.5076 against the rung's replicate-derived ceiling of 0.941 is 53.9% of achievable.** Best on
the rung is not close to the ceiling.

**5. One rung, one organism, one substrate, 10 seeds.** Nothing here speaks to any other rung, and
this must never be read against the fragment rung's 0.7928 ceiling — the census `ceiling_rule`
forbids that comparison.

**6. Four Kaggle runs, three of them failures, and only one failure was the substrate's fault.** A
tile-load contract bug (the dry run stubbed that exact function with the wrong return shape — the
one boundary it faked was the one it got wrong), a dead validity gate, and one genuine transient
NCBI `ContentTooShortError` now covered by a bounded retry.

## Consequence for the ceiling criterion

`ceiling_verdict` **stays `blocked:untried-arm`** — 10 measured, **1 planned**, 1 killed. `best`
moves to this arm. The single remaining blocker is
`substrate-change-constructed-eukaryotic-or-human-panel-UNRUN`, and it is **not drive-blocked**:
measured 2026-11-01, its substrate is open CC-BY-4.0 on Zenodo record 10558183, four tables totalling
~74 MB, fetchable in-kernel exactly as this arm fetched GSE144621. So the criterion is closable
without the detached D: volume. That arm carries its own hazard, recorded already: Table S6 publishes
per-fold incumbents on those exact folds but **never states the metric**.
