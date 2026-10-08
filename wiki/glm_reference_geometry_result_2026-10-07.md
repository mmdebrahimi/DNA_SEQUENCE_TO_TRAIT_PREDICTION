# The designed grid is PERFECTLY separable from real genomic sequence — and "genomic" is not one distribution either

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_reference_geometry_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_reference_geometry.py --seeds 10 --n 2500` (CPU, ~12 min)
**Verdict:** `DISTINCT` · validity gate PASSED on all three sets

This is the triage that an adversarial review of the GLM decomposition said should run **first**, ahead of
both build families, because it tests a premise both of them rest on. It did, and it changed two things.

## The result

Held-out discriminator AUROC between reference sets, all exactly 150 bp so length cannot be the signal,
both sides subsampled to the same n so prevalence cannot be, 10 seeds, median reported.

| pair | k-mer | k-mer + positional | sd |
|---|---|---|---|
| **grid vs tile** | **1.0000** | **1.0000** | 0.0000 |
| **grid vs natural** | **1.0000** | **1.0000** | 0.0000 |
| tile vs natural | 0.6664 | **0.8338** | ~0.01 |

**Validity gate (decision D1, and it is load-bearing):** every set was first tested against *itself* —
split in half and try to discriminate the halves. All three came out at chance: grid 0.497, tile 0.501,
natural 0.494 (medians; max across seeds 0.516). A set separable from its own halves would mean the
discriminator fabricates signal and nothing could be graded. It does not, so the pairwise numbers stand.

## Finding 1 — the premise HOLDS, decisively, and the grid number inherits a caveat

G-A rests on "fixed features work on the designed grid and fail on real genomic sequence". If the two
substrates had turned out barely distinguishable, that distribution-shift story would have been wrong and
G-A would have needed re-reading before implementation. **They are perfectly distinguishable — AUROC 1.0000
with zero variance across 10 seeds, under both feature views.** The framing survives.

But a perfect separation is not a neutral confirmation. It means **the oracle's 0.5914 was measured on
sequence that a 4-mer histogram can identify as artificial with certainty.** The grid number was therefore
never comparable to a genomic number — which was already suspected and is now measured.

## Finding 2 — the mechanism is combinatorial recombination, NOT a construct scaffold

The obvious explanation for a perfect AUROC is a constant scaffold — shared cloning flanks on every designed
sequence. **That hypothesis was tested and REFUTED:**

| | designed grid | genomic tiles |
|---|---|---|
| longest common prefix | **0 bp** | 0 bp |
| longest common suffix | **0 bp** | 0 bp |
| invariant positions (of 150) | **3** | 0 |
| mean per-position entropy | **1.513 bits** | 1.987 bits (max 2.0) |
| median longest identical block between two random members | **13 bp** | **3 bp** |

No scaffold. Instead: the grid is **recombined from a small element vocabulary**, so any two grid sequences
share a median 13 bp identical block with each other while sharing no fixed position. That signature — long
shared blocks with zero common prefix — is what a fully-crossed design produces, and it is what the
discriminator is detecting. The measurement now ships inside the artifact (`sequence_diversity`) rather than
living in a shell, and a synthetic test pins both directions: a combinatorial set must show the long-block
signature, and a set that genuinely *does* have a scaffold must report it.

**Consequence for G-A that was not in the plan:** with only ~3 × 8 × 8 × 8 × 8 element variants and 13 bp
shared blocks, a model scoring well on the grid under leave-element-out can be exploiting the four element
axes it *did* see while the held-out axis contributes little. The 0.5914 is not necessarily learned promoter
grammar. This does not change G-A's bar (which is set on genomic tiles) but it does mean the grid arm should
be read as a *consistency check*, never as evidence of transferable grammar.

## Finding 3 — "genomic" is not ONE distribution, which refines G-E's question

The two genomic sets are themselves separable: **tile vs natural 0.6664 under k-mer and 0.8338 once
positional features are added.** The jump tells you where the difference lives — positional structure, not
bulk composition. That is mechanistically sensible: peak tiles are selected from called expression peaks and
so are enriched for active promoters at characteristic offsets, while upstream-of-CDS windows are every
CDS's upstream region regardless of activity.

**So G-E's question is not "grid or genomic?" but a three-point geometry.** An oracle validated on peak
tiles is *not* thereby valid on arbitrary genomic sequence — there is a measurable 0.83 gap between those
two genomic sets. That is a real sharpening of the loop-soundness question and it was invisible before this
run.

## What this changes

1. **G-A proceeds, premise confirmed.** Its genomic-tile bar stands; its grid arm is demoted to a
   consistency check with the element-vocabulary caveat recorded.
2. **G-E's framing is upgraded** from a binary to a three-point geometry, and the reference geometry it
   needs now exists — generated sequence can be placed against all three sets rather than two.
3. **A concrete prediction for G-D/G-E:** if generated sequence separates from genomic at ~1.0 the way the
   designed grid does, it is artificial in the same measurable sense; if it lands near the tile/natural
   cluster it is genome-like. The instrument for that call is now built and validated.

## Honest limits

- **A k-mer / positional discriminator sees local composition and coarse positional structure.** Two sets it
  cannot separate could still differ in long-range organisation. A 1.0000 is therefore a strong positive
  claim, while the 0.6664 is a *lower bound* on how different tile and natural are.
- The natural set is upstream-of-CDS windows — a **proxy** for promoters, not mapped transcription start
  sites.
- Peak tiles are genomic-but-**selected** (from called peaks), not a uniform genomic sample.
- The band thresholds (0.80 / 0.65) are **asserted**; the self-control ceiling that validates the instrument
  is **derived**.
- AUROC 1.0000 is a saturated measurement: it cannot distinguish "very different" from "infinitely
  different", so the mechanism table above is what carries the interpretation, not the AUROC itself.
