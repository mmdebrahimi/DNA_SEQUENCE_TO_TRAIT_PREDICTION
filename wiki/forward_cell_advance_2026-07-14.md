# Forward cell — panel extension + Regime-A router + AlphaMissense scope wall (--advance, 2026-07-14)

Third `--advance` on the forward "edit E. coli → predict phenotype" cell, through the three named remaining
next-steps in best-judgement VOI order.

## Step A — the forward cell generalizes across E. coli proteins

BLOSUM62 deterministic forward numbers on the cached ProteinGym E. coli DMS panel (all WT-mismatch=0,
polarity confirmed):

| protein (assay) | phenotype | n | Spearman(BLOSUM, DMS) |
|---|---|---:|---:|
| TEM-1 β-lactamase (Stiffler) | ampicillin fitness | 4,996 | 0.347 |
| CcdB toxin (Tripathi) | function | 1,663 | 0.248 |
| IF1 (Kelsic) | fitness | 1,367 | 0.182 |
| MlaC (MacRae) | fitness | 4,007 | 0.182 |
| DHFR/folA (Nguyen) | growth | 2,916 | 0.152 |
| EnvZ (Ghose) | fitness | 1,121 | 0.101 |

6 E. coli proteins, real numbers — the cell is not TEM-1-specific. **The ESM2 learned lift generalizes too:**
CcdB ESM2 **0.512 vs BLOSUM 0.248 (+0.264)**, mirroring TEM-1's 0.347→0.732 (+0.385) — both roughly double.

## Step B — Regime-A forward router (the capstone)

`dna_decode/forward/router.py::predict_edit` is one entry point that auto-classifies an edit into a G2P
regime and routes it to the RIGHT predictor (`forward_router_demo.py`, real predictors on every regime):

| edit | regime | predictor | result |
|---|---|---|---|
| **rpoB S450L / rifampicin** | A (determinant) | **real WHO TB catalogue** (21 grade-1/2 RIF substitutions) | **R** |
| rpoB A286V / rifampicin | A | WHO catalogue | S (not catalogued) |
| **blaTEM-1 E210K** (ampicillin) | B (molecular) | **cached ESM2 DMS predictor** | preserved (+2.69) |
| blaTEM-1 P181R (ampicillin) | B | ESM2 | damaging (−13.78) |
| acetate-growth-rate edit | C (organismal) | — | **ABSTAIN** (closed negative) |

The router **never sends a resistance edit to the likelihood predictor** (Regime A → catalogue; the
resistance-conservativeness finding says ESM fails on the antagonistic direction) nor an organismal trait to
any predictor (Regime C → abstain). This unifies the two working regimes + honest abstention behind one call.
The WT-coordinate gate even caught a real numbering error mid-build (Ambler "M69L" — position 69 is `T` in
the ProteinGym sequence — refused loudly, as designed).

## Step C — AlphaMissense method: an honest SCOPE WALL (not built)

AlphaMissense is **human-proteome-only** (Cheng et al. 2023). The cached AM data has **0 coverage** for
TEM-1 (UniProt P62593) or any bacterial protein (0 hits across every `am_uniprots*.txt` + `am_pg.tsv`, which
is all human UniProts). So AlphaMissense **cannot** be a method for the bacterial E. coli forward cell — this
is an external scope boundary, not a code gap. The ESM2 upgrade already provides the learned lift for
bacteria; AM would only apply if/when the forward cell extends to human/eukaryotic proteins (out of the
current E. coli scope). Documented, not forced.

## Status

Forward cell now: **6 E. coli proteins** (deterministic) + **learned ESM2 lift confirmed on 2** (TEM-1 +0.385,
CcdB +0.264) + a **regime-routing capstone** with real predictors on Regime A (WHO catalogue) and B (ESM DMS)
and honest abstention on C. AlphaMissense ruled out for bacteria (human-only). 19 forward tests. Frozen
decoder surface byte-unchanged (`verify_lock OK`); `dna_decode/forward` NON-frozen. This is a genuine plateau
for the E. coli forward cell — remaining moves need a new substrate (eukaryotic proteins for AM; a real
E. coli QRDR SNP catalogue for a native Regime-A E. coli cell) that is a fresh direction, not an increment.
