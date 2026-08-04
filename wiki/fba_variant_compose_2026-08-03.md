# Compose `forward` + `fba`: point-mutation -> cell-level trait (2026-08-03)

**What shipped:** `dna-decode fba --gene <G> --mutation <M> --protein-seq <S>` chains the two decoders to
turn a **single missense edit** into a downstream **cell-level trait** — closing the gap between the FBA
cell (which takes a full gene KNOCKOUT) and a real genotype edit (a point mutation).

```
missense edit  --forward-->  does it break the enzyme? (LOF)  --fba-->  cell-level trait
```

## The three honest paths (verified live on real gltA = b0720, an essential citrate-synthase gene)

| missense | forward (blosum62) | LOF call | FBA action | cell-level trait |
|---|---|---|---|---|
| `gltA D3W` | damaging (raw −4.0) | **LOF** | knockout | **NON-VIABLE (essential gene lost)** |
| `gltA I30L` | preserved (raw +2.0) | TOLERATED | wildtype | viable (no metabolic change) |
| `gltA D13T` | uncertain (raw −1.0) | UNCERTAIN | **conditional** | **reported both ways** — if LOF: NON-VIABLE / if tolerated: viable |

`dna-decode fba --gene pgi --mutation D2W ...` (a NON-essential gene) → LOF → *viable, altered flux* (a
damaging edit at a dispensable gene doesn't kill the cell). The cell-trait tracks BOTH the edit's effect
on the protein AND the gene's role in the network.

## Honest rails (load-bearing)

- **The chain inherits two validations, and adds one heuristic.** Upstream (missense→LOF) carries
  `forward`'s DMS validation (ProteinGym / MaveDB, a validated RANK). Downstream (LOF→cell trait) carries
  `fba`'s Keio validation (accuracy 0.954). What is **NEW and unvalidated** is the binarization of a
  ranker into a LOF call — so `lof_call` uses **forward's own method-aware threshold** (BLOSUM62 ≤ −2 /
  ESM2 ≤ −5 / hybrid ≤ 0.33), labelled a **heuristic, NOT a calibrated LOF probability**.
- **`uncertain` is never forced.** When forward is uncertain/abstains, the tool reports the **conditional
  both-ways** (if-LOF vs if-tolerated) rather than fabricating a binary — this is the anti-theater rail.
- **Scope:** METABOLIC traits only; the missense must be in an iML1515 metabolic gene. NOT clinical.
- **`--forward-method esm2/prosst/gemme/hybrid`** upgrades the upstream LOF call to the stronger DMS-validated
  scorers (needs `dna-decode[forward]`); BLOSUM62 is the offline default.

## Where this sits on the edit→cell-trait ladder

This is the rung above the FBA gene-KO cell: **point mutation** (the user's original "minor edit to a
genotype"), not just full knockouts. Tests: `tests/test_fba_compose.py` (pure decision logic + 4 real-model
smokes). Frozen AMR surface byte-unchanged (READ-only composition of two existing cells).
