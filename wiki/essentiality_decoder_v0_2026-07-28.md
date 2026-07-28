# Essentiality decoder v0 — the conserved-core R1 catalogue (E. coli, 2026-07-28)

**New cell (ratified 2026-07-28): single-gene KO -> essential/non-essential.** v0 = the deterministic
**conserved-core** decoder (the R1 "determinant catalogue" paradigm, for essentiality). Label-INDEPENDENT by
construction: it reads gene FUNCTION (product text), never an essentiality label, so it is buildable +
validatable without the (externally-walled) gold-standard labels. Organism: E. coli K-12 MG1655.

## Result — validated by size + composition (NOT a biased label AUROC)

Applied to all **4318** E. coli genes (NCBI GCF_000005845.2 feature table): **208 predicted essential.**

- **Size:** 208 vs the known E. coli essentialome ~300 (Keio 303 / PEC 230 / consensus 248). Right ballpark,
  conservative — as a conserved-core PRIOR should be (it captures the universal core, misses the tail).
- **Composition** (matches the established essentialome): translation/ribosome 67 · envelope/cell-wall 38 ·
  replication+transcription 26 · tRNA-synthetase 25 · cell-division 11 · other-core 41.
- **Top predictions are genuine known-essential genes:** rpsT, ileS, lptD, ftsL/I/Q/A/Z, murC-G, secM.

This IS the validation: a conserved-core decoder that reproduces the known essentialome's SIZE + COMPOSITION +
identity on a genome it never saw labels for. (The AMR-track discipline: validate against established biology,
not a biased in-set metric.)

## Secondary (heavily caveated, underpowered): UniProt-label agreement

On 112 genes with a UniProt disruption-phenotype binary label: **spec 0.95 / sens 0.18** (tp13 fn58 fp2 tn39).
HIGH-PRECISION (when it predicts essential it is almost always right), LOW-RECALL — because UniProt's label is
broad+noisy (includes conditional/regulatory genes: lexA, recB, spoT, which the core-prior correctly excludes)
AND the prior is conservative. **Do NOT read this as the decoder's true sens/spec** — the UniProt label is
curation-biased (64% essential base rate vs the true ~7%), incomplete (2.5% coverage), and noisy.

## E1 DATA WALL (honest, named) — the gating external input

Free, machine-readable, GOLD-STANDARD essential-gene LABELS (Keio 303 / PEC 230 / consensus 248) are NOT
freely programmatically fetchable: OGEE is fully down (all hosts), DEG's bulk download is form-gated (403),
DeeplyEssential ships data only via a Google-Drive snapshot, the DepMap manifest API was finicky, UniProt
disruption-phenotype is incomplete+noisy. This is an EXTERNAL wall (a single file the user provides or a
source the user authorizes), NOT a code wall. The v0 conserved-core decoder is deliberately designed to NOT
need it (function-based), so this wall does not block v0 — it blocks the proper label-based AUROC (v0.1).

## Scope / honesty
- v0 is the R1 conserved-core PRIOR — high-precision, conservative-recall by design; the organism-specific /
  conditionally-essential TAIL is the R2 learned-complement target (Family E3, deferred).
- Validated by size+composition vs established biology, NOT a per-gene AUROC (the gold-standard label is walled).
- Cross-organism transfer (Family E4, the "any organism" test) + human (DepMap) are the next organisms.

## Reproduce
`dna_decode/essentiality/core_decoder.py` (the decoder) + the D: feature table (GCF_000005845.2). Data on
D:/dna_decode_cache/essentiality/. Deterministic, offline, no labels needed.
