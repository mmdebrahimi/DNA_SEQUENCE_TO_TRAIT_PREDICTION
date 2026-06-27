# Functional-alphabet within-lineage probe — 5-drug synthesis (2026-06-27)

The complete cross-drug result of the non-neural functional-alphabet probe: does a curated functional-
determinant alphabet, or a base-level k-mer (raw-sequence) alphabet, carry E. coli AMR signal WITHIN a
lineage (the de-confounded test that the closed 0-for-4 embedding bet failed)? This decides whether a
learned "DNA-LLM" has any headroom over the shipped deterministic decoder. CPU-only, no GPU.

## The 5-drug table (each: 20 shared MLST lineages, 120 strains 60R/60S, 174-180 within-lineage pairs)

| drug | mechanism regime | functional WL | k-mer WL | gap | perm p | verdict |
|---|---|---|---|---|---|---|
| ciprofloxacin | QRDR target-site (point) | 1.000 | 0.721 | +0.28 | 0.057 | TIES |
| tetracycline | efflux + ribosomal-protection | 0.963 | 0.540 | +0.42 | 0.000 | BEATS_KMER |
| ceftriaxone | beta-lactamase (acquired) | 0.939 | 0.606 | +0.33 | 0.000 | BEATS_KMER |
| meropenem | carbapenemase (acquired) | 0.928 | 0.644 | +0.28 | 0.000 | BEATS_KMER |
| gentamicin | aminoglycoside-modifying enzymes | 0.594 | 0.667 | -0.07 | 0.865 | FAILS |

## The conclusion: NO DNA-LLM headroom on any of the 5 drugs

The headroom signature a learned model would need is **k-mer (raw sequence) carrying within-lineage signal
the determinant alphabet misses** — i.e. k-mer WL high, or k-mer >> functional. **It never appears:**

- **4/5 drugs (cipro, tet, cef, mero):** the curated determinant alphabet separates R/S within-lineage at
  0.93-1.00 while raw-sequence k-mer is at/near chance (0.54-0.72). The signal IS the curated determinants
  = the shipped deterministic decoder. A DNA-LLM (which learns from the raw-sequence side, the weak one)
  would have to rediscover those determinants to match — becoming a worse determinant caller. BEATS_KMER /
  TIES both mean "decoder is the right tool", not "learning has an edge".
- **1/5 (gentamicin) FAILS — and it's the honest non-confirmation, not a headroom opening:** the determinant
  alphabet is at CHANCE within-lineage (0.594) — BUT so is k-mer (0.667). Within-lineage gentamicin
  resistance is hard for BOTH alphabets. That points to weak/absent within-lineage sequence signal (or the
  noisier relaxed-MIC BV-BRC label, or aminoglycoside-modifying-enzyme diversity), NOT to raw-sequence
  beating determinants. A learned model needs sequence signal to exploit; k-mer~chance says little is there.

## What this settles
The "DNA-as-language / DNA-LLM" path is a **confirmed soft-negative across 5 drugs spanning 5 mechanism
regimes** (target-site point, efflux/ribosomal, two acquired beta-lactam classes, aminoglycoside-modifying).
The cheap CPU probe did its full job: it gated the expensive GPU build and returned "no headroom over the
shipped decoder" everywhere. Re-confirms the standing project thesis from a fresh angle: **the binding
constraint is labels, not compute.** No GPU spend is justified.

The gentamicin FAILS is the one result worth a follow-up footnote (not a reopen): is it label noise
(relaxed BV-BRC MIC) or genuine mechanism diversity? A strict-MIC gent cohort would disambiguate — but it
does not change the headline (no alphabet, learned or curated, separates gent within-lineage here).

## Artifacts
Per-drug packets `wiki/functional_alphabet_probe_{ciprofloxacin_n147,tetracycline_n118,gentamicin_n120,
ceftriaxone_n120,meropenem_n120}_*.{md,json}`; cohorts `data/processed/shared_lineage_*_cohort.parquet`
(+ D: backup tarball, see `wiki/dna_llm_cohorts_durability_2026-06-27.md`); code
`dna_decode/eval/functional_tokens.py` + `scripts/functional_alphabet_probe.py`; prior closeout
`wiki/functional_alphabet_probe_closeout_2026-06-26.md`.
