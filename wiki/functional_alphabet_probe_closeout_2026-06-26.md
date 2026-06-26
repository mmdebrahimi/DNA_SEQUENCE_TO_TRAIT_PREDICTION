# Functional-alphabet probe — closeout + DNA-LLM path decision (2026-06-26)

## Result (N=147 cipro, the de-confounded within-lineage gate)
- **Verdict: TIES** (functional within-lineage 1.000 vs k-mer 0.721; gap +0.279; paired in-MLST permutation **p = 0.0565**, just over 0.05; 6 shared lineages / 43 pairs).
- Smoke (N=40): UNDERPOWERED by construction (MLST-balanced → 40 unique MLSTs, 0 shared R/S lineages).

## What it means (the honest reading — this is the load-bearing part)
1. **The alphabet DOES beat k-mer within-lineage** (1.000 vs 0.721). Directionally, a function-level alphabet discriminates R/S *within a lineage* far better than base-level k-mers — k-mers leak lineage, function tokens don't. So hypothesis #2 (alphabet matters) is **directionally supported**.
2. **BUT the winning functional alphabet IS the curated determinant set.** The tokens that drive the 1.000 are the QRDR alleles (`gyrA_S83L`, `parC_*`) + acquired `qnr` mechanisms — i.e. **exactly the determinants `amr_rules.py` (the shipped deterministic decoder) already uses.** A within-lineage concordance of 1.0 is not "a learned model discovered mechanism"; it is "the curated determinant token generalizes across lineages" — which is *precisely why the deterministic decoder works*. The win is the decoder, re-expressed.
3. **The result is underpowered.** 6 shared lineages / 43 pairs gives a coarse permutation null (95% band [−0.35, +0.31]); p = 0.0565 is suggestive, not established. Establishing significance needs **more R+S-sharing lineages** — a cohort/label problem, not a compute one.

## Decision on the DNA-LLM path: DO NOT proceed to neural #3 / DO NOT spend GPU
- The neural contrastive step (#3) is justified *only* if a function-level alphabet shows **headroom over the curated determinants**. This probe shows the opposite: the function alphabet's within-lineage signal **is** the curated determinants — there is **no demonstrated headroom** for a learned model on concentrated-mechanism cipro. A transformer to "discover" `gyrA_S83L` adds nothing; AMRFinder already calls it and the deterministic decoder already uses it.
- This **re-confirms the project thesis**: compute was never the binding constraint. The binding constraints are **labels** (more shared-lineage strains to even power the metric) and **mechanism concentration** (where signal is concentrated + already curated, a learned LLM has no edge).
- Consistent with the closed 0-for-4 embedding negative: embeddings learned lineage; here the *only* thing that beats lineage is the hand-curated determinant alphabet — i.e. the deterministic decoder.

## What WOULD change this (the non-foreclosed levers, both data not compute)
- **More shared-lineage data** to push p < 0.05 and power the metric (cohort expansion — external/data wall).
- A **DISTRIBUTED-mechanism drug** (not concentrated like cipro QRDR) where curated determinants are incomplete and a learned representation might have real headroom — e.g. the tet failure mode (mobile elements). That is the only regime where the DNA-LLM bet has a plausible edge, and it is a *new-cohort + label* question, not a GPU question.

## UPDATE 2026-06-26 PM — the tet test RAN, and it FALSIFIED the headroom prediction
Section 19 predicted tet (a distributed-mechanism drug) *might* show headroom because its determinants are
"incomplete". A powered tet shared-lineage cohort was built (118 strains, 20 shared lineages, 174 pairs) and
the **functional-determinant within-lineage concordance came back 0.963 (p<0.0001)** — NOT the predicted
fail/tie. So across BOTH mechanism regimes tested:

## FINAL — D: reconnected, the formal k-mer comparison ran (2026-06-26, both arms, both drugs)

| drug | regime | functional WL | k-mer WL | gap | verdict |
|---|---|---|---|---|---|
| ciprofloxacin | concentrated (QRDR target-site) | 1.000 | 0.721 | +0.28 | TIES (p=0.057) |
| tetracycline | distributed (efflux/ribosomal/mobile) | **0.963** | **0.540** | **+0.42** | **BEATS_KMER (p=0.0)** |

(tet: 20 shared lineages, 174 within-lineage pairs — properly powered, ~4× the cipro probe.)

### The honest interpretation — BEATS_KMER STRENGTHENS "no headroom", it does not open one
The functional alphabet is the **curated determinant set the deterministic decoder already uses** (tetA/tetB
efflux, tetM/tetO ribosomal-protection, acrB/acrR). So this comparison is **decoder-alphabet vs raw-sequence
(k-mer) alphabet**:
- On tet the **base-level k-mer alphabet is at CHANCE within-lineage (0.540)** — raw sequence carries no
  within-lineage tet-resistance signal once lineage is conditioned out (tet resistance = acquired-gene
  presence; the top-N k-mer vocab is the core genome, which misses it within a lineage).
- The **determinant alphabet carries it cleanly (0.963)**.
- A learned **DNA-LLM operates on the raw-sequence side** (the chance-level one). To match the decoder it
  would have to REDISCOVER the curated determinants — i.e. become a worse determinant caller. So
  BEATS_KMER = "the determinant decoder is the right tool; raw-sequence learning has nothing to add here."
- cipro (TIES) and tet (BEATS_KMER) differ only in *how badly* raw k-mer loses within-lineage (weak signal
  from QRDR SNPs vs chance for acquired genes) — both land at **no demonstrated headroom for a learned model
  over the curated deterministic decoder.** The distributed-mechanism escape hatch (section 19) is closed.

## Status
The cheap CPU probe did its full job: it gated the expensive GPU build and returned a formal "no headroom
over the shipped decoder" verdict on BOTH cipro (TIES) and tet (BEATS_KMER, determinant-side). **The
DNA-LLM-via-alphabet path is a confirmed soft-negative across two mechanism regimes — no GPU spend.**
Re-confirms the standing thesis: the binding constraint is labels, not compute. Artifacts:
`wiki/functional_alphabet_probe_{smoke,n147,tet_n118}_2026-06-26.{md,json}` (the tet_partial packet is
superseded by tet_n118); code `dna_decode/eval/functional_tokens.py` + `scripts/functional_alphabet_probe.py`.
