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

| drug | regime | functional within-lineage concordance |
|---|---|---|
| ciprofloxacin | concentrated (QRDR target-site) | **1.000** |
| tetracycline | distributed (efflux/ribosomal/mobile) | **0.963** |

The curated determinant alphabet separates R/S within-lineage on BOTH — even the "distributed" drug where
headroom was expected. **There is no demonstrated headroom for a learned model on either.** The distributed-
mechanism escape hatch (section 19) is now itself a soft-negative. (Caveat: the tet number is the functional
arm alone — the k-mer comparator is blocked on a disconnected D: — but both readings of a k-mer result point
to no learned-model headroom; see `wiki/functional_alphabet_probe_tet_partial_2026-06-26.md`.)

## Status
The cheap CPU probe did its job: it gated the expensive build and returned a "no headroom over the shipped decoder" signal on BOTH cipro and tet. **The DNA-LLM-via-functional-alphabet path is parked (soft-negative across two mechanism regimes), no GPU spend.** Artifacts: `wiki/functional_alphabet_probe_{smoke,n147,tet_partial}_2026-06-26.{md,json}`; code `dna_decode/eval/functional_tokens.py` + `scripts/functional_alphabet_probe.py` (commit 753a1a9).
