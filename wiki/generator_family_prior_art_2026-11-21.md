# The generator family is the wrong one, and the right one is published, on E. coli, with wet-lab numbers

**Date:** 2026-11-21 · **Type:** literature scan, NOT a measurement · **Trigger:** the user asked
whether the surviving 0.1186 gap to the 3-mer Markov null is the better target than refining a delta
that failed to resolve, and whether "AlphaGenome 2 seems to have done it".

**Status of this memo:** it records what the published record says and what that implies. It contains
no new measurement of our own system. Nothing here is validated here.

## The short answer to the question asked

**AlphaGenome 2 is real and has not done it.** It was released **11 November 2026** — genuinely after
my knowledge cutoff, so the user was right that it exists and I did not know. But DeepMind's own
announcement says it *"predicts how DNA sequences and their variations impact specific human cell
types"* and **never describes generating DNA**. The post does not mention synthetic-DNA design at
all; its closest claim is helping researchers *"prioritize which experiments to do next"*.
Independent coverage is explicit that it is *"a conditional model, which is different from a
generative model that creates new sequences"*.

**It is the half we already have.** Our oracle was measured converting rho 0.59 into **79% of
achievable selection value at k=10** (`wiki/glm_oracle_selection_result_2026-10-07.md`). Scoring was
never our bottleneck. AlphaGenome 2 is a better scorer, for **human and mouse**, **API-only with
weights unreleased**, non-commercial research only, commercial access in private preview.

**v2 moves FURTHER from us, not closer.** `wiki/alphagenome_lessons_actionable_2026-10-07.md` already
established the v1 split — the architecture does not transfer (1 Mb context, U-Net, 2D contact
branch, and eleven heads predicting splicing, histone marks and chromatin contacts that a bacterium
does not have, against a 100-300 bp bacterial promoter) while the training METHOD does. v2's headline
change is that it became **conditional on one human tissue / cell type / state at a time**, because
the original *"spread the model's capacity thin"*. Conditioning on human cell identity is further
from a bacterial promoter, not nearer.

Reported benchmarks, all DeepMind's own, against "the original AlphaGenome fine-tuned on single-cell
data and other specialized tools": on TenK10K (28 immune cell types) **30% vs 16%** causal variants
identified and **43% vs 23%** directional gene effects. Trained on **107 million** single-cell
profiles across **68 tissues** (one secondary source says "more than 100 million" — the figures
disagree across outlets).

## The literature that HAS done it is a different one, and it is on our organism

### Adversarial: Wang et al., WGAN-GP synthetic E. coli promoters

The closest published work to our problem, and it is close: **E. coli K-12 MG1655**, trained on
**14,098 experimentally identified promoters** (mostly sigma-70), generating **50 nt** upstream of
the TSS.

- **38 of 83 tested sequences (45.8%) were functional in the wet lab** (sfGFP, DH5-alpha, high-copy
  plasmid, BH FDR < 5%, three biological replicates). Three synthetic promoters matched or exceeded
  all six positive controls.
- Genuinely novel: mean CLUSTAL similarity to the training set **62.04%**, max **72%**, and no
  significant BLAST hit against the E. coli genome.
- Of the top ten 6-mers, WGAN-GP shared **five** with natural promoters, versus one each for DCGAN
  and a PSSM baseline. The -10/-35 spacer clustered at **16-18 nt**.
- Honest limits in that paper: only 45.8% of SELECTED sequences worked, so most failed; the
  expression predictor correlated at Spearman **0.55**; one host, one reporter, one vector, five
  random controls.

### Diffusion: the successors

- **DiPro** — unconditional diffusion for E. coli promoters, which explicitly aims to *"faithfully
  capture native promoter architectural rules without relying on supervised recognizers or
  discriminators"*, with functional validation. (Primary text not retrieved: PMC served a
  cookie-enablement page and europepmc returned 403, so these details come from search snippets.)
- **Multinomial diffusion (MDM)** for E. coli promoters, reporting it outperforms existing
  generative approaches in silico.
- A cyanobacterial comparison found diffusion-generated promoters more similar to natural ones than
  VAE output, *"particularly in the key -35 and -10 regions"*, while the VAE *"often generates
  additional and spurious peaks"*.

## Why our generator loses, as the literature names it

**Absorb & Escape** (NeurIPS 2024) states the two failure modes as complementary:

- autoregressive models **"fail to capture the global property of DNA sequences"**
- diffusion models **"tend to produce errors at the base pair level"**

Its toy HMM test makes it concrete: DiscDiff produced 0 illegal start tokens against HyenaDNA's 812,
but 110,192 incorrect transitions against HyenaDNA's 3,586. Its hybrid starts from a diffusion sample
and lets an AR model repair local errors.

A separate design paper puts the mechanism plainly: *"left-to-right causal modeling ... is suboptimal
for DNA where regulatory relationships are inherently bidirectional."*

**Our generator is a left-to-right autoregressive LM**, used zero-shot and then LoRA-adapted. That is
the family the literature identifies as weak on global structure.

**But the picture is not one-sided, and saying so matters.** In the Dirichlet Flow Matching paper the
autoregressive model was the **strongest baseline** in the promoter-design experiments. D3LM reports
diffusion LMs *"substantially outperform the previous best result of 29.16 from autoregressive
models"*. Most of these comparisons were run by the proponents of the winning method, and the search
did not find an independent head-to-head. So "diffusion often wins on these benchmarks" is supported;
"autoregressive models fundamentally cannot do this" is not.

## Three findings that bear on OUR measurement specifically

**1. Nobody publishes the metric we are using.** The flagship adversarial paper *"never operationalizes
'indistinguishable from natural' as a measured criterion, such as discriminator accuracy"* — it judges
quality by motif logos and k-mer statistics. Absorb & Escape likewise *"does not train a real-vs-generated
classifier"*, using S-FID (Frechet distance in Sei feature space; 57.4 for DiscDiff vs 504.0 for DDSM
at 256 bp), motif-distribution correlation, BLASTN similarity and Enformer-based SSE.

Consequence: **our held-out discriminator AUROC is a stricter and less common bar than this field's
standard, and our 0.7787 has no published comparator.** Our "loses to a 3-mer Markov chain" framing
may be harsher than how this literature would judge the same generator. That does not make the gap
unreal — it makes it unanchored.

**2. External corroboration of our over-training suspicion.** WGAN-GP trained for 160 epochs but
**reported best results at about 12**, and the paper does not explain how that point was chosen. Our
own `TRAJECTORY_UNRESOLVED` reading showed step means rising +0.0016 / +0.0187 / +0.0255 with positive
counts 6 / 8 / 9 across 10 seeds — suggestive of progressive over-training, and deliberately not
upgraded to a finding. Independent evidence that early stopping matters enormously on this exact task.

**3. Our ARCHITECTURE is already the field's architecture.** WGAN-GP selected which sequences to test
with *"a CNN expression predictor, not the discriminator"*. regLM used *"supervised regression models
... to discard generated sequences whose predicted activity did not match the prompt"* (only 1.1%
discarded). Generator plus a separate oracle for selection is exactly what we have, and our oracle is
measured adequate. **What differs is the generator family, not the loop.**

`regLM` (Genome Research 2024, HyenaDNA-based, *"the first use of language modeling for DNA in a
generative setting"*) is the closest published analogue to our approach, and it worked — on yeast
promoters and human enhancers, with generated sequences at least 25 bp edit distance from training.

## What this implies, and what it does not

**Implies:** the highest-value move is a **generator-family change** — discrete diffusion or an
adversarial generator trained on E. coli promoters, keeping our oracle — not scaling the LM (measured
dead: ESM2 peaks at 650M, 3B and 15B regress) and not more LoRA data (that proposal was withdrawn
after its own review).

**Does NOT imply** that our 0.1186 gap would close under diffusion. No paper reports our metric, so
nobody has shown a diffusion generator beating a 3-mer Markov chain under a held-out
k-mer-plus-positional discriminator. That is an open empirical question, and the honest expectation is
that the metric is strict enough to be unflattering to any of these methods.

**Does NOT touch** the surviving measured claim: the generator is worse than its own seed-matched
Markov null at 10 of 10 seeds, mean gap +0.1186.

## Honest limits of this memo

- A literature scan (6 searches, 4 primary-source fetches), not a measurement. Nothing here was
  reproduced.
- Several AlphaGenome 2 figures come from secondary coverage published days after release, and the
  outlets disagree (107M vs "more than 100M" cell profiles). DeepMind's own post is the one primary
  source used.
- DiPro's primary text was NOT retrieved (PMC cookie wall, europepmc 403); its details are from
  search snippets and should be verified before being relied on.
- The AR-vs-diffusion comparisons are almost all authored by proponents of the winning method; no
  independent head-to-head was found.
- Wet-lab fractions (45.8%) are for sequences that an oracle had already SELECTED, not for raw
  generator output, and are from one host / one reporter / one vector.
