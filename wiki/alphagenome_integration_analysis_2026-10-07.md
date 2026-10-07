# AlphaGenome: what it is, what we can integrate, and what its LOGIC gives us

**Date:** 2026-10-07
**Primary source read:** the Kohli interview transcript (`youtu.be/Wkaw03p3BrM`, 28,293 chars, 777 segments,
fetched in full). Pushmeet Kohli leads this at DeepMind, so the architecture/availability claims below are
**first-party statements**, not press summaries — which is why they are quoted rather than paraphrased.
**Verification in flight:** three agents are checking the paper/architecture/licence, the bacterial MPRA
datasets, and the sequence-to-function lineage. **Dataset accessions in §5 are NOT yet verified.**

> **§0 — PRIMARY-SOURCE VERIFICATION LANDED. Read this before anything below it.** Several figures in §1
> are now pinned, two of my own statements are corrected, and one licence clause changes the integration
> options. Corrections are in **§7**, which supersedes the body on every point it touches.

---

## 1. What AlphaGenome actually is (first-party)

- A **supervised sequence-to-function** model: DNA in, **molecular phenotypes** out (gene expression,
  splicing, and other regulatory readouts). Not a generative model.
- **1 Mb input context AT base-pair output resolution**, and Kohli's claim is that doing *both at once* is
  the novelty: its predecessor **Enformer** ("informer" in the transcript) had *both* a smaller window *and*
  coarser resolution — "we were not looking at what's happening at the base pair level... we were looking at
  a large number of base pairs taken together".
- **Human + mouse.** And the multi-species point is explicit and interesting: training on mouse too makes
  the human predictions better, because "by forcing the model to explain the same effects in the context of
  the mouse genome, we are asking the model to think broader and think about general principles" — otherwise
  "it can memorize... mechanisms that only work in the case of the human genome".
- **Matched or beat the strongest comparison on 25 of 26 variant-effect benchmarks.**
- **Effect-size correlation ρ ≈ 0.5** for mutations affecting gene activity (interviewer cites the paper,
  Kohli confirms). Worth holding onto as calibration: that is *state of the art at DeepMind scale*, and our
  own forward cell sits at 0.35–0.76 on protein variants. This problem class is simply hard.
- **Weights are downloadable**, free for academic/non-commercial use, and **fine-tunable** — "for them to use
  it for their own research, fine-tune it, and the model weights are available."
- **AlphaGenome Atlas**: precomputes all **~9 billion** possible SNVs, petabyte-scale, "**30 times the size
  of the AlphaFold database**". Carries a calibrated **AVI** (variant impact) score — high means scrutinise,
  low means "high confidence it will not have a dramatic impact".

**Two claims in the circulated ChatGPT summary that I believe are WRONG and have flagged for verification:**
"published Oct 2024" and a "Nature 2024" citation. AlphaGenome is a **2025** release; an Oct-2024 Nature
citation most likely belongs to Enformer or AlphaMissense. Do not propagate those until the agent reports.

---

## 2. Can we integrate AlphaGenome into the bacterial GLM? **No — and the reason is structural**

It is **human/mouse**, and most of what it predicts **does not exist in a prokaryote**: splicing,
nucleosome/chromatin accessibility, and the 3D-contact effects that motivate its 1 Mb window. Pointing it at
an *E. coli* promoter is an out-of-distribution query against machinery the organism lacks. There is no
prokaryotic AlphaGenome (being confirmed).

**Where it could integrate for real: the HUMAN clinical track** (`pgx` / `clinvar` / `hla`). The Atlas is a
free, precomputed effect prediction for every human SNV. **But the honest framing matters:** AlphaGenome is a
*model*, not a *label*. Using it as ground truth for our human cells would be the circular-label failure this
repo has gates for (G1). It is legitimately a **comparator** or a **doubt-layer signal**, never a label.

---

## 3. What transfers as LOGIC — five things, and one is immediately buildable

**(a) A rich representation with a TRIVIAL head.** The best architectural detail in the interview. The
interviewer notes the paper puts "a tiny super simple layer of lasso regularization" on top of the model —
"a Formula 1 engine with a little bicycle bell on top". Kohli's answer is the design principle: *"most of the
work should be done by the model in terms of the feature representation and then the final predictor on top
should take a very simple, almost trivial form."* Directly applicable to our expression oracle: learn one
promoter representation, keep the per-task head linear.

**(b) Multi-task heads over many measured tracks**, not a single scalar. Our falsifier currently emits ONE
discriminator AUROC. The AlphaGenome shape is: predict many real functional readouts at once and let variant
effects fall out as deltas.

**(c) Context × resolution — and the bacterial version of this problem is already free.** AlphaGenome's hard
engineering win was base-pair resolution *at* 1 Mb. A bacterial promoter is ~100–300 bp and an operon a few
kb, so **we get the easy regime by default**: a 150 bp window at 1 bp resolution is trivially affordable.
The expensive part of their contribution is the part we do not need.

**(d) Multi-species training improving single-species performance** — stated by Kohli for human+mouse, and
**directly testable for bacteria**: train on *E. coli* + *B. subtilis* + *Pseudomonas* promoters, evaluate on
*E. coli*. This repo already has the regime machinery to score that honestly, and the claim is falsifiable.

**(e) Generative model proposes, sequence-to-function model scores.** Kohli says this unprompted about DNA
language models: *"you would have seen DNA language models like Evo come up and people using them to sample
different genomes... AlphaGenome gives you a lens by which you can study the properties of those genomes."*
**That is exactly the propose→score→select architecture already designed here**, framed the same way by the
people who built the scorer. Independent convergence on the design — and a reason to cite rather than claim
novelty for the shape.

---

## 4. The strongest external corroboration of this repo's own regime finding

Kohli, twice, unprompted — on what AlphaGenome does **not** do:

> "how do those molecular phenotypes then affect the property of the whole organism... **that's something
> that will still need to be connected and research will need to be done on that**."

and in next steps: *"connecting these molecular phenotypes... and showing their connections and impact on
properties at the organism level."*

So the frontier lab with a petabyte of predictions and 25-of-26 benchmark wins says the **molecular→organism
link is open research**. That is precisely the boundary `dna_decode/eval/regime.py` encodes from our own
measurements: constructed→**molecular** works, constructed→**organism** is partially open, natural-population
→organism zero-shot is closed. We are not behind on that axis; nobody has it.

---

## 5. Where this lands on our actual blocker — and it closes a loop

Our Family 3 (**expression oracle**) has been blocked on one question: *is there a free, measured bacterial
regulatory-sequence→expression dataset?* Three research agents died on a session rate limit before answering.

AlphaGenome's logic says exactly what that oracle should be: **supervised sequence→measured-function, rich
representation, trivial head, multi-task.** AlphaGenome itself works *because* it has enormous measured
functional data (ENCODE/GTEx-class). The bacterial analogue of that data is **MPRA / sort-seq promoter
libraries**, and the circulated research points at specific ones.

**Those accessions are UNVERIFIED and under active check** — the summary's AlphaGenome dates are already
suspect, so its dataset citations get the same scrutiny before anything is downloaded or planned around. The
leads being verified: an *E. coli* lacZYA promoter MPRA (~8,269 variants, GEO), a σ70 promoter MPRA
(~5,000), and an RBS library (~5,000). **If even one resolves to a real, downloadable
sequence→measured-value table, Family 3 is unblocked and the expression oracle becomes buildable** — which
in turn gives the generative loop the scoring head it currently lacks.

---

## 6. Honest limits of this analysis

- Based on **one interview transcript**, not the paper. Architecture specifics (layer types, parameter count,
  exact output heads, loss functions) are *not* established here and are being verified.
- "Weights are available" is a first-party verbal claim; the **licence text and the actual download path**
  need confirming before any dependency is taken. Non-commercial terms would matter for anything shipped.
- ρ ≈ 0.5 is second-hand within the transcript (interviewer states it, Kohli assents) — treat as approximate
  until read in the paper.
- No AlphaGenome component has been run here. Nothing in this document is a measurement.

---

## 7. VERIFIED against primary sources (2026-10-07) — supersedes the body where it conflicts

Every row was reached by fetching the primary page. Anything marked UNVERIFIED must not be repeated as fact.

### 7.1 Citation — the circulated summary was WRONG, as suspected

| | |
|---|---|
| "Published Oct 2024 / Nature 2024" | **CONTRADICTED.** Nothing supports 2024. |
| Announced | **25 June 2025**, DeepMind blog |
| Preprint | bioRxiv **27 June 2025**, DOI `10.1101/2025.06.25.661532`, Avsec et al. |
| **Journal version** | **Nature, 28 Jan 2026, vol 649, pp 1206–1218**, DOI `10.1038/s41586-025-10014-0`, first author **Žiga Avsec** |

The "2024 Nature" attribution most plausibly belongs to **Enformer** (Nat Methods 2021) or **AlphaMissense**
(Science 2023) — neither of which is a 2024 Nature paper either.

### 7.2 Architecture — now pinned

- **U-Net style, conv AND transformer** (not either/or): encoder downsamples, transformers with inter-device
  communication in the middle, decoder upsamples. Convolutions carry local patterns; transformers carry
  long-range dependencies such as enhancer–promoter interactions.
- **1 Mb input.** Output resolution is **1 bp AND 128 bp for 1-D tracks, 2,048 bp for contact maps** — the
  abstract says *"up to single base pair resolution"*. **My §1 "base-pair resolution" was too blanket**: true
  of most 1-D tracks, not of everything.
- **11 output types / 5,930 human tracks / 1,128 mouse tracks.**
- **8 TPUv3** with sequence parallelism, ~4 hours, half the compute of Enformer.
- **Parameter count: UNVERIFIED** — absent from every primary page (Orbax checkpoints expose no total).
- **Polyadenylation is NOT among the 11 heads** — that is a Borzoi head and a likely source of confusion.

### 7.3 CORRECTION to my own §3(d) — the multi-species claim is interview-only

Kohli states that training on mouse improves the human predictions. **The paper contains no ablation
supporting it.** Joint training is confirmed; *the benefit is not*. So §3(d) stands as a **hypothesis worth
testing in bacteria** — which is what I proposed — but must **not** be cited as an established external
result. Downgraded.

### 7.4 THE LICENCE CLAUSE THAT CHANGES THE INTEGRATION OPTIONS

Weights **are** downloadable (since ~14 Jan 2026: HF `google/alphagenome-all-folds` + four CV folds, and
Kaggle), **gated**, under bespoke **non-commercial AlphaGenome Model Terms** — *not* Apache (that covers the
client code only). And, verbatim from the FAQ:

> *"Outputs generated by AlphaGenome and other information provided in AlphaGenome Atlas are for
> non-commercial use only and **should not be used for the training of other machine learning models**."*

**That forbids using AlphaGenome or Atlas output as a label, a distillation target, or training signal of any
kind.** I had already ruled it out on *circularity* grounds (a model is not a label — gate G1); there is now
an independent *licence* reason. Read-only **comparator** use remains permissible; training on it does not.

Also: *"must not be used for clinical decision-making"*, and local inference wants **at least an H100** —
decisively not this hardware.

### 7.5 Atlas and AVI — confirmed, and AVI is MORE than I said

- Atlas preprint: Cheng et al., **medRxiv 20 Sept 2026**, DOI `10.64898/2026.09.16.26363192`.
- **~9 billion SNVs** and **1 PB** both VERIFIED, plus **>100 million indels** (gnomAD / UK Biobank /
  All of Us); ~27,000 scalar predictions per variant.
- **"30× the AlphaFold database" is UNVERIFIED** — Kohli says it in the interview; it appears on no primary
  page. Attribute it to the interview or drop it. **1 PB is the citable figure.**
- **AVI is not just AlphaGenome + AlphaMissense.** It also folds in VEP protein loss-of-function annotations
  (termination / stop-lost / start-lost), **two conservation scores** (PhastCons 470-way, Zoonomia Cactus
  241-way) and indel indicators. My §1 understated it.
- Access: programmatic API + interactive UI. **Bulk download / UCSC track hub: UNVERIFIED either way.**

### 7.6 Benchmarks — the number moved between versions

**24 of 26** in the bioRxiv preprint → **25 of 26** in Nature 2026 (plus 22 of 24 genome-track tasks). My §1
said 25/26, correct **for Nature** — but the two must not be mixed, and the itemised 26-benchmark list is
**UNVERIFIED** (abstracts give counts only).

### 7.7 Prokaryotes — my structural conclusion is CONFIRMED, now from the FAQ verbatim

> *"The model is trained and evaluated on human and mouse DNA. Its performance on DNA from other species has
> not been determined."* … *"Prediction quality is likely to degrade as evolutionary distance from these two
> species increases."*

**hg38 and mm10 only. No prokaryotic version exists and no plan is stated** — verified as an absence across
the blog, both GitHub READMEs, the docs/FAQ, the science page and both preprints: not one mention of
bacteria, prokaryotes or cross-kingdom extension.

And the reason is stronger than "untested": the entire output vocabulary is eukaryotic — splicing, histone
modifications, chromatin contact maps. **A bacterium has no spliceosome and no nucleosomes**, so there is no
"run it at reduced accuracy" path; the heads predict quantities that do not exist in the organism.

**Net: AlphaGenome is unusable as a direct tool for the bacterial GLM, and its licence additionally forbids
using its outputs to train ours. The LOGIC in §3 is the whole of what transfers — and §3(e), the
generator-proposes / scorer-scores split, is the part DeepMind itself frames the same way.**
