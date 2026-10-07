# AlphaGenome: what it is, what we can integrate, and what its LOGIC gives us

**Date:** 2026-10-07
**Primary source read:** the Kohli interview transcript (`youtu.be/Wkaw03p3BrM`, 28,293 chars, 777 segments,
fetched in full). Pushmeet Kohli leads this at DeepMind, so the architecture/availability claims below are
**first-party statements**, not press summaries — which is why they are quoted rather than paraphrased.
**Verification in flight:** three agents are checking the paper/architecture/licence, the bacterial MPRA
datasets, and the sequence-to-function lineage. **Dataset accessions in §5 are NOT yet verified.**

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
