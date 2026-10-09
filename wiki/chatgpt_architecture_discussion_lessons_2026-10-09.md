# Lessons from the ChatGPT architecture discussion, screened against what this repo has MEASURED

**Date:** 2026-10-09 · **Source archived verbatim:**
`wiki/refs/chatgpt_alphagenome_architecture_discussion_2026-10-09.md` · **Screened against:**
`dna_decode/eval/regime.py`, the GLM umbrella, and the committed artifacts named inline.

This is the screen, not the source. The source is a **secondary** source (a model's proposal), and its
sibling summary already contained two factually wrong citations — so every item below is marked by what
*this repo has measured*, not by how persuasive the proposal is.

---

## 1. The discussion IS the plan's origin — confirmed from the ledgers, verbatim

| ChatGPT priority | Became | Originating goal (verbatim from the ledger) |
|---|---|---|
| **#1** multi-resolution + hierarchical tokens | **G-A** | *"Replace the expression oracle's fixed sequence features … with a LEARNED multi-resolution encoder plus a trivial head, and measure whether it transfers to real genomic sequence."* |
| **#3** (2nd half) cell-state conditioning | **G-C** | *"Build the bacterial analogue of **AlphaGenome's cell-state conditioning**: one oracle conditioned on growth medium, tested against two separate per-medium oracles."* |
| **#2** regulatory graph + sparse distal attention | — | **never built** |

So two of the three stated priorities were turned into project families, and **both reported NEGATIVE on
2026-10-08**. That is the headline — and the next section is why it does not mean what it looks like.

---

## 2. BOTH negatives tested the WEAKEST FORM of the proposal. Read this before re-proposing either.

This is the most important thing in this document, in both directions: it is the honest caveat on the
negatives **and** the reason not to simply re-run them bigger.

| | ChatGPT proposed | What we actually tested | Gap |
|---|---|---|---|
| **#1** | a pyramid **1 bp → 16 bp → 128 bp → 1 kb → 10 kb** with information travelling between levels, over a **1 Mb** context | a **13,217-parameter** conv at 3 kernel widths (6/20/75) over a **150 bp** tile | our substrate spans **~1 decade**; the proposal spans **four**. There is no 1 kb or 10 kb level to build on a 150 bp sequence. |
| **#3** | a **cell-state embedding modulating a learned representation** | a **linear interaction term** (`seq × medium`) over 3-mer frequencies | the artifact's own honest limit says it: *"not the same as a model whose **representation** is modulated by condition. A null here does not foreclose that."* |

**But the substrate bounds how much any of this can buy, and that is also measured.** On the same tiles, a
**single GC number** reaches **0.3112** against a replicate-derived ceiling of **0.941** — i.e. **33% of
achievable** — and the best learned arm reached **34%**. Capacity has already *hurt* here: 6-mers (4,096
features) scored **worst of every arm at 0.2022**
(`wiki/glm_genomewide_oracle_2026-10-07.json`).

**So the synthesis is not "re-run it bigger":** most of ChatGPT's architecture proposals are about
**extending range** and **adding hierarchy**, and our substrate is 150 bp with very little hierarchy to
find. Re-running a larger multi-resolution encoder on *this* substrate has a small ceiling. The honest open
question is whether the **substrate** is the binding constraint, not the architecture.

---

## 3. The sixteen proposals, screened

Five verdict classes. The two that matter most are **ALREADY SHIPPED** (we are further along than the
proposal assumes) and **MEASURED NEGATIVE** (do not re-propose without new information).

### STRUCTURALLY INAPPLICABLE to bacteria — ~6 of 16

Not "hard", **absent**. A bacterium has no spliceosome, no nucleosomes, no TADs, no CTCF, one chromosome
copy, and no cell types in the eukaryotic sense.

| # | Proposal | Why it cannot apply |
|---|---|---|
| 3 | hierarchical attention through **TAD / chromatin domain** levels | no TADs; the L4–L5 levels do not exist |
| 4 | sparse **distal** attention over 1–10 Mb | bacterial cis-regulation is overwhelmingly promoter-**proximal**; a 300 kb enhancer–promoter bridge is not a bacterial mechanism |
| 5 | regulatory graph with **enhancer / CTCF** nodes | no enhancers in that sense, no CTCF |
| 6 | 3D structure as latent variable | nucleoid supercoiling is real and does affect transcription, but there is no Hi-C-scale TAD analogue to bottleneck through |
| 10 | diploid haplotypes | bacteria are **haploid**. Relevant only to this repo's **human clinical track** (pgx / clinvar / hla) |
| 16 | mixture-of-experts over splicing / neuronal / developmental experts | those domains do not exist here; and at 13k parameters MoE is premature by orders of magnitude |

**Reusable:** roughly a third of a eukaryote-derived architecture proposal is inapplicable to a prokaryote
*by biology*, not by effort. Screening first is how you avoid spending a month on an absent mechanism.

### ALREADY SHIPPED here — 3 of 16, and this is the strongest convergence

| # | Proposal | Status in this repo |
|---|---|---|
| **8** | use **perturbation** data, because observational tracks give correlation not causation | **This is already the project's measured winning regime.** `eval/regime.py`: natural-population zero-shot is a **closed negative** (0-for-5, de-confounded), while **constructed** variation **WORKS** — DMS ρ 0.35–0.76, FBA/Keio MCC 0.70–0.74. Our substrate *is* an MPRA (GSE144621). ChatGPT reasoned to this from theory; we measured it. |
| **9** | train on `ref + alt → Δ`, making variant effect first-class | **Shipped** as `dna-forward` with a genome-edit path (ρ **0.761** on TEM-1 ampicillin fitness, 1,715 variants). Independently, the project's own R8 rail says the same thing: *a differential readout on the same background cancels a confound present in both arms.* |
| **15** | add **uncertainty**; the model should be able to say "I have not seen enough biology like this" | **Shipped** as the **L2 doubt layer** (`dna_decode/eval/doubt.py`) — and it is unusual: the project's own survey notes most of the field, including Mykrobe, *predicts SUSCEPTIBLE when no known mutation is found*. The doubt block structurally **cannot contain a call** (`assert_no_call` runs on its own output). |

**So the three proposals a model reasons to from first principles are the three this project had already
built and validated.** That is real external corroboration of the architecture choices — worth more than
any of the untried suggestions.

### ALREADY MEASURED NEGATIVE here — 4 of 16

Do not re-propose these without new information; each has a committed artifact.

| # | Proposal | What we measured |
|---|---|---|
| **12** | self-supervised DNA pretraining, **then** functional supervision | The pretrained half is **measured worse than trivial at zero-shot generation**: GENERator-1.2B-prokaryote scores distinguishability **0.7706** against a 3-mer Markov chain's **0.6679** (lower is better; `wiki/glm_generator_falsifier_REAL_2026-10-07.json` — note the `_REAL_` suffix, a different artifact from the same-named non-REAL one). A 1.2B model more detectable than a Markov chain. **The "then functional supervision" half is exactly G-D and is the un-pulled lever** — so this proposal is half-refuted, half-untested. |
| **13** | feed **evolutionary conservation** explicitly | **Measured structurally blind for our phenotype.** Resistance-mutation sites sit at ESM entropy percentile **0.494 — exactly average**, refuting the textbook "the site is poorly conserved" account; and **BLOSUM62, which has never seen an HIV sequence, ranks the same DRMs as well as ESM2** (4.0 vs 4.5 of 19). Johnston 2024 (*PNAS*, TrpB) independently: the most-fit variants carry a substitution *"nearly absent in natural sequences — a result that conservation-based predictions would not capture."* **Conservation-keyed methods are blind to the best variants.** |
| **11** | model variant **combinations** (epistasis) | **Measured to lose.** On HIV RT, pairwise and non-linear models both lost to a linear one-hot model **4/4, two with CI below zero** — the linear model is the measured **ceiling, not a floor**. Corroborated by an independent 2026-07-11 epistasis negative. *Scope: one virus, one protein family.* |
| **1** | multi-resolution tokens | **Tested as G-A → `WEAK`, not distinguishable from zero** (6/10 seeds, p=0.754, CI [−0.0099, +0.0349]). But see §2 — this tested a 1-decade version of a 4-decade proposal. |

### GENUINELY OPEN and worth the conversation — 3 of 16

| # | Proposal | Why it is live |
|---|---|---|
| **7** | cell-state conditioning **of the representation** | G-C's negative explicitly does **not** foreclose it, in the artifact's own words. The bacterial analogue of cell state is growth condition, and the substrate genuinely carries the signal: noise-matched cross-condition agreement **0.4063** vs within-condition **0.4677**, disattenuated true cross-condition **0.7887** — **~38% of real variance is condition-specific**. Something is there that a linear interaction term did not capture. |
| **2** | **learned biological tokens** rather than spatial windows | The bacterial vocabulary is different but real: promoter, operator, RBS, terminator, operon, regulon, σ-factor binding site. Discovering those as tokens is not an absent mechanism — it is an untried one. |
| **14** | mechanistic **intermediate** supervision | Partly present already: the repo has an FBA metabolic layer that *is* a mechanistic intermediate, and `constraints/` encodes axiomatic laws. The proposal's version — latent intermediates supervised wherever data exists — is untried. |

---

## 4. The two sharpest transferable ideas, and neither is an architecture

**(a) "Attention weights are not biological causality."** Quoted from the source, and it generalises far
past transformers: *a high attention score means the model found this region useful when constructing the
representation* — not that the enhancer controls the promoter. This repo has the same failure mode under a
different name: a feature's importance is not a mechanism, and the project's gates G1–G10 exist because
tool-derived and sampling-defined labels look like measurements. **Add it to the standing vocabulary: an
internal weight is evidence about the MODEL, never about the biology.**

**(b) "I would not build your model as a much bigger AlphaGenome."** The strongest framing in the
discussion, and it matches this repo's own north-star seam exactly. AlphaGenome stops at
`DNA → molecular phenotype`; the goal is `genotype → phenotype`, which needs protein effects → cellular
state → tissue → physiology → development → environment. The authors themselves say molecular predictions
do not capture complex traits because **development and environment lie outside scope**.

This is the same thing `eval/regime.py` encodes and the same thing lesson 6 ("prerequisite, not
sufficient") records: **the expression oracle predicts expression, which is a PREREQUISITE for a trait.**
ρ 0.59 and 79%-of-selection-value are not trait claims. The discussion arrives at the project's own honesty
rail from the outside, which is the best possible corroboration of it.

---

## 5. What this changes about the plan

1. **G-D (fine-tune the generator) remains #1, and this discussion strengthens it.** Proposal #12's
   sequencing is *pretraining → functional supervision → perturbational supervision*. We have measured the
   pretraining half failing at zero-shot and have never run the functional-supervision half. That is
   exactly G-D.
2. **Do NOT re-run G-A bigger on the same substrate.** The range the proposal needs (1 kb–10 kb levels)
   does not exist in a 150 bp tile, and GC already takes 33% of a 0.941 ceiling. If G-A is revisited, the
   substrate must change, not the parameter count.
3. **G-C's live successor is representation-level conditioning**, not a wider pooling grid — and the
   measured 38% condition-specific variance says the target is real.
4. **Nothing here unblocks G-B.** Self-distillation is blocked on a learned representation existing, and
   proposal #9's "train on ref+alt Δ" is already shipped by a different route.
5. **New candidate worth a plan: learned biological tokens** (proposal #2) in the bacterial vocabulary —
   the one structurally-applicable, never-tried item in the top of the list.

## 6. Honest limits on this screen

- The source is a **model's proposal**, not a paper. Its sibling summary had two wrong citations. Nothing
  above treats its claims about AlphaGenome's internals as verified — only the *ideas* are screened.
- The "structurally inapplicable" verdicts rest on standard prokaryotic biology, not on a measurement in
  this repo. They are high-confidence but they are **asserted**, not derived here.
- The epistasis negative (#11) is scoped to **HIV RT**; it is not a general claim about epistasis.
- The conservation negative (#13) is scoped to **pathogen target-site resistance** — Friedman 2013
  supplies a counterexample class (EGFR/ALK resistance mutations tend to be *radical*), so it must not be
  stated as a general law.
- §2's "33% of achievable" is specific to **peak tiles at 150 bp**; it is not a statement about bacterial
  promoters in general.
