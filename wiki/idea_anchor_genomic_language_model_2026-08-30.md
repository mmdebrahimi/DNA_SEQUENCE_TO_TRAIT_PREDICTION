# Idea Anchor — "the genome as a language, and a gene LLM that speaks it"

> **STATUS: DRAFT, NOT ANCHORED.** `/idea-anchor` is user-confirmed by standing directive. This is
> Soraya's drafted framing for ratification or redirect — it is not a decision, not a plan, and nothing
> downstream should treat it as anchored.
>
> Drafted 2026-08-30 against the repo's own measured record, not from first principles.

---

## 1. Formal Rephrase

**As stated:** treat genes/expression as a language and proteins as its letters; train a large language
model on that language; use it to predict how the language expresses genes, proteins and phenotype.

**Tightened:**

> Build a self-supervised sequence model over biological sequence at some chosen token level, such that
> the representations it learns support prediction of downstream biological state — molecular
> (protein function), cellular (expression, essentiality), or organismal (phenotype) — across sequences
> it has not seen.

Two things are worth separating out of the original wording, because they are different projects:

| | |
|---|---|
| **(A) A representation bet** | a model trained on sequence learns structure that transfers to phenotype prediction |
| **(B) A generative bet** | the same model can be *run backwards* to compose sequence achieving a target phenotype |

Your phrasing ("predict how this language **can be used to** express…") contains both. (B) is design; (A)
is prediction. They have different evidence, different failure modes, and in this repo **(B) already
exists and works in a bounded regime while (A) has failed five times under de-confounding.**

**On "proteins as the letters" — the wording is worth fixing, because the fix is the interesting part.**
Under the conventional analogy the letters are nucleotides and proteins are closer to *words*. But taken
literally, your version describes something genuinely different from what has failed here: a model whose
**tokens are whole genes/proteins**, and whose "text" is a genome — not a nucleotide-level model. That is
a different model class, and it is the one this repo has *not* tested. Keep the instinct; move the label.

---

## 2. Fundamental Clarifications (3, with drafted answers for ratification)

**Q1 — What is the token?** *(the fork; everything else follows from it)*
Nucleotide · codon · protein-domain · **whole gene/protein**.

> **Draft:** whole gene/protein. It is the level your wording actually points at, it is the level this
> repo's *working* decoders operate at (gene presence/absence, determinant catalogs, pathotype
> resolvers), and it is the level whose LM version has **not** been falsified here. Nucleotide-level is
> where the 0-for-5 record lives. **Risk of the alternative:** rebuilding the arm that already failed.

**Q2 — What is the training objective?**
Likelihood/masked-reconstruction (*fluency*) · or something function-supervised (*meaning*).

> **Draft:** this is the real question hiding inside the idea, and the answer must not be "likelihood".
> The repo has measured that **fluency ≠ function**, mechanistically, not just empirically (see §4).
> A likelihood objective answers "is this sequence plausible", which is not the question phenotype asks.
> **Risk of the alternative:** a sixth well-built model that reproduces the same negative.

**Q3 — What is the deliverable's regime: natural populations, or constructed variation?**

> **Draft:** constructed variation. Every de-confounded success in this repo is constructed (TEM-1 edits,
> yeast segregant cross, FBA knockouts); every de-confounded failure is a natural population. The
> discriminating variable is **population design, not organism complexity** — this is recorded as a
> corrected finding after I got it wrong three times.
> **Risk of the alternative:** the model learns population structure and reports it as signal.

---

## 3. Current Assumptions Embedded in the Idea

- **That scale is the missing ingredient.** *Measured false in this domain:* ESM2 peaks at **650M**;
  3B and 15B **regress** (median Spearman 0.484 / 0.467 / 0.438 over ProteinGym's 217 assays, reproduced
  in-house at 0.490 over 217/217).
- **That a language model's notion of "plausible" tracks biological function.** *Measured false where it
  matters most* (§4).
- **That one model can span molecular → cellular → organismal.** No evidence here supports a single model
  crossing those scales; the repo's regime map says the scales behave differently.
- **That the bottleneck is modelling.** The repo's own banked conclusion is that the binding constraint is
  **labels, not models** — and separately, that the live headroom is **modality (structure, MSA), not
  parameters**.
- **That "predict phenotype" is one task.** It is at least three: molecular effect, cellular state,
  organismal trait. Only the first currently works.
- **That this would be novel.** Gene-level genomic language models exist in the literature.
  **Unverified — check the prior art before assuming a clear field.**

---

## 4. Blunt Opinion

**The intuition is right and the proposal as worded is the version that has already been falsified here.**

Five de-confounded failures of zero-shot sequence embeddings on natural populations. Not one of them
failed for want of scale or engineering — the best-designed of them (Arabidopsis flowering time, canonical
model, real GPU, 3 seeds, n=1003) produced a **negative within-group r²** while structure-only scored
0.48. The models learned **population structure** and reported it as biology.

**The sharpest result is the one you should build the idea around, because it is a mechanism, not a
shrug.** On drug resistance, ESM2's AUROC is **0.454 — below chance** — while a curated catalog scores
0.926. Four probes located why: resistance is reached through chemically **conservative** substitutions at
ordinarily-conserved sites, so every exchangeability model calls the resistant residue benign. The control
that makes it a diagnosis: **BLOSUM62 — a 1992 substitution matrix that has never seen an HIV sequence —
makes the *same* error, slightly worse.** The blindness is a property of the **phenotype**, not of model
capacity. It is not fixable by scale, by conservation filtering, or by a better likelihood model.

**So: your idea is a training-objective question wearing a scale question's clothes.** "Fluency ≠
function" is the finding. A model that predicts the next token well has learned what evolution has
already accepted; phenotype often lives exactly where evolution has *not* been sampled. That reframing is
the valuable part of your instinct and it is genuinely open.

**Three further hard edges:**

1. **"Proteins as letters" is a level error in the standard analogy — but your literal version is the
   less-tested and more promising one.** Say "tokens are genes/proteins; a genome is the text" and the
   idea stops colliding with the failed arm.
2. **The generative half (B) already exists here and its honest limit is instructive.** `dna-decode
   inverse` proposes edits toward a target effect — and it **ranks, it does not dose**. Magnitude needs a
   calibrator fit on the target protein's own measurements, which would make the tool unnecessary.
   Expect the same wall.
3. **The one lever measured to work is modality, not parameters.** A naive rank-average of *orthogonal*
   modalities (ESM2 + GEMME + ProSST) beats ESM2 alone on **90.5%** of proteins — *even though GEMME
   alone loses to ESM2*. Sequence-only ceiling 0.458; +structure 0.507; +MSA 0.518. If you want a better
   model, that is the measured direction.

**What I would not fund:** a nucleotide-level foundation model trained on likelihood, evaluated on
natural-population phenotype. That is the 0-for-5 arm.
**What is genuinely open:** a **gene/protein-token** model over **genomic context**, trained on an
objective that is **not pure likelihood**, evaluated on **constructed variation**. Nothing in this repo
falsifies that, and the three constraints are each earned from a measurement.

---

## 5. Recommended Next Step

**Stay in conversation and answer Q1–Q3 first.** Q1 in particular is a genuine fork — nucleotide-level
and gene-token-level are different projects with different evidence, and probing the repo before that is
settled would ground the wrong one.

Then `/probe` — specifically against the closed-negative artifacts, so the design is forced to state how
it differs from each: `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md`,
`wiki/hiv_esm_vs_catalog_2026-07-09.md`, `wiki/organism_gp_regime_correction_2026-08-29.md`,
`wiki/forward_modality_hybrid_2026-07-17.md`.

No paste-ready block emitted — the next step is a decision, not a command.
