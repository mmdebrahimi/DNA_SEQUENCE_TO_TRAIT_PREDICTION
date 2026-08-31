# Prior-art deep dive: is the "conservative-substitution / BLOSUM62 control" finding novel?

The prior-art check on 2026-08-31 flagged one possible original micro-finding and explicitly refused to
claim it: *"resistance via chemically conservative substitutions at averagely-conserved sites, with
BLOSUM62 reproducing the error (which rules out memorisation) — NOT established, search absence is not
proof."* This is the targeted check.

**Verdict: the finding is a COMPOSITE. Three of its four components are prior art. One appears
unreported. And one of our claims needs a caveat it did not have — a published counterexample exists.**

---

## What we actually claim (pinned from the artifacts before searching)

From `wiki/hiv_esm_vs_catalog_2026-07-09.md` + `scripts/resistance_conservativeness_probe.py`:

| probe | claim | measured |
|---|---|---|
| **A** positional tolerance | are DRM *sites* just unconserved? | **No** — entropy percentile **0.494**, exactly average |
| **B** mutant specificity | does ESM rank the *resistant* residue as likely? | **Yes** — median **4.5 / 19** (null 10); 55% in top-5 |
| **D** BLOSUM62 control | is B ESM-specific, or generic exchangeability? | **Not ESM-specific** — BLOSUM62 ranks them **4.0 / 19**, *better*; ESM wins only 12/38 |
| generalisation | does it hold beyond HIV? | HIV RT **4.5** (n=38) · SARS-CoV-2 Mpro **6.5** (n=84) · fungal ERG11/FKS1 **4.0** (n=17) · **pooled 4.5** (n=139) |
| our own caveat | is it selection, or the genetic code? | codon-accessible re-rank **2.5** vs accessible-null **3.5**, P=0.614, n=22 — **underpowered** |

---

## PRIOR ART — three components are already published

### 1. "Resistance mutations are biochemically conservative" — PUBLISHED, **and contradicted elsewhere**

**Friedman R. (2013), *PLOS ONE* 8:e82059** — *Drug Resistance Missense Mutations in Cancer Are Subject
to Evolutionary Constraints.* Analysing tyrosine-kinase-inhibitor resistance:

- Of 43 Abl1 resistance mutations, **none was novel** — the same substitution type appears in
  evolutionarily related proteins at that position.
- **"Abl1 mutations tend to be biochemically conservative, whereas EGFR and ALK mutations tend to be
  radical."**
- For Abl1, in all but two cases the mutant residue is **less conserved** in the conserved-domain
  database.

**This is the most important finding of the whole check, and it cuts both ways.** The conservativeness
observation is not ours. **And Friedman supplies an explicit counterexample class: EGFR and ALK resistance
mutations are RADICAL.** So "resistance mutations are chemically conservative" is *not* a general law, and
our write-ups should stop implying one.

What survives as ours: the observation holds across **three unrelated pathogen target-site systems**
(HIV RT, SARS-CoV-2 Mpro, fungal ERG11/FKS1, pooled n=139) measured with one instrument. That is an
extension into a different domain, not a discovery — and it now has a named boundary.

*Note a mechanistic divergence:* Friedman's framing is **evolutionary accessibility** (the mutant residue
occurs in homologues at that position, and is usually *less* conserved). Ours is **chemical similarity**
at a site of *average* conservation. Related, not identical.

### 2. Drug-agnostic zero-shot PLM baselines on resistance — PUBLISHED, in our exact framing

*Contextualizing Biological Language Models across Modalities via Logit-Space Contrastive Alignment*
(arXiv 2606.18703) runs **ESM-1v, ESM-2, EVE and Tranception** zero-shot via mutant-vs-wild-type
pseudo-log-likelihood as **drug-agnostic** baselines, explicitly to test *"whether general protein fitness
or evolutionary plausibility alone explains drug-resistance effects."* That is our experimental framing,
already in print.

*Caveat: I retrieved their SETUP, not their numbers. Whether they report a below-chance result comparable
to our AUROC 0.454 is unverified.*

### 3. BLOSUM62 as a baseline — LONG-STANDING

**VESPA** (Marquet et al.) defines *"BLOSUM62bin"* as an explicit naïve baseline (negative BLOSUM62 ⇒
"effect") and also uses BLOSUM62 as a **feature** alongside conservation and pLM reconstruction
probabilities. ESM1b's genome-wide paper uses BlastP (BLOSUM62-scored) as a baseline. A 2026 preprint
extends the pattern — PSSM baselines are statistically indistinguishable from pLM zero-shot scoring.

So "BLOSUM62 as a comparator" is standard. **The specific USE is what differs — see below.**

---

## APPARENTLY NOT PRIOR ART — what survives

### A. BLOSUM62 as a **memorisation control**
Searched directly. **Not found.** The nearest hit inverts the logic: a generative-protein-model paper uses
BLOSUM62-weighted *"soft accuracy"* to argue its model learned **function-preserving invariances** — the
matrix as evidence *for* learned biology, not as a control ruling *out* memorisation.

Our inferential move — *BLOSUM62 has never seen an HIV sequence, it makes the same error slightly worse,
therefore the PLM's behaviour is generic exchangeability rather than memorised circulating variants* —
appears unreported.

### B. The refutation of the standard explanation — the sharpest piece
The **published** mechanism for variant-effect-predictor false negatives is that *"disease-causing alleles
residing in poorly or non-conserved regions will be false-negatively classified as neutral."*

**We measured that and rejected it.** DRM positions sit at entropy percentile **0.494** — exactly average
conservation. So in this system the textbook explanation is *wrong*, and the mechanism is the **chemistry
of the substitution**, not the **conservation of the site**.

Refuting the standard explanation with a measurement is a real, specific, falsifiable contribution.

### C. The *direction* of the failure is different from the published one
The literature's framing (e.g. *From Likelihood to Fitness*) is that pLM scores are effectively **unsigned
deviation magnitude** — an activating mutation and a destabilising one both score "unlikely", so the model
cannot tell direction.

**Ours is a different failure.** The model rates the resistant residue **LIKELY** — rank 4.5/19, 55% in
its top-5. It is not confused about direction; it actively calls the mutation *normal*. Distinct mechanism,
distinct remedy.

### D. A model-free, cross-pathogen instrument
BLOSUM-only, no model, no GPU, seconds, three unrelated pathogens, n=139, mid-ranks throughout. Not found
in this form.

---

## What this check COSTS us, honestly

1. **Drop any implication that "resistance mutations are conservative" is general.** Friedman's EGFR/ALK
   radical counterexample must be cited wherever we state it. Our claim is scoped to pathogen target-site
   resistance in three systems.
2. **Stop implying the drug-agnostic zero-shot framing is ours.** It is in print.
3. **The deep mechanism is classical, not ours.** Our own artifact already concedes codon accessibility
   does most of the work (2.5 vs 3.5, P=0.614, n=22, underpowered) — and the **genetic-code
   error-minimisation** literature (single-nt neighbours are chemically similar *by construction*) is
   decades old and robust across physicochemical metrics. We correctly identified the confound; we did not
   discover it.

## The honest one-sentence version

> Across three unrelated pathogen systems, likelihood-based scorers rate resistance substitutions as
> *likely* rather than merely undecidable, and this is **not** explained by the standard "the site is
> poorly conserved" account (entropy percentile 0.494) — a model-free BLOSUM62 control reproduces the same
> error, indicating generic amino-acid exchangeability rather than memorisation, with codon accessibility
> a large and unresolved share of the cause.

Everything in that sentence is measured. It is a **diagnosis assembled from known parts**, with two
genuinely unreported moves (the memorisation control; the refutation of the conservation explanation).

## Coverage limitations of this check

- **Two searches were blocked** by content safeguards on ordinary bioinformatics phrasing (conservation at
  resistance positions; Grantham-distance analysis of resistance substitutions). Those axes are
  **under-searched**, and a Grantham/Miyata-distance analysis of resistance substitutions may well exist.
- Friedman 2013's full text was not retrievable (cookie wall); findings are from the abstract and
  secondary reporting.
- The arXiv drug-agnostic-baseline paper's *results* were not retrieved, only its method.
- Search absence remains not proof of absence.

## Sources

- [Friedman R. (2013) Drug Resistance Missense Mutations in Cancer Are Subject to Evolutionary Constraints — *PLOS ONE* 8:e82059](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0082059)
- [Contextualizing Biological Language Models across Modalities via Logit-Space Contrastive Alignment — arXiv 2606.18703](https://arxiv.org/pdf/2606.18703)
- [Embeddings from protein language models predict conservation and variant effects (VESPA)](https://pubmed.ncbi.nlm.nih.gov/34967936/)
- [Genome-wide prediction of disease variant effects with a deep protein language model (ESM1b) — *Nature Genetics* 2023](https://www.nature.com/articles/s41588-023-01465-0)
- [Simple baselines rival protein language models in mutation-dense design of function tasks — bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.05.01.722313v2.full)
- [From Likelihood to Fitness — bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.05.20.655154v1)
- [How good are pathogenicity predictors in detecting benign variants? — *PLOS Comp Biol*](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1006481)
- [Updated benchmarking of variant effect predictors using deep mutational scanning](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10407742/)
- [Exceptional error minimization in putative primordial genetic codes — *Biology Direct*](https://biologydirect.biomedcentral.com/articles/10.1186/1745-6150-4-44)
- [Substitution scoring matrices for proteins — an overview, *Protein Science* 2020](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3954)
