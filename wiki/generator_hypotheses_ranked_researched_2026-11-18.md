# The eight generator hypotheses — RANKED on merit and researched against the literature

<!-- doc-status: LIVE as of 2026-11-18. Supersedes the RANKING in
     wiki/generator_hypotheses_2026-11-18.md (which ordered by EFFORT, not merit) and CORRECTS two
     false claims in it. The corrections and the C1/C2/C3 review findings in that file still stand. -->

## READ THIS FIRST — two claims in the previous memo are FALSE, and one hypothesis is already built

### F1 — H8 is NOT unbuilt. `dna_decode/fba/design.py` already IS it.

The previous memo said of H8 (*let a metabolic objective be the trait specification*): *"the chain from
metabolic objective through gene set to sequence edit has never been composed in this repo."* **Wrong on
both halves.**

- **`dna_decode/fba/design.py`** is growth-coupled strain design, in its own words: *"given a PRODUCT,
  which gene knockouts make producing it NECESSARY for the cell to grow?"*
- **`dna_decode/fba/compose.py`** closes the junction I called unbuilt: *"missense edit --forward--> does
  it break the enzyme? (LOF) --fba--> cell-level trait."*

**So H8 leaves the hypothesis set entirely and becomes the highest-value ACTION in it** — see Tier 0.

### F2 — H1's mechanism is not novel either; it already ships, one layer over

`design.py`'s docstring: *"`cobra.design`'s OptKnock was removed from cobrapy (absent in 0.31), so the
search here is an explicit bounded enumeration rather than a MILP: exhaustive single knockouts, then
pairs and triples over the most promising singles."*

**Bounded enumeration plus oracle ranking as a design generator is already working in this repo**, in the
metabolic layer. That is simultaneously evidence H1's mechanism WORKS and evidence it is less novel than
pitched. Combined with C3 (generators already ship in `glm/`), H1's relaxed constraint — *"generation
requires a learned generative model"* — is a strawman twice over.

### What survived from the previous memo

C1 (constraints are all CODING rules, vacuous on promoters), C2 (all three symmetries invalid here) and
C3 (generators already ship) **all stand** and are reinforced below by independent literature.

## The ranking criteria, fixed before scoring

Effort ascending was the previous order and it is **not merit** — it put the two refuted picks on top.
Scored here on five criteria, in this precedence:

1. **Not refuted** — by the 2026-11-18 review OR by published evidence.
2. **Decisive measurement is cheap, local, and needs no new labels.**
3. **Attacks the GENERATOR** (the measured bottleneck) rather than the scorer, which the project's own
   2026-10-07 result already calls adequate.
4. **Not already built** — a built thing is a run-it item, not a hypothesis.
5. **Worth OUR time** — an unpublished question beats re-deriving a published method.

## TIER 0 — not a hypothesis: a built capability that has never been run

**H8 · growth-coupled strain design.** `fba/design.py` exists, is CLI-adjacent under `dna-fba`, and has
**no artifact in `wiki/`** — it has never produced a result. `compose.py` has (`fba_composed_selector_2026-08-22.json`).

**And its implementation solves OptKnock's documented failure mode.** The literature records that
*"the predictions made by OptKnock will allow for increased production, they will not necessarily make
production growth-coupled, as alternative pathways can be used instead"*, and that *"standard OptKnock
may identify a productive growth-optimal solution while another equally growth-optimal solution produces
nothing."* `design.py`'s two-sided LP is exactly the guard: `min_flux > 0` means **every** flux
distribution that grows at that rate also secretes the product. That is the right design, written before
anyone asked.

**Honest limit from the literature: there is no systematic hit rate for this method.** The sources
document individual successes — the 2005 E. coli lactate validation gave 25–73% production gains with
measured yields (0.49, 0.45 g/g) close to predictions (0.36–0.51, 0.50 g/g) — but also that *"the in
silico and previously experimentally validated designs are only partially the same"*, and that validation
typically needs adaptive laboratory evolution, *"which implies that many designs don't pan out."* So a
run here produces **bench hypotheses**, which is exactly what `design.py`'s own scope note already says.

**Action: run it.** Cheapest real step in the whole set, no quota, and it serves the north star directly.

## TIER 1 — the three strongest hypotheses

### 1. H5 · predict a nearer-to-sequence target and compose upward — PROMOTED from months-tier

My own `why_not_trivial` said the composition model *"does not exist here"*. It exists in **both** places.

- **In the repo:** `fba/compose.py` is a working oracle-composition chain.
- **In the literature, quantitatively and parameter-free:** a synthetic constitutive-promoter library in
  E. coli found a ribosome-allocation model explains the fitness cost of expression and *"quantitatively
  predicts a linear relationship between microbial growth rate and heterologous protein expression with
  no free parameters."* A resource-allocation model builds a *"resource recruitment strength"* coefficient
  from promoter and RBS strengths. And a biophysical model maps promoter sequence to expression trained
  on random mutant libraries.

So the step this project treats as the hard unbuilt part — promoter activity → organism-level consequence
— has a published law with no free parameters. **This is also the hypothesis the repo's own R8 rail says
is systematically under-asked** ("is it my METHOD or my TARGET?"). Caveat from the same literature: when
an exogenous gene imposes a large burden *"the separate effects of promoter and RBS strength are partly
lost"*, so the linear law has a validity range.

### 2. H2 · bound every proposal by the oracle's measured decay radius

**The phenomenon is established and the exact measurement is NOT published** — the best outcome a
hypothesis can have.

- Established: *"the influence of a mutant protein against its wild type decreases as a function of its
  edit distance to the wild type."* A gradient-boosted model shows a sharp decay pattern — at K=2 it
  extrapolates three mutational regimes beyond training, at K=4 only one, and **at K=5 it fails
  completely**. Design work routinely evaluates at 5 and 10 mutations as "extrapolation distance".
- Not published: no paper found quantifies language-model accuracy as a clean function of mutational
  distance.

**The warning it must answer** comes from this repo's own prior-art synthesis: *"a regime in which no
offline method can avoid over-optimistic extrapolation."* So an offline radius may bound the damage
without eliminating it, and the honest deliverable is a decay curve plus a stated admissible radius —
not a guarantee.

### 3. H7 · cross-paradigm oracle disagreement as the exploration reward — PROMOTED

**A 2026 preprint across ten predictors and more than 13,000 proteins explicitly proposes using predictor
DISAGREEMENT to choose which proteins to assay next.** That is this hypothesis, endorsed in print. Our
three oracles are about as paradigm-distinct as available — sequence-learned, stoichiometric, and
correct-by-construction — and query-by-committee theory says disagreement estimates information gain,
with the signal strongest when members are structurally different.

**Caveat that matters here:** for regression, *"when the model class is misspecified, the active strategy
may be worse than passive learning"* — and an FBA model is explicitly approximate. Also unaddressed in
the literature: calibrating disagreement between models with different output scales and uncertainty
types, which is precisely our case.

## TIER 2 — viable once re-scoped

### 4. H1 · constraint enumeration as the generator — re-aimed at CODING targets

Refuted as pitched (C1, C3, F2), but the re-aimed form is supported: in binder design *"hallucination
methods have so far outperformed generative models... at the cost of greater compute cost"*, and
*"oracle-guided search tends to win on success rate when compute is available, while learned generative
models win on sampling efficiency and novelty."* Feasibility is the scoping fact: *"exhaustive
enumeration over a full sequence space is generally infeasible"* — but all single substitutions of a
150 bp window is 450 candidates, trivially enumerable. **Must beat `MarkovGenerator`**, which the repo
names as the null.

### 5. H6 · generate the experiment, not the edit — this is ALDE, already published

Not novel: active-learning-assisted directed evolution is published (Nature Communications 2025), uses
batch Bayesian optimization with uncertainty quantification, and reports that active learning *"can
significantly outperform the average performance of DE and random sampling."* **Adopt rather than
invent.** Two findings from it are directly actionable and both cut against the obvious design:
*"frequentist uncertainty quantification works more consistently than typical Bayesian approaches"*, and
*"incorporating deep learning does not always boost performance."* Its payoff still ends at a bench
assay, which brushes the no-new-labels non-goal.

## TIER 3 — weakened by evidence

### 6. H4 · distil on symmetry, not on predicted value

C2 is **confirmed and sharpened** by the literature's invariance-versus-equivariance distinction: labels
are *"RC invariant at the sequence level (e.g., promoter classification), or RC equivariant at the profile
level"*, and *"a gene's promoter comes from the '+' or '-' strand, a label that is by definition
dependent on orientation, so strand-specific labels need care."* Our readout is a **scalar activity for a
directional insert**, which is neither case — the reverse complement of a promoter tile is a different
construct. Methods exist (reverse-complement consistency regularisation, tested on Nucleotide
Transformer / HyenaDNA / DNABERT-2 backbones; parameter-sharing architectures that keep forward and RC
filters separate precisely to allow strand-specific prediction) but they do not rescue this target.
**The orthogonality logic still stands; it needs a symmetry that is real here, and reverse-complement is
not it.**

### 7. H3 · cross-substrate agreement as the selection filter — DEMOTED HARD

**Measured negative, from the same preprint that promotes H7:** across ten predictors and 13,000+
proteins, *"inter-VEP concordance showed no relationship with agreement to experimental MAVE data"* —
which, in that work's words, cuts against using predictor consensus as a proxy for accuracy. A summary
of the evidence: for a precision-focused confidence filter, consensus among predictors is not a reliable
signal by itself.

**This converges with the repo's own measurement** that ensembling similar predictors gave no paired
lift. **The inversion is the finding: agreement as a filter fails, disagreement as experiment-selection
works, and they are the same paper.** H3 and H7 are near-mirror hypotheses and the literature separates
them decisively.

## What changed versus the effort ordering, and why

| hypothesis | effort rank | merit rank | what moved it |
|---|---|---|---|
| H8 | 8 | **Tier 0** | already built, never run; implementation solves OptKnock's known defect |
| H5 | 5 | **1** | composition law is published AND parameter-free; `compose.py` is precedent |
| H2 | 2 | 2 | phenomenon established, exact measurement unpublished |
| H7 | 7 | **3** | published endorsement of disagreement-driven experiment choice |
| H1 | 1 | 4 | refuted as pitched; re-aimed form supported by binder-design practice |
| H6 | 6 | 5 | it is ALDE — adopt, do not invent |
| H4 | 4 | 6 | invariance/equivariance distinction confirms C2 |
| H3 | 3 | **7** | measured negative across 10 predictors / 13,000 proteins |

**Four of eight moved by two or more places.** Effort ranking had put the two refuted picks first and the
already-built capability last.

## Recommended order of work

1. **Run `fba/design.py`** — built, never executed, no quota, serves the north star. Record the artifact.
2. **Measure the oracle decay curve (H2)** — DMS data on disk, pure compute, and the result is publishable
   because nobody has plotted it for a language model.
3. **Test the published expression→growth law against our own substrate (H5)** — a parameter-free law is a
   falsifiable prediction, not a model to fit.
4. **Then** H7, then re-scoped H1. H6 as method adoption. H4 only after a real symmetry is measured.

## Honest limits of this research pass

- **Search-result level, not full-text.** These are abstracts and excerpts; no paper was read in full.
  Every quoted figure should be re-checked against the primary source before it anchors a bar — the
  2005 lactate numbers and the K=2/K=4/K=5 decay pattern especially.
- **No systematic hit rate exists for growth-coupled design**, so Tier 0's value is "produce a bench
  hypothesis", never "produce a validated strain".
- **One search direction came back thin**: nothing addresses calibrating disagreement between models with
  different output scales and uncertainty types, which is exactly H7's implementation problem.
- The ranking criteria were fixed before scoring, but **the scoring is my judgement**, not a computation.

## Sources

- [Neural network extrapolation to distant regions of the protein fitness landscape](https://www.nature.com/articles/s41467-024-50712-3)
- [Investigating the determinants of performance in machine learning for protein fitness prediction](https://www.biorxiv.org/content/10.1101/2020.09.30.319780v4.full)
- [Protein Language Model Fitness Is a Matter of Preference](https://www.biorxiv.org/content/10.1101/2024.10.03.616542v1.full)
- [Biophysics-based protein language models for protein engineering](https://www.nature.com/articles/s41592-025-02776-2)
- [Towards a Better Understanding of Reverse-Complement Equivariance for Deep Learning Models in Regulatory Genomics](https://www.biorxiv.org/content/10.1101/2020.11.04.368803v2.full)
- [Reverse-Complement Consistency for DNA Language Models](https://openreview.net/forum?id=DSBb2g26Bi)
- [Reverse-Complement Equivariant Networks for DNA Sequences](https://www.biorxiv.org/content/10.1101/2021.06.03.446953v1.full)
- [Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling](https://arxiv.org/pdf/2403.03234)
- [Fast and Ultra-Capable Protein Design: Advancing the Frontier Through Atomistic SE(3)-Equivariance with Genie 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC13174375/)
- [Systematic comparison of Generative AI-Protein Models](https://www.biorxiv.org/content/10.1101/2025.03.23.644844.full.pdf)
- [Active learning-assisted directed evolution](https://www.nature.com/articles/s41467-025-55987-8)
- [Disagreement among variant effect predictors guides experimental prioritization of target proteins](https://www.biorxiv.org/content/10.64898/2026.03.18.712765v1)
- [Why variant effect predictors and multiplexed assays agree and disagree](https://pmc.ncbi.nlm.nih.gov/articles/PMC13620131)
- [In Silico Constraint-Based Strain Optimization Methods: the Quest for Optimal Cell Factories](https://pmc.ncbi.nlm.nih.gov/articles/PMC4711187/)
- [OptCouple: Joint simulation of gene knockouts, insertions and medium modifications](https://www.sciencedirect.com/science/article/pii/S2214030118300373)
- [The Interrelationship between Promoter Strength, Gene Expression, and Growth Rate](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0109105)
- [RBS and Promoter Strengths Determine the Cell-Growth-Dependent Protein Mass Fractions](https://pubs.acs.org/doi/10.1021/acssynbio.1c00131)
- [Predicting bacterial promoter function and evolution from random sequences](https://elifesciences.org/articles/64543)
- [Active Learning for Regression based on Query by Committee](https://dl.ifip.org/db/conf/ideal/ideal2007/BurbidgeRK07.pdf)
