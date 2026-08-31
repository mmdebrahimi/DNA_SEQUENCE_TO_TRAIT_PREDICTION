# Methods worth borrowing, given that likelihood scoring fails on resistance

Follow-on to `wiki/prior_art_conservative_resistance_blosum_2026-08-31.md`. That memo established the
mechanism; this one asks what to do about it. **Everything below is a literature CANDIDATE — nothing here
has been tested in this repo.** Screen each against gates G1–G10 in `wiki/negative_results_map_2026-06-13.md`
before spending effort.

## The mechanism, restated as a design constraint

Resistance substitutions are **chemically conservative** at sites of **average conservation** (entropy
percentile 0.494). Any scorer whose question is *"how surprising is this sequence?"* answers *"not very"*
and misses them. BLOSUM62 — no model, no training — reproduces the same error.

**So the requirement is: a signal that is not a function of sequence plausibility.** Four such signals
exist in the literature, and they differ in how well they fit.

---

## Ranked candidates

### 1. ΔΔG of **binding** (not stability) — Rosetta `flex_ddG`, FoldX ★★★ best fit

Resistance is *reduced drug binding*. That is a **physical quantity a likelihood model never computes** —
maximally orthogonal to the failing signal.

- **Aldeghi et al., *ACS Central Science* 2018** — 134 protein–ligand mutations; combined Rosetta + MD
  reached **RMSE 1.2 kcal/mol** (0.8 on well-behaved systems). The authors explicitly frame it as
  *"opening the prospect of ... predicting drug resistance."*
- `flex_ddG`, built for protein–**protein** interfaces, was found to work for protein–**ligand** too. Its
  trick is backbone sampling ("backrub" ensembles) — earlier alanine-scanning methods froze the backbone,
  which is fine for →Ala but not for the larger side chains resistance actually uses.
- **StaB-ddG** (2025) claims FoldX-level accuracy at ~1000× speed.

**Why this is the strongest candidate for us specifically:** the catalog's measured product gap is
*novel/unseen variants* (53 resistant isolates carrying no catalog DRM). ΔΔG_bind scores **any**
substitution at a target site with no catalog entry required. And it is **a different modality** —
physics, not evolution — which is precisely the condition under which this repo has already measured
ensembling to pay (`ESM2+GEMME+ProSST` +0.056, wins 90.5%).

**Costs:** Rosetta is free academically but CPU-heavy; needs a structure of the target–drug complex.
Benchmark caveat found: ΔΔG benchmark sets over-represent →Ala and under-represent charge↔charge.

### 2. Positive / diversifying selection scans — HyPhy `MEME`, `FUBAR`, `FEL` ★★★ conceptually exact

Likelihood models capture **purifying** selection (what nature has kept). Resistance sites are under
**positive/diversifying** selection (what nature is actively exploring under drug pressure). **These are
opposite signals** — which is exactly why a site can sit at average conservation and still be a resistance
hotspot.

- `MEME` detects **episodic** positive selection at individual sites — the only HyPhy method covering both
  pervasive and episodic; recommended when power matters.
- `FUBAR` is the fast Bayesian option for large datasets.
- Free, open-source, scriptable headlessly: `hyphy meme --alignment data.fas --tree tree.nwk`. **No GPU.**

**Two caveats found, and the second is disqualifying for part of our surface:**
- Most selection-scan literature is **retrospective** (detecting past fixations), not predictive.
- **Codon-level scans systematically miss HGT-mediated resistance.** Most bacterial AMR in this repo is
  *acquired gene* resistance (blaCTX-M, rmt, sul/dfr), not target-site point mutation. So this applies to
  **HIV RT/PR/IN, SARS-CoV-2 Mpro, fungal ERG11/FKS1, TB rpoB/katG — and NOT to the acquired-gene AMR
  cells.**

### 3. Drug-**conditioned** models — ConPLex, MolTrans, Boltz-2 ★★

This is the literal build of the "missing environment input" claim, and it already exists.

**ConPLex** (Singh, Sledzieski, Bryson, Cowen & Berger, *PNAS* 2023) — protein-anchored **contrastive
co-embedding** placing proteins and drugs in a shared space, trained to contrast true binders from decoys.
Validated wet-lab: 12 of 19 predicted kinase–drug interactions confirmed, including a 1.3 nM binder.
Code: `github.com/samsledje/ConPLex`.

**Adaptation gap:** built for *"does this drug bind this protein"*, not *"does this mutation reduce
binding"*. Using it for resistance means scoring wild-type vs mutant and taking a difference — plausible,
untested.

### 4. Epistasis-aware Potts / DCA — EVcouplings, plmDCA, bmDCA ★★ with a caveat

Direct-coupling models add pairwise terms `J_ij` on top of per-site fields, capturing that a mutation's
effect depends on background. Two findings make this relevant to us specifically:

- **HIV drug-target proteins show strong long-range epistasis**, captured by a Potts covariation model and
  **not** by an independent-site model.
- **TEM-1 β-lactamase context-dependence** (Figliuzzi et al. 2016) — our own working molecular cell.

**The caveat that lowers its rank:** a Potts model is still a *likelihood* model. It may inherit the same
"looks normal" failure, just with couplings. **It fixes independence, not the plausibility framing.**

### 5. Inverse folding conditioned on the **complex** — ESM-IF, ProteinMPNN, Boltzmann-Aligned IF ★★

Useful nuance found: *"conditioning on the full complex and using SKEMPI-style ΔΔG_bind calibration
outperforms naive native-vs-mutant likelihood ratios."*

**But the ProteinGym stratification argues against it here:** *"sequence-only models tend to outperform
structure-incorporating models on most properties **except stability**."* Inverse folding is a **stability**
specialist, and resistance is **binding**. Rank 5, not 1 — go to explicit ΔΔG_bind (#1) instead.

---

## Two datasets worth acquiring, not just methods

**PEAR (β-lactamase, prospective/retrospective split)** — ~23,000 *E. coli* strains, each carrying a
**unique single-copy variant** of `blaCTX-M-14`, with relative growth measured under cefotaxime or
ceftazidime. Separate `PEAR^P` (prospective) and `PEAR^R` (retrospective) models; PEAR^P-predicted
mutations were **enriched among clinical isolates**.

**This is constructed variation at scale, on an AMR target, with measured phenotype** — the exact regime
this repo has established as the one that works, on the exact drug class where our catalog operates.
Highest-value acquisition target found in this search.

**SARS-CoV-2 prospective forecasting** (*Science Translational Medicine*) — validated against Omicron
**before emergence**, with scores rising during emergence. Combines positive-selection signal with
frequency dynamics. Directly relevant to our **prospective-lock arm**, which currently has a harness and
almost no accrued data.

---

## What I would actually do first

**A single orthogonality test, not a build.** On the HIV NNRTI set where we already have measured
PhenoSense fold-change and a catalog at AUC 0.962:

1. compute ΔΔG_bind for the catalog-**negative** subset (the 53 resistant isolates the catalog misses);
2. score against the same labels;
3. report against the catalog and against ESM2 — with the repo's usual rails (paired deltas, win counts,
   not differenced medians).

If physics recovers catalog-negative cases that likelihood cannot, that is a **real, orthogonal,
deployable complement** to the frozen decoder — and it targets a *measured* product gap rather than a
hypothesised one. If it fails, it fails cheaply and the negative is publishable in the repo's own idiom.

## Honest limits of this search

- **Three searches were blocked** by content safeguards on ordinary bioinformatics phrasing (conservation
  at resistance positions; Grantham-distance analysis; mutation-frequency-trajectory surveillance). Those
  axes are **under-searched**.
- No candidate here has been tested by us. Method-fit rankings are **my judgement against the measured
  mechanism**, not empirical results.
- ΔΔG_bind requires a structure of the target–drug complex; availability per cell is unchecked.

## Sources

- [Aldeghi et al., Accurate Estimation of Ligand Binding Affinity Changes upon Protein Mutation — *ACS Cent. Sci.* 2018](https://pubs.acs.org/acscii/article/4/12/1708/733730/Accurate-Estimation-of-Ligand-Binding-Affinity)
- [Flex ddG — *J. Phys. Chem. B* 2018](https://pubs.acs.org/doi/10.1021/acs.jpcb.7b11367)
- [HyPhy selection-method tutorials (MEME / FEL / FUBAR / BUSTED)](https://github.com/veg/hyphy-tutorials/blob/master/docs/selection/README.md) · [hyphy.org](http://hyphy.org/)
- [ConPLex — *PNAS* 2023](https://www.pnas.org/doi/10.1073/pnas.2220778120) · [code](https://github.com/samsledje/ConPLex)
- [EVmutation / GEMME background — GEMME, *MBE* 2019](https://academic.oup.com/mbe/article/36/11/2604/5548199)
- [Potts Hamiltonian models of protein co-variation — *Curr Opin Struct Biol* 2016](https://www.sciencedirect.com/science/article/abs/pii/S0959440X16301841)
- [Boltzmann-Aligned Inverse Folding for PPI mutational effects — arXiv 2410.09543](https://arxiv.org/pdf/2410.09543)
- [Predicting the mutational drivers of future SARS-CoV-2 variants of concern — *Sci. Transl. Med.*](https://www.science.org/doi/10.1126/scitranslmed.abk3445)
- [Prediction of Antibiotic Resistance Evolution by Growth Measurement of All Proximal Mutants of Beta-Lactamase (PEAR)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9087888/)
- [A computational method for predicting stepwise accumulation of resistance mutations (flex-ddG-parameterised trajectories)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10807863/)
