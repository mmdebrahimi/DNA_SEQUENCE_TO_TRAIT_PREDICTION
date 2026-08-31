# Prior-art check: the field is NOT clear, and it has independently reached every one of this repo's findings

The draft idea-anchor (`wiki/idea_anchor_genomic_language_model_2026-08-30.md`) flagged one claim as
**unverified**: *"gene-level genomic language models exist in the literature — check the prior art before
assuming a clear field."* This is that check. Web research, 2026-08-31.

**Verdict: the field is crowded, the specific model class already exists with open weights, and four
separate findings this repo derived independently are already published.** The anchor needs revision, and
the honest remaining niche is much narrower — but it is still real, and it is not where I put it.

---

## 1. The user's exact idea is built, published, and downloadable

| model | tokens | scale | what it did |
|---|---|---|---|
| **gLM** (Hwang et al., *Nat Comms* 2024) | **whole genes**, as ESM2 protein embeddings | 19 layers, contigs ≤30 genes | learned operons + co-regulation; contextualised protein embeddings |
| **gLM2 / OMG** (Cornman et al., Tatta Bio, 2024) | **mixed** — amino acids for coding, nucleotides for intergenic | 150M + **650M**, 3.1 Tbp / 3.3B proteins | beats ESM2 on most protein tasks; learns protein–protein interfaces via coevolution |
| **GenSyntax** (2025) | gene **product descriptors** as text | 49,250 RefSeq genomes | plasmid host ID, gene function, genome assembly, **essentiality** |
| **Geneformer / scGPT / scBERT / scFoundation / UCE** | **genes**, ranked or binned by expression | 30M–95M cells | cell type, perturbation response, network inference |

**gLM is the anchor's proposal, published two years ago.** Genes as tokens, protein embeddings as the
token vectors, masked prediction over genomic context. **gLM2 goes further** and is on HuggingFace
(`tattabio/gLM2_650M`) — so the "build it" step is partly a download.

Two details worth noting: gLM2's best checkpoint is **650M**, the same size this repo independently
measured as the sweet spot; and its headline win is **protein-protein interfaces**, which Evo and ESM2
miss — a *molecular* result, not a phenotypic one.

**The gap I could not close:** no paper surfaced using **orthogroup IDs** (OrthoFinder/eggNOG) directly as
a discrete vocabulary. GenSyntax's product descriptors and gLM's embeddings bracket it on either side.
*Absence from search is not proof of absence.*

---

## 2. Four of this repo's own findings are already in the literature

This is the uncomfortable half, and it is worth stating plainly.

### 2a. "Scale is dead — 650M peaks" — CONFIRMED, published
An iScience benchmark finds **300M–650M optimal** for variant-effect prediction, with larger ESM2 models
regressing — attributed to overfitting that inflates the reference amino acid's probability under
masked-margin scoring. An independent pathogenicity benchmark reports the same non-monotonicity
(ESM-1 AUCs 0.725 → 0.769 → 0.874 → **0.856**).

**Our ProteinGym re-derivation (0.484 / 0.467 / 0.438) reproduces a published result.** Good news for our
method; not a novel finding.

### 2b. "Fluency ≠ function" — TRUE, NAMED, and someone has a fix
This is the reframing I called the valuable half of the idea. **It is an active research program:**

- **Hou et al., *Nature Computational Science* (2026)** — a "Goldilocks" effect: model size, training set
  and stochastic factors bias `p(sequence)` away from fitness; at the extremes likelihoods become
  *uninformative*.
- **Gordon, Lu & Abbeel, ICLR 2025 — "Protein Language Model Fitness is a Matter of Preference"** — the
  implicit pretraining *preference* for a sequence predicts fitness-prediction capability; both
  over-preferred and under-preferred wild-types hurt. Causally linked to training data via influence
  functions.
- **Pugh et al., bioRxiv 2025 — "From Likelihood to Fitness"** — states the diagnosis almost exactly as
  the anchor does ("likelihood also reflects phylogenetic structure and sampling biases, especially as
  model capacity increases") and ships **Likelihood-Fitness Bridging**, a no-retraining inference-time
  correction. **Public code: `github.com/DiasFrazerGroup/lfb`.**

**Consequence for the anchor:** the insight stands and is correct. It is **not** ours, and there is an
off-the-shelf mitigation to beat before claiming anything.

### 2c. "Don't train on natural populations" — CONFIRMED at scale, in our exact domain
**PLOS Biology 2025 — "Biased sampling driven by bacterial population structure confounds machine
learning prediction of antimicrobial resistance."** 24,000+ genomes, five pathogens, **6,740 models**.
Findings: models perform poorly when resistance is confounded with phylogeny; **increasing training
sample size fails to rescue performance**; little overlap in predictive features across clades; most
predictions depend on phylogenetic background.

**That is our 0-for-5, published, at 100× our scale, on our organisms.** The "more data won't fix it"
clause is stronger than anything we measured.

### 2d. "A curated catalog beats the learned model" — CONFIRMED
A published AMR benchmark finds ML excels **only** on closely related strains, while the rule-based
**ResFinder is superior for divergent genomes**. And a 2020 pangenome-regression paper (*mBio*) states
outright that accuracy is similar between newer ML and simpler approaches, and that **the prediction
model matters far less than the dataset** — which is this repo's "labels, not models" conclusion,
published six years earlier.

---

## 3. The critique literature is broad, independent, and converges

| study | finding |
|---|---|
| **Ahlmann-Eltze, Huber & Anders, *Nature Methods* 2025** | 5 foundation models + 2 DL models vs *deliberately simple* baselines for perturbation prediction — **none outperformed**, including a "no change" baseline |
| **Kedzierska et al., *Genome Biology* 2025** | zero-shot scGPT/Geneformer don't consistently beat HVG selection or scVI; scGPT on PBMC performed **close to a randomly initialised model** |
| **DART-Eval** (regulatory DNA) | simpler ab-initio supervised models match or exceed fine-tuned DNALMs; DNALMs **particularly poor on counterfactual variant prediction**; zero-shot embeddings lose to k-means on motif counts |
| **BEND** | no consistent advantage over supervised baselines; no LM dominates; none approaches AUGUSTUS on gene finding |
| **"Fundamental limitations of gLMs for realistic sequence generation"** (2026) | Evo 2 generations fail to preserve long-range organisation, k-mer composition, TFBS architecture |
| **Evo 2** (40B, 9.3T nt) — own reported limits | underperforms alignment+structure methods on DMS; human gene-essentiality AUROC **0.66** |

**DART-Eval's methodological criticism is the one to internalise:** prior benchmarks rely on
*"oversimplified or flawed baselines that exaggerate the relative benefits of DNALMs."*

---

## 4. What this does to the anchor

| draft-anchor claim | verdict |
|---|---|
| "gene/protein-token is the untested model class" | **PARTLY FALSE** — gLM, gLM2, GenSyntax exist. What is untested is gene-token → **organismal phenotype** |
| "fluency ≠ function is the reframing" | **TRUE, NOT NOVEL** — named problem, ≥3 papers, public fix (LFB) |
| "scale is dead" | **CONFIRMED**, independently published |
| "natural populations are the trap" | **CONFIRMED** at 24,000 genomes; *and more data does not rescue it* |
| "curated catalog beats learned" | **CONFIRMED** (ResFinder > ML on divergent genomes) |
| **"the field may be clear"** | **FALSE — it is crowded** |

### The niche that actually survives

Every published gene-token model targets **molecular or cellular** endpoints — co-regulation, operons,
interfaces, function, gene essentiality, cell type. **None targets organismal phenotype**, and the
organismal attempts that exist (single-cell FMs, DNALMs, whole-genome AMR ML) are exactly the ones the
critique literature says lose to simple baselines.

So the surviving combination is narrow and specific:

> **gene/protein-token representation + an objective that is not pure likelihood (or LFB-corrected) +
> evaluated on constructed variation, against a curated-catalog baseline.**

Nothing found contradicts it. Every component constraint is now backed by *published* evidence rather
than only our own.

### The uncomfortable strategic read

The field has independently arrived at all four of our hard-won conclusions. That means our **findings**
are sound and our **methods** are good — but it also means the obvious versions of the idea are taken, and
several groups are further along on compute.

**What this repo has that the critique literature says the field lacks is the evaluation discipline
itself**: de-confounding by construction, margin-preserving nulls, denominator audits, curated-catalog
comparators, refusing to report a number from a concentrated cohort. DART-Eval, Kedzierska and
Ahlmann-Eltze are all, in effect, papers *about the absence of that discipline*. That is a more defensible
differentiator than another model.

### One possible genuine micro-finding

I did **not** surface a paper stating the specific mechanism this repo measured — that resistance
mutations are reached via chemically *conservative* substitutions at *averagely-conserved* sites, with
**BLOSUM62 reproducing the error** (which rules out memorisation). The nearest published statement is that
PLMs implicitly encode substitution conservativeness. The BLOSUM62 control may be genuinely novel.
**RESOLVED 2026-08-31 by a targeted check — see `wiki/prior_art_conservative_resistance_blosum_2026-08-31.md`.**
The finding is a COMPOSITE: three of four components are prior art (conservativeness published by Friedman
2013 for Abl1 — **with EGFR/ALK as a published RADICAL counterexample**; the drug-agnostic zero-shot PLM
framing is in print; BLOSUM62-as-baseline is standard). **Two moves appear unreported:** BLOSUM62 as a
*memorisation control*, and the *refutation* of the standard "the site is poorly conserved" explanation
for VEP false negatives (measured entropy percentile 0.494 = exactly average).

---

## Sources

- [Genomic language model predicts protein co-regulation and function — *Nature Communications* 2024](https://www.nature.com/articles/s41467-024-46947-9)
- [The OMG dataset / gLM2 — bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.08.14.607850v1) · [weights](https://huggingface.co/tattabio/gLM2_150M) · [code](https://github.com/TattaBio/gLM2)
- [Are genomic language models all you need? — *Bioinformatics* 2024](https://academic.oup.com/bioinformatics/article/40/9/btae529/7745814)
- [Decoding Prokaryotic Whole Genomes with a Product-Contextualized LLM (GenSyntax) — bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.03.692003.full.pdf)
- [scGPT — bioRxiv](https://www.biorxiv.org/content/10.1101/2023.04.30.538439v2.full) · [Heimdall tokenization framework](https://www.biorxiv.org/content/10.1101/2025.11.09.687403.full.pdf)
- [Zero-shot evaluation reveals limitations of single-cell foundation models — *Genome Biology* 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03574-x)
- [Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines — *Nature Methods* 2025](https://www.nature.com/articles/s41592-025-02772-6)
- [DART-Eval](https://arxiv.org/html/2412.05430) · [BEND](https://arxiv.org/html/2311.12570v3) · [Benchmarking genomic language models — *Nature Methods*](https://www.nature.com/articles/s41592-025-02829-6)
- [Biased sampling driven by bacterial population structure confounds ML prediction of AMR — *PLOS Biology* 2025](https://journals.plos.org/plosbiology/article?id=10.1371%2Fjournal.pbio.3003539)
- [Whole-genome phenotype prediction with ML: open problems in bacterial genomics — *Bioinformatics* 2025](https://academic.oup.com/bioinformatics/article/41/7/btaf206/8171528)
- [Improved prediction of bacterial genotype-phenotype associations using pangenome-spanning regressions — *mBio* 2020](https://journals.asm.org/doi/10.1128/mbio.01344-20)
- [Understanding language model scaling for protein fitness prediction — *Nature Computational Science*](https://www.nature.com/articles/s43588-026-01010-z)
- [Protein Language Model Fitness is a Matter of Preference — ICLR 2025](https://openreview.net/forum?id=UvPdpa4LuV)
- [From Likelihood to Fitness (LFB) — bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.05.20.655154v1) · [code](https://github.com/DiasFrazerGroup/lfb)
- [Efficient inference, training and fine-tuning of protein language models — *iScience* 2025](https://www.cell.com/iscience/fulltext/S2589-0042(25)01756-0)
- [Fundamental limitations of genomic language models for realistic sequence generation — 2026](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12871140/)
- [Evo 2 — Arc Institute](https://arcinstitute.org/tools/evo)
- [Assessing computational predictions of AMR phenotypes from microbial genomes — *Briefings in Bioinformatics* 2024](https://academic.oup.com/bib/article/25/3/bbae206/7665136)
