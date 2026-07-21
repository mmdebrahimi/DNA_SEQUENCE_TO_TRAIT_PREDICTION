# Plain-language glossary — the jargon in the findings analysis

Companion to `findings_deep_analysis_2026-07-21.md`. Every technical term explained simply, with an
analogy where one helps. Grouped by topic.

---

## The big picture

**Genotype → phenotype.** Genotype = the DNA sequence (the *code*). Phenotype = the observable trait — what
the organism actually *does* (e.g. "this bacterium survives ciprofloxacin"). Decoding = predicting the trait
from the DNA. **Analogy:** predicting a program's runtime behavior by reading its source code.

**Determinant.** One specific, known DNA change that *causes* a trait — e.g. "mutation X in gene gyrA makes
the bug resistant to cipro." The causal switch.

**Curated catalog / determinant-scan.** A hand-verified lookup table of determinants. "Determinant-scan" =
read the DNA, look up each change in the table, decide the trait. **Analogy:** a signature-based antivirus —
scan the file, match known bad patterns, flag it. This is the approach that *works* in this project.

**Resistant (R) / Susceptible (S).** R = the drug does NOT kill the bug. S = the drug DOES kill it. The whole
game is predicting R vs S from the DNA.

---

## The biology mechanisms (the "what makes it resistant" names)

These are just the names of specific resistance machinery — one-liners:

- **QRDR** — the exact spot in a gene where a mutation makes a bug shrug off fluoroquinolone drugs
  (ciprofloxacin, levofloxacin). Clean, single-mutation, easy to catalog.
- **β-lactamase / ESBL** — an enzyme the bug produces that *chops up* penicillin-family drugs (called
  β-lactams: penicillin, cephalosporins). **ESBL** = an "extended-spectrum" version that chops up more of
  them. If the gene is present, the bug is resistant.
- **van clusters (vanA / vanB)** — a gene package that rebuilds the bug's cell wall so vancomycin can't grab
  on. vanA also defeats teicoplanin; vanB doesn't — that's how the two drugs are told apart.
- **penA mosaic** — a patchwork ("mosaic") version of a gene that lowers a gonorrhea bug's susceptibility to
  cephalosporins like cefixime.
- **ERG11 / FKS1** — the fungal target genes: ERG11 mutations resist azole antifungals (fluconazole), FKS1
  mutations resist echinocandins (micafungin).

**The pattern:** each of these is a compact, known cause. That's *why* the catalog approach works here — there's
a real "if this gene/mutation, then this trait" rule to look up.

---

## The two words that govern how to read every number

**Clonal / clone / clonal inflation.** Bacteria reproduce by cloning — making near-identical copies of
themselves. If a single resistant bug spreads through a hospital, you might collect 60 samples that are really
57 copies of *one* original + 3 genuinely different bugs. **Clonal inflation** = counting those 60
near-duplicates as 60 independent data points, which fakes high accuracy. **Analogy:** measuring a poll's
accuracy by asking the *same person* 60 times and reporting "60 respondents." Everywhere you see a raw
sens/spec number, it's ~2–2.5× inflated by this — you have to look at the clone-collapsed number instead.

**Lineage.** A family of clones — a branch of the family tree, all descended from one recent ancestor, so
nearly identical. **Analogy:** a git branch where every commit is one tiny change from the last. The honest
fix for clonal inflation is "one vote per lineage," not "one vote per sample."

**Effective-N.** The *real* number of independent data points after collapsing clones. "60 isolates → 3
effective" means you really only have 3 independent bugs' worth of evidence. A "SCORED" result resting on
effective-N = 2 is much weaker than one resting on effective-N = 15, even if the raw number looks the same.

---

## The three "regimes" (the project's core organizing idea)

The finding that a trait falls into one of three categories, and each needs a *different* approach:

1. **Curated-catalog regime** → a trait caused by a *known handful of genes*. **Deterministic rules win.** (All
   the successful cells.)
2. **Organism-polygenic regime** → "poly" = many, "genic" = genes. A trait controlled by *many* genes
   scattered across the whole genome, with no single switch. **Nothing works well** here — there's no compact
   rule to look up, and AI models just learn who's-related-to-whom instead of the actual cause. **Analogy:** a
   bug that only appears from the interaction of 50 different config flags — no single line to point at.
3. **Molecular-property regime** → predicting a *protein's* function from its sequence. **Learned AI wins —
   but only if the trait is "fitness-aligned"** (see below).

---

## Why AI *fails* on drug resistance specifically

**Fitness-aligned vs antagonistically-selected.** This is the subtle, important one.

- **Fitness-aligned** = the trait tracks how *healthy/functional* the protein is. Mutations that break the
  protein hurt the organism, and AI models trained on evolution are good at spotting "this looks broken."
- **Antagonistically-selected** (drug resistance) = the *opposite*. Resistance mutations are ones evolution
  normally selects *against* — they slightly harm the bug in normal life. Only under drug pressure do they
  become useful. So to an evolution-trained AI they look "mildly damaging but tolerable" = *normal* → and it
  calls them benign. **Exactly backwards.** That's why off-the-shelf AI scores *below a coin flip* on
  resistance. (Even a 1992 amino-acid similarity table beats it — the failure is fundamental, not fixable by
  a bigger model.)

**Convergent evolution vs clonal (for AI training).**
- **Convergent** = the same trait arises *independently* many separate times (many unrelated bugs each evolve
  resistance on their own). Great for learning — genuinely independent examples. HIV is like this.
- **Clonal** = the trait spread by *one* lineage copying itself. Looks like many examples but it's really one.
  A trained model can't learn a general rule from it. TB is like this — which is why the trained model works
  on HIV (0.81) but collapses to a coin-flip (0.51) on TB.

---

## The AI / model terms

**Zero-shot vs supervised.**
- **Zero-shot** = use a big general AI *off the shelf*, no training on your task — just ask it. (This is what
  *failed*: "0-for-N" refers to zero-shot only.)
- **Supervised** = a model you *trained* on labeled examples of your specific task. (This *works*, in the
  convergent regime.)
- The shipped tool is a **hybrid**: the deterministic catalog + a supervised model to cover the catalog's gaps.

**Foundation model / embedding.** A foundation model = a huge general AI trained on tons of DNA/protein
sequences (like an LLM, but for biology). An **embedding** = the vector of numbers it produces to represent a
sequence — its numeric "fingerprint." The bet was: feed DNA in, use the fingerprint to predict traits. It
failed here — the fingerprint captured *ancestry* (who's related to whom), not the resistance *mechanism*.

**ESM2 / ProteinGym / DMS.**
- **ESM2** = a specific well-known protein foundation model.
- **DMS (deep mutational scanning)** = a lab technique that measures the effect of *thousands* of individual
  mutations in one experiment — a complete "what does each single change do" map. The gold-standard label for
  the protein work.
- **ProteinGym** = a standard benchmark collection of those DMS experiments.

---

## The validation / statistics terms

**Sensitivity / specificity.** Two different kinds of "accuracy":
- **Sensitivity** = of the *truly resistant* bugs, what fraction did we correctly catch? (catching threats)
- **Specificity** = of the *truly susceptible* bugs, what fraction did we correctly clear? (not false-alarming)
- Both matter. High sensitivity + low specificity = you cry wolf constantly. Low sensitivity = you miss real
  threats. A good cell needs both high.

**Powered / underpowered / SCORED_ENDORSED.**
- **Powered** = you have enough samples of *both* classes (enough R *and* enough S) for the number to mean
  something.
- **Underpowered** = too few of one class to trust it (e.g. "sens 1.0" when there were only 3 resistant bugs).
- **SCORED_ENDORSED** = powered *and* accurate enough to stand behind.

**AUC / AUROC.** A single score from 0.5 to 1.0 for how well a predictor separates two classes. 1.0 = perfect,
0.5 = a coin flip. 0.96 = very good. (Used for the HIV cells.)

**Spearman.** A correlation score (−1 to 1) for how well two *rankings* agree — used for the protein work
("did the model rank the variants in the same order the lab experiment did?").

**Wilson CI (confidence interval).** An honest error-bar on a percentage, especially when you have few
samples. "sens 0.5 [0.15–0.85]" means: best guess 50%, but it could genuinely be anywhere from 15% to 85%. A
*wide* bar = weak evidence. Always read the bar, not just the point.

---

## The data-integrity terms (why the results are trustworthy)

**Provenance-disjoint.** Provenance = where the data came from (which lab, study, country). Disjoint = zero
overlap. So the tool is tested *only* on data from completely separate sources than it was built from — you
never grade yourself on data related to your training set. Prevents accidental cheating.

**Circular label.** If the "answer key" you grade against was itself produced by a tool *like yours*, you're
only checking whether you *agree with the other tool* — not with reality. **Analogy:** grading your exam
against a copy of your own answers. A big chunk of available data is disqualified for this reason.

**In-distribution vs independent.**
- **In-distribution** = tested on data that overlaps / resembles what you built from. Weaker proof.
- **Independent** = tested on genuinely separate, unseen data. The real proof. (HIV's independent label is
  what makes it the flagship result.)

**MIC / broth microdilution.** MIC = Minimum Inhibitory Concentration = the smallest drug dose that stops the
bug growing, *measured in a real lab* (grow the bug in tubes with increasing drug, see where it stops). This
is the trusted wet-lab "ground truth" — above a threshold = R, below = S. The whole project is bottlenecked on
finding *free, lab-measured* MIC data like this.

**Reads / assembly / SRA / shallow reads.**
- **Reads** = the raw short DNA fragments a sequencer produces — like the torn pieces of a shredded book.
- **Assembly** = stitching the reads back into the full genome (reassembling the book).
- **SRA** = the public archive where raw reads are stored.
- **Shallow reads** = too few fragments to reassemble the book → gaps → a resistance gene might fall in a gap
  and get missed. (This is what blocked the Enterococcus susceptible-class this session.)

**AMRFinder.** A standard existing tool that scans a genome and lists which known resistance genes/mutations
it found. Our decoder uses its output as the raw "which determinants are present" list, then applies the
curated rules on top.

---

## One-line summary of the whole thing

We built a **signature-scanner for drug resistance** (catalog of known DNA→trait rules), it works wherever a
trait has a *known compact cause* and a *real lab-measured answer key* to check against, and the hard limit is
**finding that free lab-measured answer key** — not the AI, not the compute.
