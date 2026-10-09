# SOURCE ARCHIVE — ChatGPT discussion: AlphaGenome's transformer + how to improve the architecture

**Archived verbatim 2026-10-09, at the user's instruction, because a plan was built on it and it was LOST.**

> **Why this file exists.** This discussion is the direct origin of GLM families **G-A** (its priority #1,
> multi-resolution tokens + trivial head) and **G-C** (its priority #3, cell-state conditioning) — the
> ledgers' verbatim originating goals name both. It was pasted into a session, the session was compacted,
> and the content was gone. A later session could still see the *plan* but not the *reasoning the plan came
> from*, and the 2026-10-07 distillation had deliberately kept only the training-method half, discarding
> exactly the architecture half reproduced below. **A pasted source that a plan depends on must be archived
> to the repo before it can be compacted away.** Screening + lessons:
> `wiki/chatgpt_architecture_discussion_lessons_2026-10-09.md`.

**Provenance:** ChatGPT (model/version not recorded by the user). It is a *secondary* source — a model's
explanation and proposal, not a paper. The sibling archive note at
`wiki/alphagenome_lessons_actionable_2026-10-07.md` records that an earlier circulated ChatGPT summary of
AlphaGenome contained **two factually wrong claims** (a "published Oct 2024" / "Nature 2024" citation;
AlphaGenome was announced 25 Jun 2025 and published in *Nature* 28 Jan 2026, vol 649 pp 1206–1218). Treat
the content below as a well-reasoned **hypothesis set**, not as established fact.

---

## Part 1 — the transformer stage in AlphaGenome (as explained)

At this point the CNN has already compressed the ~1 Mb DNA sequence into roughly **7,800 regional feature
vectors**. So the transformer does not mainly ask *what nucleotide is here?* or *is there a motif here?* —
the CNN mostly handled that. The transformer asks:

> **How should information from one genomic region change the interpretation of another genomic region far
> away?**

**What goes into the transformer.** Each ~128-bp region is represented by a feature vector containing
learned local information: motif patterns · promoter/enhancer-like signals · splice-related signals ·
chromatin-related sequence features · combinations of local motifs. So the transformer sees
`Region 1 → vector`, `Region 2 → vector`, … `Region ~7,800 → vector`.

**The key mechanism: self-attention.** For each region the transformer can examine all other regions and
ask *which of these are relevant to interpreting me?* If region 2,000 looks like an enhancer and region
5,400 looks like a promoter, it can learn a strong relationship between them even across hundreds of
kilobases. An unrelated region receives very little attention.

**Query / Key / Value intuition.** Each regional representation is transformed into three learned
representations:

- **Query** — what information am I looking for?
- **Key** — what kind of information do I contain?
- **Value** — what information should I contribute if I am relevant?

A region's Query is compared with the Keys of all regions; strong matches get higher attention weight;
their Values are combined to update that region. Roughly: *region A asks a question → searches all regions
→ collects useful information → updates its representation.*

**Why multiple heads.** Different heads can specialise in different relationships — enhancer↔promoter,
CTCF↔CTCF, exon↔splice-site context, promoter↔regulatory-domain structure, relationships between nearby
regulatory elements. The model is **not told** "head 3 should learn enhancer–promoter interactions"; these
specialisations **emerge during training** if they improve prediction.

**One layer is not enough.** A stack of transformer blocks, each doing self-attention → combine distant
information → feed-forward → next layer. With depth, representations can become progressively more
abstract: early layers *region A looks relevant to region B*; middle layers *A + B + C form a regulatory
neighborhood*; deeper layers something closer to *these combined signals likely affect transcription and
chromatin structure around this gene*. These are **conceptual interpretations of latent representations,
not explicit symbolic statements**.

**Why it matters.** With `Enhancer ——— 300,000 bp ——— Gene promoter`, the CNN can identify enhancer-like
and promoter-like sequence but is poor at connecting things hundreds of thousands of bases apart. The
transformer provides that bridge. The division of labour:

| Component | Main role |
|---|---|
| CNN encoder | Understand local DNA structure |
| Transformer | Understand long-range relationships |
| Decoder | Return global understanding to fine resolution |

**The most important conceptual point.** The transformer is not simply "looking farther" — it creates a
**context-dependent representation**. Before: *"this looks enhancer-like."* After: *"this enhancer-like
region is relevant to this particular promoter, inside this regulatory neighborhood, given all the other
surrounding genomic information."*

---

## Part 2 — how Query/Key/Value mathematically create long-range interaction

From each region vector the model creates Query, Key and Value as **learned linear transformations of the
same input vector**.

1. **Query** = what this region is looking for (numerically, not in language) — e.g. for a promoter,
   *"which other regions are relevant to interpreting this promoter?"*
2. **Key** = what each region advertises about itself — an enhancer-like region has one kind of key, a CTCF
   site another, an inactive region another.
3. **Similarity score** = usually a dot product of `Query(A)` with `Key(B)`. Larger ⇒ more relevant.
   Illustrative raw relevances: nearby neutral region 0.2 · enhancer-like 4.7 · CTCF site 2.1 · unrelated
   −0.8. Softmax normalises these into attention weights, e.g. B 0.62 · C 0.21 · D 0.10 · others 0.07.
4. **Value** = the information actually transferred. The key decides *whether* a region is relevant; the
   value is *what gets passed*. So
   `Updated Region A = 0.62 × Value(B) + 0.21 × Value(C) + 0.10 × Value(D) + …`, plus the original
   representation through residual connections. **This is the actual mechanism by which distant genomic
   information gets moved around.**

**A genomic example.** Region A = promoter, B = enhancer 250 kb away, C = CTCF boundary, D = random
intergenic. The promoter's query matches the enhancer and CTCF keys strongly, the random region weakly — so
the promoter representation is updated using the enhancer and boundary, becoming *"promoter-like sequence
influenced by this distal enhancer, with this surrounding regulatory architecture."*

**The scaling factor.** `Attention(Q, K, V) = softmax(QKᵀ / √d) V`. Dividing by √d prevents dot products
from getting too large as dimensionality increases; without it softmax becomes too extreme and training
destabilises.

**Why it works at long genomic distance.** Attention does not care that two regions are 300,000 bases apart
the way a local convolution does. If both are present in the input, A can attend directly to B.

**Multi-head attention** repeats this in parallel with different learned projections, then combines — so the
model learns several notions of "relevance" rather than one.

**Important caveat (quoted, and load-bearing).** *"Attention weights are not guaranteed to correspond
directly to biological causality. A high attention score does not automatically mean 'this enhancer
biologically controls this promoter.' It means 'the model found this region useful when constructing the
representation.'"*

**Shortest summary.** For every region: Query = what do I need · Key = what do I contain · Query × Key =
how relevant are you to me · softmax = relevance → weights · Value = transfer the useful information ·
repeat across heads and layers.

---

## Part 3 — "How would you improve this architecture?"

**Stated main weakness:** AlphaGenome has a large receptive field, but that does not mean it has solved
long-range biological reasoning. DeepMind explicitly says interactions beyond **~100 kb** remain difficult,
tissue/cell specificity needs improvement, and the model is **not designed for personal-genome or
complex-trait prediction**.

**Proposed shape:**

```
DNA / haplotypes → local sequence encoder → multi-resolution genomic tokens
→ regulatory-element discovery layer → hierarchical + sparse long-range transformer
→ explicit genomic interaction graph → cell-state-conditioned biological representation
→ molecular heads + mechanistic intermediate heads → high-resolution decoder
→ molecular phenotype / eventually organism phenotype
```

### The sixteen proposed changes

**1. Don't force everything through fixed 128-bp tokens.** Biology does not divide into uniform chunks:
TF motif ~6–20 bp · splice signal a few bp · nucleosome ~147 bp · promoter hundreds bp · enhancer
hundreds–thousands bp · gene kb–Mb · TAD hundreds kb–Mb. Instead maintain **multiple resolutions
simultaneously** (1 bp → 16 bp → 128 bp → 1 kb → 10 kb) and let information travel between levels —
analogous to an image pyramid.

**2. Create biological tokens, not just spatial tokens.** Introduce learned tokens for candidate biological
entities — promoter, enhancer, splice-site, exon, gene, CTCF-site, regulatory-domain — **discovered by the
model rather than taken from fixed annotations**. So instead of `region 4,981 ↔ region 6,732`, the model
reasons `enhancer E37 ↔ promoter P12`, much closer to the structure of the underlying problem.

**3. Hierarchical attention instead of one flat transformer.** ~7,800 regions all talking to one another is
a flat representation, but the genome has hierarchy: L1 motifs↔motifs · L2 motifs→enhancer · L3
enhancers↔promoter↔gene · L4 genes↔TAD↔chromatin domain · L5 domain↔domain. Information moves
base → motif → enhancer → gene → domain and back down, drastically reducing the need for every location to
attend to every other.

**4. Sparse long-range attention.** Most DNA does not need to communicate strongly with everything else, so
full attention wastes capacity. Give each region **local attention** (±10–50 kb) plus **candidate distal
attention** — select a small number of possibly relevant distant regions (e.g. a promoter picks an enhancer
80 kb upstream, one 260 kb downstream, a CTCF boundary) and do expensive attention only over those:
`7,800 regions → retrieve top ~20 → deep attention among those 20`. Could scale 1 Mb → 10 Mb →
chromosome-scale without quadratic explosion. Matters because *merely giving the model 1 Mb does not
guarantee it learns the important distal interactions* — DeepMind notes >100 kb effects remain challenging.

**5. Add an explicit genomic interaction graph** (described as the largest structural departure). DNA
sequence is linear but **functional genome organisation is closer to a graph**. Nodes: enhancers,
promoters, genes, splice sites, CTCF sites, domains. Edges: enhancer→promoter, promoter→gene, CTCF↔CTCF,
exon→exon, element→domain. Then run a graph transformer / GNN on top, giving **sequence model + regulatory
graph** instead of sequence alone.

**6. Make 3D genome structure a central latent variable**, not merely another output:
`DNA → predicted chromatin organisation → regulatory interaction reasoning → expression`. Biologically, 3D
organisation determines which enhancer can reach which promoter — this provides a **mechanistic
bottleneck**.

**7. Explicit cell-state conditioning.** The DNA is identical in neuron, hepatocyte, muscle, B cell, yet
expression differs enormously — therefore **DNA alone cannot completely determine cellular behaviour**.
Condition explicitly: `DNA + cell-state embedding → molecular phenotype`, where cell state could include TF
abundance, chromatin accessibility, methylation state, transcriptional state, developmental stage, tissue
identity. DeepMind itself identifies improved cell/tissue-specific modelling as an important future
direction.

**8. Use actual perturbations much more aggressively.** AlphaGenome learns mostly from *observational*
assays, which creates the classic problem: **correlation ≠ causation**. If enhancer A and gene B are
correlated, the model can learn the relationship without learning that A *causes* B. Incorporate CRISPR
knockout, CRISPRi, CRISPRa, MPRA, saturation mutagenesis, Perturb-seq, eQTL, sQTL — so training includes
`normal enhancer → expression 1.0` versus `destroy enhancer → expression 0.3`. Pushes the network toward
**causal** regulatory effects rather than statistical association.

**9. Train directly on variants.** Instead of primarily `sequence → track`, train jointly on
`reference + alternative → Δ molecular state`, making variant effects a **first-class task**: Δ TF binding,
Δ accessibility, Δ chromatin interaction, Δ expression, Δ splicing. This teaches the model an explicit
concept of **intervention**. (Notes AlphaGenome's student model already sees mutationally perturbed
sequences during distillation.)

**10. Model diploid genomes.** Humans have maternal + paternal chromosomes and variants interact.
`Haplotype A + Haplotype B → diploid regulatory state` allows representing compound heterozygosity, cis
interactions, phase, allele-specific expression, dominance, recessivity, multiple nearby variants —
necessary for serious personal-genome modelling.

**11. Model combinations of variants.** Real genomes contain millions of variants; asking "what does this
one SNP do?" answers only part of the problem. Need interactions where A alone ≈ negligible, B alone ≈
negligible, but **A + B → substantial effect** — i.e. **epistasis**. Train explicitly on multi-variant
perturbations.

**12. Self-supervised genomic pretraining before functional supervision.** Recent benchmarking suggests
general DNA foundation models sometimes capture sequence-level constraint/pathogenicity better than
specialised track-prediction models, while AlphaGenome dominates context-dependent regulatory QTL tasks.
Combine: massive evolutionary/self-supervised DNA pretraining (thousands of humans, primates, mammals,
vertebrates, broader species) → functional genomic supervision → perturbational/variant supervision. So
**evolution → molecular biology → causal perturbation**, each stage teaching something different.

**13. Incorporate evolutionary information explicitly.** Evolution has run an experiment for hundreds of
millions of years; near-perfect conservation across mammals is extremely informative. Feed phylogenetic
conservation, orthologous sequences, population frequencies, evolutionary constraint — as auxiliary inputs
or pretraining objectives — giving the model another source of information about *what changes biology
tolerates*.

**14. Introduce mechanistic intermediate supervision.** Rather than `DNA → expression`, prefer
`DNA → TF binding → chromatin accessibility → enhancer activity → 3D interaction → promoter activity →
transcription → RNA processing`. **Not** forcing biology into a rigid hand-written pipeline — rather create
latent intermediate representations and supervise them wherever data exists. (AlphaGenome already benefits
strongly from multimodal training; its ablations show the value of combining modalities.)

**15. Add uncertainty.** Instead of `variant effect = 0.71`, produce `predicted effect = 0.71,
confidence = high`, or a distribution. The system should be able to say *"I have not seen enough biology
resembling this situation"* — especially for clinical/genotype applications.

**16. Mixture-of-experts.** Don't force one network to be equally expert at everything: a shared genomic
backbone plus specialised experts (splicing, enhancer, promoter, 3D genome, immune-cell, neuronal,
developmental) with a router. Addresses the known tension between general multimodal models and specialised
models — the AlphaGenome literature explicitly recognises specialised models can remain superior on some
individual tasks.

### The ten-stage architecture it would actually build

1. DNA — diploid/haplotype-aware raw DNA
2. Local encoder — motifs, splice signals, local regulatory grammar
3. Multi-resolution representation — 1 bp + 16 bp + 128 bp + 1 kb + 10 kb simultaneously
4. Regulatory-element discovery — latent enhancer / promoter / gene / CTCF / exon tokens
5. Hierarchical sparse transformer — reason locally, retrieve candidate distant interactions
6. Regulatory graph — enhancer ↔ promoter ↔ gene ↔ chromatin domain
7. Cell-state conditioning — cell type + TF state + developmental state + epigenetic state
8. Mechanistic latent biology — TF binding → accessibility → chromatin architecture → transcription → splicing
9. Decoder — return global reasoning to nucleotide resolution
10. Counterfactual variant engine — `genome A → state A` vs `genome B → state B` ⇒ **Δ biological state**

### And for the actual goal (the strongest framing in the discussion)

AlphaGenome approximately stops at `DNA → molecular phenotype`. The user's target is
`genotype → phenotype`, which requires several additional reasoning levels:

> genotype → molecular regulatory effects → protein effects → cellular state → cell–cell interactions →
> tissue state → organ state → physiology → development over time → environmental interaction →
> observable phenotype

AlphaGenome's authors explicitly acknowledge that molecular predictions do not capture complex traits and
disease because **development and environment lie outside its scope**. So: *"I would not build your model as
simply 'a much bigger AlphaGenome.' I would use an AlphaGenome-like system as the **bottom genomic layer of
a hierarchical biological world model**."*

### The three highest-priority changes, if resources were limited

1. **Multi-resolution + hierarchical tokens**, rather than one fixed 128-bp transformer representation.
2. **Explicit regulatory graph + sparse distal attention**, so the model learns enhancer/promoter/gene
   relationships instead of arbitrary window-to-window correlations.
3. **Causal perturbation + cell-state conditioning**, because sequence alone and observational tracks
   cannot get you to robust biological causality.

> *"Those three changes, in my view, are more important than simply making the transformer larger."*
