# Track C step 2 is aimed at the wrong target (2026-08-12)

Step 2 of Track C was specified as **"protein-function prediction over dark-matter genes → candidate
missing reactions → measure the FBA accuracy delta."** Before building it, I analysed the error set it is
meant to fix. It cannot fix it, for reasons that are measured rather than argued.

The premise established on 2026-08-11 still holds — **gaps do explain the false negatives** (FN
gap-adjacency 66.3% vs TN 41.2%, p = 2.1×10⁻⁶). What does not follow is the jump from *"gaps explain the
errors"* to *"therefore predict protein function."*

## Finding 1 — the false-negative set contains zero dark matter, and cannot contain any

**Structurally:** FBA gene-essentiality is computed over the *model's own* gene set — the 905 genes that
have a GPR in iMM904. A dark-matter gene is by definition one whose function is unknown, so it has no GPR
and cannot appear in that set. The FN set therefore cannot contain dark matter, ever.

**Empirically**, all 57 gap-caused false negatives are named, fully-annotated enzymes:

| group | n | examples |
|---|---|---|
| **aminoacyl-tRNA synthetases** | **21** | Isoleucyl-, Glycyl-, Seryl-, Lysyl-, Arginyl-, Phenylalanyl-tRNA synthetase … |
| heme biosynthesis | 5 | 5-aminolevulinate synthase, porphobilinogen synthase, uroporphyrinogen decarboxylase, coproporphyrinogen oxidase, heme O monooxygenase |
| cofactor activation | 6 | FMN adenylyltransferase, riboflavin kinase, pantothenate kinase, GTP cyclohydrolase I, biotin–acetyl-CoA carboxylase ligase, biotin uptake |
| lipid / membrane | 8 | serine C-palmitoyltransferase ×2, 3-dehydrosphinganine reductase, PI 3-kinase ×2, PI 4-kinase ×2, dolichol kinase |
| cell wall / glycosylation | 5 | chitin synthase, UDP-GlcNAc diphosphorylase, GlcN6P synthase, phosphoacetylglucosamine mutase, GDP-mannose antiport |
| other named enzymes | 12 | fatty acid synthase, deoxyhypusine synthase, N-myristoyltransferase, peptide α-N-acetyltransferase … |

Not one is an unknown. A function predictor has nothing to contribute here.

## Finding 2 — the gaps are missing *demand*, not missing enzymes

The dead-end metabolites blocking these genes are orphaned **products**, and the pattern is unmistakable:

```
trnaphe_c   no_producer      phetrna_c   no_consumer
trnaala_c   no_producer      alatrna_c   no_consumer
trnacys_c   no_producer      cystrna_c   no_consumer      ... 21 of these
```

iMM904 contains the charging reaction (amino acid + tRNA → aminoacyl-tRNA) but **translation is not
represented**, so nothing consumes the charged tRNA and nothing regenerates the free tRNA. Same shape for
heme (no cytochrome assembly), sphingolipids and PIs (membrane structure not in biomass), and the
cofactors.

The missing consumers are **non-metabolic processes** — translation, respiratory-complex assembly,
membrane biogenesis. A metabolic-reaction function predictor cannot supply them *even in principle*.

## Finding 3 — adding the demand doesn't help. Measured, exactly zero effect

A label-blind rule — add a demand reaction for **every** `no_consumer` dead-end, chosen by model structure
alone and never by which genes are false negatives — was applied and the whole gene set re-scored:

| | TP | FP | FN | TN | MCC |
|---|---|---|---|---|---|
| baseline | 34 | 13 | 101 | 757 | 0.3773 |
| **+106 demand reactions** | 34 | 13 | 101 | 757 | **0.3773** |

Identical to four decimal places. **Unblocking a reaction does not make its gene essential.** A demand
sink lets the product drain away; it does not make the cell *require* it. Only membership in the biomass
objective does that.

## Finding 4 — the fix that would work is circular, so it must not be used

Adding charged tRNAs, heme and sphingolipids to the biomass function *would* flip these genes to essential
— because we would be adding to biomass precisely the products whose producers we already know are
essential. That is fitting the test set, and the resulting MCC gain would mean nothing.

**Conclusion: the yeast essentiality metric cannot be honestly improved by gap-filling or by function
prediction of any kind.** Its remaining error is a biomass-scope property of a 2009 reconstruction. The
ceiling if every gap-caused FN flipped is MCC 0.734 (from 0.377), and that ceiling is not honestly
reachable by us.

## Where this leaves Track C

**Do not build step 2 as specified.** Three redirects, ranked.

### 1. Test the reconstruction, not the annotation *(cheap, non-circular, in flight)*

iMM904 is from 2009. **yeast-GEM v9.1.0 shipped 2026-07-01** and is actively maintained — free, 8 MB.
Its biomass comes from measured cell composition, not from our labels, so scoring it against the same SGD
gold standard is non-circular. If a modern reconstruction closes the gap that no gap-filling can, then the
FBA cell's cross-organism weakness is **model vintage**, and the fix is a model registry refresh rather
than an ML build. *Status: downloaded and scoring; the sweep is slow on a 4,000-reaction model.*

### 2. Point gap-filling at the metric where it already demonstrably works

Gap-filling has a **measured win** on carbon-source growth, not essentiality: iML1515 predicted 0.000 /h on
sucrose, the structural diagnostic localised `suc6p_c` as a dead end, and adding one donor reaction (`FFSD`)
restored growth to 1.7798 /h. That is the shape Track C wants — and it is *design-relevant*, which
essentiality is not.

**But it needs a two-sided gold standard first.** The current carbon validation is `SCORED_RECALL`: 21
mapped positive sources, recall 1.0, **no negatives**. Gap-filling can only ever *increase* growth
capability, so a recall-only benchmark rewards over-filling and cannot validate it — a model that grows on
everything scores 1.0. The prerequisite is a substrate set with genuine non-growth cases (Biolog PM
plates), and acquiring that is the real first step here.

### 3. Genuine dark-matter prediction, aimed at growth rather than essentiality

This is the only framing in which function prediction earns its place: a dark-matter gene may encode a
**transporter or catabolic enzyme** that explains growth on a substrate the model says is impossible. The
gene is *not in the model*, so predicting its function genuinely adds a reaction — unlike the essentiality
case, where every candidate gene is already present. Metric: growth-phenotype accuracy on the two-sided
substrate set from (2). Gated on that gold standard existing.

## The reusable lesson

**"Gaps explain the errors" does not imply "predict the missing functions."** Both halves were true here —
the gaps are real and they do concentrate on the errors — but the missing pieces were *known biochemistry
with no downstream consumer*, not unknown biochemistry. The check that separated them cost about twenty
minutes: list the genes and read their names.

This is the third time in two days that pricing the cheaper explanation first redirected the work
(medium-vs-structural, model-class-vs-data, and now annotation-vs-representation).

## Reproduce

The three probes are in the session scratchpad; the durable ones are folded into
`scripts/fba_gap_premise_check.py`. The FN anatomy is derivable from that script's
`wiki/fba_gap_premise_yeast_2026-08-11.json` plus the model.
