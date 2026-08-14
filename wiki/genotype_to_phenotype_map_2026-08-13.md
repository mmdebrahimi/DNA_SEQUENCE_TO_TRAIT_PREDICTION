# Where we are on "change one thing in E. coli → predict the phenotype change" (2026-08-13)

A step-back map of the whole quest, assembled from the committed artifacts rather than recall, plus a
**go/no-go on the next epoch's candidate experiment** (verdict: NO-GO, blocked on data).

## The short answer

For the narrow reading of the question — **change a regulatory part, predict the expression change** —
we are essentially there, and it is measured. The gap is not expression prediction. It is everything
*downstream* of expression.

## Each link in the chain

The design epoch plan (`wiki/design_epoch_plan_2026-08-07.md`) framed this as three questions.

| question | status | number |
|---|---|---|
| **Q1 — will the part work?** (amino acid → protein function) | **works** | ESM2-650M median Spearman **0.49** over 217 ProteinGym assays; TEM-1 0.73 |
| **Q2 — will the host express it?** (promoter/RBS → protein level) | **works, with error bars** | novel RBS within **~1.8×**, novel promoter within **~2.8×** of measured |
| **Q3 — can the host be rewired for yield?** (edits → flux → product) | **the wall** | growth-coupled succinate design: **0 designs found** |

### Q2 detail (`wiki/kosuri_expression_2026-08-11.md`)

Kosuri 2013, 12,563 promoter × RBS combinations with measured DNA/RNA/protein:

| what is being asked | R² |
|---|---|
| new *pairing* of two known parts | **0.893** |
| novel RBS + characterised promoter panel | 0.776 |
| novel promoter + characterised RBS panel | 0.497 |
| novel RBS, sequence alone, library confound removed | 0.495 |
| novel promoter, sequence alone, library confound removed | 0.352 |

Features are mechanistic (Shine-Dalgarno core + spacing; σ70 −35/−10 boxes + spacer), not learned. It
**failed its own pre-registered bar** (0.82 on an element split) — but that bar was mis-specified: 0.82
is a combination-level in-sample number, and an element-strength model has no strength for an unseen
element.

**The asymmetry worth remembering:** the promoter explains *more* protein variance (~54% vs ~30%) yet is
*less* predictable from its letters. Offered as interpretation, not result: translation initiation is
dominated by one short well-understood motif; promoter strength additionally depends on UP elements,
discriminator, spacer geometry, TSS selection and supercoiling, which these features do not represent.

## The shape the evidence makes

Sorted by causal distance from the edit:

| change → phenotype | chain | works? |
|---|---|---|
| promoter/RBS → expression of that gene | 1 step | **yes** |
| amino acid → protein function | 1 step | **yes** (ranking) |
| gene deletion → growth, one medium | medium | **yes**, MCC 0.652 |
| gene deletion → growth *pattern across media* | needs the network right | **no** |
| edits → maximise product yield | full network | **no** (0 designs) |
| whole genome → organism-level trait | longest | **closed negative, 0-for-5** |

**Prediction works where the causal chain is short and local. It fails as soon as the answer routes
through the metabolic network.** That is not a slogan — the failure mechanism is measured:
**76.9% of the model's missed conditional-essentiality calls are deletions where the model thinks nothing
happens** (`wiki/fba_constant_gene_diagnostic_2026-08-13.md`). The network carries routes the real cell
does not use.

All three repair levers are closed:

| lever | verdict |
|---|---|
| add reactions (gap-fill) | 154/5,425 flips, exact-set −1 |
| retune the threshold | ≤11% of the deficit is readout-recoverable |
| constrain routes (pFBA) | per-cell −0.0529, **+554 false positives**, exact-set +0 |

## The structural gap

**Nothing connects the layers.** We predict a construct's expression, and separately predict growth from
gene deletions. There is no path from *"I changed this promoter"* → *"expression changes 2×"* →
*"flux changes"* → *"phenotype changes"*. Each cell is validated in isolation and the **joint has never
been tested**.

## The candidate bridge experiment — assessed and NO-GO

The obvious next move, and the one that converges with where the FBA work landed: the last advance
concluded that the only untried idea is a **selective** restriction (target the specific redundant routes
that make a gene flat, rather than pFBA's blanket sweep). **Expression data is the selection principle.**
So the strategic bridge and the technical dead-end point at the same experiment:

> Constrain FBA with condition-matched measured expression (E-Flux / GIMME / iMAT style) and test whether
> it fixes the conditional-essentiality failure that a blanket parsimony sweep made worse.

**Assessed before building. Two independent reasons not to build it now:**

**1. Unfavourable literature prior.** Machado & Herrgård 2014 systematically evaluated GIMME, iMAT,
E-Flux, MADE, GX-FBA and Lee-12 on *E. coli* and *S. cerevisiae*: *"none of the methods outperforms the
others for all cases"*, and for many conditions **plain FBA with growth-maximisation and parsimony
criteria was as good or better** than expression-integrating methods.

That prior bites here in a specific way: the comparator the literature says expression-integration fails
to beat is **pFBA** — and on this substrate pFBA restriction is *worse than doing nothing*
(−0.0529 per-cell). "Expression-integration ≈ pFBA" would therefore also land below baseline.

*Scope caveat, and it is real:* that paper benchmarked **flux** prediction, not gene essentiality, and it
did not decompose by the flat-gene stratum that defines this failure. So the prior is unfavourable and
**indirect**, not decisive on its own.

**2. The condition-matching requirement, which is the harder blocker.** The 25 conditions here are
specific Biolog-style carbon sources — D-galacturonic acid, L-fucose, glycolate, D-sorbitol,
α-ketoglutarate. The experiment requires expression measured **under those same sources**; applying
glucose-grown expression to a galactose FBA is meaningless. Standard *E. coli* transcriptomic compendia
cover conventional lab conditions, not an exotic 25-carbon-source panel.

> **Labelled as an unfalsified hypothesis, not a finding:** I did not exhaustively search for a
> condition-matched dataset. What is verified is that no expression data exists locally (`feba.db` has
> 93 tables — fitness, cofitness, specific phenotypes — and no transcriptomics). Someone should check
> whether a carbon-source-panel transcriptome exists before this is treated as settled.

> ### ⚠ FALSIFIED THE SAME DAY — the data blocker does not hold
>
> One search overturned it. **PRECISE-NP881** (Nucleic Acids Research / PNAS, SBRG) is an 881-condition
> *E. coli* transcriptome compendium built around exactly this axis: **346 RNA-seq profiles generated
> during growth on 43 individual carbon sources**, each added to M9 minimal medium, combined with 535
> public MG1655 profiles from PRECISE-1K. It exists *because* PRECISE-1K lacked a carbon-source panel —
> three carbon iModulons (CRP-3, dmlA, SgcABCEQX) only emerged once those substrate conditions were
> added.
>
> **So the data-availability leg of the NO-GO is dead**, and this is the case I predicted would convert
> an external wall into a code wall: E-Flux itself is roughly a day's work.
>
> **What is verified vs not:**
> - **Verified:** a 43-carbon-source *E. coli* K-12 RNA-seq compendium exists, purpose-built for this axis.
> - **NOT verified:** the concrete overlap with my 25 Biolog sources, and access/licence terms. Overlap is
>   *likely* high — glucose, galactose, gluconate, glycerol, acetate, xylose, succinate, malate, pyruvate,
>   lactate, mannitol, sorbitol, ribose, maltose, mannose, fucose, NAG, glucosamine, glucuronate and
>   galacturonate are all standard carbon-panel members — but **that is an expectation, not a count.**
>   The decisive check is the substrate list in the paper's Dataset S1 or the per-sample carbon-source
>   annotations in `SBRG/precise1k`'s metadata table.
>
> **The other leg of the NO-GO still stands, and it is now the only one:** the Machado & Herrgård prior
> is unfavourable (plain FBA + parsimony as good or better than expression integration) — though
> indirect, since it benchmarked flux rather than essentiality and did not decompose by the flat-gene
> stratum.
>
> **Revised verdict: NO-GO → CONDITIONAL GO**, gated on counting the overlap. If ≥10 of the 25 sources
> match, the experiment is buildable now and the only remaining risk is the literature prior — which is
> exactly the kind of risk worth spending a day to test, because the failure mode is informative either
> way.
>
> **Process note worth keeping:** the blocker was recorded as an *unfalsified hypothesis* rather than a
> finding, and checking it flipped the verdict in a single search. Had it been written as "no such data
> exists", a whole direction would have been wrongly closed.

**Why not substitute the fitness data itself as the selection signal:** RB-TnSeq fitness *does* say which
genes the cell uses per condition, and it is sitting right there — but constraining the model with the
same labels it is scored against is circular, and this project has a documented history of exactly that
trap. Expression is the right input **because** it is independent of the fitness labels.

## What would unblock it

- A carbon-source-panel *E. coli* transcriptome overlapping ≥10 of the 25 sources. That converts an
  external wall into a code wall — the E-Flux implementation itself is a day's work.
- Or: re-scope to the conditions where matched expression *does* exist (glucose / acetate / glycerol /
  LB), accepting a much smaller condition panel and correspondingly weaker power.

## Three genuinely different next investments

1. **Deepen Q2** — codon usage (named in the plan, completely untouched), native chromosomal context, a
   non-GFP reporter. Extends the one layer that works.
2. **Bridge Q2→Q3** — the experiment above. Highest value, currently data-blocked.
3. **Accept Q3's wall** — 0 growth-coupled designs plus a comprehensively closed conditional-essentiality
   question suggest constraint-based modelling may be the wrong instrument for the design direction.

This is a strategic call about what the project is for, not a technical one.

## Provenance

Assembled 2026-08-13 from: `wiki/kosuri_expression_2026-08-11.md`,
`wiki/design_epoch_plan_2026-08-07.md`, `wiki/fba_constant_gene_diagnostic_2026-08-13.md`,
`wiki/fba_regulatory_carbon_test_2026-08-13.md`, `wiki/fba_label_threshold_sweep_2026-08-13.md`,
`wiki/fba_infeasibility_finding_2026-08-13.md`, `wiki/fba_keio_validation_2026-08-03.json`,
`wiki/fba_strain_design_succ_2026-08-07.md`, `wiki/proteingym_esm2_650m_full_2026-07-09.md`.

Literature: Machado D, Herrgård M (2014) *Systematic Evaluation of Methods for Integration of
Transcriptomic Data into Constraint-Based Models of Metabolism.* PLOS Comput Biol 10(4):e1003580.
