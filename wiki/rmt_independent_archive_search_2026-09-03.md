# Finding an independent archive: three candidates eliminated on structure, one that worked

The named external wall — test the gentamicin `rmt` rescue's specificity somewhere other than NCBI-PD.
The search is recorded because **three of the four obvious candidates fail for reasons that generalise**,
and those reasons are cheaper to reuse than to rediscover.

**Outcome: BV-BRC worked.** See [the result](gentamicin_rmt_bvbrc_2026-09-03.md) — 67 susceptible carriers,
organism-stratified, controls survived.

---

## The requirement, stated before searching

`rmt` sits at roughly 1% prevalence, so a usable archive needs **all** of:

1. isolate-level, **measured** gentamicin phenotype (not aggregate, not model-predicted);
2. a linked genotype or a fetchable genome;
3. **broad ascertainment** — isolates not selected on aminoglycoside resistance;
4. not a feeder into NCBI-PD;
5. large enough that ~1% prevalence yields a usable carrier count.

Requirement 3 is the one that does the eliminating, and it was not obvious in advance.

## Eliminated, with the reason

**`rmt`-enriched surveillance studies — fail (3), fatally.** The corpora richest in `rmt` are assembled by
screening on high-level aminoglycoside resistance. Real inclusion criteria: *"MIC > 256 mg/L to both
amikacin and gentamicin, then multiplex PCR"* (Bulgarian genomic surveillance, 150 RMTase carriers from
10,731 isolates); *"isolates growing on Mueller-Hinton with 256 mg/L amikacin"* (UK prospective study, 79
carriers); *"37 of 44 amikacin-resistant isolates"* (Wenzhou, 680 *E. coli*). **Every one conditions on the
outcome whose exceptions we are hunting**, so each contains zero susceptible carriers by construction. The
corpora most enriched for the determinant are precisely the ones structurally incapable of answering the
question.

**NARMS — fails (4).** The obvious independent national archive with isolate-level MICs and routine WGS.
But NARMS *"genomic sequences are uploaded to the NCBI Pathogen Detection web portal on a weekly basis"* —
it is a **feeder** to the archive already exhausted, not an alternative to it. Checking data-flow direction
costs one search and saves a full fetch-and-analyse cycle.

**ATLAS / Vivli — fails (2).** ~634,000 isolates, MICs against a median of 11 antibiotics, 20 years, 89
countries: the largest broad-panel MIC resource anywhere. But its genotype content is presence/absence for
**a limited number of β-lactamase genes only** — no `rmt`. Phenotype-rich and genotype-poor, so the join
cannot be made at all. (Access is also request-gated, but that never became the binding constraint.)

## The field-level pattern, which is our own label wall from outside

The resources split cleanly and do not intersect:

| | phenotype | genotype |
|---|---|---|
| ATLAS / Vivli | ~634k isolates, full MIC panels | β-lactamase presence only |
| AllTheBacteria / EnteroBase | — | 660k–1.1M assembled genomes |

A survey of >6,500 *A. baumannii* strains found only ~10% of public datasets carry the metadata needed to
apply CLSI/EUCAST breakpoints at all. The literature's own summary of the gap: for a genuinely unselected,
thousands-scale, genome-plus-MIC-panel dataset, **nothing public currently fits**.

That is this project's label wall, described independently by the field. It is worth recording that the
wall is structural rather than a failure of our searching.

## What worked, and why

**BV-BRC** is different on **both** axes, which is what makes it a second opinion rather than the same data
re-served:

- **Phenotype** — BioSample antibiograms *plus* ~300 hand-curated publications, per-row `pmid`,
  `laboratory_typing_method`, `testing_standard`. **49,553 measured gentamicin records**
  (`evidence = Laboratory Method`), against 20,816 labelled isolates in the whole PD sweep. The
  publication-curated portion is content PD's `AST_phenotypes` field does not hold.
- **Genotype** — `sp_gene` is **CARD/BLAT**, not AMRFinder: a different caller with different thresholds
  from both our pipeline and PD's.

Measured overlap: **7 of 48** resolvable carriers shared with the PD sweep → **162 new**.

## Also considered

- **ResFinder 4.0 validation set** — a convenience sample across six species with ECOFF-interpreted MICs;
  per-isolate accessions live in the *source* papers rather than the validation paper, so assembly is
  multi-hop. Not needed once BV-BRC worked.
- **"One Day in Denmark"** — 488 isolates, broth microdilution, genuinely unselected (all Danish clinical
  labs, one day). Ideal ascertainment, but ~1% prevalence gives ~5 carriers: underpowered.
- **"Two Weeks in the World"** (>3,000 unselected diagnostic isolates) and **MRSN** (3,878 *K. pneumoniae*,
  63 facilities) — both plausible; MRSN releases full AST only for a 100-isolate diversity panel. Left
  unexplored.

## Honest limits

- A scan, not a systematic review: ~8 searches and several fetches in one session, biased toward bacterial
  AMR.
- The eliminations are **structural**, so they hold regardless of how hard one looks — but "ATLAS has no
  `rmt` genotype" is a statement about the public release, not about what Pfizer holds.
- BV-BRC and PD both ultimately draw on public assemblies; the independence claim rests on the **measured**
  overlap, not on the archives being disjoint in principle.
