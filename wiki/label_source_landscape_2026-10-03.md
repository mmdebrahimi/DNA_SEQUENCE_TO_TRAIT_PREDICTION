# The current problem, the ideal label sources by type, and which are actually free (2026-10-03)

## Part 1 — What the problem actually is, derived from the ledgers

Read from the four active families' open candidates, not from memory:

| family | open candidate | blocker |
|---|---|---|
| `label-acquisition` | re-run the accrual sweep | **ran today → accruing zero.** Rate-limited by NCBI-PD ingestion, not by us |
| `learned-narrow` | "a 4th gene means a new LABEL, not a new estimator" | **label**, stated outright |
| `doubt-layer` | "Re-screen when any NEW independent label set lands" | **label**, stated outright |
| `catalog-curation` | restrict the gentamicin rescue to E. coli | **user authority** |

**Four of four converge on the same thing, and it is not a modelling problem.** The record already says
"the binding constraint is LABELS not models" — what this adds is the *shape* of label each blocked move
needs, because they are not the same label:

- **(A)** a **source-diverse bacterial AMR** cohort — unblocks today's Klebsiella meropenem specificity
  question, and would let `source_diverse_validate` reach non-E. coli organisms at all (only E. coli
  currently clears its bar).
- **(B)** a **non-HIV viral target-site** phenotype, free, quantitative, with a real negative class —
  unblocks `learned-narrow`'s one-VIRUS limit *and* `doubt-layer`'s screen, which are the same gap.
- **(C)** a **free fungal isolate-level phenotype** — the fungal cells are `NO_FREE_PHENOTYPE_SOURCE`.
- **(D)** post-lock accrual rate — **external**, unfixable by acquisition. Name it and stop pulling on it.

So (A)-(C) are acquisition questions and (D) is not. That is the actual decomposition.

## Part 2 — The immediate solve was already on disk

**`data/raw/ar_isolate_bank/kleb_details.json` holds 157 Klebsiella isolates with a CDC `Meropenem` INT
call: 104 R / 39 S / 14 I, of which 143 R/S carry a BioSample.** Measured, not inferred, offline.

Against today's prospective run:

| | prospective (2026-10-03) | AR Bank (on disk) |
|---|---|---|
| scorable | 16 (7R/9S) | **143 (104R/39S)** |
| panels/sources | 1 batch | **14, largest share 0.154** |
| susceptible class | 9, one batch | **39 across 10 panels, largest 0.308** |

It clears the repo's own ≥5-source / ≤0.60 bar with room, and the **susceptible** class — the one a
specificity claim rests on — is itself diverse. MICs are real measured values (`<=0.12`, `>8`).

This matters because the cell now has three disagreeing specificities:

| measurement | N | sens | spec | sources |
|---|---|---|---|---|
| calibration (the frozen rule's own `validated` string) | 30 | 1.000 | 0.733 | — |
| provenance-disjoint (report card) | 60 | 0.467 | **0.900** | 15, largest 37% |
| prospective (today) | 16 | 1.000 | **0.000** | **1** |

AR Bank would be the adjudicating fourth reading, and it is free, offline, and already fetched.

### The honest catch, which halves the opportunity

**134 of the 143 already appear in sibling AR Bank cohorts** (ceftriaxone 128, ciprofloxacin 105,
gentamicin 44). Only **9 are strictly disjoint**.

These are the *same isolates on different drugs*. The meropenem **label has never been used** and the
frozen meropenem rule was never fit on them — so this is not label leakage. But
`dna_decode/eval/cohort_manifest.py` is **accession-based and drug-blind** by design (`prior_accessions`
excludes only the exact-self cohort), so applied strictly it cuts 143 → 9 and the opportunity dies.

**My judgement, stated for ratification rather than asked:** per-drug independence is the right unit for
a per-drug rule, so the full 143 is scoreable with the genome-overlap **disclosed**, and the 9-isolate
strictly-disjoint set is the fail-closed floor that will be underpowered and should say so. Report both.
The drug-blind manifest convention is a real limitation worth recording either way.

**A second finding en route: zero AR Bank cohorts appear in the leakage manifest** (43 cohorts,
`incomplete=False`) even though the AR Bank arm clearly ran — 40 wiki artifacts. That is a **registry
gap**, not evidence of non-use, and it means the manifest cannot currently answer "was this isolate used
before?" for any AR Bank isolate.

## Part 3 — Ideal label sources BY TYPE

What each type must satisfy, written against the repo's own G1–G10 gates so a candidate can be screened
before any labour is spent (`uv run python scripts/screen_candidate_gates.py`).

| # | Type | Ideal shape | Gates that bite hardest |
|---|---|---|---|
| T1 | **Bacterial AMR, isolate-level** | measured MIC/AST + downloadable assembly, ≥20/class, ≥5 independent sources, ≤0.60 largest share | G2 study==class · G4 surveillance domination · G5 assembly attrition · G8 clonal collapse |
| T2 | **Viral target-site, quantitative** | per-isolate fold-change (not a rule-interpreted call), real negative class, non-HIV | **G1 circular label** (the catalog must not be the label's source) |
| T3 | **Fungal AMR, isolate-level** | AFST MIC + WGS, both classes populated, multi-clade | G6 censoring at the breakpoint · G4 |
| T4 | **Molecular typing, wet-lab gold standard** | agglutination / Quellung / full-locus typing, per isolate | G1 (must not be an in-silico call) |
| T5 | **Constructed-variation organism-level** | segregant / mutant panel — the regime that *works* (yeast 12/12, r 0.46–0.80) | none of G1–G8; this regime's constraint is availability |
| T6 | **Molecular forward (DMS)** | per-variant continuous fitness | G6 **assay degeneracy** (the CcdB trap: 79.3% tied at a ceiling) |
| T7 | **Human clinical** | consented genotype + measured phenotype | G3 sampling-defined · consent/DUA |
| T8 | **Conditional essentiality** | knockout × condition with **condition breadth** | not a label gate — the measured bottleneck is condition coverage |

The ranking that matters: **T1, T2, T3 are the live gaps.** T4 is partly served (Salmonella
agglutination, pneumo Quellung). T5/T6 are served (Bloom, ProteinGym). T7 is served for calling
(GeT-RM). T8 is screened and closed — "Public K-12" makes glucose dominance *worse*, not better.

## Part 4 — Deep search: what is actually free

### T2 — non-HIV viral target-site → **HerpesDRG is the strongest lead in this whole scan**

[HerpesDRG](http://cmv-resistance.ucl.ac.uk/herpesdrg/) ([paper](https://link.springer.com/article/10.1186/s12859-024-05885-5)):
a curated herpesvirus antiviral-resistance database carrying **EC50 fold-changes relative to wild-type**,
released **MIT-licensed on GitHub**. CMV is the one human herpesvirus with consensus guidelines binding
fold-change to R/S.

**Why this is pointed:** the project ships an **HCMV cell with 5 CLI-routable drugs** registered
`KNOWLEDGE_BASELINE` / abstention **`NO_FREE_PHENOTYPE`**. A free, quantitative, MIT-licensed
fold-change resource is exactly what that field says does not exist.

**Two caveats that decide its tier, and must be resolved before any number is quoted:**
1. It is **per-MUTATION, literature-curated — not isolate-level.** The HCMV contract's own note says
   HCMV phenotyping is per-mutation recombinant marker transfer. So this is the **CoV-RDB shape**, which
   this project has used before and labelled `..._IN_DISTRIBUTION_KNOWLEDGE_BASELINE` — not the
   independent isolate-level win HIV is.
2. **G1 risk is live.** If our HCMV catalog and HerpesDRG both descend from the same Chou/consensus
   literature, comparing them is a catalog-vs-catalog audit (`FAITHFUL_TO_TOOL`), not validation. **This
   must be screened, not assumed** — it is the one gate that would make the whole thing circular.

HBV / HCV / influenza: **no equivalent isolate-level repository.** Fold-changes are scattered across
individual phenotyping papers; [geno2pheno\[HCV\]](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0155869)
is a prediction *server*, not a downloadable phenotype table — and a server's output is a tool call, so
it trips G1 by construction.

### T3 — fungal → **two free sources, and they dwarf what we have**

1. [**CDC FKS1 echinocandin benchmark**](https://journals.asm.org/doi/full/10.1128/spectrum.03147-24) —
   **100 WGS *C. auris* isolates, 53 echinocandin-S / 47 echinocandin-R**, AFST MIC-categorised,
   purpose-built as a *FKS1* validation benchmark. Our existing AR Bank micafungin cohort is **32
   (6R/26S)** — this is ~8× the resistant class, and it targets the FKS1 cell directly.
2. [**Archived global *C. auris* population**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10821368/) —
   **387 isolates with both WGS and MIC across five clades** (152 I / 12 II / 119 III / 99 IV / 5 V).
   Our fungal G1 cohort is N=24 over clades I+III, so this adds both scale and clade breadth.

Honest limit: **no single downloadable resource pairs MIC+WGS across *Candida* and *Aspergillus***. The
assembly path is SRA reads joined to supplementary MIC tables. *Aspergillus* stays sparse —
[24 isolates](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12654709/) is the best single study found.
[ATLAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC9769544/) is the largest MIC surveillance but ships
**no genomes**, so it fails G5 outright.

### T1 — bacterial AMR

The AR Bank is the answer and it is already here (Part 2). Beyond it, the gap is not *more* AMR data but
**source-diverse** AMR data; NCBI-PD supplies breadth but its per-cell source concentration is measured
and often fails our own bar (3 of 10 SCORED cells rest on one BioProject).

## Part 5 — Ranked recommendation

1. **Score the frozen meropenem rule on the AR Bank Klebsiella cohort.** Free, offline, 143 isolates,
   14 panels. Adjudicates spec 0.900 vs 0.000. Disclose the sibling-drug genome overlap; report the
   9-isolate strict floor beside it. ~33 of the needed AMRFinder runs are already cached.
2. **Screen HerpesDRG against G1 before anything else.** If its curation is independent of ours, it is a
   quantitative phenotype for a cell that currently claims none. If not, it is a tool-tier audit — still
   worth doing, but it must not be called validation.
3. **Pull the CDC FKS1 100-isolate benchmark.** The cleanest labelled fungal set found, aimed at a cell
   we already ship, and ~8× our resistant class.
4. **Register every AR Bank cohort in `cohort_manifest`.** Cheap, and right now the leakage registry
   cannot answer a reuse question about 40 artifacts' worth of isolates.
5. **Stop pulling on post-lock accrual.** It is external. Today's sweep reproduced the eligible set
   exactly; the clock advances on NCBI-PD's schedule, not ours.

## Honest limits of this memo

- Parts 1–2 are **measured** from files on disk. Part 4 is a **web scan** — every source named needs its
  own G1–G10 screen before labour, and `screen_candidate_gates.py` refuses a verdict without a human
  judgement string for G1 and G3 precisely because those cannot be automated.
- I have **not** fetched or verified the contents of HerpesDRG or either fungal dataset. Record counts
  are as published; the repo's own standing rule is that a cited accession can be empty or wrong until
  resolved.
- The HerpesDRG/HCMV circularity question is **open**, and it is the difference between a validated cell
  and a restated one.
