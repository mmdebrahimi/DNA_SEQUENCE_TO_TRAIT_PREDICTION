# What the rest of the field is building, and the five methods worth stealing

A deliberate scan of tools that do what this project does — genotype → phenotype decoding with an honest
trust surface — read for **borrowable methodology**, not for a survey. One borrowed method is applied to a
live open question at the end and it materially strengthens a committed claim.

**Headline for orientation:** an ESCMID working group of ~200 curators is building the *same layer we
built* (AMRrules), and an independent 78-dataset benchmark **reproduces our central finding** that
catalogs beat ML on divergent strains. We are not alone, we are not redundant, and two of our design
choices are ahead of the deployed tools while three of theirs are ahead of ours.

---

## 1. The closest sibling: AMRrules / ESGEM-AMR

[AMRrules](https://github.com/AMRverse/AMRrules) is organism-specific interpretation of AMR **genotypes**,
explicitly modelled on how EUCAST/CLSI standardise interpretation of AST **phenotypes**. Users run
AMRFinderPlus, then AMRrules turns the determinant list into an S/I/R call. That is, precisely, our
`amr_rules.DRUG_RULE` layer — built by an ESCMID working group (convenors Kat Holt, Natacha Couto, Jane
Hawkey), partnered with the EUCAST WGS subcommittee, ~200 members in organism subgroups.

**Read the rule table, not the description.** Each rule is one row with **24 fields**
(`rules/Escherichia_coli.tsv`, fetched 2026-09-03):

```
ruleID · txid · organism · gene · nodeID · protein accession · HMM accession ·
nucleotide accession · ARO accession · mutation · variation type · gene context ·
drug · drug class · phenotype · clinical category · breakpoint · breakpoint standard ·
breakpoint condition · PMID · evidence code · evidence grade · evidence limitations ·
rule curation note
```

Controlled vocabularies actually in use across the 114 E. coli rules:

| field | values (count) |
|---|---|
| `gene context` | acquired (62) · **core (52)** |
| `phenotype` | nonwildtype (78) · wildtype (36) |
| `clinical category` | R (71) · S (38) · I (5) |
| `evidence code` | ECO ontology terms — `ECO:0001103 natural variation mutant evidence` (85), `ECO:0005027 genetic transformation evidence` (24), `ECO:0001113 point mutation phenotypic evidence` (3) |
| `evidence grade` | high (103) · moderate (8) · low (3) |
| `breakpoint standard` | `EUCAST v15.0 (2025)` (62) · `EUCAST v16.0 (2026)` (21) · **`ECOFF (May 2025)` (20)** · `ECOFF (January 2024)` (11) |
| `evidence limitations` | `-` (103) · "lacks evidence across diverse allelic backgrounds" (5) · "statistical geno/pheno evidence but no experimental evidence" (3) · "lacks evidence for this drug" (3) |

**Scope, measured not assumed.** Total ~300 rules across the six organism files checked, and coverage is
early: E. coli is **β-lactams only** — no ciprofloxacin, no gentamicin, no tetracycline rules at all.
Klebsiella has **6** rules. **No `rmt`/`armA`/`npmA` rule exists in any of the six organisms.** So AMRrules
does not supersede this project on the drugs we actually cover, and it neither confirms nor refutes our
gentamicin v2 `rmt` rescue — the community catalog has not reached that determinant class yet.

## 2. The benchmark that independently reproduces our central finding

Hu et al., *Briefings in Bioinformatics* 25(3) bbae206 (2024) — 78 species-antibiotic datasets, four ML
methods (Kover, PhenotypeSeeker, Seq2Geno2Pheno, Aytan-Aktug) vs rule-based ResFinder, under **random,
phylogeny-aware, and homology-aware** fold splits.

- Random splits: F1-macro ≥ 0.9 on **79%** of species-antibiotic pairs. Phylogeny-aware: **60%**.
- **ResFinder — the catalog — is best on 44% of datasets under phylogeny-aware splitting and 50% under
  homology-aware**, versus 25% under random.

That is our regime map, derived independently and at far greater scale: **ML wins in-distribution,
the deterministic catalog wins as you move away from the training population.** It is the strongest
external corroboration this project has for shipping the catalog and closing the embedding track.

## 3. The five borrowables, ranked

### (a) SOLO grading — steal this first, and it is already applied below

The [WHO TB mutation catalogue](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(21)00301-3/fulltext)
grades a variant **only on the isolates where it appears alone**, then scores PPV + OR + FDR-corrected p,
and assigns one of five confidence groups. Grade 1 ("associated with resistance") requires **≥5 solo
occurrences and a PPV 95% CI lower bound ≥ 0.25**.

This is the published, adopted answer to the confound that **wrecked our NNRTI curation**: the
multivariate OLS deleted canonical `Y181C` because it co-occurs with `K103N`, which absorbed its
coefficient. SOLO handles the same confound by *exclusion* rather than *regression*, so it cannot silently
reassign a determinant's effect to its co-traveller. Counterpoint worth carrying: a *Nature Communications*
2025 analysis found multivariable penalised regression **beats** SOLO on sensitivity for TB, recovering
compensatory and hypersusceptibility variants solo-only grading misses. SOLO is the *conservative*
estimator — it buys an uncontaminated denominator by discarding data.

### (b) `gene context: core | acquired` as a first-class rule field

Our own recorded lesson — *intrinsic genes break broad AMR class rules* (Acinetobacter OXA-51, Klebsiella
OqxAB drive specificity to 0) — is a field in their schema, populated 52 core / 62 acquired for E. coli.
We rediscovered by failure what they encode by design.

### (c) ECOFF as the anchor, with the clinical breakpoint secondary

The EUCAST WGS subcommittee's stated principle: anchor genotype-phenotype comparison on **wild-type MIC
distributions / ECOFFs**, clinical breakpoints second. The logic is exact — acquiring a determinant moves
an isolate out of the wild-type distribution *even when the shift is not clinically meaningful*, so the
ECOFF is the cutoff that corresponds to "has a mechanism", which is what a genotype predicts. ResFinder
4.0's residual errors cluster exactly there: most discordant E. coli sat **one dilution above the ECOFF**.
Our `mic_tiers.py` carries CLSI/EUCAST **clinical breakpoints only**, and AMRrules versions its breakpoint
source per rule (including `ECOFF (May 2025)`) — we version nothing per entry.

### (d) VME / ME and susceptible-class metrics as the headline

Clinical microbiology reports **very major errors** (reference R → test S: the dangerous direction) and
**major errors** (S → R), against CLSI M52 / FDA criteria (CA ≥ 90%, VME and ME < 3%; FDA's real criterion
is CI-based, not the commonly-cited flat 3%). Hu et al. likewise report **F1-negative, precision-negative,
recall-negative** as the "clinically relevant" metrics. We report sensitivity/specificity — which contains
the same information but does not put the dangerous error in the headline. Our entire doubt layer exists
to qualify exactly the false-susceptible call that VME names.

### (e) Per-rule provenance is representable after all

The catalog-curation family recorded that per-entry sourcing is **"not representable"** because our
catalogs are bare `set[str]` with nowhere to hang a source. AMRrules is the existence proof that it is
representable, at production scale: `PMID`, `evidence code` (ECO ontology), `evidence grade`,
`evidence limitations`, `rule curation note`, per row. That conclusion of ours was true about our *current
schema*, not about the problem.

## 4. Where we are ahead — worth knowing, not just borrowing

**Mykrobe has no concept of "unknown significance": it predicts susceptible when no known resistance
mutation is found.** So does most of the field. Best-practice reviews for targeted-NGS TB now *recommend*
returning **indeterminate** for uncertain-tier variants or inadequate coverage, and explicitly warn that
susceptibility must not be inferred from missing data — which is our L2 doubt layer's founding argument,
already shipped and now firing on real calls (`V179F` → `DOUBT [strong]`, call stays S).

CARD RGI's Perfect/Strict/Loose tiering is the field's closest analogue to our evidence tiers, with a
cautionary tale attached: by default RGI **silently promotes** Loose hits ≥95% identity to "Strict" unless
`--exclude_nudge` is set, so a published "Strict" call may sit below the model's own threshold. That is the
same class of defect as our own placeholder and self-to-self traps — a tier that quietly means something
other than it says.

## 5. SOLO applied: it strengthens the `rmt` claim

Our committed hunt reported a **pooled** PPV of 146/206 = 0.709 for `rmt` → gentamicin-R. Pooled is the
wrong denominator by WHO's standard, because **47 of the 146 resistant carriers also carry `aac(3)`**, the
classic gentamicin-modifying enzyme — those isolates cannot tell you what `rmt` did.

`scripts/solo_ppv.py`, run against the committed artifact:

| stratum | R | S | PPV | 95% CI |
|---|---|---|---|---|
| pooled (as reported) | 146 | 60 | 0.709 | [0.643, 0.767] |
| **SOLO — no `aac(3)` co-carriage** | 99 | 21 | 0.825 | [0.747, 0.883] |
| co-carriage with `aac(3)` | 47 | 39 | 0.547 | [0.442, 0.647] |
| **SOLO, artifact project excluded** | **99** | **0** | **1.000** | **[0.963, 1.000]** |

**The honest evidence for the rule is stronger than what we published, and better controlled:** 99
isolates carrying `rmt` with no co-carried `aac(3)`, outside the label-artifact submission, every one
resistant — clearing WHO's grade-1 bar (≥5 solo, CI lower ≥ 0.25) by a wide margin. Note the middle row:
the co-carriage stratum sits at 0.547, which is where the pooled number's dilution came from.

This does **not** touch the open specificity question. Zero susceptible solo carriers outside the artifact
project is still an absence of counter-examples, not a bound.

## 6. Screened against our own gates

Every borrow above is a **method**, not a label source, so the label wall is untouched by all of it.
Running the candidates through `scripts/screen_candidate_gates.py`'s framework: AMRrules rule tables are
**curated interpretations, not measurements** — using them as a validation label would trip **G1
(circular label)** exactly as HBV did. They are a catalog *source* and a schema *model*; they are not a
phenotype. The one genuinely new label-adjacent lead is **EUCAST's public MIC distribution / ECOFF data**,
which is aggregate rather than isolate-level and so cannot serve as a per-isolate label either.

## 7. Honest limits

- **Scan, not a systematic review.** Nine searches plus five fetches across one session; selected for
  borrowable method, biased toward bacterial AMR because that is where our frozen surface lives. The
  viral/fungal/human tracks were not scanned.
- **AMRrules coverage was measured on six organism files on 2026-09-03** and it is under active
  development — the β-lactam-only E. coli finding will date quickly, and it would be wrong to cite it as
  a standing property of the project.
- **The AMRgen preprint was read through a summarizer**, not in full; its solo-PPV and Firth-regression
  claims are second-hand here.
- **Nothing was adopted in this pass** except the SOLO computation in §5. (b)–(e) are proposals with a
  named cost, not decisions — several would touch the frozen surface and none is authorized.
- **The Hu et al. corroboration is about the same *phenomenon*, not the same experiment**: their
  phylogeny-aware split is a training-time control, our lineage collapse is a reporting-time correction.
  They agree in direction; they are not interchangeable evidence.

## Sources

- [AMRverse/AMRrules](https://github.com/AMRverse/AMRrules) · rule tables read directly from `rules/*.tsv`
- [Hu et al. 2024, *Brief Bioinform* 25(3) bbae206](https://academic.oup.com/bib/article/25/3/bbae206/7665136) — phylogeny-aware AMR benchmark
- [WHO TB mutation catalogue, *Lancet Microbe* 2022](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(21)00301-3/fulltext) — SOLO grading
- [Ellington et al., EUCAST WGS subcommittee, *CMI* 2017](https://www.sciencedirect.com/science/article/pii/S1198743X16305687) — ECOFF-anchored principle
- [ResFinder 4.0, *JAC* 2020](https://academic.oup.com/jac/article/75/12/3491/5890997) — discordance clusters at the ECOFF
- [hAMRonization / PHA4GE spec](https://www.biorxiv.org/content/10.1101/2024.03.07.583950v1.full) — cross-tool output harmonization
- [AMRgen preprint, bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.05.01.722195v1.full) — solo-PPV + Firth regression + ECOFF tooling
- [CARD RGI docs](https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst) — Perfect/Strict/Loose and the nudge
- [Mykrobe, *Wellcome Open Res* 2019](https://wellcomeopenresearch.org/articles/4-191) — quality filters, no unknown-significance category
