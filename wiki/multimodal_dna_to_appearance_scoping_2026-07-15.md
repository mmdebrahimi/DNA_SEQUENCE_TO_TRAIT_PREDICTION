# Scoping the multimodal "DNA → appearance" moonshot (2026-07-15)

The project's stated long-term north star is a multimodal genotype→phenotype platform expanding to
**image-paired phenotype** ("DNA → animal/appearance image"). This memo scopes whether that frontier is
reachable, against the project's hard-won constraint: **learned genomic representations capture lineage, not
mechanism** (embeddings 0-for-5 across the kingdom boundary), so any path that *requires* a learned
representation is walled by prior evidence, not just by data. A `/research` probe (executed by hand), no build.

## The split the pre-bar framing forced

"DNA → appearance" is two fundamentally different problems, and they have opposite verdicts:

### (1) DNA → IMAGE (pixel generation) — **WALLED** (doubly)

- **Architecture wall.** Generating pixels from a genome is inherently a *learned-representation* problem —
  there is no curated "genome → pixel" catalog to apply deterministically. That is the exact regime the
  project's 0-for-5 embedding record already closed. Worse, the forensic-genomics literature finds **facial
  SHAPE prediction AUC < 0.6** even with dedicated panels (VISAGE) — structural/morphological *shape* does not
  decode from genotype at useful accuracy, independent of model. So the pixel path is walled by *both* the
  project's own negative and the field's.
- **Data wall.** No free, at-scale, genome↔face-image paired dataset exists (human faces are privacy-walled).
  The genotype+image datasets that DO exist are **microbial colony morphology** (Scientific Data 2026: 19
  species/151 strains; Pseudomonas PLOS 2023: 69 strains; yeast deletion-library images via TAMMiCol) — and
  those are the *reverse* direction (identify strain FROM image) via learned CV, not genome→predict-image.

  **Verdict: do NOT pursue pixel generation.** It is the moonshot's un-reachable half — walled by evidence,
  not just effort. Scaling compute won't fix a signal-vs-structure + no-catalog problem.

### (2) DNA → VISIBLE TRAIT (discrete phenotype call) — **REACHABLE + SCORABLE**

The determinism-compatible reframe: don't generate the image — **decode the catalogued visible traits a
genome implies** (eye / hair / skin pigmentation). This is textbook visible-trait genetics, and it is the
*exact* paradigm the project has already proven (a small curated marker set + published coefficients → a
phenotype call, identical in shape to the human PGx cells that already read defining-SNPs from a VCF).

- **Free deterministic model:** **HIrisPlex-S** — predicts eye/hair/skin colour from 41 SNPs via *logistic
  regression with published, open-access β-coefficients* (Forensic Sci Int: Genetics 2018 + supplementary +
  free online tool `hirisplex.erasmusmc.nl`). Eye AUC ~0.9+ (blue/brown), hair strong, **skin lower (0.80
  intermediate / 0.999 white)**. It is a curated-coefficient model — NOT a learned genomic embedding — so it
  sits in the project's *winning* regime (curated-catalog), not the closed one.
- **Free INDEPENDENT validation cohort:** **openSNP** — 6,401 open genotype files × 668 self-reported
  phenotypes incl. eye/hair colour (a commonly-used labelled eye-colour subset: 806 people, 404 blue-green /
  402 brown). Downloadable (opensnp.org / `rsnps` R pkg / the QC'ed `ofrei/opensnp` archive). This is a free,
  observed-phenotype, genotype-paired cohort → the cell can be **SCORED** (sens/spec/AUC on held-out people),
  not merely re-apply published coefficients.
- **Infrastructure already exists:** the `dna-decode pgx` decoder already reads defining-SNPs from a GRCh38
  VCF → applies a coefficient/rule → emits a trait call. A pigmentation cell is the same pipeline with the
  HIrisPlex-S coefficient table.

  **Verdict: this is the reachable form of the "DNA → appearance" north star** — deterministic, free labels,
  free coefficients, reuses existing infra, clears every negative-results-map gate (free + independent +
  observed + not-tool-derived + deterministic).

## Honest caveats (load-bearing — the disciplines that bit AMR apply here)

1. **openSNP labels are self-reported + DTC-chip + USA-population-skewed** → noisy labels + SNP missingness +
   population confound. The "suspect the label / high-sens-low-spec" + per-record-provenance disciplines apply;
   expect to validate on the cleanest binary subset (blue vs brown) first and report the label-noise ceiling.
2. **Validate the wrapper vs the underlying tool.** Applying HIrisPlex-S coefficients IS HIrisPlex-S — the
   value is the *decoder-tool integration + an independent-openSNP number*, not a new model. Headline the
   independent-cohort delta, and never claim novelty over the published model (mirrors the PGx PHENOTYPE
   faithful-to-CPIC posture + the validate-wrapper-vs-underlying-tool lesson).
3. **Application-context sensitivity.** Forensic DNA phenotyping carries governance baggage (ENFSI restricts
   it to severe crime / unidentified remains; public unease). The **clean framing for this project is a
   visible-trait *genetics* decoder** (the same benign SNP→trait genetics as the PGx cells / a textbook
   eye-colour example) — explicitly NOT a forensic/surveillance tool. Ship it with that scope rail, or don't
   ship it. This is a genuine authority call, not a technical one.

## Recommendation

- **Kill the pixel-generation moonshot** — walled by the project's own 0-for-5 + facial-shape AUC<0.6 + no
  free paired data. It is not a compute problem.
- **The reachable north-star step is a visible-trait pigmentation decoder cell** (`dna-decode pigment`, or a
  visible-traits module): HIrisPlex-S deterministic coefficients from a VCF, scored on the openSNP eye-colour
  subset. It is buildable now, no money, no learned embedding, methodologically identical to the PGx cells,
  and it is the *first honest instance* of "DNA → appearance" in the project's proven paradigm.
- **It is a build DECISION, not an auto-spillover:** it lands in DNA-11's human-SNP-from-VCF lane (coordinate)
  and carries the application-context authority call above. Surfaced for the user to choose.

Sources: HIrisPlex-S — [FSI:Genetics 2018 developmental validation](https://www.sciencedirect.com/science/article/abs/pii/S1872497318302205) · free tool `hirisplex.erasmusmc.nl`; VISAGE facial-shape AUC<0.6 — [Predicting Physical Appearance review, PMC12841266](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12841266/); openSNP — [opensnp.org](https://opensnp.org/) + [rsnps / rOpenSci](https://docs.ropensci.org/rsnps/) + [ofrei/opensnp](https://github.com/ofrei/opensnp); genotype+colony-image datasets — [Scientific Data 2026](https://www.nature.com/articles/s41597-026-07095-5), [Pseudomonas PLOS Comp Biol 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011699).
