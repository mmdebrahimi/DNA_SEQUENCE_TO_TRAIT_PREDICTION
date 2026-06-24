# Klebsiella K-antigen capsule typing — label-gate feasibility (2026-06-24)

> Executes recs 1+2 of the non-AMR pivot (`plans/Non_AMR_Phenotype_Pivot_Assessment_2026-06-24.md`): run
> the LABEL GATE for the first GREEN candidate BEFORE any build. Verdict: **GO** — both gates clear (rare for
> a non-AMR trait). The build is a real new-data cell -> ratify-worthy, not an auto-increment.

## The label gate (the project's binding test)

| Gate | Result | Evidence |
|---|---|---|
| **1. Curated determinant catalog?** | **YES** | **Kaptive** — free + open K-locus (KL) reference DB, the standard genomic capsule typer. Kaptive 2.0 added 16 novel KL types from >17,000 KpSC genomes; **>160 KL types** defined. (`github.com/katholt/Kaptive`; Microbial Genomics `mgen.0.000800`). Exact shape of the shipped `serotype`/`plasmid` callers (blastn vs a curated locus DB). |
| **2. Free, independent, MEASURED, isolate-level validation label?** | **YES** (newly available) | **KlebNET-GSP Technical Report 2025 (Zenodo 15742130): 731 isolates with matched SEROLOGICAL + genomic K-type calls, 84.5% concordant.** Plus the original Kaptive paper: 86 genomes with matched serological typing. The **serological K-type is a wet-lab MEASUREMENT** (the Kauffmann ~77-K-type scheme), independent of the genomic Kaptive call -> a non-circular genotype->measured-phenotype validation, exactly the AMR-grade label shape. |

**Both gates clear** -> K-antigen capsule typing is a **GREEN-VALIDATABLE** cell (not merely faithful-to-tool
like the existing typing callers). This is the cleanest non-AMR validatable candidate the project has found.

## Honest caveats (must ship with any build)
- **Serological non-typeability: 10-70% of isolates** are serologically non-typeable (novel capsule or
  non-capsulated) -> the 731-isolate validation set is BIASED toward typeable isolates; spec on
  non-typeables is a known blind spot.
- **>160 genomic KL types vs ~77 serological K-types** -> the genotype->serotype map is many-to-fewer +
  INCOMPLETE; many KL loci have no serological equivalent. A v0 must report KL-type (genomic) + the
  serological K-type ONLY where the mapping is established.
- **The validation compares genomic-prediction vs serology** (84.5% concordant) -> our caller, if
  faithful-to-Kaptive, would inherit ~that ceiling; the honest number is "genotype vs measured serology,"
  reported with the non-typeable caveat. NOT an embedding/learned claim (that stays a closed negative).

## Build complexity (why this is a real cell, not a quick increment)
Two implementation shapes:
- **(a) Faithful-to-Kaptive (locus-level).** WRAP Kaptive (a NEW external dependency + its KL DB) -> highest
  accuracy; mirrors the `vf_runner`/`resfinder` faithful-to-tool pattern. Cost: a new dependency on a
  disk-tight host + the KL reference DB.
- **(b) wzi-single-gene v0 (lighter).** `wzi` allele -> K-type prediction (Kp BIGSdb) is a SINGLE-gene blastn,
  the EXACT shape of the shipped `serotype` caller (wzx/wzy/fliC). No Kaptive dep (blastn already installed).
  Lower accuracy than full locus typing; the honest "smallest credible slice."

Either is a real new-DATA cell (a new reference DB to fetch + a validation cohort of genomes to download +
the 731-isolate serological label set). Comparable in effort to the AMR-Portal independent-validation arm,
NOT a 20-minute increment.

## Verdict: **GO** — but the build is RATIFY-WORTHY (new dependency + cohort), per rec 1's "before any build"
The label gate (rec 1's explicit precondition) CLEARED with a strong result. The build itself is a genuine
new-external-data cell:
- new reference DB (Kaptive KL or wzi-allele) to fetch,
- a validation cohort (the 731-isolate paired set + their genomes) to download (disk-tight host -> stage on D:),
- a faithful-to-tool caller + CLI + offline-safe degrade + a reference-control test (serotype pattern),
- a validation run scoring genotype vs measured serology, reported with the non-typeable caveat.

This is `/idea-anchor`-class (a new dependency-accepting, multi-step, ~AMR-Portal-effort cell), not an
auto-increment. **Recommended: greenlight a focused `--until-mvp` on the K-antigen cell** (I'd default to
shape (b) wzi-v0 first — lightest, serotype-pattern, no new heavy dep — then optionally (a) Kaptive-wrap for
accuracy), with the cohort staged on D: per the disk discipline.

## Scoped build plan (when greenlit)
1. Pick shape: (b) wzi-v0 (default) or (a) Kaptive-wrap (accuracy). [user/judgement]
2. Fetch the reference DB (wzi allele set + wzi->K-type map, OR Kaptive KL DB). [download, stage on D:]
3. `dna_decode/ktype/` caller (blastn vs the DB; faithful-to-tool; offline-safe), mirroring `serotype`.
4. CLI `dna-ktype` / `dna-decode ktype` + a reference-control test (a known-K-type genome).
5. Validation: score vs the KlebNET-GSP 731-isolate serological set (fetch genomes -> call -> concordance),
   report genotype-vs-measured-serology + the non-typeable caveat + a NEW namespace-separate report card
   (`ktype_report_card`, NOT keyed into the frozen AMR/canonical_cell_key card).
6. Trust-surface: a `ktype` tier (faithful-to-tool + measured-serology-validated), NOT an AMR cell.

## Sources
- Kaptive 2.0: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000800 ; code+DB https://github.com/katholt/Kaptive
- KlebNET-GSP 2025 paired serological+genomic validation (731 isolates, 84.5%): https://zenodo.org/records/15742130
- Original Kaptive (K-locus DB, 86 serologically-typed genomes): https://pmc.ncbi.nlm.nih.gov/articles/PMC5359410/

---

## Rec 3 — the phenotype question + the productization fork (resolved)
- **The productization fork is ALREADY RESOLVED:** Path A (packaging gate) shipped this session (commit
  `4b50ae1`) -- the wheel ships the trust cards; a `pip install`ed `dna-amr`/`dna-decode` serves correct
  badges from the artifact. The deterministic decoder is now a genuinely pip-installable, honest tool.
- **The only remaining productization item is PyPI publish** -- a publish/outward step gated on USER
  authority (the wheel is publish-READY; recommend a TestPyPI dry-run first).
- **The phenotype question is NOT banked away as "no" -- it has a concrete GREEN GO** (K-antigen above). The
  decoder is trait-general; K-antigen is the first non-AMR cell that clears the FULL validation gate (free
  measured label), making it a genuinely worthwhile build when ratified.

## Bottom line (all 3 recs resolved)
1. **Build a GREEN cell:** label gate RUN -> **GO** for Klebsiella K-antigen capsule typing (both gates clear;
   the rare free-measured-label win for a non-AMR trait). Build is ratify-worthy (new dep + cohort); scoped
   plan above; surfaced for greenlight.
2. **/research the label availability:** DONE -> the free measured validation set EXISTS (KlebNET-GSP 731
   isolates) -- this is what flips K-antigen from faithful-to-tool to GREEN-VALIDATABLE.
3. **Bank + productization:** the fork is RESOLVED (Path A done); PyPI is the only open item (user-authority).
