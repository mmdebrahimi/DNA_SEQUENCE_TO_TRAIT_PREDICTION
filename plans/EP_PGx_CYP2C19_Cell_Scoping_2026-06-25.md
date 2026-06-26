# EP-PGx — human pharmacogenomics CYP2C19 cell (decompose + plan, 2026-06-25)

The first **human** decoder cell + the honest form of the "higher species" jump (catalog-tractable trait in a
eukaryote; the complex-trait/embedding path stays a closed negative). Decomposed + feasibility-gated here;
build is a scoped next step.

## Decompose — the cell (deterministic-pattern fit, new organism)
`genome/VCF → detect defining SNPs → assign CYP2C19 star-alleles → diplotype → CPIC metabolizer phenotype`.
Same proven variant→catalog→call pattern as AMR point-mutations / the TB WHO catalogue, on human DNA.

## Phase 1 — feasibility probe: **GO** (2026-06-25)
All three free-data dependencies confirmed (the cleanest "higher organism" cell — has a free INDEPENDENT panel):
| Dependency | Status | Source |
|---|---|---|
| Star-allele definitions | ✅ free download + REST API, version-pinnable | PharmVar `pharmvar.org/genes` (CYP2C19) |
| Diplotype→phenotype table | ✅ free standardized files | CPIC/PharmGKB `pharmgkb.org/page/cyp2c19RefMaterials` |
| **Independent validation panel** | ✅✅ GeT-RM consensus diplotypes ⋈ public 1000 Genomes VCFs by Coriell NA/HG id | Gaedigk 2022 *J Mol Diagn*; PharmCAT validated **59/59 CYP2C19 concordant** |

- **Input contract:** phased VCF, GRCh38. Variant-based caller (rsID/position match).
- **v0 scope:** CYP2C19 (SNP-defined). Core alleles *2 (rs4244285), *3 (rs4986893), *17 (rs12248560) + the
  common set; CYP2D6 DEFERRED (copy-number + hybrid genes — notoriously hard).

## Honest tier (load-bearing)
- **Star-allele CALLING = INDEPENDENTLY validatable** vs the GeT-RM consensus panel (a real measured/consensus
  label, free) → the genuine number.
- **Phenotype (PM/IM/NM/RM/UM) = FAITHFUL-TO-CPIC** (assigned from diplotype via CPIC's table, not measured) —
  like the serotype/ktype cells. A truly-measured phenotype (probe-drug PK) is not free at scale.
- Reference baseline = **PharmCAT** (the gold-standard tool); our cell is a deterministic PharmCAT-style caller,
  validated vs GeT-RM. caller_is_independent_baseline = False for the phenotype step.

## Known caveats (pin before building)
- **`*38` is the true reference allele** (no variants vs GRCh38), NOT `*1` (`*1` has rs3758581 chr10:94842866 A>G).
- Some allele-defining variants are MISSING in the NYGC 30× VCFs (*16/*24/*30) — use Phase-3 10× OR handle absence.
- NA19122 is a known *2/*35-ambiguity sample (phasing) — expected, not a defect.
- Homozygous-reference confidence varies by VCF source — handle no-call explicitly (never silent ref-by-absence).

## Build plan (Phase 2–4 — the scoped next step)
- **P2 catalog ingestion:** PharmVar CYP2C19 allele defs (rsID→star-allele) + CPIC diplotype→phenotype table →
  a gitignored DB + a build script (reproducible from the pinned versions).
- **P2 caller:** `dna_decode/pgx/` — VCF → per-allele defining-variant match → star-allele set → diplotype
  (two haplotypes) → CPIC phenotype. Offline-safe degrade; explicit no-call.
- **P3 validate:** run on GeT-RM Coriell samples (their public 1000 Genomes VCFs) → diplotype vs consensus →
  the independent calling number. Report card (calling independent / phenotype faithful-to-CPIC).
- **P4 ship:** `dna-pgx` console script + `dna-decode pgx` dispatcher + pyproject + CHANGELOG + tests.

## Success bar / falsifier
- **Bar:** caller reproduces ≥1 GeT-RM consensus CYP2C19 diplotype correctly (MVP), → cohort concordance on the
  GeT-RM-overlapping 1000 Genomes samples (the real number; PharmCAT's ceiling is ~59/59).
- **Falsifier:** if our calls diverge materially from GeT-RM consensus on the common *2/*3/*17 set, the caller
  is wrong (not a label problem — GeT-RM is consensus truth).

## Why this is the right "higher species" move
- It's the honest catalog-tractable form of the eukaryote jump (the complex-trait/embedding path is a closed
  0-for-4 negative — NOT reopened here).
- Uniquely, it has a FREE INDEPENDENT validation panel (GeT-RM) — better-grounded than most bacterial cells.
- FROZEN E. coli AMR surface untouched (new `dna_decode/pgx/` package; non-frozen).

## Status: BUILT (P2–4 shipped 2026-06-25, v0.6.0)

> **UPDATE 2026-06-25 — the cell is BUILT.** `dna_decode/pgx/` ships the full VCF→star→diplotype→CPIC-phenotype
> caller + `dna-pgx` console script + `dna-decode pgx` dispatch + 23 tests (all pass). Coordinates grounded
> (GRCh38 chr10: *2=94781859, *3=94780653, *17=94761900); CPIC function + diplotype→phenotype table encoded
> (Caudle 2020). Honest tier wired into the record (`calling_is_independent_baseline=True` /
> `phenotype_is_independent_baseline=False`). Edge handling verified live: phased/unphased trans assumption,
> multiallelic, no-call, assumed-reference-at-absent (all flagged, never silent). FROZEN AMR surface
> byte-unchanged (leak guard green). Released as v0.6.0; CHANGELOG updated.
>
> **P3 PARTIAL DONE 2026-06-25 — faithful-to-PharmCAT number landed; GeT-RM independent number is a tooling wall.**
> `scripts/pgx_cyp2c19_validate.py` validates the caller against the reference tool's OWN real VCF fixtures
> (PharmCAT, vendored under `tests/data/pgx_cyp2c19/`): **core diplotype 6/6 + phenotype 6/6**
> (`wiki/pgx_cyp2c19_report_card.{md,json}`). Honest tier = FAITHFUL-TO-PHARMCAT (in-distribution; the
> "validate-wrapper-vs-tool" discipline), NOT independent. Blind spots surfaced: *35→*1/*1, *4b→*1/*17
> (shares the *17 SNP — clinically meaningful). 14 harness tests.
>
> **P3 REAL-1000G RUN DONE 2026-06-25 via Docker (the "tooling wall" was dissolved — Docker IS on this host).**
> `scripts/pgx_1000g_population.py` ran the caller on the REAL 1000 Genomes 30x panel (n=3202), region
> fetched via Docker bcftools remote-querying the public VCF (htslib via the `quay.io/biocontainers/bcftools`
> image — no native htslib needed). Results (`wiki/pgx_1000g_population_2026-06-25.{md,json}`):
> biologically-sane population phenotype distribution (NM 38.3% / IM 32.4% / RM 19.5% / PM 6.7% / UM 3.2%);
> blind-spot exposure QUANTIFIED on real data (*35→*1 mis-call 44/3202=1.37%, *4-family 5/3202=0.16%);
> and a GROUNDED GeT-RM data point — NA19122 (consensus *2/*35) → v0 calls *1/*2, genotypes confirm the
> *35 haplotype is real + invisible to v0 (blind spot confirmed end-to-end on the real genome).
>
> **REMAINING (the full GeT-RM consensus concordance % — wall reclassified: data-access, NOT tooling):**
> the per-sample GeT-RM consensus labels live in paper supplements (Gaedigk 2022 Table S; ursaPGx S1), not
> a clean public TSV. `pgx_cyp2c19_validate.py --source getrm --expected-tsv` consumes them once extracted.
> v0.1 refinements: the sentinel layer (rs28399504/*4 + rs12769205/*35 -> withhold, per the brainstorm),
> *1-vs-*38 (rs3758581=chr10:94842866, verified). PyPI publish of 0.6.0 is a separate (gated) step.
