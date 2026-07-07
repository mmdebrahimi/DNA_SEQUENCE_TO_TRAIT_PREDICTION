# PGx decoder on REAL PGP-UK humans — first real-people deterministic decode (2026-07-05)

**Milestone (F4, teed up by the mosfaer/Soraya-V3-5 session's handoff → this DNA-11 session):** the
deterministic PGx decoder now runs **end-to-end on real, named, open-consent human individuals** from a
FREE, non-application-gated cohort — no UK Biobank, no committee, no GPU, no money.

## What ran
- **Cohort:** PGP-UK (Personal Genome Project UK), ENA `PRJEB17529` — open-consent, fully-public individual
  WGS VCFs. **GRCh37/hg19**, called (GATK) `*.pass.recode.vcf.gz`, ~90–100 MB each.
- **Individuals decoded (N=3):** `FR07961000`, `FR07961006`, `FR07961007` (pilot deep-WGS batch `ERZ389…`).
- **Decoder:** the shipped `dna_decode.pgx` cells (8 genes) via a new **targeted position-liftover** —
  `scripts/pgx_decode_pgp_uk.py`. rsIDs are build-stable, so instead of lifting a 100 MB whole-genome VCF,
  it extracts only the ~22 PGx sites at their GRCh37 positions (resolved once via Ensembl GRCh37 REST,
  hardcoded + provenance-stamped), relabels to the GRCh38 catalog, and runs the normal callers. No liftOver
  binary, no CrossMap, no chain file, no dep-install.

## Result — the decoder discriminates real inter-individual variation

| sample | CYP2C19 | CYP2C9 | CYP3A5 | TPMT | CYP2B6 | VKORC1 | SLCO1B1 |
|---|---|---|---|---|---|---|---|
| FR07961000 | *1/*2 **IM** | *1/*1 NM | *1/*6 IM | *1/*1 NM | *1/*6 IM | G/G (normal) | T/T (normal) |
| FR07961006 | *1/*2 **IM** | *1/*2 **IM** | *1/*1 NM | *1/*1 NM | *1/*6 IM | A/A (**high sens**) | T/T (normal) |
| FR07961007 | *1/*17 **RM** | *1/*1 NM | *1/*1 NM | *1/*1 NM | *1/*6 IM | A/A (**high sens**) | T/T (normal) |

**The decoder reads real biology, not plumbing:** CYP2C19 spans Intermediate → **Rapid** Metabolizer
(FR…007 carries the *17 increased-function allele); CYP2C9 catches a *2 carrier (FR…006 IM); VKORC1
discriminates normal vs **high warfarin sensitivity** (G/G vs A/A). Every call is a valid CPIC diplotype +
metabolizer phenotype. All 3 CYP2B6 *1/*6 calls were **raw-genotype spot-verified** (each genuinely `0/1` at
rs3745274) — real, not an artifact (uncommon-but-genuine in a small cohort).

## Honest tier (load-bearing — do NOT inflate)
This is a **DEPLOYMENT / robustness demonstration**, NOT a new accuracy number:
- The decoder runs on an **independent real-world cohort** (PGP-UK, a *different* sequencing pipeline and
  genome build than the 1000G/GeT-RM benchmark) and produces sane, discriminating calls with the same honest
  abstention. That is the added value: it works outside the benchmark it was tuned/validated on.
- PGP-UK ships **no GeT-RM/CPIC truth labels**, so this is **NOT** an accuracy-vs-truth measurement. The
  accuracy validation lives in the GeT-RM concordance cells (`wiki/pgx_report_card.md`: CYP2C19 72/72,
  CYP2C9 73/73, CYP2C8 82/82, CYP3A5 8/8, TPMT 85/85, CYP2B6 62/62).
- **GRCh37→38 assumption:** a called WGS VCF omits reference sites, so a PGx site absent-as-variant on a
  chromosome the VCF covers is encoded homozygous-reference (dna-pgx flags it `assumed_reference`). Sites on
  a chromosome the VCF never touches are omitted (→ honest `no_input`, never a fabricated ref call).
- **NOT a clinical tool** (every cell carries this rail).

## Why it matters
UK Biobank died (externally walled, `dd4fa92`); PGP-UK is the **sole surviving open route to
individual-level human genomes** without an application/committee. This turns "validated on a benchmark" into
"runs on real people from an independent source" — the reachable real-human frontier.

## Artifacts
- Runner: `scripts/pgx_decode_pgp_uk.py` (reusable; `--vcf <PGP-UK GRCh37 vcf> --sample-id <id>`).
- Per-individual provenance JSON: `wiki/pgp_uk_pgx_results/FR0796100{0,6,7}_pgx.json`.
- Cohort + gotcha handoff (mosfaer): `wiki/pgp_uk_realization_handoff.md` (git `00cd99c`).
- Big VCFs live on `D:/dna_decode_cache/pgp_uk/` (gitignored-class; regenerable from the ENA URLs above).

**Non-duplication note (R4):** PGP-UK was already used for the eye-colour + ABO *trait* cells (2026-06-30);
this is the first **PGx** decode on PGP-UK — the specific gap the mosfaer session handed to DNA-11. Additive,
not a redo.
