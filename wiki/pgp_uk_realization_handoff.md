# HANDOFF → DNA-11 session: run the deterministic decoder on real PGP-UK people (2026-07-04)

> Teed up by the parallel (mosfaer) session. **Non-duplication checked:** DNA-11's last 30 commits +
> uncommitted work show ZERO PGP/personalgenomes/real-people activity (it's on the J2 Kaggle ESM runner +
> PGx report cards). So this is a genuinely new milestone, not a redo. Everything below is **verified
> real-surface** 2026-07-04, including one gotcha that would otherwise cost you a debugging cycle.

## The opportunity
**PGP-UK** = a FREE, open-consent, NON-application-gated, individual-level human cohort with **called VCFs
already deposited** in the open ENA archive (no dbGaP/EGA, no committee). The deterministic **PGx decoder**
(`dna_decode/pgx/cli.py`, `dna-pgx`) consumes a VCF **directly** → this is real-people validation of the
decoder without waiting on the UK Biobank application.

## THE GOTCHA (verified, load-bearing) — assembly, but liftOver is CONDITIONAL
- All 99 PGP-UK VCFs are **GRCh37 / hg19** — ENA `ASSEMBLY = GCA_000001405.14` (GRCh38 starts at `.15`),
  verified on both analysis batches (`ERZ389532`, `ERZ1065952`).
- **The PGx caller matches by (chrom,pos) FIRST, then FALLS BACK to rsID** — verified in
  `dna_decode/pgx/caller.py:101–129` (`by_pos.get((chrom,pos))` → else `by_rsid.get(vcf_ID_col)`). The
  defining CYP2C19/CYP2C9 sites are all standard dbSNP rsIDs (rs4244285, rs4986893, rs12248560, …).
- **⇒ liftOver is CONDITIONAL, likely UNNECESSARY:**
  - **If the VCF `ID` column carries rsIDs** (common for GATK/dbSNP-annotated VCFs) → the rsID fallback
    matches despite the GRCh37 coords → **run `dna-pgx` DIRECTLY, no liftOver.**
  - **Only if the `ID` column is `.`** → the position path misses on GRCh37 → **liftOver GRCh37→GRCh38 first**
    (free: UCSC `liftOver`+`hg19ToHg38.over.chain`, `CrossMap.py vcf`, or `bcftools +liftover`).
- **2-second check at download time:** `zcat FR07961000.pass.recode.vcf.gz | grep -v '^#' | head`
  — is column 3 an `rs…` id? **Yes → skip liftOver.** `.` → liftOver.

## Verified decoder command (from `dna_decode/pgx/cli.py`)
```
dna-pgx <sample.vcf> --gene cyp2c19            # CYP2C19 diplotype + CPIC metabolizer phenotype (default)
dna-pgx <sample.vcf> --gene cyp2c9             # CYP2C9 activity-score
dna-pgx <sample.vcf> --gene vkorc1             # VKORC1 warfarin sensitivity (single-SNP)
dna-pgx <cohort.vcf> --sample <COL> --json-only --out prov.json
```
Phased VCF preferred (unphased → `phase_ambiguous` flag, still calls); pure-stdlib parse, no Docker.

## Download (free FTP, no gate) — index already captured
- Manifest + pointers on D: from the mosfaer session: `D:/dna_decode_cache/pgp_uk/pgp_uk_manifest.tsv` (via `scripts/capture_pgp_uk.py`).
- **99 called VCFs**, e.g.:
  - pilot (11, deep multi-omics): `ftp.sra.ebi.ac.uk/vol1/analysis/ERZ389/ERZ389532/FR07961000.pass.recode.vcf.gz`
  - broader WGS (~88): `ftp.sra.ebi.ac.uk/vol1/analysis/ERZ106/ERZ1065952/ERS1176551_noHLA.vcf.gz`
- Full list: ENA `result=analysis` for `PRJEB17529`. Put downloads on **D:** (per the heavy-artifacts rule), not C:.

## Suggested minimal run (done-criterion)
1. Download **one** PGP-UK VCF (~MB) to D:.
2. Check the `ID` column (`zcat …vcf.gz | grep -v '^#' | head`): rsIDs present → **skip liftOver**; `.` → liftOver GRCh37→GRCh38.
3. `dna-pgx <vcf> --gene cyp2c19 --json-only` → a real CPIC metabolizer call on a real person.
4. Done = a non-INDETERMINATE diplotype/phenotype on ≥1 PGP-UK individual → first real-people deterministic
   decode on a free non-gated cohort. Scale to the cohort + a report card as you see fit.

## ClinVar (Mendelian) extension — small gap, not a wall
`dna_decode/data/clinvar.py::ClinVarDecoder.call(chrom,pos,ref,alt)` is a **per-variant API** (no VCF CLI) →
needs a ~10-line VCF-iterate wrapper (+ the same GRCh37→38 liftOver). Cheap follow-on once PGx lands.

## Provenance
Accessions, VCF paths, and the full verification trail are in `wiki/independent_verification_audit_2026-07-04.md`
on the `mosfaer` branch. Assembly + command + non-dup all verified real-surface 2026-07-04.
