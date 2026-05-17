# Cipro AST/MIC audit — N=38 cohort (2026-05-17)

**Purpose:** rejoin raw BV-BRC AST MIC values onto the cohort to surface borderline/no-MIC strains whose binary R/S label is structurally noisy.
**Source:** `C:\Users\Farshad\Downloads\BVBRC_genome_amr.csv`
**Breakpoints:** CLSI 2024 (S ≤0.5 / R ≥2.0), EUCAST 14.0 (S ≤0.25 / R ≥1.0) ug/mL.

## Verdict: **NOISY**
- Decisive-R (HIGH_R only): **7** / 20
- Decisive-S (HIGH_S only): **0** / 20
- Total decisive subset: **7** strains
- Clean-R fraction: **35.0%**

## Tier distribution

| tier | R strains | S strains | description |
|---|---:|---:|---|
| HIGH_R | 7 | 0 | MIC ≥8.0 (4× CLSI-R) — decisive R |
| HIGH_S | 0 | 0 | MIC ≤0.125 (1/4× CLSI-S) — decisive S |
| DECISIVE | 0 | 1 | CLSI/EUCAST agree, not borderline |
| BORDERLINE | 3 | 12 | MIC in [0.25, 4.0] — within 2× breakpoint |
| AMBIGUOUS | 0 | 3 | CLSI and EUCAST disagree on call |
| CONFLICT | 1 | 0 | multiple AST rows disagree R vs S |
| NO_MIC | 9 | 4 | no MIC value parsed |

## Per-strain audit

| strain_id | accession | label | tier | median MIC | n_rows | n_mic | mlst |
|---|---|---|---|---:|---:|---:|---|
| 1328432.3 | GCA_000492655.1 | S | AMBIGUOUS | 0.500 | 1 | 1 | MLST.ecoli_achtman_4.2569 |
| 562.50295 | GCA_004568615.1 | S | AMBIGUOUS | 1.000 | 2 | 2 | MLST.Escherichia_coli_1.131 |
| 562.7717 | GCA_001285145.1 | S | AMBIGUOUS | 0.500 | 1 | 1 | MLST.ecoli_achtman_4.2346 |
| 562.16325 | GCA_002056065.1 | S | BORDERLINE | 0.250 | 2 | 2 | MLST.ecoli_achtman_4.1408 |
| 562.16326 | GCA_002056145.1 | S | BORDERLINE | 0.250 | 2 | 2 | MLST.ecoli_achtman_4.1809 |
| 562.28389 | GCA_002948655.1 | S | BORDERLINE | 0.250 | 2 | 1 | MLST.ecoli_achtman_4.1276 |
| 562.50287 | GCA_004567805.1 | S | BORDERLINE | 0.250 | 2 | 1 | MLST.ecoli_achtman_4.442 |
| 562.52722 | GCA_009650035.1 | S | BORDERLINE | 0.250 | 2 | 1 | MLST.ecoli_achtman_4.29 |
| 562.7575 | GCA_001277595.1 | S | BORDERLINE | 0.250 | 2 | 1 | MLST.ecoli_achtman_4.337 |
| 562.7627 | GCA_001283345.2 | S | BORDERLINE | 0.250 | 2 | 2 | MLST.Escherichia_coli_1.4554 |
| 562.7690 | GCA_001284605.1 | S | BORDERLINE | 0.250 | 3 | 2 | MLST.ecoli_achtman_4.2144 |
| 562.7695 | GCA_001284705.1 | S | BORDERLINE | 0.250 | 1 | 1 | MLST.ecoli_achtman_4.2089 |
| 562.7710 | GCA_001285005.1 | S | BORDERLINE | 0.250 | 3 | 1 | MLST.ecoli_achtman_4.3 |
| 562.7784 | GCA_001286485.1 | S | BORDERLINE | 0.250 | 3 | 2 | MLST.ecoli_achtman_4.2678 |
| 562.7789 | GCA_001286605.1 | S | BORDERLINE | 0.250 | 2 | 1 | MLST.ecoli_achtman_4.123 |
| 1045010.61 | GCA_008727135.1 | S | DECISIVE | 0.133 | 3 | 2 | MLST.ecoli_achtman_4.11 |
| 562.28565 | GCA_003000595.1 | S | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.38 |
| 562.45853 | GCA_004567665.1 | S | NO_MIC | - | 1 | 0 | MLST.ecoli_achtman_4.315 |
| 562.50301 | GCA_004569855.1 | S | NO_MIC | - | 1 | 0 | MLST.ecoli_achtman_4.1317 |
| 562.7641 | GCA_001283625.2 | S | NO_MIC | - | 1 | 0 | MLST.Escherichia_coli_1.5543 |
| 1328434.3 | GCA_000522325.1 | R | BORDERLINE | 4.000 | 1 | 1 | MLST.ecoli_achtman_4.405 |
| 562.7572 | GCA_001277535.1 | R | BORDERLINE | 4.000 | 3 | 1 | MLST.ecoli_achtman_4.301 |
| 562.7699 | GCA_001284785.1 | R | BORDERLINE | 4.000 | 2 | 1 | MLST.ecoli_achtman_4.335 |
| 562.45851 | GCA_004567345.1 | R | CONFLICT | - | 2 | 1 | MLST.ecoli_achtman_4.372 |
| 1328433.3 | GCA_000522345.1 | R | HIGH_R | 8.000 | 1 | 1 | MLST.ecoli_achtman_4.131 |
| 562.12960 | GCF_001874845.1 | R | HIGH_R | 16.000 | 2 | 1 | MLST.ecoli_achtman_4.101 |
| 562.17621 | GCA_002192295.1 | R | HIGH_R | 16.000 | 1 | 1 | MLST.ecoli_achtman_4.167 |
| 562.17721 | GCA_002201835.1 | R | HIGH_R | 12.000 | 3 | 2 | MLST.ecoli_achtman_4.1284 |
| 562.28805 | GCA_003073955.1 | R | HIGH_R | 8.000 | 1 | 1 | MLST.ecoli_achtman_4.410 |
| 562.30362 | GCA_003204155.1 | R | HIGH_R | 16.000 | 2 | 1 | MLST.ecoli_achtman_4.4 |
| 562.50250 | GCA_004566955.1 | R | HIGH_R | 8.000 | 1 | 1 | MLST.ecoli_achtman_4.44 |
| 562.13502 | GCF_001747365.1 | R | NO_MIC | - | 1 | 0 | MLST.ecoli_achtman_4.10 |
| 562.22426 | GCA_002180135.1 | R | NO_MIC | - | 0 | 0 | MLST.Escherichia_coli_1.167,MLST.Escherichia_coli_2.2 |
| 562.28563 | GCA_002999075.1 | R | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.156 |
| 562.45848 | GCA_004567065.1 | R | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.349 |
| 562.50237 | GCA_004566685.1 | R | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.224 |
| 562.50245 | GCA_004566865.1 | R | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.1431 |
| 562.50252 | GCA_004566985.1 | R | NO_MIC | - | 1 | 0 | MLST.ecoli_achtman_4.393 |
| 562.50304 | GCA_004569995.1 | R | NO_MIC | - | 1 | 0 | MLST.ecoli_achtman_4.1193 |
| 562.59220 | GCA_902807115.1 | R | NO_MIC | - | 0 | 0 | MLST.ecoli_achtman_4.354 |

## Decisive subset (recommended for clean re-run)

- **Decisive-R (HIGH_R):** ['1328433.3', '562.12960', '562.17621', '562.17721', '562.28805', '562.30362', '562.50250']
- **Decisive-S (HIGH_S):** []
- **N = 7** (vs cohort N=38)

## How to use

1. If verdict is CLEAN: cohort labels are trustworthy; downstream NT/k-mer FAIL is a model issue, not label noise.
2. If verdict is MIXED or NOISY: re-run Stage 1 on the decisive subset; if AUROC ≥0.10 better, label noise was a real confounder.
3. Use the AMBIGUOUS / CONFLICT / NO_MIC strains as a flag list to exclude or re-label.
4. Feed `decisive_R_ids` + `decisive_S_ids` to a curated-baseline experiment for a cleaner ground truth.

_JSON sidecar: `wiki\cipro_mic_audit_2026-05-17.json`_