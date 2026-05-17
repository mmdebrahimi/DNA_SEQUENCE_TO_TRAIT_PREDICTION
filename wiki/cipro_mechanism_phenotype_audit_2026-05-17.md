# Cipro mechanism x phenotype audit — N=40 cohort (2026-05-17)

**Purpose:** the load-bearing combined audit table from Codex's cross-question integration plan. Joins AMRFinder mechanism classification (Experiment 1) x MIC noise tiering (Experiment 2) to decide whether the N=38 cipro labels are clean enough for any genomic predictor to fit.

## Verdict: **NOISE_DOMINATES**
- Cohort signal quality (CLEAN / total): **0.17**
- Clean strains: 7R + 0S = **7**
- Suspect strains: 0R + 6S = **6** (likely mislabeled or atypical biology)
- Opaque-R strains (HIGH_R + no primary mechanism): **0** — tool/parser miss vs novel biology
- Total opacity flagged: **0** (informs whether to suspect labels OR AMRFinder)
- Noisy strains: 13R + 13S = **26** (label is structurally unreliable)

## Pre-curated-baseline gate: **SUSPEND_CONDITION_4**
- clean_count=7 too low AND opacity_count=0 not the bottleneck — the N=38 cohort is structurally unusable for PIVOT TRIGGER condition 4. Next experiments: (a) cohort expansion to N=150 with strict MIC filter, OR (b) per-gene NT windows on the small clean set as a diagnostic.

Definitions:
- **Primary cipro mechanism** = QRDR target alteration (gyrA/parC/parE) OR plasmid quinolone protection/modification (qnr* / aac6-Ib-cr). Textbook cipro-conferring.
- **Co-resistance modifiers** = efflux + regulatory + porin_loss. Real biology but do not confer cipro-R alone; reported separately, do not drive noise classification.
- **Mechanism opacity flag** = HIGH_R + no primary mechanism. AMRFinder may have a catalog gap; distinct from label noise.

## Noise class distribution

| class | count | meaning |
|---|---:|---|
| CLEAN_R_primary_mechanism | 7 | HIGH_R MIC + QRDR or qnr/aac present |
| OPAQUE_R_co_resistance_only | 0 | HIGH_R + only efflux/regulatory/porin found — tool gap likely |
| OPAQUE_R_no_mechanism | 0 | HIGH_R but no AMRFinder hits at all |
| NOISY_R_borderline | 4 | borderline/ambiguous MIC — label may be wrong |
| NOISY_R_no_mic | 9 | no MIC in BV-BRC — label is opaque |
| CLEAN_S_no_primary_mechanism | 0 | HIGH_S MIC + no QRDR / no qnr — clean susceptible |
| SUSPECT_S_silent_primary_mechanism | 0 | HIGH_S MIC but primary mechanism present — silent / non-functional |
| SUSPECT_S_borderline_primary_mechanism | 6 | borderline MIC + primary mech — likely mislabeled to S |
| NOISY_S_borderline | 10 | borderline MIC with no primary mech — label opaque |
| NOISY_S_no_mic | 3 | no MIC in BV-BRC — label opaque |

## Per-strain merged table

| strain_id | accession | label | mic_tier | med MIC | primary? | mechs | co-res | noise_class | opacity | mlst |
|---|---|---|---|---:|---|---|---|---|---|---|
| 1328432.3 | GCA_000492655.1 | S | AMBIGUOUS | 0.500 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.2569 |
| 562.28389 | GCA_002948655.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.1276 |
| 562.50287 | GCA_004567805.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.442 |
| 562.52722 | GCA_009650035.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.29 |
| 562.7575 | GCA_001277595.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.337 |
| 562.7627 | GCA_001283345.2 | S | BORDERLINE | 0.250 | no | regulatory | regulatory | NOISY_S_borderline |   | MLST.Escherichia_coli_1.4554 |
| 562.7690 | GCA_001284605.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.2144 |
| 562.7710 | GCA_001285005.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.3 |
| 562.7784 | GCA_001286485.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.2678 |
| 562.7789 | GCA_001286605.1 | S | BORDERLINE | 0.250 | no | - | - | NOISY_S_borderline |   | MLST.ecoli_achtman_4.123 |
| 562.28565 | GCA_003000595.1 | S | NO_MIC | - | no | - | - | NOISY_S_no_mic |   | MLST.ecoli_achtman_4.38 |
| 562.45853 | GCA_004567665.1 | S | NO_MIC | - | no | - | - | NOISY_S_no_mic |   | MLST.ecoli_achtman_4.315 |
| 562.7641 | GCA_001283625.2 | S | NO_MIC | - | no | - | - | NOISY_S_no_mic |   | MLST.Escherichia_coli_1.5543 |
| 1045010.61 | GCA_008727135.1 | S | DECISIVE | 0.133 | no | - | - | OTHER_S |   | MLST.ecoli_achtman_4.11 |
| 562.16325 | GCA_002056065.1 | S | BORDERLINE | 0.250 | yes | plasmid_protect_modify | - | SUSPECT_S_borderline_primary_mechanism |   | MLST.ecoli_achtman_4.1408 |
| 562.16326 | GCA_002056145.1 | S | BORDERLINE | 0.250 | yes | plasmid_protect_modify | - | SUSPECT_S_borderline_primary_mechanism |   | MLST.ecoli_achtman_4.1809 |
| 562.50295 | GCA_004568615.1 | S | AMBIGUOUS | 1.000 | yes | QRDR_target_alteration | - | SUSPECT_S_borderline_primary_mechanism |   | MLST.Escherichia_coli_1.131 |
| 562.50301 | GCA_004569855.1 | S | NO_MIC | - | yes | QRDR_target_alteration | - | SUSPECT_S_borderline_primary_mechanism |   | MLST.ecoli_achtman_4.1317 |
| 562.7695 | GCA_001284705.1 | S | BORDERLINE | 0.250 | yes | QRDR_target_alteration | - | SUSPECT_S_borderline_primary_mechanism |   | MLST.ecoli_achtman_4.2089 |
| 562.7717 | GCA_001285145.1 | S | AMBIGUOUS | 0.500 | yes | QRDR_target_alteration,regulatory | regulatory | SUSPECT_S_borderline_primary_mechanism |   | MLST.ecoli_achtman_4.2346 |
| 1328433.3 | GCA_000522345.1 | R | HIGH_R | 8.000 | yes | QRDR_target_alteration,regulatory | regulatory | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.131 |
| 562.12960 | GCF_001874845.1 | R | HIGH_R | 16.000 | yes | QRDR_target_alteration,plasmid_protect_modify,regulatory | regulatory | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.101 |
| 562.17621 | GCA_002192295.1 | R | HIGH_R | 16.000 | yes | QRDR_target_alteration,regulatory | regulatory | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.167 |
| 562.17721 | GCA_002201835.1 | R | HIGH_R | 12.000 | yes | QRDR_target_alteration,plasmid_protect_modify | - | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.1284 |
| 562.28805 | GCA_003073955.1 | R | HIGH_R | 8.000 | yes | QRDR_target_alteration,plasmid_protect_modify | - | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.410 |
| 562.30362 | GCA_003204155.1 | R | HIGH_R | 16.000 | yes | QRDR_target_alteration,regulatory | regulatory | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.4 |
| 562.50250 | GCA_004566955.1 | R | HIGH_R | 8.000 | yes | QRDR_target_alteration | - | CLEAN_R_primary_mechanism |   | MLST.ecoli_achtman_4.44 |
| 1328434.3 | GCA_000522325.1 | R | BORDERLINE | 4.000 | yes | QRDR_target_alteration | - | NOISY_R_borderline |   | MLST.ecoli_achtman_4.405 |
| 562.45851 | GCA_004567345.1 | R | CONFLICT | - | no | - | - | NOISY_R_borderline |   | MLST.ecoli_achtman_4.372 |
| 562.7572 | GCA_001277535.1 | R | BORDERLINE | 4.000 | yes | QRDR_target_alteration | - | NOISY_R_borderline |   | MLST.ecoli_achtman_4.301 |
| 562.7699 | GCA_001284785.1 | R | BORDERLINE | 4.000 | yes | QRDR_target_alteration | - | NOISY_R_borderline |   | MLST.ecoli_achtman_4.335 |
| 562.13502 | GCF_001747365.1 | R | NO_MIC | - | yes | QRDR_target_alteration | - | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.10 |
| 562.22426 | GCA_002180135.1 | R | NO_MIC | - | yes | QRDR_target_alteration,plasmid_protect_modify,regulatory | regulatory | NOISY_R_no_mic |   | MLST.Escherichia_coli_1.167,MLST.Escherichia_coli_2.2 |
| 562.28563 | GCA_002999075.1 | R | NO_MIC | - | yes | QRDR_target_alteration,plasmid_protect_modify,regulatory | regulatory | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.156 |
| 562.45848 | GCA_004567065.1 | R | NO_MIC | - | no | - | - | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.349 |
| 562.50237 | GCA_004566685.1 | R | NO_MIC | - | yes | QRDR_target_alteration,plasmid_protect_modify | - | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.224 |
| 562.50245 | GCA_004566865.1 | R | NO_MIC | - | yes | QRDR_target_alteration,plasmid_protect_modify | - | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.1431 |
| 562.50252 | GCA_004566985.1 | R | NO_MIC | - | yes | QRDR_target_alteration,regulatory | regulatory | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.393 |
| 562.50304 | GCA_004569995.1 | R | NO_MIC | - | yes | QRDR_target_alteration,regulatory | regulatory | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.1193 |
| 562.59220 | GCA_902807115.1 | R | NO_MIC | - | yes | QRDR_target_alteration | - | NOISY_R_no_mic |   | MLST.ecoli_achtman_4.354 |

## How to use (gated)

Pre-curated-baseline gate fired: **SUSPEND_CONDITION_4**.
- clean_count=7 too low AND opacity_count=0 not the bottleneck — the N=38 cohort is structurally unusable for PIVOT TRIGGER condition 4. Next experiments: (a) cohort expansion to N=150 with strict MIC filter, OR (b) per-gene NT windows on the small clean set as a diagnostic.

Verdict-level interpretation:
- **SIGNAL_DOMINATES (>=0.70 clean):** the cohort is clean enough that NT's FAIL on cipro is a genuine model issue, not label noise. Curated baseline verdict is load-bearing.
- **MIXED (0.40-0.70 clean):** label noise is a real confounder. Curated baseline at full N is descriptive; consider a clean-subset rerun if `clean_count >= 20`.
- **NOISE_DOMINATES (<0.40 clean):** the N=38 cohort is structurally noisy. Two branches: if `opacity_count >= 5` -> AMRFinder may be missing mechanisms, debug before declaring labels unusable; otherwise -> N=150 expansion or per-gene NT diagnostic.

_JSON sidecar: `wiki\cipro_mechanism_phenotype_audit_2026-05-17.json`_