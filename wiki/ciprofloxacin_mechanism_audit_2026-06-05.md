# Ciprofloxacin AMRFinderPlus mechanism audit - N=147 cohort (2026-06-05)

**Drug:** ciprofloxacin
**Tool:** AMRFinderPlus `ncbi/amr:4.2.7-2026-03-24.1` on `Escherichia` organism mode + `--mutation_all`.
**Cohort:** `data\processed\stage2_n150_cipro_cohort.parquet` (147 pool strains for ciprofloxacin).
**Per-drug catalogs:** sourced from `dna_decode/data/mic_tiers.py`.

## Verdict: **PRIMARY_DOMINANT**
- Primary mechanism found in: **69** / 72 R strains
- Any known mechanism found in: **69** / 72 R strains
- R strains with NO known mechanism: **3** / 72
- S strains with silent mechanism hit: **20** / 75

## Mechanism distribution

| mechanism | R count | S count |
|---|---:|---:|
| QRDR_target_alteration | 63 | 12 |
| efflux | 0 | 1 |
| plasmid_protect_modify | 27 | 7 |
| regulatory | 10 | 3 |

## Per-strain mechanism table

| strain_id | accession | label | status | primary mech | mechanisms | mlst |
|---|---|---|---|---|---|---|
| 1045010.61 | GCA_008727135.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.11 |
| 1045010.62 | GCA_008727155.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.11 |
| 562.109859 | GCA_025200555.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.6617 |
| 562.111034 | GCA_025200615.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.398 |
| 562.111035 | GCA_025200675.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.12782 |
| 562.12962 | GCF_900093145.1 | S | OK | plasmid_protect_modify | efflux,plasmid_protect_modify | MLST.klebsiella.48 |
| 562.16325 | GCA_002056065.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.1408 |
| 562.16326 | GCA_002056145.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.1809 |
| 562.34086 | GCA_003571685.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.62 |
| 562.45846 | GCA_004566475.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.7522 |
| 562.50224 | GCA_004566395.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.38 |
| 562.50230 | GCA_004566545.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.69 |
| 562.50238 | GCA_004566715.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.501 |
| 562.50248 | GCA_004566905.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.70 |
| 562.50287 | GCA_004567805.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.442 |
| 562.50295 | GCA_004568615.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.Escherichia_coli_1.131 |
| 562.50298 | GCA_004569045.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50299 | GCA_004569675.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.457 |
| 562.50303 | GCA_004569925.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.38 |
| 562.52722 | GCA_009650035.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.29 |
| 562.7565 | GCA_001277395.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.517 |
| 562.7567 | GCA_001277435.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5488 |
| 562.7568 | GCA_001277455.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5339 |
| 562.7570 | GCA_001277495.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.40 |
| 562.7575 | GCA_001277595.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.337 |
| 562.7576 | GCA_001277615.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5241 |
| 562.7577 | GCA_001277635.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.301 |
| 562.7578 | GCA_001277655.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.382 |
| 562.7581 | GCA_001277715.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.121 |
| 562.7583 | GCA_001277755.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.10 |
| 562.7584 | GCA_001277775.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.583 |
| 562.7618 | GCA_001283165.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.10 |
| 562.7619 | GCA_001283185.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.29 |
| 562.7623 | GCA_001283265.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.328 |
| 562.7624 | GCA_001283285.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.328 |
| 562.7627 | GCA_001283345.2 | S | OK | regulatory | regulatory | MLST.Escherichia_coli_1.4554 |
| 562.7632 | GCA_001283445.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.5234 |
| 562.7634 | GCA_001283485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.301 |
| 562.7641 | GCA_001283625.2 | S | OK | NO_MECHANISM |  | MLST.Escherichia_coli_1.5543 |
| 562.7642 | GCA_001283645.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.4550 |
| 562.7645 | GCA_001283705.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.34 |
| 562.7646 | GCA_001283725.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5426 |
| 562.7649 | GCA_001283785.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5353 |
| 562.7652 | GCA_001283845.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.121 |
| 562.7656 | GCA_001283925.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5618 |
| 562.7665 | GCA_001284105.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.382 |
| 562.7669 | GCA_001284185.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.337 |
| 562.7671 | GCA_001284225.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.40 |
| 562.7673 | GCA_001284265.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5553 |
| 562.7678 | GCA_001284365.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.378 |
| 562.7679 | GCA_001284385.1 | S | OK | regulatory | regulatory | MLST.ecoli_achtman_4.589 |
| 562.7683 | GCA_001284465.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.378 |
| 562.7684 | GCA_001284485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.17 |
| 562.7686 | GCA_001284525.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.381 |
| 562.7687 | GCA_001284545.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.302 |
| 562.7690 | GCA_001284605.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2144 |
| 562.7695 | GCA_001284705.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.2089 |
| 562.7696 | GCA_001284725.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.803 |
| 562.7700 | GCA_001284805.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.34 |
| 562.7701 | GCA_001284825.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.335 |
| 562.7702 | GCA_001284845.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2089 |
| 562.7710 | GCA_001285005.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.3 |
| 562.7712 | GCA_001285045.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.4590 |
| 562.7724 | GCA_001285285.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.5406 |
| 562.7730 | GCA_001285405.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5485 |
| 562.7733 | GCA_001285465.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5330 |
| 562.7734 | GCA_001285485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.381 |
| 562.7742 | GCA_001285645.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.5381 |
| 562.7761 | GCA_001286025.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5582 |
| 562.7776 | GCA_001286345.2 | S | OK | regulatory | regulatory | MLST.ecoli_achtman_4.4554 |
| 562.7779 | GCA_001286385.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5251 |
| 562.7783 | GCA_001286465.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5445 |
| 562.7784 | GCA_001286485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2678 |
| 562.7785 | GCA_001286505.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5480 |
| 562.7790 | GCA_001286625.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.335 |
| 1438684.3 | GCA_000692695.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 1438686.3 | GCA_000692735.1 | R | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.5441 |
| 562.101974 | GCA_022694385.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.10 |
| 562.101975 | GCA_022694405.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.10 |
| 562.102124 | GCA_022758785.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.10 |
| 562.109860 | GCA_025200635.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.155 |
| 562.111036 | GCA_025200635.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.155 |
| 562.112391 | GCA_025782825.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.361 |
| 562.115101 | GCA_026420965.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.Escherichia_coli_2.398 |
| 562.115102 | GCA_026421005.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.224 |
| 562.115105 | GCA_026420945.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.683 |
| 562.115106 | GCA_026420925.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.Escherichia_coli_2.398 |
| 562.115107 | GCA_026421045.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.4380 |
| 562.16327 | GCA_002056635.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.16466 | GCA_002142695.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.617 |
| 562.17620 | GCA_002192275.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.617 |
| 562.17721 | GCA_002201835.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.1284 |
| 562.17722 | GCA_002202175.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.19631 | GCA_002180215.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.167 |
| 562.19632 | GCA_002180195.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.20499 | GCA_002180135.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.167 |
| 562.28566 | GCA_003000615.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.167 |
| 562.28805 | GCA_003073955.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.410 |
| 562.29003 | GCA_003075555.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.29006 | GCA_003073875.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.156 |
| 562.30349 | GCA_003194265.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.648 |
| 562.33281 | GCA_003324505.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.167 |
| 562.34085 | GCA_003571665.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.45845 | GCA_004564175.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.405 |
| 562.45849 | GCA_004567165.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.38 |
| 562.45850 | GCA_004567315.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.38 |
| 562.45851 | GCA_004567345.1 | R | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.372 |
| 562.45910 | GCA_004569165.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.131 |
| 562.50226 | GCA_004566415.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.50233 | GCA_004566595.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.410 |
| 562.50235 | GCA_004566665.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.131 |
| 562.50239 | GCA_004566755.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.7358 |
| 562.50241 | GCA_004566775.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50242 | GCA_004566785.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50244 | GCA_004566835.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.50247 | GCA_004566895.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50249 | GCA_004566915.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50251 | GCA_004566965.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50257 | GCA_004567115.1 | R | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.547 |
| 562.50259 | GCA_004567175.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.457 |
| 562.50261 | GCA_004567225.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.50266 | GCA_004567285.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.8671 |
| 562.50269 | GCA_004567375.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50270 | GCA_004567405.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.131 |
| 562.50273 | GCA_004567495.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50275 | GCA_004567535.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.50277 | GCA_004567585.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50281 | GCA_004567675.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.88 |
| 562.50283 | GCA_004567715.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50285 | GCA_004567745.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50289 | GCA_004567895.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.131 |
| 562.50291 | GCA_004568145.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50293 | GCA_004568225.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.131 |
| 562.50305 | GCA_004570035.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.744 |
| 562.59218 | GCA_902807085.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.410 |
| 562.59222 | GCA_902807135.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.5474 |
| 562.59224 | GCA_902807205.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.410 |
| 562.59226 | GCA_902807265.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.224 |
| 562.59229 | GCA_902807295.1 | R | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.5474 |
| 562.7572 | GCA_001277535.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.301 |
| 562.7573 | GCA_001277555.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.443 |
| 562.7660 | GCA_001284005.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.443 |
| 562.7693 | GCA_001284665.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.443 |
| 562.7699 | GCA_001284785.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.335 |
| 562.7703 | GCA_001284865.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.443 |
| 562.7771 | GCA_001286225.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.443 |
| 562.7775 | GCA_001286305.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.301 |

## Verdict interpretation

- **PRIMARY_DOMINANT:** ≥ 70% of R strains carry a primary mechanism for this drug. NT/classifier should find this signal; if it does not, model failure is implicated.
- **MIXED_MECHANISMS:** 50-69% of R strains carry primary OR co-resistance mechanisms. Per-gene attribution is unlikely to converge on one locus.
- **MOSTLY_UNKNOWN:** < 50% of R strains have any known mechanism. AMRFinder alone may miss novel/regulatory mechanisms; pair with Bakta gene-presence + curated baseline.

_JSON sidecar: `wiki\ciprofloxacin_mechanism_audit_2026-06-05.json`_