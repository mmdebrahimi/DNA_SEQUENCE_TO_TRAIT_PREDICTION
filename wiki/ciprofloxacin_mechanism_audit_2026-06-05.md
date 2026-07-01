# Ciprofloxacin AMRFinderPlus mechanism audit - N=50 cohort (2026-06-05)

**Drug:** ciprofloxacin
**Tool:** AMRFinderPlus `ncbi/amr:4.2.7-2026-03-24.1` on `Escherichia` organism mode + `--mutation_all`.
**Cohort:** `data\processed\gate_b_cohort.parquet` (50 pool strains for ciprofloxacin).
**Per-drug catalogs:** sourced from `dna_decode/data/mic_tiers.py`.

## Verdict: **PRIMARY_DOMINANT**
- Primary mechanism found in: **18** / 20 R strains
- Any known mechanism found in: **18** / 20 R strains
- R strains with NO known mechanism: **2** / 20
- S strains with silent mechanism hit: **10** / 30

## Mechanism distribution

| mechanism | R count | S count |
|---|---:|---:|
| QRDR_target_alteration | 18 | 6 |
| plasmid_protect_modify | 7 | 3 |
| regulatory | 8 | 2 |

## Per-strain mechanism table

| strain_id | accession | label | status | primary mech | mechanisms | mlst |
|---|---|---|---|---|---|---|
| 1045010.61 | GCA_008727135.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.11 |
| 1328432.3 | GCA_000492655.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2569 |
| 562.16325 | GCA_002056065.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.1408 |
| 562.16326 | GCA_002056145.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.1809 |
| 562.28389 | GCA_002948655.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.1276 |
| 562.28565 | GCA_003000595.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.38 |
| 562.45853 | GCA_004567665.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.315 |
| 562.50287 | GCA_004567805.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.442 |
| 562.50295 | GCA_004568615.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.Escherichia_coli_1.131 |
| 562.50301 | GCA_004569855.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.1317 |
| 562.52722 | GCA_009650035.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.29 |
| 562.7570 | GCA_001277495.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.40 |
| 562.7575 | GCA_001277595.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.337 |
| 562.7578 | GCA_001277655.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.382 |
| 562.7580 | GCA_001277695.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.4424 |
| 562.7581 | GCA_001277715.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.121 |
| 562.7623 | GCA_001283265.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.328 |
| 562.7625 | GCA_001283305.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.378 |
| 562.7627 | GCA_001283345.2 | S | OK | regulatory | regulatory | MLST.Escherichia_coli_1.4554 |
| 562.7641 | GCA_001283625.2 | S | OK | NO_MECHANISM |  | MLST.Escherichia_coli_1.5543 |
| 562.7645 | GCA_001283705.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.34 |
| 562.7684 | GCA_001284485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.17 |
| 562.7686 | GCA_001284525.1 | S | OK | plasmid_protect_modify | plasmid_protect_modify | MLST.ecoli_achtman_4.381 |
| 562.7687 | GCA_001284545.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.302 |
| 562.7690 | GCA_001284605.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2144 |
| 562.7695 | GCA_001284705.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.2089 |
| 562.7710 | GCA_001285005.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.3 |
| 562.7717 | GCA_001285145.1 | S | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.2346 |
| 562.7784 | GCA_001286485.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2678 |
| 562.7789 | GCA_001286605.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.123 |
| 1328433.3 | GCA_000522345.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.131 |
| 1328434.3 | GCA_000522325.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.405 |
| 562.12960 | GCF_001874845.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.101 |
| 562.13502 | GCF_001747365.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.10 |
| 562.17621 | GCA_002192295.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.167 |
| 562.17721 | GCA_002201835.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.1284 |
| 562.22426 | GCA_002180135.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.Escherichia_coli_1.167,MLST.Escherichia_coli_2.2 |
| 562.28563 | GCA_002999075.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify,regulatory | MLST.ecoli_achtman_4.156 |
| 562.28805 | GCA_003073955.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.410 |
| 562.30362 | GCA_003204155.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.4 |
| 562.45848 | GCA_004567065.1 | R | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.349 |
| 562.45851 | GCA_004567345.1 | R | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.372 |
| 562.50237 | GCA_004566685.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.224 |
| 562.50245 | GCA_004566865.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,plasmid_protect_modify | MLST.ecoli_achtman_4.1431 |
| 562.50250 | GCA_004566955.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.44 |
| 562.50252 | GCA_004566985.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.393 |
| 562.50304 | GCA_004569995.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,regulatory | MLST.ecoli_achtman_4.1193 |
| 562.59220 | GCA_902807115.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.354 |
| 562.7572 | GCA_001277535.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.301 |
| 562.7699 | GCA_001284785.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration | MLST.ecoli_achtman_4.335 |

## Verdict interpretation

- **PRIMARY_DOMINANT:** ≥ 70% of R strains carry a primary mechanism for this drug. NT/classifier should find this signal; if it does not, model failure is implicated.
- **MIXED_MECHANISMS:** 50-69% of R strains carry primary OR co-resistance mechanisms. Per-gene attribution is unlikely to converge on one locus.
- **MOSTLY_UNKNOWN:** < 50% of R strains have any known mechanism. AMRFinder alone may miss novel/regulatory mechanisms; pair with Bakta gene-presence + curated baseline.

_JSON sidecar: `wiki\ciprofloxacin_mechanism_audit_2026-06-05.json`_