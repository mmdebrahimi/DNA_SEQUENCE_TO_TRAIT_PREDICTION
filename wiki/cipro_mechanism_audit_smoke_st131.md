# Cipro AMRFinderPlus mechanism audit — N=1 cohort (2026-05-17)

**Purpose:** classify the actual cipro-resistance mechanism in each R strain (and check for silent mechanisms in S strains) before judging whether NT's preflight INCONCLUSIVE_MISS reflects TRUE biology vs a model failure.
**Tool:** AMRFinderPlus ncbi/amr:4.2.7-2026-03-24.1 on `Escherichia` organism mode + `--mutation_all`.
**Cohort:** `data\processed\gate_b_n40_cipro_cohort.parquet` (40 ciprofloxacin-labeled strains).

## Verdict: **QRDR_DOMINANT**
- QRDR mechanism found in: **1** / 1 R strains
- Plasmid (qnr / aac6-Ib-cr) found in: **0** / 1 R strains
- R strains with NO known cipro mechanism: **0** / 1
- S strains with silent mechanism hit: **0** / 0

## Mechanism distribution (R set)

| mechanism | R count | S count |
|---|---:|---:|
| QRDR_target_alteration | 1 | 0 |
| efflux | 1 | 0 |
| regulatory | 1 | 0 |

## Per-strain mechanism table

| strain_id | accession | label | status | primary mech | mechanisms | mlst |
|---|---|---|---|---|---|---|
| 1328433.3 | GCA_000522345.1 | R | OK | QRDR_target_alteration | QRDR_target_alteration,efflux,regulatory | MLST.ecoli_achtman_4.131 |

## How to use

- If verdict is QRDR_DOMINANT: NT's preflight INCONCLUSIVE_MISS is a MODEL failure (the signal is there in 70%+ of R strains; NT just isn't finding it).
- If verdict is MIXED_MECHANISMS: cipro resistance is heterogeneous; per-gene attribution should NOT expect to converge on one locus. Use per-strain known-mechanism overlap as ground truth for attribution.
- If verdict is MOSTLY_UNKNOWN: the biology is unclear; AMRFinderPlus alone may miss novel/regulatory mechanisms. Need Bakta gene-presence + downstream curated baseline.
- Silent-S count of 0 suggests label noise upper bound — these strains may be mislabeled or have functional-but-not-clinical-MIC resistance.

_JSON sidecar: `wiki\cipro_mechanism_audit_smoke_st131.json`_