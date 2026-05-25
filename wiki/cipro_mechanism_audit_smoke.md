# Cipro AMRFinderPlus mechanism audit — N=2 cohort (2026-05-17)

**Purpose:** classify the actual cipro-resistance mechanism in each R strain (and check for silent mechanisms in S strains) before judging whether NT's preflight INCONCLUSIVE_MISS reflects TRUE biology vs a model failure.
**Tool:** AMRFinderPlus ncbi/amr:4.2.7-2026-03-24.1 on `Escherichia` organism mode + `--mutation_all`.
**Cohort:** `data\processed\gate_b_n40_cipro_cohort.parquet` (40 ciprofloxacin-labeled strains).

## Verdict: **MOSTLY_UNKNOWN**
- QRDR mechanism found in: **0** / 0 R strains
- Plasmid (qnr / aac6-Ib-cr) found in: **0** / 0 R strains
- R strains with NO known cipro mechanism: **0** / 0
- S strains with silent mechanism hit: **0** / 2

## Mechanism distribution (R set)

| mechanism | R count | S count |
|---|---:|---:|

## Per-strain mechanism table

| strain_id | accession | label | status | primary mech | mechanisms | mlst |
|---|---|---|---|---|---|---|
| 1045010.61 | GCA_008727135.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.11 |
| 1328432.3 | GCA_000492655.1 | S | OK | NO_MECHANISM |  | MLST.ecoli_achtman_4.2569 |

## How to use

- If verdict is QRDR_DOMINANT: NT's preflight INCONCLUSIVE_MISS is a MODEL failure (the signal is there in 70%+ of R strains; NT just isn't finding it).
- If verdict is MIXED_MECHANISMS: cipro resistance is heterogeneous; per-gene attribution should NOT expect to converge on one locus. Use per-strain known-mechanism overlap as ground truth for attribution.
- If verdict is MOSTLY_UNKNOWN: the biology is unclear; AMRFinderPlus alone may miss novel/regulatory mechanisms. Need Bakta gene-presence + downstream curated baseline.
- Silent-S count of 0 suggests label noise upper bound — these strains may be mislabeled or have functional-but-not-clinical-MIC resistance.

_JSON sidecar: `wiki\cipro_mechanism_audit_smoke.json`_