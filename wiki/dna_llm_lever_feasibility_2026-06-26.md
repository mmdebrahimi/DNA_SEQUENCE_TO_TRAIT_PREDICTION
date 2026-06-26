# DNA-LLM forward-lever feasibility census (2026-06-26)

_Read-only ceiling from the raw BV-BRC tables on disk; relaxed binary R/S phenotype. The within-lineage metric needs MLSTs carrying BOTH R and S with a downloadable assembly._

| drug | AST R/S | downloadable R/S | shared R+S lineages (ceiling) | within-lineage pairs (ceiling) | strains in shared | current probe |
|---|---|---|---|---|---|---|
| ciprofloxacin | 11014/28879 | 6124/15378 | 317 | 722845 | 12191 | 6 lin / 43 pairs (N=147) |
| tetracycline | 17667/20500 | 10770/10595 | 685 | 1187909 | 17677 | — |

_Ceiling = the most a re-selected cohort could achieve IF every downloadable strain in a shared lineage were fetched + AMRFinder'd. The gap to the current probe is the new-genome download + AMRFinder cost (Docker), NOT a compute/GPU question._