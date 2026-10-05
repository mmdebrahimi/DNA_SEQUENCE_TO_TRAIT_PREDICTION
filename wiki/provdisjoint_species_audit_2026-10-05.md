# Species composition of the provenance-disjoint validation cohorts (2026-10-05)

Read-only audit. **It discloses; it does not re-score.** The 10
`provenance_disjoint_validation_*.json` artifacts are frozen units of the reproducibility
freeze and are untouched here.

**Species come from each assembly's own GenBank `ORGANISM` line** — the submitter's assertion,
not a wet-lab identification. Strong evidence about what a cohort *contains*; not an
independent species call.

**No cohort is decomposed unless its committed confusion matrix reproduces exactly** from the
cached AMRFinder runs through the deployed `call_resistance`. A cohort with incomplete cached
runs reports composition and **no** outcome cross-tab.

| cohort | scored as | composition | off-species | reconciled | FN off-species |
|---|---|---|---|---|---|
| `campylobacter_provdisjoint_ciprofloxacin` | `Campylobacter` | MATCHES_SCORED_ORGANISM | 0 | yes | 0/0 |
| `escherichia_coli_shigella_provdisjoint_ceftriaxone` | `Escherichia` | MATCHES_SCORED_ORGANISM | 0 | yes | 0/1 |
| `escherichia_coli_shigella_provdisjoint_ciprofloxacin` | `Escherichia` | MATCHES_SCORED_ORGANISM | 0 | yes | 0/2 |
| `escherichia_coli_shigella_provdisjoint_gentamicin` | `Escherichia` | MATCHES_SCORED_ORGANISM | 0 | yes | 0/3 |
| `escherichia_coli_shigella_provdisjoint_tetracycline` | `Escherichia` | MATCHES_SCORED_ORGANISM | 0 | yes | 0/2 |
| `klebsiella_provdisjoint_ceftriaxone` | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | no | — |
| `klebsiella_provdisjoint_ciprofloxacin` | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | yes | 0/1 |
| `klebsiella_provdisjoint_gentamicin` | `Klebsiella_pneumoniae` | MIXED_SPECIES | 1 | no | — |
| `klebsiella_provdisjoint_meropenem` | `Klebsiella_pneumoniae` | MIXED_SPECIES | 23 | yes | 16/16 |
| `klebsiella_provdisjoint_tetracycline` | `Klebsiella_pneumoniae` | MIXED_SPECIES | 2 | no | — |

## campylobacter_provdisjoint_ciprofloxacin  (`ciprofloxacin`, scored as `Campylobacter`)

- verdict: **MATCHES_SCORED_ORGANISM**
  - 31x `Campylobacter jejuni`  (in_organism_vocab)
  - 9x `Campylobacter coli`  (in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 20, 'fp': 0, 'tn': 20, 'fn': 0, 'n_scored': 40}
- false negatives: **0**, of which **0 are off-species** and 0 are the expected species

## escherichia_coli_shigella_provdisjoint_ceftriaxone  (`ceftriaxone`, scored as `Escherichia`)

- verdict: **MATCHES_SCORED_ORGANISM**
  - 60x `Escherichia coli`  (in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 29, 'fp': 1, 'tn': 29, 'fn': 1, 'n_scored': 60}
- false negatives: **1**, of which **0 are off-species** and 1 are the expected species
  - FN: 1x `Escherichia coli`

## escherichia_coli_shigella_provdisjoint_ciprofloxacin  (`ciprofloxacin`, scored as `Escherichia`)

- verdict: **MATCHES_SCORED_ORGANISM**
  - 60x `Escherichia coli`  (in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 28, 'fp': 9, 'tn': 21, 'fn': 2, 'n_scored': 60}
- false negatives: **2**, of which **0 are off-species** and 2 are the expected species
  - FN: 2x `Escherichia coli`

## escherichia_coli_shigella_provdisjoint_gentamicin  (`gentamicin`, scored as `Escherichia`)

- verdict: **MATCHES_SCORED_ORGANISM**
  - 60x `Escherichia coli`  (in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 27, 'fp': 0, 'tn': 30, 'fn': 3, 'n_scored': 60}
- false negatives: **3**, of which **0 are off-species** and 3 are the expected species
  - FN: 3x `Escherichia coli`

## escherichia_coli_shigella_provdisjoint_tetracycline  (`tetracycline`, scored as `Escherichia`)

- verdict: **MATCHES_SCORED_ORGANISM**
  - 60x `Escherichia coli`  (in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 28, 'fp': 2, 'tn': 28, 'fn': 2, 'n_scored': 60}
- false negatives: **2**, of which **0 are off-species** and 2 are the expected species
  - FN: 2x `Escherichia coli`

## klebsiella_provdisjoint_ceftriaxone  (`ceftriaxone`, scored as `Klebsiella_pneumoniae`)

- verdict: **MIXED_SPECIES**
  - 58x `Klebsiella pneumoniae`  (in_organism_vocab)
  - 2x `Klebsiella aerogenes`  (not_in_organism_vocab)
- outcome cross-tab **withheld**: cached AMRFinder runs incomplete (54/60); attributing outcomes on a partial run set would compute a rate over a silently-shrunken denominator. Composition above is unaffected -- it needs only the assemblies, which are complete.

## klebsiella_provdisjoint_ciprofloxacin  (`ciprofloxacin`, scored as `Klebsiella_pneumoniae`)

- verdict: **MIXED_SPECIES**
  - 58x `Klebsiella pneumoniae`  (in_organism_vocab)
  - 2x `Klebsiella aerogenes`  (not_in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 29, 'fp': 1, 'tn': 29, 'fn': 1, 'n_scored': 60}
- false negatives: **1**, of which **0 are off-species** and 1 are the expected species
  - FN: 1x `Klebsiella pneumoniae`

## klebsiella_provdisjoint_gentamicin  (`gentamicin`, scored as `Klebsiella_pneumoniae`)

- verdict: **MIXED_SPECIES**
  - 59x `Klebsiella pneumoniae`  (in_organism_vocab)
  - 1x `Klebsiella aerogenes`  (not_in_organism_vocab)
- outcome cross-tab **withheld**: cached AMRFinder runs incomplete (3/60); attributing outcomes on a partial run set would compute a rate over a silently-shrunken denominator. Composition above is unaffected -- it needs only the assemblies, which are complete.

## klebsiella_provdisjoint_meropenem  (`meropenem`, scored as `Klebsiella_pneumoniae`)

- verdict: **MIXED_SPECIES**
  - 37x `Klebsiella pneumoniae`  (in_organism_vocab)
  - 23x `Klebsiella aerogenes`  (not_in_organism_vocab)
- committed confusion reproduced exactly: {'tp': 14, 'fp': 3, 'tn': 27, 'fn': 16, 'n_scored': 60}
- false negatives: **16**, of which **16 are off-species** and 0 are the expected species
  - FN: 16x `Klebsiella aerogenes`

## klebsiella_provdisjoint_tetracycline  (`tetracycline`, scored as `Klebsiella_pneumoniae`)

- verdict: **MIXED_SPECIES**
  - 58x `Klebsiella pneumoniae`  (in_organism_vocab)
  - 2x `Klebsiella aerogenes`  (not_in_organism_vocab)
- outcome cross-tab **withheld**: cached AMRFinder runs incomplete (33/60); attributing outcomes on a partial run set would compute a rate over a silently-shrunken denominator. Composition above is unaffected -- it needs only the assemblies, which are complete.

