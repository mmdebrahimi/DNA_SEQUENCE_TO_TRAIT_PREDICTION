# Does the SHIPPED (blosum62) rank inverse generalize? — ProteinGym N=200 (2026-07-17)

**No, not usefully. On 200 scored ProteinGym assays the wheel-only default is materially better than guessing on only 27 (13.5%), and often WORSE.**

Tonight's inverse-design generalization rested on FOUR proteins and reported the rank inverse beating the no-oracle null 4/4. That 4/4 was **ESM** (the learned method, GPU-only). The shipped `dna-decode inverse` CLI DEFAULT is `blosum62` (pure-python, no GPU) -- so this asks the honest deployability question at scale: across the whole ProteinGym substitution benchmark, does the wheel-only default actually work?

## Result

| | |
|---|---:|
| assays scored | 200 |
| beats the no-oracle null **materially** (>=25%) | **27 (13.5%)** |
| beats it at **all** (>0%) | 118 (59.0%) |
| margin vs null (median) | +7.1% |
| margin vs null (q25 .. q75) | -14.4% .. +19.3% |
| margin vs null (min .. max) | -81.2% .. +46.9% |

The margin is NEGATIVE across the whole lower half below the median-ish: q25 is -14.4%, i.e. on a quarter of proteins the blosum62 inverse is >14% WORSE than picking a variant at random. 'Any positive edge' at 59% is barely a coin flip.

## It is not a clade artifact -- the weakness is general

| taxon | n | material wins | rate |
|---|---:|---:|---:|
| Human | 88 | 12 | 14% |
| Prokaryote | 44 | 3 | 7% |
| Eukaryote | 38 | 6 | 16% |
| Virus | 30 | 6 | 20% |

## Cross-check: tonight's own 4 proteins agree

Re-reading tonight's deployable run through the null (not the ESM-vs-blosum lens it was reported in): blosum62-vs-null was material on only **1 of the 4** (TEM-1 +47%), marginal on PTEN (+23.5%), and WORSE than guessing on RL40A (-42%) and SR43C (-51%). The 4/4 headline was ESM. So N=200 does not contradict tonight -- it reveals what tonight's headline hid: the default is weak.

## What this corrects (the load-bearing part)

The `dna-decode inverse` cell earns its keep ONLY with the learned method (ESM), which needs a GPU per protein. **The shipped wheel-only default (blosum62) is materially useful on ~1 in 7 proteins and is frequently worse than random.** So the cell's evidence must NOT let the shipped default inherit the learned method's 4/4 -- the 'per-protein check required' rail is now quantified: run the falsifier on YOUR protein first, because the default helps materially only ~13.5% of the time.

## Scope

- blosum62 ONLY -- the shipped default. ESM (the learned ceiling) needs a GPU per protein this host lacks; tonight's 4-protein sweep showed the learned oracle beats blosum on only 3/4, so this is a floor, not the ceiling.
- RANK/percentile inverse only (the deployable form). It ranks; it does not dose.
- censored assays EXCLUDED (they flatter the metric); coordinate-mismatched variants dropped.
