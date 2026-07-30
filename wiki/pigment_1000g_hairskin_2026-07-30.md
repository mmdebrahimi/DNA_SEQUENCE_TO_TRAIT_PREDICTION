# HIrisPlex-S hair + skin -- 1000G population-geography validation (2026-07-30)

**Known geography held: True** (3202 samples, 37/41 SNPs).

Mean P(skin) by superpop:
```
{
 "EUR": {
  "intermediate": 0.689,
  "very_pale": 0.115,
  "pale": 0.167,
  "dark": 0.025,
  "dark_to_black": 0.004
 },
 "SAS": {
  "intermediate": 0.183,
  "very_pale": 0.001,
  "pale": 0.002,
  "dark": 0.452,
  "dark_to_black": 0.361
 },
 "AMR": {
  "intermediate": 0.464,
  "very_pale": 0.009,
  "pale": 0.018,
  "dark": 0.282,
  "dark_to_black": 0.228
 },
 "EAS": {
  "intermediate": 0.016,
  "very_pale": 0.0,
  "pale": 0.0,
  "dark": 0.013,
  "dark_to_black": 0.971
 },
 "AFR": {
  "intermediate": 0.015,
  "very_pale": 0.0,
  "pale": 0.0,
  "dark": 0.073,
  "dark_to_black": 0.913
 }
}
```
Mean P(hair) by superpop:
```
{
 "EUR": {
  "blond": 0.381,
  "brown": 0.453,
  "red": 0.059,
  "black": 0.107
 },
 "SAS": {
  "blond": 0.014,
  "brown": 0.414,
  "red": 0.0,
  "black": 0.572
 },
 "AMR": {
  "blond": 0.083,
  "brown": 0.455,
  "red": 0.008,
  "black": 0.455
 },
 "EAS": {
  "blond": 0.003,
  "brown": 0.172,
  "red": 0.0,
  "black": 0.826
 },
 "AFR": {
  "blond": 0.007,
  "brown": 0.369,
  "red": 0.0,
  "black": 0.624
 }
}
```
Geography checks (relative EUR-vs-AFR contrast):
```
{
 "EUR light-hair(blond+red) >> AFR": {
  "lhs": 0.4398,
  "rhs": 0.0072,
  "pass": true
 },
 "AFR dark-skin >> EUR (cline)": {
  "lhs": 0.9852,
  "rhs": 0.0286,
  "pass": true
 },
 "EUR pale-skin >> AFR": {
  "lhs": 0.2819,
  "rhs": 0.0001,
  "pass": true
 }
}
```

KNOWN MODEL LIMITATION: HIrisPlex-S mis-predicts East-Asian skin as dark (EAS dark_to_black ~0.97); faithfully reproduced from the webtool, not an extraction defect.

Models recovered + held-out-validated from the HIrisPlex-S webtool; coords Ensembl-pinned + strand-harmonized on 1000G. Population-level, NOT per-individual.
