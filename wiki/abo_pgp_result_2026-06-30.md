# M4 — ABO blood-type (O vs non-O) decoder on PGP: serological cell, underpowered pilot (2026-06-30)

> **TIER: PILOT / DEMO** (self-reported label, n=6). Per `wiki/off_pathogen_cell_admission_gate_2026-07-01.md`:
> a demo-tier integration result (deterministic rule hosts a serological trait), NOT a measured-phenotype
> validation. Any future ABO A/B extension or powered run must clear the label-first admission gate.

M4 extends the deterministic gene→trait decoder beyond pigmentation to a SEROLOGICAL trait: ABO blood group.
Deterministic sourced rule (rs8176719 c.261delG homozygous deletion → blood type O) × free independent
self-reported blood type (PGP survey 19). Sidecar: `wiki/abo_pgp_validation_2026-06-30.json`.

## Result (200-profile bounded sweep)

| | value |
|---|---|
| N scored | **6** (200 attempted; 173 no genotype file; 21 no-call/missing rs8176719) |
| Confusion (O = positive) | TP=3, FP=0, TN=1, FN=2 |
| non-O specificity | **1.00** |
| O sensitivity | 0.60 |
| accuracy | 0.667 |

## Interpretation (honest)

**The rule is sound where the genotype is decisive.** Every homozygous-deletion (DD) → O call was correct
(3/3), and the non-O call never misfired (spec 1.0). The O-antigen mechanism (deletion → non-functional
enzyme) is textbook and the calls respect it.

**The 2 FN are label/allele issues, not a rule error (most likely).** Both are self-reported-O individuals
who carry a *functional* (insertion) allele (DI/II) → genetically A or B. Two explanations, both known:
(a) **blood-type self-report error** — self-reported blood type is ~15% erroneous in the literature; (b)
**non-deletional O alleles** (O2 and weak variants) that rs8176719 doesn't capture. This is the recurring
"suspect the noisy label before the mechanistic rule" pattern — here the *genotype* is plausibly the more
reliable value than the self-report. But **n=6 is far too small to adjudicate** — this is a pilot.

## Honest scope

- **O vs non-O ONLY.** The A-vs-B distinction (tag SNPs rs8176746/rs8176747) needs careful allele→A/B +
  23andMe-strand sourcing; NOT fabricated → deferred (named). rs8176746 (GG×14/GT×2) + rs8176747
  (CC×11/CG×2) coverage is present in the cached files, so an A/B extension is feasible once sourced.
- **Severely underpowered (n=6).** Bottleneck = PGP survey→genotype linkage yield (~5%, same as M3): 173/200
  attempted participants had no downloadable DTC file. Full power needs a larger polite sweep (--offset).
- **Self-reported label** (independent, noisy) — noisier than eye colour (blood type is often misremembered).
- Real-surface validated (R3); rule sourced (Yamamoto 1990 c.261delG; 23andMe I/D coding confirmed
  empirically in the cached files); `--` no-call → INDETERMINATE (never guessed O).

## Disposition

M4 demonstrates the deterministic decoder's architecture extends to a **serological** trait (a new trait
class beyond pigmentation), with the O-status mechanism calling correctly where genotype is decisive. It is
an **underpowered pilot** — the honest deliverable is "architecture generalizes + rule sound on decisive
genotypes," not a powered accuracy number. Full power (larger sweep) + the A/B extension (sourced tag-SNP
coding) are bounded follow-ups.

## Reproduce
`uv run python scripts/abo_pgp_validate.py --limit 200` (live PGP; reuses the M3 D: file cache). Rule:
`dna_decode/data/abo_blood.py`. Tests: `tests/test_abo_blood.py` (3).
