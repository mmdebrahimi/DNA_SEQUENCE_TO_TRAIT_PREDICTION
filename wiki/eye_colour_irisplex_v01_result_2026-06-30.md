# Eye-colour v0.1 (IrisPlex 6-SNP) — REAL OpenSNP result + value-over-v0 (2026-06-30, M1 complete)

The deployed IrisPlex 6-SNP model (Walsh/HIrisPlex; coefficients sourced, not fabricated) scored on the same
OpenSNP archive as v0. This is the "does the 6-SNP deployed model EARN its complexity over the 1-SNP
baseline" test — the wrapper-vs-underlying-tool rail applied to our own v0→v0.1 step. Machine sidecar:
`wiki/eye_colour_irisplex_v01_validation_2026-06-30.json`.

## Headline

| | N (binary) | accuracy | brown-sens | blue-spec |
|---|---|---|---|---|
| **v0** — rs12913832 alone (paired, complete-case) | 263 | 0.996 | 0.989 | 1.00 |
| **v0.1** — IrisPlex 6-SNP | **391** | **0.985** | 0.986 | 0.983 |

**v0.1 makes 49% MORE confident calls (263 → 391) for a ~1 pp accuracy cost.** That is the correct
engineering trade: the extra 128 calls are exactly the heterozygotes v0 abstained on.

## The rescue (why v0.1 exists)

Of **128** users v0 abstains on (rs12913832 AG → "intermediate"), v0.1 correctly resolves **123 → recall
0.961**. Label→v0.1-prediction breakdown:

| v0-abstained label | →brown | →blue |
|---|---|---|
| brown (n=123) | **121** | 2 |
| blue (n=5) | 3 | 2 |

The AG heterozygotes are overwhelmingly brown (123/128) — the established dominance of the brown allele —
and the 6-SNP model captures it (brown β 5.41 ≫ intermediate β 3.16). v0's abstention wasn't wrong, it was
*incomplete*; v0.1 completes it.

## Honest accounting

- **Complete-case:** 391 scored; **154 excluded** for missing ≥1 of the 6 SNPs on their chip (v0 needed only
  rs12913832, so it tolerated more chips). Coverage cost of the 6-SNP requirement — reported, not hidden.
- **Trade is real, not free:** overall accuracy dips 0.996→0.985 (FP 0→3, FN 1→3). Still excellent; the dip
  buys the +49% coverage. If maximal per-call accuracy matters more than coverage, v0 on homozygotes is
  slightly better; if calling more people matters, v0.1 wins.
- **Same caveats as v0:** self-reported label (near-independent, noisy); ancestry-confounded (see
  `wiki/eye_colour_ancestry_confound_2026-06-30.md` — EUR/EAS blue-allele ratio 318×; within-European
  re-score = M2-full); rs16891982 palindromic literal-forward-strand assumption (bounded).
- Tier `IN_SAMPLE_SELF_REPORT`, namespace-separate from the AMR/HIV trust surfaces.

## Verdict

The deployed 6-SNP IrisPlex model **adds real value over the 1-SNP v0**: +49% coverage at 96% rescue recall,
~1 pp accuracy cost. The off-pathogen decoder now has a v0 (high-precision, homozygote) AND a v0.1
(high-coverage, deployed forensic model) tier, both validated on real independent labels.

## Reproduce
`uv run python scripts/eye_colour_irisplex_validate.py` (reads the D: OpenSNP zip; streams only the phenotype
CSV + each user's 6 SNPs). Model: `dna_decode/data/eye_colour_irisplex.py`. Tests: `tests/test_eye_colour_irisplex*.py` (18).
