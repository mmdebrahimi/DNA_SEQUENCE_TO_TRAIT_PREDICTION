# Eye-colour M3 — PGP cross-cohort replication (clean signal, modest N, 2026-06-30)

M3 replicates the eye-colour decoder on a SECOND, independent cohort: the Personal Genome Project
(PGP-Harvard) — the ethically-clean OpenSNP successor (open-consent, CC0 public-domain). Same v0
(rs12913832) + v0.1 (IrisPlex 6-SNP) callers, PGP's "Basic Phenotypes 2015" survey for the eye-colour label,
each participant's donated 23andMe/AncestryDNA file for the genotype. Sidecar:
`wiki/eye_colour_pgp_validation_2026-06-30.json`.

## Result (300-profile bounded sweep)

| Model | N (binary) | correct | accuracy |
|---|---|---|---|
| v0 — rs12913832 | 11 | 11/11 | **1.00** |
| v0.1 — IrisPlex 6-SNP | 15 | 15/15 | **1.00** |

**Zero errors across 26 calls on an independent cohort** — the eye-colour rule replicates. (v0.1 scores more
than v0 because it resolves heterozygotes, exactly as on OpenSNP.)

## Honest verdict: PILOT_UNDERPOWERED (pre-set gate held, not moved)

The pre-registered gate was `REPLICATES iff n≥20 AND acc≥0.90`. Observed v0 n=11 < 20 → the code returns
**PILOT_UNDERPOWERED**. I am NOT relaxing that bar post-hoc (R2 discipline). So the honest statement is: *a
clean, zero-error replication SIGNAL at modest N*, not a fully-powered replication.

**Why N is modest — the real finding:** PGP survey respondents ≠ genome donors. Of 300 eye-labelled
blue/brown participants attempted, **246 (82%) had no downloadable DTC genotype file**, and of the ~54 that
did, ~38 were WGS/Complete-Genomics formats without a parseable rs12913832 → **~5% effective linkage yield**.
The 1065-eye-labelled / 574-blue-brown survey is large, but the genotype-linked subset is small.

## Path to full power (bounded, deferred — not blocked)

A fully-powered PGP replication needs a larger sweep (~1000+ profile attempts at ~5% yield → n≥50). Not run
now: it's ~1000 requests to PGP's servers (politeness) for marginal value over the already-clean zero-error
signal. If desired, add an `--offset` to `eye_colour_pgp_validate.py` (attempt profiles N..M without
re-fetching earlier ones — genotype files already cache to D:) and accumulate. Left as a bounded follow-up.

## Honest scope
- Self-reported label (2nd independent cohort; non-circular, noisy) — same tier as OpenSNP.
- PGP is also Euro-majority → same ancestry caveat (M2 disentangled it on OpenSNP; a PGP within-ancestry
  re-score would need PGP ancestry fields — deferred).
- Real-surface validated (R3): ran on the live PGP API, produced the real artifact; harness parsers pinned
  by 3 offline tests.

## Disposition

The eye-colour decoder now has evidence on TWO independent cohorts: OpenSNP (v0 0.993 / v0.1 0.985 at large
N) + PGP (v0 11/11, v0.1 15/15, zero errors at modest N). The generalization claim — deterministic curated
rule + free independent labels, off-pathogen — holds across cohorts. Full-power PGP is a bounded follow-up.

## Reproduce
`uv run python scripts/eye_colour_pgp_validate.py --limit 300` (live PGP; rate-limited; caches to D:).
Tests: `tests/test_eye_colour_pgp.py` (3 offline parser pins).
