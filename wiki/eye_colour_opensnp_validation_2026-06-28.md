# Eye-colour decoder — REAL OpenSNP validation result (2026-06-28)

> **TIER: PILOT / DEMO** (self-reported label). Per the admission gate
> `wiki/off_pathogen_cell_admission_gate_2026-07-01.md`, self-report cells are demo-tier: they show the
> decoder can host deterministic non-microbial rules with honest abstention — NOT a measured-phenotype
> validation. "Generalizes" is reserved for measured labels; this is the engineering/integration claim.

The flagship off-pathogen cell, validated on real data. The deterministic gene→phenotype decoder, proven
across the AMR/typing/viral surface, **generalizes to a human visible trait** — scored against FREE,
INDEPENDENT, INDIVIDUAL-LEVEL labels (OpenSNP self-reported eye colour). Machine sidecar:
`wiki/eye_colour_opensnp_validation_2026-06-28.json`.

## Headline

**rs12913832 single-locus v0 (strand-agnostic): accuracy 0.993, brown-sensitivity 0.98, blue-specificity
1.0 on N=298 binary-scored users** (TP=100 brown, TN=196 blue, FP=0, FN=2). Source: the Internet Archive
OpenSNP 2017-12-08 dump (~20 GB; live site deleted 2025-04). User-ratified human-data run.

## The load-bearing scope caveat (this is NOT a 99% "eye-colour predictor")

The number is on the **DECISIVE subset only** — rs12913832 *homozygotes*. The rule calls homozygous
GG/CC→blue and AA/TT→brown; **heterozygotes (AG) abstain as "intermediate"** by design. The strata make
this explicit:

| self-reported | →blue | →brown | →intermediate (abstain) | indet |
|---|---|---|---|---|
| **blue** (n=203) | 196 | 0 | 6 | 1 |
| **brown** (n=249) | 2 | 100 | 146 | 1 |

So of 452 blue/brown-labelled users with a callable genotype, **298 (66%) get a confident binary call and
the rule is 99.3% right on them**; the other **152 are heterozygotes the single-locus rule abstains on**.
This is the SAME abstention discipline as the AMR decoder (confident where the determinant is decisive,
abstains where it isn't) — high accuracy is bought with honest abstention, not by over-calling.

Why so clean: rs12913832 (HERC2/OCA2) is the single best-established human visible-trait variant (~74% of
blue-vs-brown variance alone), and at the *homozygous* genotype it is near-deterministic. A ~99% number on
homozygotes is the EXPECTED textbook result, not a label-noise artifact (it is high in BOTH directions —
not the high-sens/low-spec shape that signals a bad label).

## The two honest failure modes (the v0.1 levers)

1. **146 brown heterozygotes abstain.** Among AG genotypes, self-reported colour is overwhelmingly brown
   (146 brown vs 6 blue) — brown is dominant. A v0.1 that calls AG→brown would recover ~146 brown calls.
   The principled v0.1 is the **IrisPlex 6-SNP model** (Walsh 2011 — adds rs1800407/rs12896399/rs16891982/
   rs1393350/rs12203592); coefficients must be SOURCED, not fabricated (deferred, named).
2. **FN=2 / FP=0:** two brown-labelled users are GG/CC (predicted blue). Real het-of-other-loci override
   or label noise — not a systematic error (blue-spec is a perfect 1.0).

## Coverage accounting (honest denominators)

- 923 users self-report an eye colour (the real CSV has 6 "eye"-containing columns — **R3 caught that the
  parser must pick exact "Eye color" (idx 7), NOT the decoy "Eye pigmentation" (idx 4)**; regression-pinned).
- 377 excluded as non-binary (green/hazel/other — a real third class, reported not hidden).
- 93 had a genotype file but no parseable rs12913832 (chip/format gap — e.g. 23andme-exome-vcf; a parser
  coverage item, NOT a miscall — missing → excluded, never guessed).
- 298 binary-scored (the headline N).

## Significance

This is the **first off-pathogen cell with a real independent number**. It demonstrates the project's
central claim at a new boundary: a deterministic curated rule + a free independent measured/observed label
= a validatable decoder cell — now shown for a HUMAN VISIBLE TRAIT, not just AMR. The label is
self-reported (near-independent, non-circular, noisy), so this is a strong proof-of-generalization, not a
forensic-grade calibration. **Honest tier:** `IN_SAMPLE_SELF_REPORT_DECISIVE_SUBSET` — namespace-separate
from the AMR/HIV trust surfaces (do NOT fold into the frozen AMR report card).

## Caveats (carried in the JSON)
- self-reported label (not a lab assay); ancestry-confounded (rs12913832 European-calibrated → a
  within-ancestry split using 1000G/HGDP allele freqs is the v0.1 confound control);
  intermediate/green/hazel excluded from the binary; 2017 dump; single-locus v0.

## Reproduce
`uv run python scripts/eye_colour_opensnp_ingest.py` (reads the committed-absent 20 GB archive zip on D:;
streams only the phenotype CSV + each eye-colour user's genotype member — no full extraction, no network).
Rule: `dna_decode/data/eye_colour.py`. Tests: `tests/test_eye_colour.py` + `tests/test_eye_colour_opensnp_ingest.py` (11).
