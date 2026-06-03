<!-- result.md for 2026-06-02-1600-ep4-v0-cli-build -->
# Result — v0 pathotype compatibility-resolver CLI SHIPPED (MVP reached)

Verdict: **mvp-reached** — all 3 frozen MVP predicates live-MET.

## MVP predicate evaluation (re-checked at stop)
| # | predicate | status |
|---|---|---|
| 1 | `dna_decode/pathotype/cli.py` exists | MET |
| 2 | `python tests/test_pathotype_resolver.py` exit 0 | MET — **18/18** decision-table tests pass |
| 3 | CLI emits schema-complete provenance JSON on a real genome | MET — 2 genomes, all 7 contract keys + full marker-hit fields + DB sha256 |

## What shipped (`dna_decode/pathotype/` package)
- `markers.py` — 23 marker clusters → VirulenceFinder gene prefixes; ExPEC strong/support split; supported (ExPEC/EPEC/ETEC) vs scope-limited classes.
- `resolve.py` — pure 11-class decision table + abstention. 18 unit tests, all green.
- `detect.py` — k=15 dual-strand seed presence, anchored clusters (LEE→eae etc.), located provenance (contig:start:strand), assembly QC.
- `cli.py` / `__main__.py` — `python -m dna_decode.pathotype <fasta> [--gff3] [--db] [--out]` → provenance JSON + human summary.

## End-to-end demonstration (the north-star deliverable)
- **EPEC AIEY01 → `tEPEC_COMPATIBLE` [CONFIDENT | supported]**, driven by `eae` (alpha-1 intimin subtype, 98.9% cov, ENA|AIEY01000062:10052) + `bfpA` (100%). FASTA → pathotype call → which gene drove it, with location. Exactly the goal.
- **ExPEC JSIS01 → `AMBIGUOUS`** ("single strong ExPEC marker; UPEC needs ≥2"), papC + iron/capsule support. On-spec abstention (conservative UPEC rule), NOT a bug — abstention-first design working.

## Honest characterization (not a bug — a v0.1 calibration item)
The contract's UPEC rule (≥2 strong fimbrial/toxin markers) is conservative, so a real ExPEC with one fimbrial system + iron/capsule support abstains to AMBIGUOUS. This is the intended abstention-first behavior, but it under-calls ExPEC. Calibrating ExPEC sensitivity vs abstention rate is exactly ledger Goal 5 / H4 (abstention-quality eval). Also: detection is k-mer-seed COVERAGE (a presence proxy), not BLAST %identity — `percent_identity` is null with `method=kmer_seed_coverage_v0`; a true CGE VirulenceFinder side-by-side diff (ledger requirement) is a v0.1 add needing the VF software.
