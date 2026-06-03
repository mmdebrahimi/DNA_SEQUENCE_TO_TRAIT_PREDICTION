<!-- intent-contract.md for 2026-06-02-1600-ep4-v0-cli-build -->
# Intent Contract — 2026-06-02-1600-ep4-v0-cli-build

**Persona:** Soraya `--until-mvp` (money-only gate; attended).
**Terminal goal:** ship the ledger-locked v0 deterministic pathotype compatibility resolver CLI (FASTA → pathotype call + virulence-cluster provenance + abstention).

## MVP bar (frozen at run start; checkable predicates)
1. `file-exists` `dna_decode/pathotype/cli.py` (CLI entry).
2. `test-exit-0` `python tests/test_pathotype_resolver.py` (11-class decision table + abstention).
3. `file-exists` a provenance JSON produced by the CLI on a real genome (schema-complete).

## Build (per ledger `## v0 Output Contract`)
- `markers.py` — 23 marker clusters → VirulenceFinder gene prefixes; ExPEC strong/support split; supported vs scope-limited classes.
- `resolve.py` — pure 11-class decision table + abstention (HYBRID / AMBIGUOUS / UNCLASSIFIED / COMMENSAL_LOW_MARKER_BURDEN / AMBIGUOUS_LOW_QC).
- `detect.py` — k=15 seed presence (both strands), anchored clusters, located provenance (contig:pos), assembly QC.
- `cli.py` + `__main__.py` — FASTA-in → provenance JSON (contract schema) + human summary.
- `tests/test_pathotype_resolver.py` — 18 tests over the full surface.

Framing (ledger v5): COMPATIBILITY resolver, not predictor; supported = ExPEC/EPEC/ETEC; rest scope-limited/abstain.

## Gate
Money: none. Deletes: none. Dep-install: none (DB is data; pure-Python detection). `git push`: per just-expressed user intent this turn ("push … and then proceed").

## Stop condition
All 3 MVP predicates live-MET (re-checked) → `mvp-reached`.
