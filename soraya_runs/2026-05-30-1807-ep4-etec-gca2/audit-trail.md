# Soraya Run Audit — 2026-05-30-1807-ep4-etec-gca2

| # | action | gate | result |
|---|---|---|---|
| 1 | WebFetch PMC8085198 for supplementary file URLs | auto | DONE — MOESM4/5 = XLSX |
| 2 | Download MOESM4/5 (PMC bin 403'd curl → Springer CDN) + install openpyxl (uv) | auto / irreversible-dep-install (un-gated) | DONE — real XLSX (PK), openpyxl 3.1.5 |
| 3 | Parse MOESM5 → 8 reference chromosome accessions | auto | DONE — E925→LR883050 etc. |
| 4 | Parse MOESM4 → 558-genome ETEC GCA collection; write CSVs; update memo+ledger | auto-edit | DONE |

## Outcome: PASS
ETEC arm fully resolved by-accession: 8 complete reference chromosomes (LR) + 558 GCA-accession collection genomes. Closed the PARTIAL from run 1135.

## Notes / deviations
- PMC /bin/ supplementary 403s direct curl → used Springer static-content CDN mirror (authoritative same files, sizes matched listing).
- uv venv lacks pip → used `uv pip install openpyxl` (dep-install; gate=irreversible, un-gated default, money-free).
- /project-state path-gated from cwd → ledger via direct Edit (surfaced).

## Stop reason: batch complete (4/4 under cap), ETEC arm closed, no further money-free in-cwd EP-4 action (remaining = workhorse caller exec + user-only Gate B).
