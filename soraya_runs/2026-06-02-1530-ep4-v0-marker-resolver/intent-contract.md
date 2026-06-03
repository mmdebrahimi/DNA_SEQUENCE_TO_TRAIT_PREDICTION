<!-- intent-contract.md for 2026-06-02-1530-ep4-v0-marker-resolver -->
# Intent Contract — 2026-06-02-1530-ep4-v0-marker-resolver

**Persona:** Soraya `--advance` (money-only gate; attended dogfood flight, run 2).
**Terminal goal (north star):** read DNA → which sequence regions map to which trait. Here: emit a pathotype-relevant call WITH provenance (which known virulence genes are present).

## Why this run
Run 1 found the learned k-mer/NT signal is confounded (study==class) + batch-prone. The confound-IMMUNE path is the ledger-locked **v0 deterministic virulence-gene-cluster resolver**: search the 24 cached genomes for SPECIFIC known marker genes (EPEC: eae/LEE + bfpA; ExPEC/UPEC: papC/sfa/afa/iutA/fyuA/hlyA/cnf1/usp/kpsMII; + EHEC stx, ETEC LT/ST, EAEC aggR for completeness). Whole-gene presence is robust to the rare-k-mer assembly artifact that drove the learned signal.

**Confound-immune biology test:** if known markers track the labels (eae/LEE in EPEC-Hazen; ExPEC adhesin/iron systems in ExPEC-Salipante), that is REAL pathotype biology — not assembler batch — because we look for specific genes, not whatever-separates-the-studies.

## Planned batch (≤ dynamic cap 8–10; this run ~4)
| # | action | gate | decision |
|---|---|---|---|
| 1 | fetch + cache VirulenceFinder E. coli DB (4942 alleles / 680 genes) | auto | DONE |
| 2 | write `scripts/pathotype_v0_marker_screen.py` (pure-Python k-mer presence, both strands) | auto | run |
| 3 | screen 24 genomes → per-genome marker profile + per-pathotype-group calls + label-tracking test | auto | run |
| 4 | result/recommendation/audit + dogfood notes + emit `/project-state` | mixed | run local; EMIT |

**Held:** `git push` (weekly-sync). **Money:** none. **Deletes:** none. **Dep-install:** none (no BLAST — pure-Python seeding; DB is data, not software).

## Stop condition
Batch complete OR money-declined (n/a) OR user-only command blocks the only path (emit + continue).
