<!-- recommendation.md for 2026-06-02-1530-ep4-v0-marker-resolver -->
# Recommendation — next moves

## Stop condition reached
Batch complete (DB fetch + screen built + run + analysis). No money gate. No deletes. `/project-state` emitted (below). Lock released.

## What this run established
The ledger-locked **v0 deterministic virulence-gene-cluster resolver is the right architecture and it works**: known markers separate EPEC vs ExPEC perfectly (eae AUROC 1.0), interpretably (per-gene provenance), and confound-immune (specific genes, not batch k-mers). The learned k-mer/NT track is decisively deprioritized on this substrate.

## Ranked next actions (VOI)
1. **Promote the screen into the v0 CLI resolver** (the project's actual deliverable). `pathotype_v0_marker_screen.py` already has the marker→family→group map + presence engine; wrap it in the ledger's `## v0 Output Contract` (FASTA-in → multilabel cluster profile → 11-class decision table + abstention → provenance JSON). Decision logic is standard: `LEE+ ∧ stx− ∧ bfp+ → typical EPEC`; `LEE+ ∧ stx− ∧ bfp− → atypical EPEC`; `LEE+ ∧ stx+ → EHEC`; `LT+ ∨ ST+ → ETEC`; `aggR+ ∧ AAF+ → EAEC`; `(pap ∨ sfa ∨ afa) ∧ (iut ∨ fyu) → ExPEC/UPEC`; else commensal/abstain. Add hybrid handling (multilabel, not forced single class).
2. **Harden detection**: the k=15 exact-seed COV≥0.80 is strict (misses divergent alleles). For v0.1 consider lowering k or adding a gapped/identity-tolerant pass; benchmark against CGE VirulenceFinder on a few genomes (the ledger's required side-by-side diff).
3. **Generalize the substrate test**: run the resolver on the ETEC (von Mentzer) + any EHEC/EAEC genomes already resolved, to exercise the full 11-class table beyond the ExPEC/EPEC pair.
4. **(Lower priority) learned track**: only revisit if a within-study / lineage-matched dataset materializes; on current data it is uninterpretable. Per-CDS NT is NOT worth it given the marker resolver already hits AUROC 1.0 interpretably.

## EMIT (user-only `/project-state` — Soraya cannot self-invoke; run when convenient)
```
/project-state ecoli-pathotype-prediction-cli-2026-05-26 --append-decision \
  "v0 marker resolver VALIDATED on 24 genomes: VirulenceFinder marker presence (pure-Python k=15 seeding) separates EPEC vs ExPEC perfectly and interpretably (eae 12/12 EPEC vs 0/12 ExPEC, AUROC 1.0; ExPEC virulence-score AUROC 0.882). Confound-immune (specific known genes, not batch k-mers), unlike the learned k-mer 0.729. Decision: promote scripts/pathotype_v0_marker_screen.py into the v0 CLI per the ledger v0 Output Contract; deprioritize the learned k-mer/NT track on this study==class substrate."
```
Consider `--update-hypothesis` on any 'deterministic resolver viable' hypothesis → mark supported.

## Held for user
- `git push origin main` — local commits ahead (now incl. this run). Push at your weekly sync; not pushed.
