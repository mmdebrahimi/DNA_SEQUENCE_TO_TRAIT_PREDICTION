<!-- recommendation.md for 2026-06-02-1500-ep4-interpret-confound -->
# Recommendation — next moves

## Stop condition reached
Planned batch (4 steps) complete. No money gate hit. One user-only command emitted (below). Lock released.

## Ranked next actions (VOI)

1. **Stop scaling the learned k-mer/NT track on this substrate.** Every bake-off number (k-mer 0.514/0.604/0.729, NT 0.38) was computed on a study==class subset → none is interpretable as pathotype biology, and the top features look like assembler batch (CTAG-motif rare k-mers). Do NOT invest in per-CDS NT or more pooling variants here.

2. **Advance the ledger-LOCKED v0 deterministic virulence-gene-cluster resolver** (the interpretable north-star deliverable). It matches KNOWN markers (eae/LEE, bfpA, stx1/2, LT/ST, AAF/aggR, papC/afa/hlyA, …) → immune to the study==class confound because it looks for specific biology, not whatever-separates-the-batches. This is also the v0 the project already designed (23 markers + 11-class decision table in the ledger `## v0 Output Contract`). Needs a marker reference set (VirulenceFinder/ETECFinder DB or a curated probe set) — a real dependency to resolve next.

3. **To ever trust a LEARNED number, break study==class:** obtain ExPEC + EPEC from a SINGLE study, or ST/lineage-matched cross-study pairs. Until then, no learned AUROC on this data means biology.

## EMIT (user-only `/project-state` — Soraya cannot self-invoke; run when convenient)
Record the confound + override as a project decision:
```
/project-state ecoli-pathotype-prediction-cli-2026-05-26 --append-decision \
  "EP-4 bake-off substrate is study==class (ExPEC=Salipante/EPEC=Hazen); the 0.729 presence/absence signal is sparse (~10 k-mers) but those are rare CTAG-motif k-mers with partial presence => assembler batch suspected, not biology. Learned-track numbers on this subset are uninterpretable as pathotype biology. Decision: stop scaling learned k-mer/NT here; advance the v0 known-marker virulence-gene-cluster resolver (confound-immune) and obtain within-study/lineage-matched genomes before trusting any learned AUROC."
```
Also worth an `--update-hypothesis` if a 'learned > classical for pathotype' hypothesis is tracked → mark refuted-on-this-substrate.

## Held for user
- `git push origin main` — 2 local commits ahead (`2d878bb` + the pending Soraya-run commit). Push when you next sync; I did not push (your ~weekly cadence).
