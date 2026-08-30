# Open threads — what a fresh session should pick up

Short by design. **Transient** state only: what is in flight, what is waiting on the user, what was just
learned that isn't durable yet. Durable findings belong in `CLAUDE.md` / `wiki/`; scope facts belong in
`scripts/project_status.py` (derived, never written).

Prune aggressively. A stale entry here is worse than an empty file.

_Last updated: 2026-08-29._

---

## Waiting on the user (authority calls — not executor tasks)

1. **v2 gentamicin lock.** The frozen AMR rule matches AMRFinder `Subclass == GENTAMICIN`, which cannot
   see 16S rRNA methyltransferases (`rmtB/E`, `npmA`) — AMRFinder files those under the generic
   `AMINOGLYCOSIDE`. Measured cost: **+0.369 sensitivity** recoverable on 131 leakage-gated disjoint
   isolates at **zero measured specificity cost** (`wiki/gentamicin_rmt_disjoint_validation_2026-08-28.md`).
   Patching the frozen surface **invalidates the prospective lock and the reproducibility freeze** — a fix
   is an unfrozen revision needing its own validation and a NEW lock. Evidence is now substantially
   stronger than when this was first raised.
   *Honest limit:* **zero S-labelled `rmt` carriers exist in any of the three datasets**, so "specificity
   unchanged" is arithmetic, not evidence. Over-calling risk is untested, not bounded.

2. **Whether a single-source cell warrants more than disclosure.** 3 of 10 SCORED AMR cells rest on one
   BioProject (`wiki/provdisjoint_source_concentration_2026-08-28.md`). Current answer is *disclose*, in a
   namespace-separate layer. Demoting them is a scope decision.

3. **The 7 unscreenable colour cells** — no existing evidence tier fits (`NO_FREE_SOURCE` is about labels;
   `NOT_CENSUSED` means never-scored). And whether curating the 40 unrecorded colour loci is worth doing
   (fabrication hazard unless every locus is OMIA/literature-sourced).

## Cheap untried levers (executor work, no authority needed)

- ~~FBA conditional switch — continuous ratio as a ranking~~ **DONE 2026-08-29, bounded PASS.**
  Within-gene AUROC 0.7308 (non-flat, n=26, p=0.001); all-genes 0.5896 because 61% are flat. Oracle
  ceiling 11/67 exact-set vs deployed 3/67, and it ranks rather than calls. The failure is silence, not
  error. `wiki/fba_within_gene_ranking_2026-08-29.md`. **Follow-on, still untried:** estimate k (a gene's
  essential-condition count) so the ranking becomes a callable rule — that is the only thing standing
  between 3/67 and 11/67.
- **FBA axis choice is a free lever** — the nitrogen axis has almost no dynamic range by construction
  (6 of 13 conditions give identical wildtype growth 0.92593). Carbon has real condition-specificity.
- **Staleness auditor** — one clean 110/110 corpus run at `TOTAL_TOKEN_BUDGET=5500` to verify the OOM
  mitigation. Still labelled unverified in `scripts/kaggle/staleness_corpus_kernel.py`.

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**. `scripts/project_status.py` is the authority — it
  caught two of my own figures wrong within an hour of writing them (46 traits → **44**; "~3x" → **4.1x**).
- `wiki/project_distillation_2026-08-29.md` says **46 CLI traits**. It is **44**. Left uncorrected as a
  worked example of exactly the drift this file exists to route around.
