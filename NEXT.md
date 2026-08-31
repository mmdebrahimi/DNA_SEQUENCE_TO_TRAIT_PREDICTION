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

3. **Whether to compress the other 17 long CLAUDE.md bullets (~9,100 words).** The file costs ~36,800
   tokens EVERY session; 18 bullets are long AND cite a resolvable memo, so their derivations could become
   pointers without losing anything. I compressed only my own (1,022 → 309). The rest is other sessions'
   institutional memory and rewriting it wholesale is a call about what belongs in the always-loaded
   surface. Measure first: `uv run python scripts/claude_md_weight.py`.
   *Two long bullets have NO external store and must stay whole — the tool already protects them.*

4. **The 7 unscreenable colour cells** — no existing evidence tier fits (`NO_FREE_SOURCE` is about labels;
   `NOT_CENSUSED` means never-scored). And whether curating the 40 unrecorded colour loci is worth doing
   (fabrication hazard unless every locus is OMIA/literature-sourced).

## The gene-LLM idea (raised 2026-08-30, prior-art-checked 2026-08-31)

Draft anchor `wiki/idea_anchor_genomic_language_model_2026-08-30.md` (**NOT ANCHORED** — user-confirmed
skill). Prior-art check `wiki/prior_art_genomic_language_models_2026-08-31.md`: **the field is crowded.**
gLM / gLM2 / GenSyntax already build gene-token models, and four of this repo's findings (650M peak,
fluency≠function, population-structure confounding, curated-rules-beat-ML) are already published.

**Cheapest decisive experiment, no training run:** score off-the-shelf `tattabio/gLM2_650M` on our
de-confounded benchmarks against the curated-catalog baseline. If a published gene-token model can't beat
a hand-written determinant catalog on constructed variation, that answers the idea for an inference pass
— and it is exactly the comparison the critique literature says nobody runs.

**Open user questions:** token level · training objective · natural vs constructed regime (drafted answers
in the anchor).

## Cheap untried levers (executor work, no authority needed)

- ~~FBA conditional switch — continuous ratio as a ranking~~ **DONE 2026-08-29, bounded PASS.**
  Within-gene AUROC 0.7308 (non-flat, n=26, p=0.001); all-genes 0.5896 because 61% are flat. Oracle
  ceiling 11/67 exact-set vs deployed 3/67, and it ranks rather than calls. The failure is silence, not
  error. **REPLICATED on the 25-source Keio carbon axis: AUROC 0.8133 (n=69, p=0.0005), flatness 68.2%.**
  `wiki/fba_within_gene_ranking_2026-08-29.md`. **The follow-on is now CLOSED, not open:** on carbon the
  deployed rule already gets 23/217 and the oracle ceiling is 27/217 (+4 genes, +1.8pp), so estimating k
  would buy almost nothing. The 4-media "3->11" was a small-axis artifact. Do not build it.
- ~~FBA axis choice is a free lever~~ **SPENT + QUANTIFIED 2026-08-29.** All three axes measured; the
  lever is now a rule: **pick an axis whose WILDTYPE growth spreads** (distinct-growth fraction / CV,
  42 LP solves, seconds — `scripts/fba_axis_dynamic_range.py`). Flatness 61.2/68.2/75.5% tracks it
  monotonically on both summaries. Nitrogen was the worst available choice. n=3, so a direction that
  survived a common yardstick, not an established relationship.
- ~~Abstaining conditional-essentiality CLI~~ — checked 2026-08-30: there is **no shipped surface** to fix
  (`dna-fba`'s `conditional` is the variant-LOF-uncertain case; `carbon` is utilization). Shipping it would
  be a NEW phenotype-claiming surface on a published package = authority, listed above.
- **Staleness auditor** — one clean 110/110 corpus run at `TOTAL_TOKEN_BUDGET=5500` to verify the OOM
  mitigation. Still labelled unverified in `scripts/kaggle/staleness_corpus_kernel.py`.

## The FBA switch cell, as it now stands (2026-08-29)

Open, and its terms are now separated. **Direction is fine** — where the model varies it points the right
way (AUROC 0.71–0.81 on three substrates, all p≤0.001). **Accuracy when it commits is fine too** — 23/33
= 70% exact-set on carbon against a chance expectation of 0.78, ~30× chance. **Coverage is the whole
problem** — the model's call is constant for 85–94% of genes, so it is silent rather than wrong, and that
silence is structural (flat stoichiometry), unreachable by any readout change, and predicted in seconds by
the axis's own dynamic range. The readout lever is closed (+1.8 pp on the best-measured axis).

**Do not quote "10.6% exact-set" for this cell** — that scores the model against a target it cannot hit
for most genes. Quote the anatomy. The remaining bottleneck is the one already measured: the conditioning
signal is not measured in the conditions the phenotype data uses (PRECISE-1K ∩ Keio carbon = 11 of 28,
621 of 1,035 samples glucose).

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**. `scripts/project_status.py` is the authority — it
  caught two of my own figures wrong within an hour of writing them (46 traits → **44**; "~3x" → **4.1x**).
- `wiki/project_distillation_2026-08-29.md` says **46 CLI traits**. It is **44**. Left uncorrected as a
  worked example of exactly the drift this file exists to route around.
