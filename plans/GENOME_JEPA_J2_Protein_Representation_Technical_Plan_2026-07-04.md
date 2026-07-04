# Technical Plan — J2: a real protein-sequence model + the DMS falsifier (+ CLIP head)

> **Plain-language summary.** We proved (2026-07-04) that a learned model can tell which protein mutations
> are damaging, by matching lab measurements at ~0.42 (using AlphaMissense as a stand-in). Now we do it "for
> real": use an actual protein-language model (ESM-2) to score the same mutations, and re-run the *same*
> test — the number to beat is **~0.48** (the score the best public models get). If we clear it, we've built
> the first genuinely working piece of the JEPA/CLIP idea, on protein data where the biology cooperates.
> **Status:** technical plan (no build yet). Anchors on `wiki/jepa_dms_learned_signal_result_2026-07-04.md`,
> `scripts/dms_learned_model_falsifier.py`, `plans/GENOME_JEPA_CLIP_Decoder_Plan_2026-07-04.md`.

## 1. Problem statement
We have a working *test harness* (`dms_learned_model_falsifier.py`) that measures how well a model's
per-mutation scores agree with wet-lab deep-mutational-scan (DMS) measurements, across ~90 proteins, with a
shuffled-data safety check. So far the "model" was AlphaMissense (a proxy). We want to plug in a **real,
self-run protein representation** and see if it reaches the field-standard ~0.48 — turning a borrowed result
into one we own end-to-end, and standing up the substrate for the CLIP head.

## 2. Goal (one measurable bar) + non-goals
- **GOAL (crisp):** using an ESM-2 protein-language model, produce our own per-mutation scores and re-run the
  falsifier → **median rank-correlation (|Spearman|) with the DMS lab data ≥ 0.45**, on the same joinable
  assays, with the shuffled control staying ~0. Stretch: ≥ 0.48 (match the ESM2-650M field number).
- **NON-GOALS:** (a) predicting *complex organism traits* — that direction is confound-blocked and out of
  scope here; (b) beating the state of the art on the public leaderboard — matching it is enough to prove the
  representation works; (c) any wet-lab / sequence-*design* claim — that carries safety obligations and is not
  in this plan.

## 3. What we already have (reuse, don't rebuild)
- **The test harness** — `dms_learned_model_falsifier.py`: joins model-scores ↔ DMS by mutation key
  (`M1A` = wild-type + position + mutant), computes Spearman + a shuffled control + a leaderboard context, and
  a pre-registered PASS bar. **J2 only swaps the score source** (AlphaMissense → our ESM-2 scores).
- **The data (cached, on D:, free):** 217 ProteinGym DMS assays (`pg_dms/`), the leaderboard
  (`pg_spearman_dms.csv`), the entry-name↔accession map. 90 proteins are already join-ready.
- **The ML stack:** `transformers` + `torch` (CPU build) are installed; ESM-2 loads via HuggingFace
  (`facebook/esm2_t33_650M_UR50D`). No new heavy framework needed.

## 4. Design — three phases (each with its own gate)

### Phase 1 — ESM-2 "zero-shot" scoring → hit ~0.48 (the crisp goal)
**What "zero-shot" means (plain):** ESM-2 was already trained on ~65M protein sequences to fill in masked
amino acids. "Zero-shot" = we *don't train it further*; we just ask it, at each mutated position, *"how
surprised are you by this amino acid vs the original?"* A mutation the model finds very unlikely tends to be
damaging. This is the standard ProteinGym "masked-marginals" score:
`score(variant) = log P(mutant aa | context) − log P(wild-type aa | context)`, summed over the variant's
positions, from ESM-2's per-position amino-acid probabilities.

**Steps:** load ESM-2 → for each joinable protein, get its reference sequence (from `pg_reference.csv`
`target_seq`) → mask each mutated position, read off the amino-acid probabilities → compute each DMS
variant's score → feed those scores into the existing falsifier in place of AlphaMissense → read the median
Spearman.

**Compute reality (the honest gate — corrected):** ESM-2 weights are **free** (MIT license); this is
*inference only*, no training.
- **ESM2-650M** (the size that scores ~0.48): weights ~2.5 GB; inference needs ~3–4 GB VRAM (fits your
  workhorse RTX 3500 Ada / 12 GB **easily**) or runs on CPU but slowly (hours across 40+ proteins).
- **Smaller ESM2** (35M/150M, ~150–600 MB) runs on the laptop CPU in minutes but scores *lower* (the
  leaderboard shows size matters a lot) — good for a pipeline smoke-test, not for hitting 0.48.
- **So the primary gate is disk + machine, NOT money:** put the 2.5 GB weights on **D:** (C: has only ~6 GB
  free — do **not** download to C:), and run the 650M pass on the **workhorse GPU** (free, fast) or laptop CPU
  (free, slow). **Cloud/money is optional** — only if you'd rather rent a GPU for speed.

### Phase 2 — fine-tune / protein-JEPA → try to *beat* 0.48 (stretch, real GPU)
Take ESM-2 and either (a) continue training it a little on the target protein families ("domain adaptation"),
or (b) train the JEPA prototype (`genome_jepa_prototype.py`) at real scale on protein sequences, then re-run
the falsifier. **This is genuine training → needs the workhorse GPU** (or cloud = money). Leakage rule: any
fine-tuning must hold out the *test* proteins entirely (train on some protein families, test on others), or
the result is a cheat. Success = beating the Phase-1 zero-shot number by a meaningful margin; if it doesn't,
zero-shot ESM-2 is the honest product and we say so.

### Phase 3 — the "CLIP head" (point 2), with an honest framing correction
CLIP (the image↔text method) lines up **two rich modalities**. Here the two things are a protein *sequence*
and a single measured *number* (its function score) — a scalar, not a rich modality — so **literal CLIP is
the wrong tool**. The right, honest version of your idea: train a small head on ESM-2's embedding so the
embedding's *geometry reflects function* — either a **regression head** (embedding → predicted DMS score) or
a **contrastive functional split** ("functional-preserving" vs "damaging" variants pulled apart). This keeps
the CLIP *spirit* (a function-aware embedding space) without misusing the CLIP *mechanism*. Same held-out-
proteins leakage rule. Gate: GPU (small, but real training).

## 5. Implementation steps (ordered, concrete)
1. `scripts/esm_zeroshot_score.py` — load ESM-2, compute masked-marginals variant scores for one protein;
   **verify** on one protein (print a handful of scores) before scaling. [free; small model on laptop OK]
2. Extend the falsifier: add a `--score-source esm` mode that reads our ESM scores instead of `am_pg.tsv`
   (keep AlphaMissense as the baseline column for a side-by-side). [free]
3. Run Phase-1 at 650M over the 90 joinable proteins **on the workhorse GPU**; write
   `wiki/esm_zeroshot_dms_result_<date>.md` with the median Spearman + per-protein table + the shuffled
   control + a column comparing ESM vs AlphaMissense. [workhorse; free]
4. Acceptance = the falsifier's PASS verdict at the ≥0.45 bar (≥0.48 stretch). If it clears, Phase 1 done.
5. (Phase 2, only if you greenlight GPU) fine-tune/JEPA with held-out test proteins; re-run; compare.
6. (Phase 3, only if Phase 1 clears) functional-embedding head; report on held-out proteins.

## 6. Testing / acceptance
- **Offline unit tests** (`tests/test_esm_zeroshot_score.py`): masked-marginals math on a tiny fixed input
  (known probabilities → known score), the variant-key parser, and a mocked-model smoke test (no weight
  download in CI). Fast, no GPU.
- **The falsifier IS the acceptance test** (real-surface-first): the run must PASS the pre-registered bar on
  real DMS data. A green unit test alone does not count — the real ESM run must produce the real number.

## 7. Risks + mitigations
- **Disk (C: at ~6 GB):** download all weights to **D:** (`HF_HOME=D:/...`), never C:. Care-check before any
  download.
- **CPU-only laptop is slow at 650M:** run the real pass on the workhorse GPU; use small ESM2 on the laptop
  only for the pipeline smoke-test.
- **Leakage in Phase 2/3:** hold out test proteins; zero-shot Phase 1 has no leakage risk (never sees the DMS).
- **CLIP-scalar mismatch:** addressed in §4 Phase 3 (use a regression/contrastive head, not literal CLIP).
- **ESM licensing/availability:** ESM-2 is MIT, on HuggingFace — free; no money, no access gate.
- **Directionality sign:** AlphaMissense↔DMS was negative (pathogenic→low fitness); ESM zero-shot sign
  depends on each assay's convention — the falsifier already uses |Spearman| + reports the sign, so this is
  handled.

## 8. Open questions for you (the decisions that are yours, not mine)
1. **Compute path for the Phase-1 650M run:** (a) workhorse GPU — free, fast, my recommendation; (b) laptop
   CPU — free but slow (subset only); (c) cloud GPU — fast but **money**. Your call. (No money needed for a/b.)
2. **Do Phase 2 (real fine-tuning/JEPA training) at all?** It needs the workhorse GPU and only matters if you
   want to *beat* off-the-shelf ESM-2, not just match it. Can defer.
3. **Model size vs disk:** 650M (2.5 GB, hits ~0.48) vs 150M (600 MB, laptop-easy but lower score). I'd stage
   150M for the smoke-test, 650M for the real number.

## 9. One-line
Phase 1 is the whole point and it's *not money-gated* — free ESM-2 weights + your workhorse GPU, reusing the
test harness we already built. Clear ≥0.45 and the JEPA/CLIP idea has its first real, self-owned win; Phases
2–3 (beat-it + function-aware embedding) are optional GPU follow-ons.
