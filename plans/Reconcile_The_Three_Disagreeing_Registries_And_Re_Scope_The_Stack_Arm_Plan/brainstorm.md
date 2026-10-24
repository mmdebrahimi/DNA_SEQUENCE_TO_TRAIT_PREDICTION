# Brainstorm — the VOI assessment that became this plan (2026-10-23)

> Note: the independent reviewer CLI was unavailable on this host (not installed; verified absent from
> PATH in both bash and PowerShell, no npm shim, no `~/.codex`), and LM Studio was unreachable. This
> review is single-model, which raises anchoring risk — so the scrutiny was pointed at the executor's
> OWN recommendation and every load-bearing claim was verified against the artifacts. It did not
> survive.

**Subject:** the claim that a `fragments_LB` learning curve was the highest-VOI next move, plus two
reported defects. **Depth:** deep · **Rounds:** 1 (reviewer unavailable) · **Stop reason:** graceful
degradation; self-critique verified against artifacts instead.

## 1. Critical issues

**C1 — The row-count argument is close to vacuous, and the refutation is inside the cited artifact.**
`[grounded]` 321,123 fragments × median 244 bp = **78.4 Mb of sequence over a 4.64 Mb genome — 16.9×
redundancy**, mean step between fragment starts **14.5 bp** against median length **244 bp**. The split
has **10 genomic blocks**, 2 held out. The existing 30,000-row subsample already spans the whole
genome; going to 320k adds **depth at the same positions, not new blocks**. The census's `size_lever`
("never been pulled") is literally true about N and nearly uninformative about available information.
The "639,580 rows at ~9% utilisation" framing treated redundant rows as independent observations.

**C2 — The cross-rung table was apples-to-oranges.** `[grounded]` It compared **0.3968 on tiles** (a
*biophysical* promoter-calculator scalar) against **0.1313 on fragments** (a *learned* k-mer ridge). On
the same split and seed — GC is **0.3000 in both artifacts**, confirming identity — the like-for-like
learned arms are **0.2769 (29.4%) on tiles vs 0.1313 (16.6%) on fragments**. The reported 42.2% vs
16.6% gap is really 29.4% vs 16.6%.

**C3 — The two rungs ask different questions, and the artifact says so in a field name.** `[grounded]`
`why_this_is_the_like_for_like_arm` states peak tiles ask *"how strong is this promoter"* while the
fragment library asks a *"rare-event detection question"*. Ranking them by "% of own ceiling" is not
like-for-like. Related: tiles permit positional-onehot (fixed 150 bp); fragments cannot without an
anchor and *"no anchor is obviously right"*. The fragment rung's binding constraint is
**representational**, not sample size.

**The correction these produce — the ladder inverts the ranking.** On peak tiles, same split:

| arm | Spearman | % of 0.941 ceiling |
|---|---|---|
| promcalc + conv, 2-channel | **0.4254** | 45.2% |
| promcalc `max_tx` (biophysical) | 0.3968 | 42.2% |
| promcalc full decomposition | 0.3925 | 41.7% |
| **GC (1 feature)** | **0.3000** | **31.9%** |
| kmer4 (256) | 0.2769 | 29.4% |
| onehot (600) | 0.2691 | 28.6% |
| kmer6 (4096) | 0.2022 | 21.5% |

**Every learned sequence feature loses to a single GC number.** The only arms that beat GC are
orthogonal-modality ones, and the 2-channel hybrid adds **+0.0286 over biophysics alone, 10/10 seeds,
CI excluding zero, with a shuffled-conv control that HURTS (−0.1232)** — so the gain is channel
*content*, not averaging.

## 2. Medium issues

**M1 — The learning curve as specified tests the wrong arm.** `[inferred]` kmer4 is 256 features at
n_train 24,006 — ~94 samples per feature, not data-starved. The arms that could be variance-limited
are kmer6 (4096) and onehot (600), and both lose to kmer4 and to GC at current N on both rungs.

**M2 — The defect-1 claim was stronger than the evidence.** `[grounded → softened]` Solid: the
census's `fragments_M9.fraction_of_ceiling` of 0.1656 is `0.1313 / 0.7928` = **LB's** ceiling, not
M9's 0.8014 — definitively mis-normalised. Not solid: "M9 is unmeasured." The artifact's `substrate`
field says "LB + M9" and the single 30,000 subsample's pool is not stated, so whether M9 rows entered
it is **undetermined**.

**M3 — A stale blocker note, found during verification.** `[grounded]` The oracle artifact's own
`post_hoc_corrections` says regenerating it *"needs the absent `D:/dna_decode_cache/mpra`
substrate."* That substrate is now on C: at `data/mpra/` (96 MB). Feasibility was never the problem
with the recommendation; the rationale was.

## 3. Deferred

- `UNDECLARED DOCUMENT VALIDITY 694/694` — no `doc-status` headers. M2 and M3 are both instances of
  what that warning predicts; fixing it repo-wide is its own project.
- The essentiality third rung — label file absent from this host.
- The missing `lmstudio-auditor.sh` in the skills install.

## 4. What's solid

- **Defect 2 verified:** both `scripts/glm_finetune_window_split.py` and the 101 KB
  `wiki/glm_finetune_window_split.json` exist, so the pending candidate row really is already shipped
  and really would be redone.
- **The census does mis-attribute a number** across two rungs with different ceilings, in a document
  feeding the mandated orientation script — aggravated by that same artifact's
  `CONDITION_CARRIES_SIGNAL` verdict (noise-matched gap −0.0614), which is exactly why LB cannot stand
  in for M9.
- The hard constraints (CPU-only, unreliable drive, frozen units) were correctly applied.

## 5. Open tradeoffs

**Is the fragment size lever worth a kill-test at all?** Not a bet either way: minutes of CPU and the
substrate is local, but C1 makes a flat curve near-certain, so it mainly buys closure.
**Should the fragment rung be re-aimed rather than re-sized?** Its one positive is that kmer4 *beats*
GC there (+0.0290) where on tiles every learned arm loses — which argues the interesting question is
the anchor-free feature representation, not N.

## 6. Recommended next step

Reorder to put the measured-positive first: the 3-or-more-channel orthogonal stack on peak tiles. Then
the `[done:…]` marker. Then, only if closure is wanted, the fragment curve framed as a kill-test at
`kmer6`/`onehot` rather than `kmer4`, with its headline stated in blocks rather than rows.

---

## POSTSCRIPT — what happened next, and why this brainstorm's own recommendation also fell

Written into this sidecar because the sequence is the lesson, not the conclusion.

While drafting the technical plan, **both of this review's top recommendations turned out to be already
executed**:

- The fragment size lever was pulled on **2026-10-18** and recorded in the arms ledger as
  `fragment-rung-321k-rows-size-lever-LB-M9`, `status: measured`, `metric: 0.185` —
  `DOES_NOT_GENERALISE` at both 20k and 150k, with the margin over GC **shrinking** with more data
  (+0.0288 → +0.0191). A real sub-finding points the other way and is worth keeping: `kmer6` DOES gain
  from size (+0.0749 at 20k → +0.1172 at 150k), confirming M1's capacity reading, but it never reaches
  kmer4's 20k value and GC rises alongside.
- The 3-channel stack was published on **2026-10-14** as `hybrid_conv_maxtx_gc`: median **0.4217** over
  10 seeds against the 2-channel **0.4254** — it did not stack, it cost −0.0037. The combo was already
  in `scripts/glm_tile_hybrid_arm.py::COMBOS`, and the artifact predates the arm's own declaration by
  five days.

So three successive recommendations in one session all pointed at finished work. **That — not the
substrate question — is what the plan addresses**, and it is why the plan contains no new measurement
arm.

**One further reversal, recorded because it was measured rather than argued.** The plan's original
Step 5 proposed two string-level registry guards. Both were run against the historical state they were
designed to catch: one fired on **nothing** (rung key `fragments_LB` vs arm `fragment-rung-…` — a
singular/plural mismatch), the other on **all four** planned arms. They were withdrawn and replaced by
a referential-integrity check. The most useful thing this review chain produced was disproving a guard
*before* it shipped green and silent through exactly the drift it was written to catch.
