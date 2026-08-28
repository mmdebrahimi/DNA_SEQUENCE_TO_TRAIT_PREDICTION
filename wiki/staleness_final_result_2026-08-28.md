# The first clean 110/110 run — and the OOM was never a token budget

Closes C1 and C2 from the review in one run: the first full-corpus pass whose prompt surface matches what
the benchmark was supposed to validate (facts present in the user turn), and the first that completes.

## Result

| | |
|---|---:|
| items scored | **110/110** (no OOM) |
| unparseable | 1 |
| recall | **2/3** |
| flags | 8 |
| false positives | 6 |
| specificity | **0.944** |
| precision | 0.25 |

**Reproduction check passes.** On the 80 items the crashed runs did complete, this run returns
**recall 2/3, FP 3 — identical** to the prior `v3 @ 6000`. The fix changed the one crashing item and
nothing else, which is what a fix should do and is worth more than the completion itself.

## The OOM diagnosis was wrong twice before it was right

1. **"Truncate the excerpt."** Worked, and cost half the recall (0.667 → 0.333). Deleting evidence.
2. **"Free the KV cache per item."** Did not work — the run died at the same item.
3. **"Cap the generation."** Did not work either — raising `TOTAL_TOKEN_BUDGET` to 5500 made no
   difference, because at 5500 it binds on 9/110.

The signature that finally gave it away: **two runs died at the *identical* item with an *identical*
3.94 GiB allocation.** That is one specific input, not accumulated pressure. Item 81 is
`amr_portal_tb_disjoint_cohort.tsv` — 6000 characters of dense accession IDs, **2% whitespace, 39 distinct
characters**. Out-of-vocabulary IDs like `SAMN03648746` shred into many tokens each, so a **character cap
is not a token cap** and my 3.3-chars/token estimate was wildly wrong for that shape.

**The fix is not a bigger budget.** A claim about a data file is a claim about its *shape* — "39,193
isolates with measured DST and a leaked flag" is answered by the row count and the header, never by reading
6000 characters of accessions. Tabular artifacts are now summarised by `tabular_digest`: **613 chars
instead of 6000, and strictly better evidence for the claim being judged.** Removing the crash is a
side effect of fixing the evidence.

Exactly one artifact in the corpus is tabular, and it is the one that killed both runs.

## The 30 items nobody had ever scored

They produced 3 flags, **all three false positives**, both new ones textbook cases the prompt explicitly
warns against:

- `tests/test_genome_map_contig_collision.py` — the claim says "3 tests" and there **are** 3. The model
  latched onto "Deferred (named, not done)" elsewhere in the same bullet, which refers to pathway/KEGG and
  hmmer/Pfam, while the contig work says "SHIPPED 2026-06-28". **Keyword proximity.**
- `wiki/colour_cell_substrate_screen_2026-08-26.md` — the claim is that a question was *left open*; the
  model read the current state as contradicting it. **An open question read as a settled assertion.**
- `dna_decode/genome_map/browser.py` — the correction-text case again. Note this **did not** flip here,
  though it flipped at 3000 chars. The correction arm of v3 is **not robust**, consistent with the
  isolation run. The finding arm is.

## Honest standing

- **v3 @ 6000 with tabular digests is the measured, complete configuration.** Prompt, excerpt length and
  runtime policy have now all been exercised in one finished run.
- Recall rests on **3 positives**; it moves in steps of 0.333 and is not statistically established. The
  better-powered number is specificity 0.944 (6 false flags across 107 negatives).
- The tail's recall is **unmeasurable** — no known positives fall in it.
- Framing is unchanged: a **triage funnel**, ~1-in-4 flags real, every flag adjudicated by hand, no
  documentation edited on the model's say-so.
