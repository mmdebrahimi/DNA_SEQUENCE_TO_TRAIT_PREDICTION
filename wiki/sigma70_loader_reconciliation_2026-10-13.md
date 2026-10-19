# A second σ70 loader reproduced the exact failure the first one was written to prevent

**2026-10-13.** Machine-readable twin: `wiki/sigma70_loader_reconciliation_2026-10-13.json`.

## DECISIONS NEEDED

None. This is a correction of a defect introduced on 2026-10-12; nothing here is an authority call.

## CORRECTIONS

Three numbers shipped wrong on 2026-10-12. All are now corrected in place, in both the markdown and
JSON twins of the census, plus a correction banner on the plan that rested on them.

| claim shipped 2026-10-12 | truth, measured | consequence of the error |
|---|---|---|
| σ70 rung has **3,821 usable** of 10,898 (35.1%) | **10,898 usable, join rate 1.00** | the rung looked like a third of its real size |
| the 65% gap is **absent sequence data** | a **hyphen-vs-underscore** separator mismatch | an open question was filed to hunt for data that was never missing |
| design is **8 × 8 × 8 × 8**, `UP_element` collapses to 1 level | **3 × 8 × 8 × 8 × 8 = 12,288 cells**, 10,898 filled; `UP_element` has **3** | one whole design axis declared constant |
| fragment ceiling **0.4583 / 0.4731** | **0.7928 / 0.8014** | reversed a strategic recommendation — see below |

## Why a trap that was already found, fixed and test-pinned happened anyway

The user asked this directly, and the answer is a mechanism, not carelessness.

**The protection was INSIDE the function that got bypassed.** `MIN_JOIN_RATE = 0.95` guards callers of
`load_pairs`. Writing a second reader of the same files is the one action that defeats every guard
attached to the first one — and it is also the action that feels like making progress, so it does not
feel like a bypass. The knowledge was perfect and sat on a path not taken.

**"It was hard to find" is false, and that was measured.** Of ten plausible search terms, **nine
locate `expression.py`**: `sigma70`, `GSE108535`, `hyphen`, `underscore`, `norm_key`,
`barcode_mapping`, `join_rate`, `naive_join`, `loader`. The only miss is `constructed` — the word
chosen for the duplicate module. **The search ran in the vocabulary of the concept being built
("constructed variation", from the regime framing) rather than of the data being loaded (the
accession).** Searching the data's nouns would have landed on the owner immediately.

### The structural fix, and what it does not cover

`tests/test_substrate_loader_ownership.py` moves the guard off the code path and onto the **data**: a
GEO accession may be named in exactly ONE module. A second reader has to name the files to open them,
so it trips the test instead of silently re-deriving a solved trap. The failure message names the
owner, the entry point to call, **and the trap that owner encodes** — a message reporting only a
collision would leave the next session to re-derive the fix, which is the whole failure mode.

**Proven non-vacuous against the real historical shape:** recreating the duplicate's two filename
constants in the package fails the guard and names the owner; removing them returns it to green.

Legitimate prose cross-references are handled by a **visible one-line exemption list**
(`ALLOWED_MENTIONS`), not a narrower regex — `genomewide.py` mentions GSE108535 six times in its
docstring (that is where the 0.7925 ceiling is recorded) and `encoder.py` mentions GSE144621 once.
A guard that fired on those would be the documented trap where guards fire on the comment explaining
what they ban.

**Two defects in the guard itself, both found by running it, neither by review** — and both are the
same error as the one being fixed, recurring inside the fix:

1. The first draft **asserted** that a bare `"GSE108535"` string was caught. It was not: the regex
   then required `GSE\d+_<filename>`, so a module naming the accession once and building filenames
   from it evaded the guard entirely. A claim about the guard's own coverage, written without running
   it.
2. The second draft anchored the pattern as `\bGSE\d+\b`, which **cannot match
   `GSE108535_sigma70...` at all** because `_` is a word character — so it caught the bare accession
   and missed every real filename, exactly backwards from what is needed.

**Remaining hole, recorded not closed:** a module that never writes the accession digits — importing
the owner's filename constant, or composing `f"GSE{series}_..."` from a number held elsewhere — still
evades. Closing it needs import-graph or AST analysis of what reaches `gzip.open`, a much more brittle
instrument. The measured case is a literal copied between modules, and that is now caught two ways.

### Is the same collision live anywhere else? Measured: no

Two registry entries are not an answer to "how do we avoid this", so the whole package was swept for
data-file literals (AST string constants whose basename carries a data extension). **48 distinct
literals; 11 are named by more than one module.** Each of the two risky patterns was then checked
rather than assumed:

- **A cell's `cli.py` + its own `runner.py`** (`wzi.fasta`, `cps_references.fasta`,
  `salmonella_antigens.fasta`) — the CLI names the file to locate it and hands the path to its own
  runner. One reader, two mentions.
- **`profile/cli.py` sharing four DBs with the cells that own them** (`enterobacteriales.fsa`,
  `virulence_ecoli.fsa`, `serotypefinder.fsa`, `resistens-overview.txt`) — this was the one that
  could have been the ResFinder locus-collapse bug again, an aggregator re-implementing four callers'
  selection logic. It **delegates**: it imports `plasmid.runner.call_replicons`,
  `serotype.runner.call_serotype`, `pointfinder.runner.call_point_mutations` and
  `resfinder.runner.call_resistance_genes`. That is exactly what the new guard's failure message tells
  a future author to do, already done here.

**`main.tsv` has SEVEN independent readers and is the one case where that is correct** — it is
per-genome AMRFinder *output*, not a fixed dataset, so each cell parsing its own genome's output is
the intended design. One-owner-per-substrate is a rule about **named external datasets**; applying it
to tool output would be wrong, which is why the registry is keyed on accessions.

**A near-miss worth recording, because the resolution is a column name.** `eval/amr_rules.py:54`
filters `main.tsv` on `Method` containing `POINT` — the frozen cipro rule's entire input — while this
repo's own notes record *zero POINT rows in `main.tsv`*. Those read as a contradiction and are not:
measured across the real cached runs, `Method` is `{EXACTX 13, POINTX 12, ALLELEX 10, POINTN 7,
INTERNAL_STOP 1, BLASTX 1}` and `Type` is `{AMR 43, STRESS 1}` — **19 POINT-in-`Method` rows, 0
`Type == POINT` rows, and 7 QRDR rows the frozen rule counts.** Both notes are true of different
columns. Checked on real output rather than the fixture, because `tests/fixtures/amr_mini/main.tsv`
is hand-written and could have been built to match the code it tests. No action: `qrdr_point` is
pinned by six test files and POINT methods appear in ten, so "cleaning up" the apparently-dead filter
— which would silently zero the frozen cipro rule — fails loudly already.

## The root cause, and why it is worse than an ordinary bug

`dna_decode/glm/expression.py` already existed. It is the canonical σ70 loader and it already:

- carried `_norm_key`, the hyphen/underscore fix, **with the mutation arrow `->` protected** so a blanket
  `replace("-", "_")` cannot rewrite `35T->A` into `35T_>A`;
- enforced `MIN_JOIN_RATE = 0.95`, so a degraded join **raises** instead of returning a third of the data;
- documented the trap in its own module docstring, in the first 25 lines;
- had `tests/test_glm_expression.py` pinning **`assert rep.naive_join == 3821`** with the comment
  *"the naive-join counterfactual is the point of the fix"*;
- had a second live assertion that every element column carries **≥3** levels — which the
  "`UP_element` collapses to 1" claim directly contradicted.

A second loader, `dna_decode/glm/constructed.py`, was written beside it and performed the naive join.
**So the correct figure (10,898), the wrong figure (3,821), the cause, the fix, and two green tests
asserting both were all already in the repository when the wrong number was published.** This is not a
case of missing information. It is a case of not looking for the thing that already existed.

It also duplicated the **split**: `constructed.py::held_out_background_split` could hold out only
`Background`, where the canonical `expression.py::leave_element_out` holds out element identities on
**any of the five** axes. The duplicate was narrower on the axis that matters most.

## What was done

- **`dna_decode/glm/constructed.py` and `tests/test_glm_constructed.py` DELETED.** Neither had ever been
  committed, so git history never carried the bad loader. No module imported `constructed` outside its own
  test, so the deletion is clean.
- **`irrelevant_feature_control` — the one original contribution — MOVED into `expression.py`**, beside the
  splits it validates, and generalised from Background-only to any element column. It keeps its
  non-vacuity design: it reports a planted leak's score alongside the irrelevant feature's, because a
  control that reads "clean" without being able to read "dirty" certifies nothing.
- **Six tests added** to `tests/test_glm_expression.py`: the control at chance, the control non-vacuous,
  the control on three different columns, a single-level column yielding **zero** folds rather than one
  degenerate fold, the real-data control, and
  `test_the_design_is_FIVE_way_crossed_and_UP_element_is_NOT_constant` — which pins
  `UP_element == 3` and `3·8·8·8·8 == 12288`, i.e. it fails on exactly the claim that shipped.
- **Census corrected** (`wiki/mpra_substrate_census_2026-10-12.{md,json}`) with a `corrections_2026_10_13`
  block naming each wrong field, its replacement, and the cause.
- **The ambiguous JSON field was removed, not just re-valued.** `replicate_ceiling` is gone from every
  rung, replaced by the two explicit fields `single_replicate_agreement` and
  `ceiling_for_a_sequence_model`. The conflation was only possible because one key stood for two
  quantities; re-valuing it would have left the trap armed.
- **Plan banner** on `plans/Promoter_Activity_Substrate_Ladder_And_The_Mechanistic_Baseline_Plan/`,
  with the body left **unedited** as the record of what was planned on a false premise. Step 5 is marked
  VOID.

## The ceiling error reversed a recommendation, and that is the costlier half

`noise_ceiling()` returns two different quantities:

- `single_replicate_agreement` — how well replicate 1 tracks replicate 2 (LB **0.4583**);
- `ceiling_for_a_sequence_model` — the Spearman-Brown step-up of that, square-rooted, because the
  published readout is the **mean** of two replicates and is therefore more reliable than either
  (LB **0.7928**).

The census read the first as the ceiling and concluded the fragment rung's headroom "falls from 0.941 to
~0.46-0.47" — a **steep** trade. It actually falls to **~0.79-0.80**: 7.3× the rows and 3.3× the median
length for roughly **85%** of the headroom, not 49%.

**That flips the advice.** A substrate offering 7.3× the data at 85% of the ceiling is attractive; at 49%
it was not. The wrong number did not merely understate a figure, it argued against the right next move.

Peak-tile values for comparison, measured on the 44,106 filtered rows: single-replicate **0.7945**,
sequence-model ceiling **0.9410**.

## Also corrected: the arms ledger named the wrong best arm

The ledger reported `conv-3width-13k-params-150bp` at **0.3234** as best, while the biophysical Promoter
Calculator had measured **0.3968** the previous day and was never registered. Now recorded with its full
evidence string. `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` is **declared** (`planned`,
medium) so it cannot quietly evaporate.

The criterion stays `blocked:untried-arm` — correctly, since three declared arms remain unrun. A
criterion with unrun declared arms must not read as plateaued.

**`measure_arm` mutates memory only; `save()` is separate.** A first attempt registered both arms and
exited without saving, and the verdict read unchanged — a silent no-op that looks exactly like a
successful write.

## Why this was NOT filed as an "arm"

An arm in this project requires a metric: `measure_arm(criterion, name, metric: float, evidence)` will not
accept anything without a number. Arms are competing attempts to predict promoter activity, ranked against
the 0.9410 ceiling. **A loader fix produces no metric and has no competitor**, so filing it as an arm would
have put a bug-fix into the ledger that decides which predictor is best, and left the criterion
`blocked:untried-arm` on a maintenance ticket. The mechanism that tracks a correction is a test plus an
artifact, which is what this is.

## Honest limits

- **The test additions are not blinded.** They were written after the failure was known, so they prove the
  current state is right; they are not evidence that the guard would have caught an unanticipated variant.
- **`measured_headline` for the σ70 rung (one-hot ~0.59 vs GC ~0.14) has no denominator.** GSE108535's two
  published tables correlate at exactly 1.0 — one is a monotone transform of the other, not a replicate
  pair — so no ceiling is computable and 0.59 cannot be read as a fraction of achievable.
- **Circularity on this rung is UNRESOLVED and is now a named open question.** The Promoter Calculator's
  346 parameters and PromoterAtlas's `urtecho2018` expression head may both have been fit on GSE108535.
  Either would be a circular comparator here. Neither has been checked.
- ~~Only the σ70 rung's loader was reconciled; the same pattern may exist elsewhere.~~ **AUDITED the
  same day and CLOSED.** Every loader-shaped function in `dna_decode/glm/` was enumerated: there is
  exactly **one loader per GEO series** — `expression.py` owns GSE108535 (σ70 designed) and
  `genomewide.py` owns GSE144621 (peak tiles + fragments), each with its own `GEO_BASE`. No other
  duplication exists. One cosmetic residue: `DEFAULT_CACHE = Path("D:/dna_decode_cache/mpra")` is
  defined identically in both modules. Same value, different series, no correctness consequence —
  recorded, not fixed, because consolidating it would couple two independent substrate loaders.
