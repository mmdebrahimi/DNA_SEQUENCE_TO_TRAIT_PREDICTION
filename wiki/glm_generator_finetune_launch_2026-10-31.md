# G-D's LoRA fine-tune is RUNNING — how to read it when it lands

<!-- doc-status: LIVE as of 2026-10-31. The run was IN FLIGHT when this was written and NO
     fine-tuned number existed. Nothing here reports a result. The bar is frozen in
     wiki/glm_generator_finetune_acceptance_bar.json; read that BEFORE the verdict artifact. -->

## What is running

`glm-generator-finetune` on a Kaggle T4, launched 2026-10-31 ~19:46Z, status RUNNING at the time of
writing. Authorised by the user that day: *"re-open the 3-day fine-tune now"*, which overrode the
HOLD they had recorded earlier the same day and the frozen plan's Step 6.

It LoRA fine-tunes `GENERator-v2-prokaryote-1.2b-base` on 3,737 MG1655 upstream-of-CDS windows and
asks the family's own question: does that bring distinguishability **below the 3-mer Markov null**?

## Read these, in this order

1. `wiki/glm_generator_finetune_acceptance_bar.json` — the FROZEN bar. All four outcomes were named
   before any number existed. **Outcome B carries an explicit prohibition on re-specifying the
   0.6679 threshold**, which is not a chosen number: it is what the dumb baseline scores.
2. The verdict artifact `glm_generator_finetune_verdict.json`, pulled into
   `wiki/fragment_null_batches/` by `scripts/fragment_null_overnight.py`.
3. **If `glm_finetune_failure.json` is there instead, the run DIED and that file is NOT a verdict.**
   It exists because Kaggle returns zero output for an ERROR kernel, so the only way a failure is
   readable at all is for the kernel to exit 0 and report it.

## The two things most likely to be misread

**THE COMPARATOR IS THE IN-RUN ZERO-SHOT ARM, NOT THE RECORDED 0.7706.** That figure was produced on
CPU; this is CUDA, and this project has already measured that the two yield different sampled
sequences for this generator. The 2026-10-07 artifact also stamps neither its shuffle seed nor
`limit_natural` — both are recoverable from the launcher (it passes no `--seed`, so argparse's
default 7 applies, and `--limit-natural 1200` is explicit), which is weaker than stamped. So the run
measures zero-shot and fine-tuned in the SAME session: same device, dtype, eval windows, seed and
CLI, with the weights as the only difference. Both arms go through the same falsifier CLI with only
`--model` differing, so the delta cannot be an artefact of the code path.

**CLEARING THE BAR WOULD BE NECESSARY, NOT SUFFICIENT.** The discriminator uses k-mer + positional
features, which measure LOCAL composition. Published audits name long-range organisation as the real
failure mode of DNA generators. A pass here would mean the local statistics improved, not that the
generated promoters work.

## What the run proves about itself before it measures anything

Verified on Kaggle by a CPU diagnostic that reached `step2-complete` (`glm-finetune-diag`):

| check | result |
|---|---|
| corpus rebuilt | 4,340 windows, fingerprint +2118/−2222, 2,708 overlapping — matches the baseline |
| split | 600 fit / 600 eval / 3,140 unused / **3,737 train** — matches the committed manifest |
| train∩eval overlap | **0**, asserted on sequence CONTENT, not index arithmetic |

Disjointness is the integrity crux: a generator reproduces what it memorised, so training on even
one eval window would look like a win. The train set is derived BY EXCLUSION from the eval half by
`scripts/glm_finetune_window_split.py`, and the kernel re-asserts it and REFUSES on any drift.

## Context that should temper expectations

The sampling sweep reported first, by design, so the two generator levers stay attributable. Its
answer: **no temperature setting beats the null** (kmer gap 0.1365 → 0.0401 → 0.0733 across
T=1.0/0.8/0.6, best at 0.8, none crossing). That **falsifies GD-H4** — sampling alone had not
already done the work — so the fine-tune is a genuine untested lever rather than a redundant one.

Pre-registered forecast: **P(bar met) = 0.2**. 3,737 windows is ~560 kb for a 1.2B model, and the
ledger itself rates "is that enough" at LOW confidence. **Do not defer to that number** — my
previous forecast in this project (~80% that the fragment null gate would still fail) was falsified
by measurement, and so was its successor.

**A miss does not separate "fine-tuning does not work" from "not enough data."** Whole-genome
continued pre-training is the ledger's candidate 2 and is already recorded as a Pending Decision — a
scope call, not an executor task.

## Cost correction

The ledger's "+3 days" was a **GTX-1070-at-fp32** estimate (consumer Pascal runs fp16 at ~1/64 of
fp32). On a T4 with a ~560 kb corpus this is a sub-session job. That revision is an estimate from
corpus size and device, not a measurement; the artifact stamps its own `wall_minutes`.

## Three defects of mine, and the guards that now exist

Three pushes were spent before this ran, all on my own errors, all avoidable:

1. **A `.py`-only payload** omitted the baseline artifact the split script needs to check the
   fingerprint. Its rc=3 refusal was CORRECT.
2. **The builder validated its own syntax and not the emitted kernel's.** The body is a raw
   triple-quoted string, so a collapsed escape is legal in the builder and an unterminated string
   literal in the artifact. It shipped twice. Now `ast.parse(body)` runs before the write.
3. **The split subprocess had no `PYTHONPATH`.** Step 1's import check could not catch it, because
   the PARENT imported everything fine.

And a wrong diagnosis worth recording: **I read "no log, no output files" as evidence of an
infrastructure failure.** It is not evidence of anything — Kaggle returns exactly that for ANY failed
kernel, including one with a plain syntax error. That misreading also means the earlier b3 batch
failure was never actually established as infrastructure; its retry simply succeeded, which is the
only fact there. All of this is pinned by `tests/test_glm_finetune_kernel_builder.py` (10 tests).
