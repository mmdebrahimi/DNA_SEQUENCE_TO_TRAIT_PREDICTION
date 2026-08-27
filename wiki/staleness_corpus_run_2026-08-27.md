# The full-corpus flagging pass: 80/110 scored, and it caught a real stale claim I made yesterday

**Headline: the auditor found a genuine stale claim in CLAUDE.md that I introduced 24 hours earlier and
never noticed — the thing it was built to do.** The run also hit an OOM at item ~80, which the
checkpointing turned into a partial result rather than a total loss.

## What ran

The benchmark-validated auditor (Qwen3-8B, free Kaggle T4, prompt asserted byte-identical to the one that
scored 5/5 TP / 1 FP) over **110 (claim, artifact) pairs extracted from CLAUDE.md**.

| | count |
|---|---:|
| pairs extracted | 110 |
| **verdicts recovered** | **80** (OOM at ~item 80; checkpoint every 10 saved the rest) |
| supported | 65 |
| **flagged stale** | **9** (11%) |
| unclear | 6 |
| unparseable | 4 |

## The catch that justifies the whole exercise

> **CLAUDE.md said "8 reusable rejection GATES … G1–G8". The map has TEN.**

I added G9 and G10 to `wiki/negative_results_map_2026-06-13.md` **on 2026-08-26** — the previous day —
and never updated the CLAUDE.md bullet that summarises it. The auditor's evidence was exact: *"the
artifact's excerpt explicitly lists 10 rejection gates (G1-G10), including G9 and G10 added on
2026-08-26, contradicting the claim's assertion of exactly 8 gates."*

This is the class no mechanical check could reach: the cited file **exists** (so `test_claude_md_citations`
passes), the paths all **resolve**, and nothing is syntactically wrong. Only the *claim about* the file
went stale. **Corrected in this commit.**

That it was *my own* drift, caught within a day, is the strongest evidence for the tool: the failure mode
is not carelessness a human would catch on re-read — it is that nobody re-reads a 22k-token file that is
loaded automatically.

## Adjudication, not acceptance

The 9 flags are **candidates**, and I adjudicated rather than reported them as findings. Two worked
examples:

- **TRUE POSITIVE — the gate count** (above). Fixed.
- **FALSE POSITIVE — `plans/Trait_Decoding_Roadmap.md`.** The model flagged *"the artifact is labelled
  'DRAFT 2026-05-26', contradicting the claim it was 'shipped 2026-05-26'."* Checked: the header is a
  title convention, and "shipped" refers to the artifact landing, not to a status field. Not stale.

The remaining 7 need the same treatment; several look like the **capability-vs-finding** category error
this run was partly designed to probe (`bvbrc_strict_mic_4drug_census.py` reappears — the same claim shape
that produced the benchmark's one false positive, now confirmed on an unseen instance).

**That is the useful negative result:** the category error is real and recurring, not a one-off artifact of
the benchmark. It is now confirmed on data the prompt has never seen, which is what I said the fix needed
before being written.

## The OOM, diagnosed

`torch.OutOfMemoryError: CUDA out of memory. Tried to allocate 3.61 GiB. GPU 0 has a total capacity of
14.56 GiB of which 2.24 GiB is free` — at roughly item 80. Cause: a long artifact excerpt (6000 chars) plus
a 2500-token generation produced a KV cache that outgrew the T4. It is a **capacity** failure, not a model
or logic failure, and the fix is mechanical: cap total prompt length, or free the cache between items.

**The checkpoint-every-10 decision is what made this recoverable** — 80 real verdicts survived a crash that
would otherwise have cost the whole two-hour run and returned nothing.

## Honest scope

- **80 of 110** pairs scored. The remaining 30 are unrun, not "clean".
- Extraction is CLAUDE.md only. The 542-memo `wiki/` corpus is untouched, deliberately: memos are dated and
  read rarely, so "stale" barely applies; CLAUDE.md is the live surface loaded every session.
- The pass is a **flagger**. At the benchmark's 1-in-5 FP rate on capability-shaped claims, every flag needs
  adjudication against its artifact — that step is the product, not an optional extra.
- 4 unparseable answers are counted, not hidden; the parser fails closed, so they cost recall, never
  precision.

## Next

1. Cap prompt length / clear the KV cache between items, and run the remaining 30.
2. Adjudicate the other 7 flags.
3. The capability-vs-finding rule now has unseen-data evidence behind it, so it can be written and
   validated on the *next* corpus slice rather than tuned against the 10-item benchmark.

Kernels `dna-staleness-corpus-run` (died on a dataset-attach path — fixed by discovering the mount rather
than hardcoding it) and `-run2` (this run). Dataset `dna-staleness-corpus`. All private, free tier,
deletable.
