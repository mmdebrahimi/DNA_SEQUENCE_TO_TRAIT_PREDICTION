# CLAUDE.md costs ~37k tokens every session, and half of it is stored somewhere else

Two runs ago I diagnosed CLAUDE.md's stale opening as the cause of a ~4x scope error and built a
derivation layer so figures could not go stale. Then I appended **1,022 words to a single bullet over
three runs** — making it the second-longest in the file, in a file whose size I had just argued was the
problem. This measures the cost and fixes my own share of it.

## The measurement

`scripts/claude_md_weight.py` (read-only, offline, never edits):

| | |
|---|---:|
| words | 18,289 |
| bullets | 78 |
| **approx tokens, loaded EVERY session** | **~36,800** |
| bullets ≥250 words **with** a resolvable external store | **18 = 9,437 words** |
| bullets ≥250 words with **no** external store | 2 = 1,065 words |

## The distinction that makes this safe

**Shorter is not better.** The body's value is hard-won gotchas that nothing else records, and deleting
those to save context would be strictly destructive. The answerable question is narrower: **does this
bullet's detail exist anywhere else?**

- A bullet citing a `wiki/` memo that is really on disk can keep the rule and the headline and point at
  the memo for the derivation.
- A bullet citing nothing **is** the only store and must stay whole, however long.

The tool refuses to call anything compressible unless a cited file actually resolves. Two long bullets
are protected by that rule right now.

## What I changed, and what I did not

**Fixed one genuinely broken citation.** `plans/Genome_Map_Virulence_Overlay_Plan/` completed and moved to
`executed_plans/`; the pointer didn't follow. A bullet promising a memo that isn't on disk is worse than
no pointer — it reads as authoritative provenance and sends a future session hunting for nothing.

**Compressed my own bullet: 1,022 → 309 words**, reclaiming 713. Nothing was lost — the derivation, the
defect narrative and the caveats live in `wiki/fba_within_gene_ranking_2026-08-29.md`,
`LESSONS_LEARNED.md` and the auto-memory. What stayed is what a future session must act on: the four
rules, the headline numbers, the "do not build the relative rule" warning, and the pointer.

**I did not touch the other 17 candidates (~9,100 words).** That is other sessions' institutional memory,
and rewriting it wholesale is a judgment call about what the project keeps in its always-loaded surface —
not a mechanical fix. The measurement is here so that call can be made with numbers.

## The checker had to be debugged before its output could be believed

Its first run reported **six broken citations. Five were its own bugs:**

| bug | effect |
|---|---|
| expanded only the FIRST brace group | `..._2026-07-1{6,7}.{md,json}` kept a literal `.{md,json}` → 3 false positives |
| treated `<date>` placeholders as literal filenames | 2 template citations reported broken |
| flattened citation groups | a satisfied `{a,b}` citation reported its non-existent siblings |

A checker whose false positives look exactly like the defect it hunts is worse than no checker. All three
are pinned by tests, along with the guard itself (`tests/test_claude_md_citations.py`, 8 tests) and a
non-vacuity check so a broken regex can't make the guard pass by finding nothing.

Reproduce: `uv run python scripts/claude_md_weight.py`
