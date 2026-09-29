# The orientation doc is the thing eating the window it exists to protect

**2026-09-29.** Measured diagnosis of the "we keep compacting and re-doing work" complaint, plus the
increment shipped against it. Two of my own claims were corrected mid-run and both corrections are here.

## The premise was half wrong, and the wrong half matters

The hypothesis was *"we need a long-term lessons ledger, which we were doing — not sure what happened
to that."* **Nothing happened to it.** `LESSONS_LEARNED.md` carries **96 entries in September**, was
written to three times that week, and stands at 243k chars. It is the healthiest store in the repo.

There are **seven** stores, not zero: git history · `LESSONS_LEARNED.md` · per-ledger `## Action Log` ·
`wiki/` memos · `CLAUDE.md` · `NEXT.md` · cross-project `memory/`. A single session routinely writes to
five of them.

**Recording is not the gap. RETRIEVAL is.** The same week, real effort went into re-deriving the
pneumoserotype no-call denominator that was stated *verbatim* in a sibling artifact nobody had cited.
Nothing was lost; it could not be found.

## What is actually consuming the window

| auto-loaded on EVERY session AND every compaction | tokens |
|---|---|
| `dna_decode/CLAUDE.md` | **~69,300** |
| `PythonProjects/CLAUDE.md` (parent) | ~5,600 |
| cross-project `MEMORY.md` | ~4,900 |
| global `~/.claude/CLAUDE.md` | ~600 |
| **before a single word of work** | **~80,400 (~40% of a 200k window)** |

`CLAUDE.md` grew **41k → 116k → 277k chars since June — 6.7×**. 53 of its 101 bullets are May–July arcs;
one single bullet is 18,500 characters; it carries ~45 `CORRECTED`/`SUPERSEDED`/`stale` markers — a
historical record of things that turned out wrong, re-read in full every time.

**So the loop is self-defeating: the document written to PREVENT rework is the main consumer of the
window whose exhaustion CAUSES it.**

The fix is already invented *in this repo, at the top of that very file* — *"every figure below is written
down and therefore goes stale; the script derives it live; when they disagree the script is right."* That
pattern was applied to **scope** (`project_status.py`) and to **data** (`data_inventory.py`). It was never
applied to the orientation document itself.

## The diet is available and measured — by a tool that already existed

`scripts/claude_md_weight.py` (already in-repo, needs `PYTHONIOENCODING=utf-8` on this host) uses a
**measured** bar, not a judgment one: a bullet is compressible only if its cited file **really exists**.

- **39 bullets / 26,718 words — 72% of the file — are PROVABLY stored elsewhere in `wiki/`.**
- **5 bullets / 2,035 words cite nothing: `CLAUDE.md` is the ONLY copy. They must stay whole.**
- Every cited path resolves today.

Realistic target **68,665 → ~22,000 tokens (saves ~47k per session)**. **My first estimate of ~10k was
optimistic by about 2×** — I proposed a judgment bar ("narrative leaves, guardrails stay") when a stricter
measured one already existed, then guessed the magnitude instead of running the tool. Both corrected by
measurement.

**The compression format is load-bearing and is NOT a bare pointer.** It must be
`headline + the guardrail verdict + pointer`, because the highest-value text in those bullets is the
`do NOT reopen` rail. A bare pointer would preserve the bytes and lose the thing that stops a future
session re-proposing a closed negative.

## The replay claim I nearly published, and the correction

The candidate-row scan looked damning: **30 rows across 7 ledgers, every one flagged replay-risk; total
`[done:]` markers repo-wide = 1.** Several rows contained hand-written 300-word defensive essays
(*"DELIVERED … do NOT run this"*, *"CLOSED BY ENUMERATION"*, *"Do NOT keep grinding"*) — the same disease
as `CLAUDE.md`, one layer down: the record grew to defend against replay instead of the row being retired.

**But "30 rows cause rework" is wrong.** The umbrella's `## Project Families` table lists only **5**
families, and `advance_ranker` keeps only families listed `ACCEPTED` there. Three ledgers
(`…-scratch`, `ecoli-pathotype-prediction-cli`, `eukaryotic-trait-decoding-cycle`) **are not in the table
at all**, so the ranker never reads them — their 15 stale rows cannot cause a code-driven replay. They are
a trap for a *reader*, which is a real but different and smaller problem.

Code-reachable stale surface: **5 rows.**

## What shipped

1. **First resume checkpoint ever written in this repo.** `resume_state.write_checkpoint` had **0** prior
   uses here despite being the mechanism designed for exactly this failure. Lives at
   `~/.claude/soraya-resume/` — outside the repo, so it survives any repo state.
2. **5 `[done:<id>]` markers journaled**, on **self-declared** terminality only (the row's own text says
   DELIVERED / CLOSED / PLATEAU / RESOLVED) — never on the advisory lexical detector, because retiring a
   live row loses real work and that is the dangerous direction. Verified: `RETIRE 0 / KEEP-LIVE 28 /
   ALREADY-JOURNALED 5`.
3. **Three terminal ledgers bannered** with `<!-- ledger-status: … -->` + why-terminal + superseded-by,
   placed after the schema marker so `/project-state` still validates (verified: marker at line 2, all six
   end-markers unique in each).
4. **`tests/test_claude_md_growth_ceiling.py` (5 tests).** A **growth CEILING, not a diet** — it does not
   claim the file shrank; it pins today's 36,861 words at a 37,600 ceiling with deliberately thin (~2%)
   headroom, so the next addition is a recorded decision instead of silent drift. Also pins that every
   cited path resolves, that the ceiling is not vacuously slack, and that the weight tool can still
   identify the single-copy bullets a future diet must not touch. **Proven non-vacuous:** adding 1,200
   words of filler fails it; `CLAUDE.md` byte-restored after.

## My own defect this run — third instance of one class in one week

My throwaway classifier sliced the Action Log on the **bare substring** `project-state:end:action-log`
instead of the full `<!-- … -->` delimiter. That substring also appears **inside Action Log rows written
earlier that week describing the marker**, so the slice truncated early and reported a clean
`journaled_ids = []` — i.e. it said the markers I had just written did not exist. Caught only because two
of my scripts disagreed about the same file.

Same family as the `2.1`-inside-`72.1` substring bug and the `0.9`-verifies-`0.939` truncation: **an
under-specified match returns a confident empty/clean answer rather than an error.**

## Still pending — named, not started

- **The diet itself** (39 bullets → headline+guardrail+pointer, `wiki/arcs/`). The dominant lever at
  ~47k tokens/session. Deliberately not half-done: a partial extraction is worse than none.
- **`scripts/where_is.py`** — the retrieval router (which of the seven stores answers a given question
  shape), derived rather than hand-listed; the hand-listed-exclusion trap has bitten this repo 5×.
- **The parent `PythonProjects/CLAUDE.md` (~5,600 tok)** is auto-loaded into *every* project on this
  machine, not just this one. Untouched here; same growth question applies.

## Honest limits

- The ceiling **stops the trend; it recovers nothing.** ~80k tokens still load per session today.
- The 28 `KEEP-LIVE` rows are kept on a *conservative* rule (no self-declared terminal ⇒ assume open).
  The lexical detector thinks ~20 of them are done. Verifying each against artifacts is real work and was
  not done — so the row tables remain partly stale by choice, in the safe direction.
- Nothing here was A/B tested against rework rate; the causal chain (auto-load size → compaction
  frequency → rework) is **argued and measured at the input end only**.
