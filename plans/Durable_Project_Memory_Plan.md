# Durable project memory — surviving compaction without lying to the next session

**Status:** F1/F2/F5 shipped 2026-08-29; F3 partly done 2026-08-30 (measurement + citation guard shipped;
the 17-bullet compression is a user call). F4b remaining.

## The problem, stated correctly

Not "context is lost on compaction." Context loss is unavoidable and already handled — `CLAUDE.md` +
`MEMORY.md` + `wiki/` reload every session. The real failure is that **the reloaded surface was wrong**,
and wrong in a way that reads as authoritative.

Three measured instances, all 2026-08-29, all in one session, all with the grounding artifact sitting on
disk the whole time:

| # | error | shape |
|---|---|---|
| 1 | quoted the AMR arm's 27 cells / 10 SCORED as the tool's evidence surface (actually 110 / 28) | **scope collapse** — a subset read as the whole |
| 2 | called organism-level g→p a "closed negative" (3rd occurrence) while a 12/12 r 0.46–0.80 positive sat in-repo | **over-compression** — scope dropped from a scoped negative |
| 3 | proposed the FBA/Keio line as unexplored — it is the deepest line here | **staleness** — memory vs. record |

The proximate cause for #1 was mechanical: **CLAUDE.md's opening paragraph** — the only project text
auto-loaded into every session — described a "Phase 1 E. coli platform predicting cipro/cef/tet." Three
months and six tracks out of date. A session's working model is downstream of that paragraph.

## The design principle

**Distil pointers to derivations, not facts.**

A written number is true when written and silently false later. A *derived* number cannot drift, because
it is computed from the thing it describes. So the orientation surface states as few figures as possible,
and every figure it does state is pinned by a test to its live source.

Three layers, each with exactly one truth-owner:

| layer | file | owner | drift risk |
|---|---|---|---|
| **L0 orientation** | `CLAUDE.md` §READ THIS FIRST (<45 lines) | points at L2 | pinned by test |
| **L1 detail** | `CLAUDE.md` body, `wiki/` | prose | high — never authoritative |
| **L2 ground truth** | `cell_registry`, `pyproject`, `cli.TRAITS`, artifacts | code | none |

`scripts/project_status.py` renders L2 on demand. It is read-only, offline, ~2s, exit 0 always — a report,
never a gate.

**It validated itself on first run** by correcting two figures written an hour earlier in
`wiki/project_distillation_2026-08-29.md`: 46 CLI traits → **44**, and "~3x understated" → **4.1x**.

## What shipped

- **F2 · `scripts/project_status.py`** — derives version, entry points, traits, the full 110-cell tier
  distribution, per-track counts, and the AMR card **labelled with its scope**, plus the regime map.
- **F1 · `CLAUDE.md` §READ THIS FIRST** — replaces the stale Phase-1 opening. Names all three errors
  above as worked examples, because a rule with a case attached survives compression better than a rule.
- **F5 · `NEXT.md`** — transient open threads + user-authority calls. Explicitly prune-not-grow.
- **F4 (partial) · `tests/test_project_orientation.py`** — 6 tests. Pins the quoted cell/SCORED counts to
  the live registry (drift → loud failure naming the line), pins the population-design correction, pins
  the exact stale sentence out of the file, and caps the block's length so orientation can't become a
  second document.

## Remaining

**F3 — PARTLY DONE 2026-08-30.** `scripts/claude_md_weight.py` + `tests/test_claude_md_citations.py`
measure the always-loaded surface (~36,800 tokens/session) and guard that every citation resolves — which
immediately caught a plan that moved to `executed_plans/` while its pointer stayed at `plans/`. The
session's own 1,022-word bullet was compressed to 309. Remaining: 17 candidate bullets (~9,100 words),
deliberately left as a user call (see `NEXT.md`). Memo `wiki/claude_md_weight_2026-08-30.md`.

**F3 (original framing) — scope tags on L1 claims.** Every count in the CLAUDE.md body should name the arm it covers. The
report-card bullet already does this (corrected 2026-08-23). Cheapest form: a convention, checked by the
staleness auditor rather than by a new test.

**F4b — point the staleness auditor at L0.** The auditor already runs a 110-item corpus at ~1-in-6 true
positives. The orientation block is ~30 lines and the highest-consequence text in the repo — running it
against L0 on every corpus pass is nearly free and targets the text that caused three errors.

## What was deliberately NOT built

- **A bigger distillation document.** `wiki/project_distillation_2026-08-29.md` already exists and is
  itself now partly stale — which is the argument. Prose scales the problem, not the fix.
- **Auto-regenerating CLAUDE.md.** Tempting and wrong: the body's value is hard-won judgment
  (gotchas, closed negatives, traps) that no generator can derive. Only the *figures* are derivable.
- **Freezing the orientation block.** It should change. The test makes change loud, not hard.

## Reproduce

```bash
uv run python scripts/project_status.py          # the live surface
uv run pytest tests/test_project_orientation.py  # the drift guards
```
