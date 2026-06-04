# Minister launch handoff — EP-4 v0.1 VF side-by-side diff (2026-06-03)

> The next DNA move, chosen by two converging Soraya analyses + user greenlight. Completes the v0
> resolver spec (its last originating-goal item). IN_PLACE mission (improve the incumbent v0 CLI, no
> family generation). Run from a dna_decode-rooted session — the minister cwd-preflight (v0.3.3) HARD-
> refuses otherwise.

## Why this move (vs the alternatives)
| move | why not now |
|---|---|
| **VF side-by-side diff** ← CHOSEN | finishes the v0 spec; in_place; highest completion value; only friction is BLAST install (an environment step for the live session) |
| ETEC 3-class expansion | opens a new front + needs full-GCA plasmid-toxin fetch (the chromosome-only footgun); lower completion value than closing v0 |
| break study==class (de-confounded labels) | the real scientific frontier, but it's the NEXT epoch — do it after v0 is cleanly closed, not instead of finishing v0 |

## Launch (from a dna_decode session ONLY)
```
cd C:\Users\Farshad\PythonProjects\dna_decode
claude                       # fresh; do NOT resume a wrong-cwd session
/soraya minister ep4-v01-vf-diff
```
The minister loads `project_state/ep4-v01-vf-diff/big-idea.md` (mission-shape: **in_place**), binds the
v0 resolver CLI as the sole frontier member, PARKS round 1 on the interrogation gate (in_place gates
EXECUTION, not promotion), then on the receipt runs `--until-mvp` to ship the diff.

## What the run will do (expect)
1. preflight (money hook hardened ✓, Agent ✓, **cwd==dna_decode** ✓, big_idea bounded ✓).
2. in_place gate: emit `/interrogate-me` for the 2 design decisions below — **with Soraya's drafted
   answers** (draft-then-ratify); you RATIFY.
3. `--until-mvp`: install BLAST+ (or wire the Docker VF path) → add a `vf_runner` (canonical VF over the
   FASTA) → add a `vf_diff` builder (resolver cluster profile vs VF per-gene calls + concordance +
   honesty flag) → wire into `cli.py` provenance JSON → write `tests/test_pathotype_vf_diff.py` →
   green both endpoints → AC9 ledger row "VF side-by-side diff shipped".

## Drafted interrogation answers (ratify or redirect)
**Q1 — canonical-VF execution path: native BLAST+ install, or Docker VF image?**
Soraya draft: **native BLAST+ install** (`blastn`), invoking `python -m virulencefinder` against it —
matches the Gate-A-proven runtime on the workhorse (CLAUDE.md / gate-A pass used `python -m
virulencefinder` + cloned VF DB + local blastn), keeps the laptop path identical, no Docker-image pin to
maintain. Risk of the alternative (Docker): another pinned image + mount plumbing for marginal benefit
when BLAST+ is a 5-min install. Authority note: this is a tooling choice, yours only if you have a
standing Docker-only policy.

**Q2 — what "agreement" means in the diff (the honesty bar):**
Soraya draft: report **per-gene presence concordance** (both-call-present / both-absent / disagree) at the
resolver's confident-coverage bar (0.80), PLUS the explicit `caller_is_independent_baseline: false` flag,
AND a one-line caveat in the section that both callers use the same VF DB so high concordance is expected
and is an audit of the fast caller, NOT independent validation. Do NOT report a single headline
"agreement %" without that caveat — it would read as independent corroboration, which it is not. This Q2
secretly carries an AUTHORITY decision (how honestly the tool represents its own non-independence) —
surfaced explicitly per draft-then-ratify.

## After this mission
v0 is COMPLETE → freeze it (ledger row already sequences this) → the frontier epoch = de-confounded
(within-study / lineage-matched) ExPEC+EPEC labels, which is where real biology-vs-batch separation
becomes measurable. That epoch is research/substrate, not modeling.

## Provenance
Chosen per Soraya analysis (this session) + the parallel Skill_Development Soraya session, both
converging; user greenlit 2026-06-03. Ledger action row 20 records the halt-recall + complete-v0-first
sequence.
