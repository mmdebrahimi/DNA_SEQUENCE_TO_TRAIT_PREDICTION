# Soraya minister — first complete real bounded run (objective observations)

**Mission:** `ep4-v01-expec-recall` · **Run:** `2026-06-03-1309-ep4-v01-expec-recall` · **Verdict:** `mission-mvp-reached`
**Soraya:** v0.3.x (canonical entry `run_minister.py`, lease auto-arm) · **Repo:** dna_decode (git, main)

This closes Soraya's "#3 — one complete real bounded minister run" item (prior 4 dogfoods parked at
blocked:user-only or were disposable).

## Per-round ledger

| Round | Φ (in-proc) | Verdict | Emitted | Promote vs Park | --until-mvp outcome |
|---|---|---|---|---|---|
| 1 (session A) | 2 → 1 | `blocked:user-only` | `/interrogate-me` ×2 | PARK (0/2 receipts) | n/a (no frontier yet) |
| 2 (session A, resume) | 2 (fresh epoch) | `mission-mvp-reached` | — | preempted by terminal-check | n/a (endpoints already MET) |

## What worked
- **Preflight clean + accurate:** money hook hardened ✓, Agent ✓, python/pytest resolve ✓, cwd=repo-root ✓,
  bounded-mission ✓, no-prior-decomposition ✓. No false blocks.
- **Code-owned sequencing held:** `run_round` evaluated mission terminal → executed/ascended → resolved via
  `resolve_stop_portfolio` deterministically. The model only supplied the 3 seams (generate / until-mvp /
  project-init) — no hand-wired driver drift (OT1 fix working as intended).
- **Gate + park + resume:** round 1 correctly parked on the 2-interrogation gate (user-only skill emitted,
  not self-invoked); receipts recorded between turns cleared it; resume reached terminal. Journal idempotent.
- **Money gate dormant correctly:** pure-Python CPU-only mission; lease armed for the resume's `run_to_stop`
  and disarmed in finally; no money/destructive classification fired.

## Frictions / observations (for the dogfood corpus)
1. **`mission_met()` fast-path preempts formal promotion.** Because the model-driven engineering was done
   in-session BEFORE the final resume, the round's first check returned terminal and the family ledger
   (`project_state/fam-per-gene-expec-scoring.md`) was never seeded via `/project-init`. Not a bug — endpoints
   are the terminal, not the promotion artifacts — but the audit trail can be misread as "skipped promotion".
   The handoff doc's expected flow (promote → project-init → until-mvp does the work) assumes the work happens
   INSIDE `until_mvp_runner`; in practice `until_mvp_runner` is a subprocess endpoint-checker (can't host model
   reasoning), so the work lands in the session and the terminal is reached early. **Candidate fix:** either (a)
   document this as the sanctioned fast-path, or (b) have the wrapper seed the family ledger at GENERATE time so
   the promotion artifact exists even when the terminal is reached early.
2. **Φ on resume reads 2 (fresh in-memory epoch), not 1.** The durable journal has the round-1 debit; the
   resume process rebuilds the epoch from disk gaps (frozen) but recomputes potential from scratch. Finiteness
   is still sound (journal is source of truth), but a reader of the resume stdout sees `Phi=2` after a debit
   already happened — mildly confusing. Worth surfacing journal-derived Φ in the resume report.
3. **`sklearn` import at module top of `scripts/pathotype_kmer_bakeoff.py` blocks importing `load_strains`.**
   Had to replicate the offline genome loader. Pre-existing repo issue, not Soraya's — but it bit the run.
4. **Interrogation-vs-data ordering:** the load-bearing fork (which 2 of 3 strains to rescue, and at what K)
   could only be decided AFTER computing per-gene counts, which came AFTER the interrogation. The minister
   has no built-in "re-gate when implementation reveals a blocking-class cost" step; I surfaced it manually
   via an extra AskUserQuestion. Candidate: a lightweight "implementation-discovery re-gate" affordance.

## Net
First real bounded minister run reached `mission-mvp-reached` with a genuine engineering deliverable
(ExPEC recall 0.75 → 0.917, precision invariant held, 25/25 resolver tests, 0 reversed). The
finiteness/gate/park/resume machinery behaved as specified. The two friction points (terminal fast-path
preempting promotion; resume-Φ display) are reporting/observability nits, not correctness failures.
