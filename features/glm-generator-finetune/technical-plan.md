# Technical plan — DEFERRED BY DESIGN (G-D)

**There is no technical plan here yet, and that is deliberate rather than an oversight.**

Two blockers, one physical and one epistemic. The 8 GB GPU is in a machine being physically relocated. And the queued sampling sweep must report FIRST: sampling and weights are two generator levers, and running them together means neither can be attributed. If the sweep already clears 0.6679, this family's premise weakens materially before any training starts.

**One slice IS startable now and is listed in the ledger as short-term action 2** — building the disjoint train/eval window split, so fine-tuning can never train on the discriminator's own natural comparison set. That is CPU work and needs no plan beyond the ledger row.

**What a plan would cost right now:** a detailed plan written against a predecessor's unknown outcome goes
stale the moment that outcome lands, and this repo has a recorded lesson about exactly that
(`feedback_saved_plan_steps_vs_amendments_divergence` — a saved plan whose premises moved became a
misleading artifact that still looked authoritative). Writing five plans when two families can start is
planning theatre; it produces paper, not progress.

**What exists instead, and is enough to start the moment the block clears:**
- the ledger's `## Goal Hierarchy → Short-term` table (concrete, file-level actions)
- the ledger's `## MVP Criteria` (checkable predicates, so `--until-mvp` can run it)
- the ledger's `## Decisions Made` (the constraints a future plan must honour)
- `feature.md`'s **Explicitly NOT in scope** section (each item a mistake already made once here)

**Unblock trigger:** the sampling sweep artifact exists AND the 8 GB host is reachable.

Write the plan then, from the predecessor's actual result.
