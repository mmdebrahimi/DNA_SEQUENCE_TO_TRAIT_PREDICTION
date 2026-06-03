# Dogfood — Soraya `--advance` on EP-4 v0 marker resolver (2026-06-02, run 2)

Run: `soraya_runs/2026-06-02-1530-ep4-v0-marker-resolver/`. Second `--advance` of the session, executed back-to-back after run 1's confound finding redirected the work. Goal: build the confound-immune v0 marker resolver.

## What worked (cleaner than run 1)
- **De-risk-dependency-FIRST paid off.** Before writing any screen code, Soraya test-fetched the VirulenceFinder DB (step 0) and confirmed all needed markers were present. Had it failed, the whole batch would have re-planned cheaply. This is the right `--advance` instinct: probe the one external dependency before building on it.
- **VOI chain across runs held.** Run 1 (confound discovery) → run 2 (build the confound-immune alternative). The persona didn't get stuck re-litigating the failed learned track; it pivoted to the locked v0 design and executed.
- **Gate behaved; no surprises.** DB fetch + script write/run all `auto`. The deliberate choice to do pure-Python k-mer seeding (NOT install BLAST/KMA) kept the run dep-install-free — a good fit for the money-only/laptop-local posture.
- **No analyst override needed.** Unlike run 1 (where the auto-verdict was wrong and the model had to override on feature identity), run 2's auto-verdict `BIOLOGY_TRACKS_LABELS` matched the substance (eae 12/0, textbook ExPEC arsenal). Healthy: the override muscle from run 1 was available but correctly not invoked.

## Minor friction
- 🟡 Same as run 1: audit-trail appended once at the end (model-discipline), not per-step. Fine at 4 steps.
- 🟡 The run produced an interpretable, high-value deliverable (eae→EPEC, AUROC 1.0) in ~one batch — arguably it could have continued straight into wrapping the CLI (action 1 of the recommendation). Soraya stopped at the batch boundary + emitted the next step instead. For a maiden-flight dogfood that's the right conservative call; for a trusted run, `--until-mvp` against a "v0 CLI emits a pathotype call" bar would have carried it further autonomously.

## Net across both runs
Soraya `--advance` did genuinely valuable, correct science across two chained runs: caught a confound that invalidated the learned track (run 1), then built + validated the confound-immune v0 resolver that the project had locked as its design (run 2). Safety-critical pieces (gate, lock, de-risk-first) behaved. The one durable lesson for the persona: an auto-metric verdict needs a model-level substance check (run 1's CTAG catch) — keep that override step explicit in the loop.
