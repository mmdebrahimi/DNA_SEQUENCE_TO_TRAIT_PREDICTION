# /innovate round 2 — the catalog is already a partial interventional world model (2026-07-13)

**What this is.** A *second, deeper* `/innovate` pass on the genome world model. Round 1 (`wiki/innovate_genome_world_model_2026-07-13.md`) reframed the **product** (federation of cells, observational→interventional lever, self-awareness flag). Round 2 hunts a concretely novel **mechanism**, and attacks the assumption round 1 left standing: that the world model is something to *build from scratch*. Engine: `falsification.py --run-ledger … --cwd .` (executes each kill-test cheap-first through the gated runner, one call). Fail-closed validator on the stamped ledger: **exit 0 — 5 candidates → 4 survived / 1 killed / 0 unfalsified.** Each kill-test is unique to its own claim (the 2026-07-13 shared-kill-test dogfound lesson). Kill-tests: `innovate_killtests/test_innovate_round2_2026_07_13.py` (out of the default `tests/` suite).

**Discrimination proof (round 2).** A fresh tempting-but-wrong control — *"gene-wide interpolation: ANY mutation ANYWHERE predicts effect"* — was **KILLED** (known-functional enrichment R-over-S = **4.203**: the signal is *localized* to functional sites, not gene-wide). This both proves the gate still discriminates and **bounds the g1r2 survivor**: generative interpolation works *at catalogued/functional sites*, not across the whole protein.

---

## The round-2 thesis (what the 4 survivors say together)

> **The curated catalog is not a lookup table waiting to be replaced by a world model — it is ALREADY a partial interventional (edit→effect) world model.** Each entry is a compressed *interventional* statement ("this substitution at this site → this effect"). The forward path is not "build a big learned model" — it is to make the existing catalog (a) *generative* at sites, (b) *transfer-driven* across organisms, (c) *interventionally-seeded* for new cells. This is a sharper, more actionable answer than round 1, and it is novel: it relocates the "world model" from a thing-to-build to a thing-that-partly-exists-and-should-be-extended.

### g8r2-cross-organism-transfer  *(G8 — the strongest survivor)*
**Move.** A target-site catalog is **protein-family-level, not organism-level**, so a validated catalog transfers to a new organism with **zero new labels** — the fastest way to add census cells, bypassing the label wall for mechanism-conserved targets.
**Kill-test (executed, STRONG).** DEAD iff fewer than 3 distinct genera share a SCORED cipro gyrA/parC-QRDR cell. Committed: **4 genera** (E. coli, Klebsiella, Campylobacter, Salmonella) all SCORED on cipro via the same QRDR → SURVIVED. This is not a precondition — it is *four validated cells demonstrating the transfer already works*. `[grounded]`

### g9r2-interventional-kb  *(G9 abduce — the deepest explanation)*
**Move.** Abduce *why* the catalog beats learned models: not because it is curated, but because each entry is a **background-INDEPENDENT interventional statement** that transfers across genetic backgrounds where an observational correlation cannot. Consequence: build the world model as a catalog-shaped interventional KB; seed every new cell from interventional experiments (DMS / marker-transfer), never observational isolate mining.
**Kill-test (executed, STRONG).** DEAD iff catalogued DRMs are background-*dependent* (observational-like) → catalog leave-study-out balanced accuracy < 0.70. Committed: **0.824** → the catalog generalizes across independent held-out studies (the interventional signature) → SURVIVED. `[grounded]`

### g1r2-generative-catalog  *(G1 reframe — serves the original north star)*
**Move.** Make the catalog **generative**: predict an *unseen* edit's effect by interpolating across a site's already-catalogued substitutions ("any non-WT substitution at a multiply-catalogued resistance site tends to confer resistance"). This directly serves *"guess the phenotype of a minor modification"* — the original north-star framing — and turns a lookup table into a predictor.
**Kill-test (executed).** DEAD iff there is no interpolation substrate → fewer than 2 DRM positions carry ≥2 catalogued substitutions. Committed: **4 positions** (106, 181, 188, 190) → SURVIVED. *Strengthened by prior committed evidence:* the position-novelty flag already validates the interpolation rule (`FLAG_RECOVERS_BLINDSPOT`, lift 3.98). **Scope (from the killed control):** works at catalogued/functional sites, *not* gene-wide. `[grounded]` `[inferred]`

### g5r2-self-supervised-catalog  *(G5 relax "needs external labels" — WEAKEST survivor, flagged)*
**Move.** Relax *"a learned world model needs external phenotype labels"* → train a generative effect-predictor **self-supervised on the catalog's own edit→effect entries** (thousands of interventional statements across organisms = a training set for "predict the effect of an edit"; no new wet-lab label).
**Kill-test (executed — PRECONDITION ONLY, honest caveat).** DEAD iff the pooled catalog is too small to be a training set (< 100 entries). Committed: **454+** (16 HIV NNRTI + 438 TB WHO grade-1/2, before the other cells) → SURVIVED. **⚠ This kill-test proves only that the training set EXISTS, not that self-supervised training WORKS** — g5r2 survived the weakest test in the set. It belongs one rung below the other three until a real train/held-out experiment scores it. `[grounded]` `[speculative]`

---

## Killed
- **control-genewide-interpolation** — gene-wide interpolation is KILLED (signal localized to functional sites, enrichment 4.203). Bounds g1r2 to catalogued/functional sites.

## Honest strength ranking of survivors
1. **g8r2-cross-organism-transfer** — 4 validated cells demonstrate it; immediately actionable (transfer a catalog to a census organism).
2. **g9r2-interventional-kb** — the deepest *explanation* (background-independence signature, bal_acc 0.824); reframes how every new cell should be seeded.
3. **g1r2-generative-catalog** — validated at sites (lift 3.98), scope-bounded by the killed control.
4. **g5r2-self-supervised-catalog** — precondition-only survival; needs a real experiment before it can be trusted.

## Handoff
The round-2 mechanism is: **treat the federated catalog as the partial interventional world model it already is, and extend it by transfer (g8r2) + interventional seeding (g9r2) + site-local generative interpolation (g1r2).** The most actionable next build is **g8r2 cross-organism catalog transfer** — take a validated target-site catalog (e.g. QRDR/cipro or the WHO-TB rpoB/katG) and transfer it to a census organism that shares the mechanism, with a leakage-clean check that the transferred catalog scores above chance with zero organism-specific labels. `/technical-plan`-ready. g5r2 is the one survivor that must be *scored*, not built on, first. Read-only over committed artifacts; frozen decoder surface byte-unchanged (verify_lock OK).
