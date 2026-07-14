# /innovate — a way forward toward the Genome world model (2026-07-13)

**What this is.** A `generate → falsify → survive` pass (the `/innovate` skill) on the project north star. 7 candidates generated under 4 forced operators (G1 reframe / G9 abduce / G8 recombine / G5 relax-constraint), each with the **cheapest executable kill-test that would DISPROVE it**, run through the gated falsification runner (`~/.claude/skills/soraya/scripts/falsification.py`) against **committed repo data**. A candidate is `survived` ONLY if an *executed, valid* kill-test failed to kill it. Fail-closed validator: **exit 0** — after a post-run hardening pass (see the g5 section), **8 candidates → 4 survived / 2 killed / 2 unfalsified** (the naive additive-rule form of g5 was killed by a deeper committed-data kill-test; the self-awareness-flag refinement survived). Machine-readable ledger: `wiki/innovate_genome_world_model_2026-07-13.json`; kill-tests: `innovate_killtests/test_innovate_killtests_2026_07_13.py` (kept OUT of the default `tests/` suite — each test *passes iff the idea is falsified*, so they are red by design and must not pollute CI).

**The pipeline discriminates (anti-theater proof).** A deliberately-planted tempting-but-wrong control — *"scale the embedding to 3B/15B to finally learn the mechanism"* — was **KILLED** by an executed test (ProteinGym median Spearman 650M=0.484 > 3B=0.467 > 15B=0.438; bigger regresses). The gate did not wave a plausible-but-false idea through. This reconfirms: **do not scale the monolith.**

---

## Survivors (executed kill-test did not kill) — and they compose into ONE reframe

The four survivors are not four separate bets; they cohere into a single strategic reframe of what the "genome world model" *is*:

### g9-interventional  *(G9 abduce — the deepest lever)*
**Move.** Build and validate the world model on **INTERVENTIONAL (edit→measure) data** — DMS + marker-transfer fold-change — not observational isolates. **Abduction:** the reason every learned whole-genome approach failed the *same way* (learned phylogeny, not mechanism, 0-for-5) is not fundamentally a label-scarcity problem — it is that **observational isolate data is phylogeny-confounded by construction**. The curated catalog wins precisely because each entry encodes *interventional* (site-directed / marker-transfer) knowledge. In an edit→measure design the data-generating process is not confounded by lineage.
**Kill-test (executed).** DEAD iff a learned scorer fails even on clean interventional data → median ProteinGym ESM2-650M Spearman < 0.10. Committed: **0.484** → SURVIVED. `[grounded]` interventional DMS works; `[grounded]` observational embeddings 0-for-5; `[inferred]` the observational↔interventional split is the deepest explanation.

### g1-federation  *(G1 reframe — what the product IS)*
**Move.** The genome world model is a **FEDERATION of de-confounded deterministic cells + calibrated abstention + the genome-map honesty layer** — a composed, auditable system — **not a single monolithic learned model**. The north star is reached by *composition and orchestration*, not by one big model that "understands the genome".
**Kill-test (executed).** DEAD iff too few independent SCORED cells exist to compose (< 3). Committed: **10 provenance-disjoint SCORED cells** → SURVIVED. (Necessary-condition test — surviving proves the substrate exists, not that composition adds value; that is the build's job.) `[grounded]` 10 cells + shipped genome-map + abstention; `[inferred]` composition-with-abstention is a different product than a monolith.

### g5-selfaware-flag  *(G5 relax "must beat the catalog" — refined after a post-run hardening)*
**Move.** Relax the assumed goal *"the learned model must BEAT the catalog"* → the learned arm's job is to **EXTEND** the catalog. **But the FORM matters, and a deeper kill-test split it in two:**
- **g5-naive-additive — KILLED.** The naive form ("add an accessory-mutation rule and it net-improves") is falsified on committed data: catalog+accessory `improves=false` (Δbal_acc **−0.006**; sensitivity 0.864→0.908 but specificity 0.783→0.727). Adding blind-spot mutations as a hard rule *trades away specificity* and does not net-improve.
- **g5-selfaware-flag — SURVIVED.** The value is a **self-awareness position-novelty FLAG** ("the catalog call may be incomplete here"), not an additive resistance rule. Committed: `FLAG_RECOVERS_BLINDSPOT`, median lift **3.98**, sens-on-blind-spot 0.604 — it raises TRUST without the specificity cost that killed the additive form.
**Kill-tests (executed).** naive: DEAD iff catalog+accessory does not net-improve → committed `improves=false` → **KILLED**. selfaware: DEAD iff the flag does not recover the blind spot → committed `FLAG_RECOVERS_BLINDSPOT` lift 3.98 → **SURVIVED**. `[grounded]` both numbers from committed HIV artifacts.

> **Post-run hardening (2026-07-13, `--advance`) — /innovate applied to itself.** The original g5 kill-test was *too shallow*: it checked only "do blind spots exist" (53 isolates) and let the naive additive-rule form survive. An R4 non-duplication scan found the extension was **already built for HIV** (`hiv_catalog_accessory_extension.py`) and its committed result is a NEGATIVE for the additive form. Adding the deeper net-improvement kill-test **KILLED** the naive form and left the honest self-awareness-flag refinement as the survivor. This is the anti-theater discipline working reflexively — a shallow "survived" that deeper committed evidence would complicate was caught and corrected before it could mislead a build decision.

### g8-residual-detector  *(G8 recombine: DEAD-END-A × ASSET-5)*
**Move.** Cross *"embeddings learn phylogeny not mechanism"* with *"we built Mash-clade de-confounding machinery"* → turn the **de-confounding machinery itself into the product**: a phylogeny-removed **residual-signal detector** that reports which genome regions carry mechanism signal *beyond* lineage. The thing that killed the embeddings (phylogeny dominance) becomes the thing the tool explicitly measures and strips.
**Kill-test (executed).** DEAD iff no residual signal survives phylogeny removal → cross-axis determinant clade-grouped AUC ≤ 0.55. Committed: **0.908** (naive 0.975) → SURVIVED. `[grounded]` real residual signal beyond lineage exists.

### The synthesis (the way forward)
> The genome world model is a **federation of de-confounded deterministic cells** (g1), whose **learned components live ONLY where interventional edit→measure data exists** (g9), tasked with **EXTENDING the catalog rather than replacing it** (g5), with **de-confounding promoted from a hidden guard to a first-class product surface** (g8). The single organizing lever is the **observational→interventional shift** — it is why the catalog wins, why embeddings lose, and where the one validated learned arm (HIV PhenoSense, ProteinGym DMS) already works. This is a genuinely different target than "a big model that reads genomes," and every piece of it survived an executed falsification against committed data.

---

## Unfalsified — next-investigations lane (named-only, NEVER "survived")

These have no *cheap* executable falsifier on committed data, so `/innovate` refuses to call them survived. They are the highest-value **next experiments** (each needs one new artifact to become falsifiable):

- **g9-objective-swap** — the catalog blind spot may be an **objective** problem, not capacity: likelihood/exchangeability scorers (ESM *and* BLOSUM62) call conservative substitutions benign by construction, so a model trained to predict **functional change** (not likelihood) on interventional DMS *at the exact conserved DRM sites* could fill the blind spot. *To falsify:* obtain a function-readout DMS overlapping the HIV NNRTI DRM positions and test whether a function-trained signal separates them where ESM/BLOSUM cannot. `[grounded]` BLOSUM62 ties ESM on DRM ranking; `[speculative]` a function objective separates them.
- **g5-sim-label** — relax *"labels must be free wet-lab"* → use **computational structure-based ddG / docking** predicted fold-change as a **weak bootstrap label** for census cells lacking a free wet-lab label (independent-but-simulated). *To falsify:* generate ddG for a committed protein and check it correlates with a real fold-change label (HIV / ProteinGym) above a usable threshold. `[grounded]` measured DMS already works; `[speculative]` computational ddG is a usable weak label.

---

## Handoff
The synthesis (g9 + g1 + g5-selfaware + g8) is the recommended reframe for the genome world model. **The post-run hardening changed the recommended first build:** g5's HIV instance is *already built* (`hiv_catalog_accessory_extension.py` + `hiv_blindspot_position_novelty.py`), and its honest form is the self-awareness FLAG (the additive-rule form is a committed negative — do not wire the accessory set into `hiv_amr.py`, it regresses balanced accuracy). So the two genuinely un-built, `/technical-plan`-ready first builds are:
- **g8-residual-detector** — package the existing `crossaxis_lineage_deconfound.py` de-confound (per-determinant clade-grouped AUC 0.908) as a reusable *product surface*: "which genome regions carry mechanism signal beyond lineage" as a first-class decoder output. Reuses shipped machinery; smallest genuine build.
- **g5-selfaware-flag GENERALIZED** — the position-novelty flag is validated for HIV NNRTI only; generalizing it to the bacterial/other cells (a per-cell "catalog may be incomplete here" self-awareness surface) is the un-built extension.

The two unfalsified ideas (g9-objective-swap, g5-sim-label) each still need one new artifact before they can be scored. Read-only over committed artifacts; frozen decoder surface byte-unchanged (verify_lock OK).
