# Phase 5 — the global genome LM: where my closure was wrong, and the decomposition

**2026-10-01.** The user's correction: *by definition we do not have wet-lab encodings for all possible life
forms, so the goal is a global genome LM that learns trends from abundant sequence and adapts to a new
organism from the few sequences we can get — a multi-shot genome LM.* This records where my 2026-10-01
"every technical gLM lever is measured and none pays" closure was **over-scoped**, what the measurements
genuinely do constrain, and the family decomposition + phased plan for the Phase 5 path.

---

## 1. Where I was wrong — three receipts, not an opinion

**(a) My negatives all live in ONE regime, and it is the smallest one.** `eval/regime.py` carries
`curated_catalog_exists = LOSES_TO_CATALOG`. Every result I cited to close the gLM was measured *inside*
that regime — zero-shot ESM2 against the HIV NNRTI catalog (0.926), capacity and modality on that catalog's
blind spot, "bacterial 5-drug: the catalog already captures the signal". **A curated determinant catalog
exists for a vanishingly small fraction of life.** Outside it the question is not "can a learned layer beat
the catalog" but "does anything work at all". My closure answered the first and was stated as if it
answered the second.

**(b) The project's own regime map has NO TRANSFER AXIS.** Regimes are indexed
`(population, endpoint, method)` — eight of them, none expressing *trained on organism A, applied to
organism B from few labels*. So the few-shot cross-organism question is **not merely untested, it is not
representable** in the taxonomy I used to declare the question closed. A grep for `few-shot|multi-shot`
across `wiki/`, `scripts/` and `plans/` returns **nothing**. The central claim of the Phase 5 thesis has
never been measured here.

**(c) The July scoping brief already said so, and I judged by the criterion it forbade.**
`plans/GENOME_WORLD_MODEL_Scoping_Brief_2026-07-04.md` states its honest prior: a genome world model
*"must be designed to **earn its keep on the polygenic residual + generative/design axis**, not to beat
mechanism scans head-on."* I beat it head-on and reported the result as a general verdict. Its adversarial
review goes further: *"generative/design thesis = MORE credible than phenotype-prediction… may be the
higher-VOI first win."*

**Corrected scope of the closure:** *within well-catalogued systems, a learned layer does not beat or
usefully exceed a curated catalog, and capacity/scale/modality do not change that.* That statement stands.
It says nothing about transfer to un-catalogued life, which is the Phase 5 objective.

---

## 2. What the measurements DO constrain — engineering input, not pushback

These are the hardest-won facts in the repo and they shape the design rather than forbid it.

| finding | consequence for a global gLM |
|---|---|
| `natural_organism_zeroshot = CLOSED_NEGATIVE` (0-for-5 de-confounded) | a gLM predicting organism-level phenotype from natural-population sequence, zero-shot, is **closed**. Phase 5 must not re-enter this cell. |
| Arabidopsis G2: the embedding learned **population structure**, not causal signal (within-group r² −0.13 vs structure-only 0.48) | **THE central technical risk.** In a model trained across all life, *phylogeny is the dominant axis of variation* — a global representation will encode taxonomy unless de-confounding is a first-class objective. |
| ESM2 scale **regresses** 650M → 3B → 15B | parameters are not the lever. "Global" must mean *breadth of modality and taxon*, not bigger. |
| orthogonal modality **WORKS**: +0.056, 90.5% of proteins (and a *weaker* model lifting a stronger one) | **supportive** — the measured lever is combining orthogonal information, which is what a multimodal genome model is. |
| constructed variation **WORKS**: yeast segregant 12/12 traits r 0.46–0.80; forward cell 0.761 / 0.352 | **supportive, and it names the de-confounding condition** — randomised ancestry is what makes learning work. Evaluate there. |
| `natural_molecular_supervised_blindspot = WORKS` (RT 0.8102 / IN 0.8923) | learned layers work on the **molecular** endpoint even in natural populations. The failures are organism-level. |

**The synthesis:** pretraining needs **no labels** (sequence is abundant — Evo2 used ~9.3 T nt / ~128 k
genomes). Few-shot adaptation needs **a few** labels per new organism, not a catalog. What the 0-for-5
establishes is that **validation is the hard part**, and that a phylogeny-confounded evaluation will
manufacture either a false positive or a false negative depending on direction.

---

## 3. Decomposition — the families

Updated from the brief's G1–G7, with the two families it was missing (**F1 transfer benchmark**,
**F5 few-shot adaptation**) — which are precisely where the user's objective lives.

| # | Family | Deliverable | Depends on | Gate |
|---|---|---|---|---|
| **F1** | **Transfer-axis benchmark + regime extension** | a TRANSFER axis in `eval/regime.py`; a k-shot protocol with held-out-**organism**/clade splits, the mandatory leakage audit (E16) and the baseline gauntlet + negative controls (E19) | — | **none — laptop, $0** |
| **F2** | Corpus assembly | multi-kingdom pretraining corpus + paired-modality eval sets, manifest-tracked | F1 (eval shape) | storage (D: 4.0 TB free) |
| **F3** | **Phylogeny de-confounding** | a validated protocol separating taxonomic from functional signal, calibrated on **known-answer** substrates (yeast segregant = works, Arabidopsis = fails) | F1 | none — $0 |
| **F4** | Backbone + continued pretraining | backbone choice (Evo2 / Caduceus / Gene42) + a domain-adapted checkpoint | F2, F3, **F9** | **COMPUTE (hard)** |
| **F5** | **Few-shot adaptation mechanism** | measured k-shot curves: in-context vs LoRA vs linear probe vs prototypical | F1, F4 | compute (small) |
| **F6** | Generative / design axis | does it propose edits satisfying **measured** molecular constraints better than existing design tools? Reuses forward/inverse + ProteinGym + PEAR | F1 | free GPU (Kaggle T4) |
| **F7** | Multimodal alignment (CLIP-style) | contrastive alignment sequence ↔ a phenotype modality | F2, F4 | paired labels (scarce) |
| **F8** | Hybrid fusion | learned layer carries the **residual**; the 115-cell deterministic surface stays primary where a catalog exists | F5/F7 + existing decoders | none |
| **F9** | Compute provisioning | a usable GPU path | — | **MONEY / AUTHORITY** |

**Critical path:** `F1 → F3 → (F6 ‖ F2) → F9 → F4 → F5 → F7 → F8`.

**Requirements flow-down (what each family owes the next):**
- F1 owes everyone the **verdict shape** — no gLM number is interpretable before it exists. It is the
  foundation precisely because the 0-for-5 was a *validation* failure as much as a model failure.
- F3 owes F4/F5 a **de-confounding protocol**; without it F4 reproduces the Arabidopsis outcome at scale.
- F9 gates F4 and everything downstream. **GTX 860M (CC 5.0, 4 GiB) cannot run these models** — CLAUDE.md
  records Evo and DNABERT-2 as not callable on this host.
- F8 owes the project **non-regression**: the frozen AMR surface and the 115-cell evidence surface must
  stay byte-intact; the learned layer is additive and fail-closed.

---

## 4. --plan — phases, and the first two are free

- **Phase 0 (now, laptop, $0) — F1.** Build the transfer benchmark + regime axis. *Falsifier:* it must be
  able to return a NEGATIVE on a model known to fail (the Arabidopsis embedding) and a POSITIVE on one
  known to work (yeast segregant). A benchmark that cannot do both is not a benchmark.
- **Phase 1 (now, $0 / free Kaggle T4) — F3 + F6.** De-confounding protocol on known-answer substrates, and
  the generative/design axis on existing cached embeddings + DMS data. **This is the cheapest decisive test
  of the whole thesis and it needs no new labels and no money.** Per the adversarial review: prove ONE axis
  before bundling.
- **Phase 2 (compute-gated) — F2 + F4 + F5.** Corpus at scale, continued pretraining, and the k-shot curves
  that are the actual Phase 5 claim. *Gate:* **F9.** Free Kaggle T4 covers small-scale; beyond that is a
  **money decision**.
- **Phase 3 (money, large) — F7 + F8, scale-up.** Only if Phase 1–2 earn it.

**Design rails carried from the brief, non-negotiable:** do NOT bundle MLM+AR+CLIP+JEPA before one
component beats simple baselines on a clean substrate (§11); E16 leakage audit is **mandatory**; E19
baseline gauntlet + negative controls (shuffled-label, geography-only, nearest-neighbour-sequence,
PCA-only, held-out-clade) required before **any** value claim.

---

## 5. The honest risk, stated once

The Phase 5 thesis is **not** refuted by this repo's negatives, and it is **not** supported by them either.
What the repo supplies is an unusually good set of **traps already mapped**: the organism-level zero-shot
cell is closed, phylogeny-confounding is the measured failure mode, scale is not the lever, and
orthogonal-modality + constructed-variation are the two levers that demonstrably pay. A Phase 5 build that
respects those four facts is a genuinely open bet. One that ignores them repeats a 0-for-5.

**Biggest single risk:** F3. If taxonomic signal cannot be separated from functional signal in a global
representation, the model will transfer *phylogenetic nearest-neighbour lookup* and score well on any
evaluation that does not hold out clades — which is exactly how a false positive would be produced here.

---

## 6. What needs the user

1. **Open the Phase 5 lane?** This supersedes a recorded closure (`Phase 5 DEFERRED`) and is a strategic
   direction decision, not a technical one.
2. **Compute (F9)** — the money gate. Free Kaggle T4 gets Phase 0–1 done; Phase 2+ needs a budget call.
3. **Ratify the family set** before ledgers are instantiated (9 live of a 25 cap; this adds up to 9 more).

Phases 0 and 1 are **$0 and decisive**, so the cheapest honest move is to run them before any budget
conversation.
