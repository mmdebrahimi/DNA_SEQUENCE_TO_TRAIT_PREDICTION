# Big idea — name a trait, get the genome edits

<!-- mission-shape: greenfield -->

**The north star, in the user's words (2026-10-04):** *name a trait, get the genome edits.* Four steps:
(1) identify the organism — **SHIPPED** as `dna-identify`; (2) work out which sections jointly make a
phenotype; (3) **act in bacteria**; (4) climb to animals. Bacteria are the **starting substrate, not the
subject** — a substrate-derived exclusion is DEFERRED, never "inapplicable".

**The product is the generative engine; the deterministic catalog was a BRIDGE** built only so there was
a number to compare against. A lookup table cannot generate, so a head-to-head AUROC against it ranks
nothing. **Every discriminative zero-shot negative in this repo is evidence about SCORING, not about
GENERATION.**

## Where this mission actually stands (2026-10-09)

**Shipped and validated.** 128 registered cells across 8 tracks on PyPI `v0.13.x`; 33
`INDEPENDENT_MEASURED`. Two oracles that would score a design: `forward` (DMS ρ 0.35–0.76) and `fba`
(Keio 0.954 / MCC 0.70–0.74), plus `constraints/`, correct by construction.

**The measured bottleneck is the GENERATOR, not the scorer.** 2026-10-07: the oracle converts ρ 0.59
into 79% of achievable selection value at k=10. 2026-10-08: G-A and G-C both ran to **valid negatives** —
a 13,217-parameter conv encoder beat one GC number by a median +0.0122 whose paired CI includes zero
(6/10 seeds), and conditioning on medium added +0.0001. **Do not quote +0.0122 as a win.**
Self-distillation stays **forbidden** until a learned representation lands (student–teacher agreement
1.0000 on a linear oracle).

**The generator-bottleneck claim stands; the implied LEVER was wrong.** On this substrate a
346-parameter BIOPHYSICAL model that was never trained on it scores rho **0.3968** — above every purely
learned arm here (conv 0.3234, pyramid 0.3167, GC 0.3112). The only tile-substrate arm whose paired CI
EXCLUDES zero is the ORTHOGONAL HYBRID at **0.4254** (10/10 seeds), and its shuffled-conv control
**HURTS** (**-0.1232**, 0/10) — so the gain is channel CONTENT, not the averaging. That is a second
independent replication of rank-average-of-orthogonal-modalities (the first: ESM2+GEMME+ProSST, +0.056
on 90.5% of proteins). **It does NOT overturn the G-A encoder-gate negative**: conv-versus-one-GC-number
still has a paired CI including zero, and that is a DIFFERENT question from whether a learned channel
pays when combined with an orthogonal non-learned one. Declared as the arm
`orthogonal-modality-stack-promcalc-plus-N-channels`; **not measured** —
`fragment-rung-321k-rows-size-lever-LB-M9` (0.185) is against a different ceiling (0.7928) and is not
comparable to any number above.

**The label wall is lower than the record said, and the user has decided to act on it.** CABBAGE /
EMBL-EBI AMR Portal: 1.7M lab-measured genome–phenotype pairs, CC-BY, 8.3 MB, **9 of 10 gates pass**, and
through this repo's own `tier_for_isolate` it clears the STRICT bar on four cells where the 2026-05-18
census concluded *no* drug clears it. **Authority decision, user, 2026-10-09: this reopens the banked
public-label AMR track.**

**Diversity is already measured per-cell, and it is uneven — which is why gap 1 names the cell it does.**
`source_concentration_of_scored_set` ships in every one of the four scored-cell artifacts: **3 of 4 read
`source_concentrated`** (Campylobacter×cipro, Salmonella×ceftriaxone, Salmonella×cipro) and only
Escherichia×gentamicin reads `source_diverse` (10 database labels, largest share 0.3534, 21 countries).
Geography is the harder limit and it does NOT clear on any cell: the three concentrated cells are
**96-98% USA** (largest-country shares 0.9637 / 0.9775 / 0.9597) and even the `source_diverse` cell is
**58.4% USA**, so every cell rests on a single-country majority. **No BioProject count removes that**,
which is the honest reason gap 1's remaining `cabbage_bioproject_diversity.json` predicate is a
diversity question and not a naming one.

**F3 is CALIBRATED — the 2026-10-09 `NOT CALIBRATED` status is superseded.** The measurement ran
(`wiki/f3_deconfound_calibration_1001g_2026-10-10.json`, `smoke: false`) and returned
**`VALID_DIFFERENTIAL_RESPONSE`**: the de-confounder responds DIFFERENTIALLY to real structure —
median paired drop **+0.0441** on the natural population (1,135 Arabidopsis accessions) against
**0.0000** on the constructed bi-parental cross, with the negative control still beating its own
grouped null. Controls were chosen BY CONSTRUCTION from each panel's `design` field, never by `pc1`.
Carrying the artifact's own limit verbatim: *"A PASS bounds the INSTRUMENT's differential response. It
does not validate any downstream organism-level prediction."* **ONE positive and ONE negative panel**,
imputed positive-control genotypes, `n_perm=100`. The live status is `eval/regime.py`'s
`DECONFOUNDER_CALIBRATION` — read that, not this paragraph.

## What must NOT be re-proposed (closed on measured evidence)

- **Zero-shot embeddings on natural populations for organism-level g→p.** 0-for-5 de-confounded.
  Independently corroborated at far larger scale (Yu/Barquist 2025, 24k genomes, 6,740 models: models key
  on lineage markers and *more data does not rescue it*). **More labels do not address this** — CABBAGE
  does not reopen it.
- **ESM2 scale.** 650M is the peak; 3B and 15B regress. Modality, not parameters.
- **Absolute cross-assay magnitude for the inverse cell.** F4, measured 2026-10-09: scale transfer is
  tractable WITHIN an assay modality and degrades ACROSS it (VEFill). The field has no method. The
  inverse cell's rank-only deployability **matches** the state of the art rather than falling short.
- **The ETD/thesis-harvest route.** EThOS offline since the 2023 British Library attack; web search does
  not index thesis full text; and the manual search found CABBAGE in one query, so search volume was
  never the bottleneck.

**Screen every candidate against `eval/regime.py` and the ten gates in `eval/rejection_gates.py` BEFORE
calling it promising.** `natural × organism × zero_shot` is the one regime that REFUSES; a SUPERVISED
natural-population proposal gets `REQUIRES_DECONFOUNDING` — a condition, not a wall — and that condition
now carries the de-confounder's own `BUILT_NOT_CALIBRATED` status.

## Mission Gaps

### gap: amr-label-substrate-validated
<!-- strategy-budget: 3 -->

The user's reopen decision, made checkable. CABBAGE's one open gate is G8 (effective lineages, which
needs the assemblies clustered), and its source diversity is unmeasured at the only granularity that
counts — the `database` column is an overlapping aggregator label, not a BioProject, and the screen flags
its own G2 as "weak". A scored cell must then exist on a cohort that is leakage-gated against the repo's
own 1,509 manifest accessions.

**Flattering-selection risk, recorded.** The third predicate names ONE of four scored cells, and it is
the one whose own `source_concentration_of_scored_set` reads `source_diverse` (10 database labels,
largest share 0.3534, 21 countries). Naming the best of four flatters. It is defensible only because the
mission prose below carries the other three — 3 of 4 read `source_concentrated`, and all four are
96-98% USA. Stricter alternatives, available as a one-line amendment: name the WEAKEST cell
(`cabbage_scored_cell_salmonella_ciprofloxacin_2026-10-11.json`), or require all four.

- file-exists wiki/cabbage_lineage_metrics.json
- file-exists wiki/cabbage_bioproject_diversity.json
- file-exists wiki/cabbage_scored_cell_escherichia_gentamicin_2026-10-11.json

### gap: deconfounder-calibrated
<!-- strategy-budget: 2 -->

The instrument exists and its response to real structure is unmeasured, so `REQUIRES_DECONFOUNDING`
currently means "de-confound with a tool of unknown sensitivity". Calibration needs the 919 MB genotype
matrix for the ridge, and the grouping must be cross-checked against this repo's own
`ibs_distance_matrix` rather than the 1001G project's converted kinship — a disagreement there is a
finding to diagnose, not a number to choose between.

- file-exists wiki/f3_deconfound_calibration_1001g_2026-10-10.json
- file-exists wiki/f3_ibs_crosscheck_1001g_2026-10-10.json

### gap: generator-beats-its-baseline-panel
<!-- strategy-budget: 4 -->

The measured bottleneck. A generator must beat the **mandatory baseline panel** — mutation-count,
BLOSUM62, one-hot ridge, site-independent, shuffled null — on biologically-motivated splits only (random
splits inflate ~3x). Comparators are other GENERATIVE methods or a random-proposal null, **never the
catalog**. The highest budget sits here deliberately: it is the one gap where the project has a measured
negative on its first attempt and the lever is genuinely unexhausted.

- file-exists wiki/glm_generator_baseline_panel.json
- file-exists wiki/glm_generator_vs_random_proposal_null.json
- test-exit-0 pytest tests/test_glm_edit.py -q

### gap: trait-to-edit-end-to-end
<!-- strategy-budget: 3 -->

The north star's step 3, reduced to one checkable instance: name ONE bacterial trait, get a ranked edit
set, scored by a validated oracle rather than by the generating model. It must carry the inverse cell's
four shipped rails — no magnitude, per-protein oracle check, regime-B only, propose-k/assay-k/keep-best —
and an explicit statement of what it does not support.

- file-exists wiki/trait_to_edit_demo_2026.json
- file-exists wiki/trait_to_edit_rails_audit_2026.json

## Honest limits of this mission definition

- **The gaps are checkable, not sufficient.** `file-exists` proves an artifact landed, not that its
  measurement was sound. Every artifact still needs its own honest-limits section, and the standing
  surfaces (`cell_registry`, the report card, `eval/regime.py`) remain the trust record.
- **Φ bounds proposals, not quality.** Σ budgets = 12, so the mission is finite by construction; that is
  a termination guarantee, not evidence that any gap will close.
- **Step 4 (climb to animals) is CARRIED by a declared, unrun arm — not absent.** It is held by
  `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN`, declared against the promoter-activity
  ceiling in `soraya_runs/arms-ledger.json` and still `planned`. "Deliberately absent" read as CLOSED;
  it is neither closed nor in flight, and naming the carrier makes it findable. The eukaryotic Path-B
  arm remains a closed negative on de-confounded EMBEDDINGS, and the corrected regime map says the
  discriminating variable is population DESIGN rather than organism complexity — which is exactly why
  the live arm is a CONSTRUCTED eukaryotic panel and not another natural-population embedding run. This
  answers `substrate-is-not-the-subject` without spending Σ budget on a duplicate gap, so Σ stays 12.
  An animal GAP would still be premature until `generator-beats-its-baseline-panel` has a result.
- **Reopening the AMR track is a scope decision the user made, and its scope is narrow:** what CABBAGE
  contradicts is the premise that no usable label substrate exists. It does not reopen the embedding
  negatives and does not change the generative north star.
