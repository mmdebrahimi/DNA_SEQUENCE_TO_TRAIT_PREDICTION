# Phase 2 Framing Brainstorm Plan

> Pre-/idea-anchor ideation pass capturing the Phase 2 anchor search space, axis critique, and 2 candidate anchor sentences after Phase 1 closeout 2026-05-17.

---

## Problem Statement

Phase 1 of dna-decode closed 2026-05-17 with an internal cross-drug architectural finding (`wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`). The synthesis stopping rule requires fresh /idea-anchor + /project-init before any Phase 2 experiment fires. This brainstorm is the pre-/idea-anchor ideation pass: **what is Phase 2 FOR? What is the project optimizing for now?**

Phase 1's anchor was "biologically interpretable E. coli AMR prediction platform; cipro + cef + tet binary R/S; foundation-model embeddings + classical baselines; LOMO-clade CV." That anchor produced infrastructure + an architectural finding + 7 residual uncertainty items. It did NOT produce a deployable classifier. The architectural finding partially invalidates the original anchor — frozen-NT-pooling is mechanism-class-bounded, not universal.

The Phase 2 anchor is genuinely open. project_state already calls Axis A1 (cohort expansion → Databricks burst) the Bellman row-1 next action, which is inertia carrying the original anchor forward instead of re-anchoring on Phase 1's actual learning.

## What Phase 1 produced (building blocks for Phase 2)

- Cross-drug architectural finding (mechanism-class-bounded NT-pooling) — concentrated-signal PASS (cipro QRDR 0.750 + cef plasmid β-lactamases 0.833); distributed mobile-element FAIL (tet 0.400 anti-predictive).
- 4-tier audit infrastructure (mechanism × MIC × opacity merge with structurally-enforced SUSPEND gate) — most generalizable artifact; 78 regression-guard tests.
- Smoke-gate runner drug-parameterized (`--drug`) + N=12 mini-cohort pattern + N=38 cipro cohort cache covering cef + tet via post-hoc filter.
- NT v2 100M frozen embedding cache at `D:/dna_decode_cache/embeddings/nt_n40_cipro.h5`.
- AMRFinderPlus + Bakta + Mash Docker toolchain installed + smoke-validated 2026-05-15.
- Stage 2 N=147 cipro cohort built (NT cache NOT populated — deferred from Databricks burst).
- Project-state ledger discipline + plan-language drift detection + 4-round brainstorm cycle + 2-layer verdict + statistical-sanity-check discipline (PC2 catch).

## Design Decisions

### D1: Axis menu structured by tier, not by topic

**Decision:** 4-tier scope (Product / Research / Horizon / Meta) covering 13 candidate axes (A1-D3).

**Rationale:** Axis labels make trade-offs visible: timeframe, artifact type, success metric, compute cost all vary by tier. Topic clustering would hide that frozen-NT-pooling work and audit-tooling work optimize for different deliverables.

### D2: Generative-ideation v2.1 payload pattern

**Decision:** Brainstorm payload structured to elicit BOTH critique AND generation in one round (axis critique + missing anchors + cross-domain analogs + sub-problem decomposition + ranking reassessment + user-intent read).

**Rationale:** Open-ended ideation use case per `/brainstorm` v2.1 pattern. Output is the framing artifact itself, not a refined implementation plan.

### D3: Stopping discipline preserved

**Decision:** Brainstorm produces framing options, NOT commitments. No Phase 2 experiment fires from this output. Fresh /idea-anchor + /project-init required before action.

**Rationale:** Phase 1 closeout synthesis stopping rule narrowly scoped reopens to (a) internal contradiction OR (b) factual/source mismatch. This brainstorm is pre-ideation, not a reopen.

## Axis Critique

**A1. Ship cipro classifier at strict-MIC N=150** — weak anchor. EP1 explicitly says no cipro Databricks burst because strict-MIC feasibility is uncertain and architecture remains a co-bottleneck; A1 is inertia unless a feasibility census first proves clean labels exist.

**A2. Ship audit infrastructure as product** — strong, but one cohort old. The mechanism × MIC × opacity SUSPEND gate is the most reusable artifact; risk is turning an emergent internal scaffold into a "toolkit" before proving it works on ≥2-3 non-cipro cohorts.

**A3. Mechanism-only production classifier (AMRFinder + acquired-gene panel)** — intellectually honest but low-upside. AMRFinder caught QRDR mechanisms where NT struggled; shipping a curated AMR panel mostly validates known biology, not the DNA-AI ambition.

**B1. Per-gene NT windows on tet-distributed-mechanism** — promising architecture repair, but premature unless tet label quality is audited. Tet mini-cohort was reused from the cipro pool; per-locus windows become meaningful only after mechanism/MIC audit says tet labels are interpretable.

**B2. 4th-mechanism-class smoke test (colistin mcr / aminoglycoside)** — best research-efficiency option. Already named as cheapest falsifier in the synthesis + future-features; fits solo/GTX constraints because it reuses the smoke runner and avoids N=150 embedding population.

**B3. Cipro pure architectural verdict at clean N=40 strict-MIC** — useful, but backward-looking. Resolves residual uncertainty; doesn't change cross-drug story unless it contradicts cipro's role in the partition.

**C1. Second organism (M. tuberculosis)** — too early. E. coli audit/model stack not yet stabilized; MTB would be a new data-governance and biology project, not a Phase 2 continuation.

**C2. MIC continuous regression head** — conceptually good, practically blocked by MIC sparsity. MIC ambiguity caused the EP1 SUSPEND gate; continuous MIC modeling should follow a MIC feasibility census, not precede it.

**C3. Pan-genome graph (Panaroo + GNN)** — attractive but overbuilt for current uncertainty. Project hasn't established clean labels or whether mechanism panels dominate; graph modeling adds architecture complexity before data substrate is trustworthy.

**C4. Multimodal DNA + image** — horizon only. CLAUDE.md marks it long-term, not a direct Phase 1 stepping stone; anchoring Phase 2 here would discard actual Phase 1 learning.

**D1. Publish EP1+EP2 architectural finding (arXiv/blog)** — not yet. Synthesis says external publication deferred; cef/tet are reused cipro-pool mini cohorts; publishable architecture claim needs BV-BRC-wide independent cohorts.

**D2. Publish audit framework methodology paper** — stronger than D1. Audit infrastructure is real, tested, methodologically general; publication-grade version needs ≥1 non-cipro replication of the SUSPEND/audit workflow.

**D3. Publish Claude Code skill ladder pattern** — useful meta-artifact, not dna-decode Phase 2. Separate process-writing project; would dilute the science/tooling anchor.

## Missing Anchors (not in original axis menu)

1. **Audit-first evidence-packet factory.** Phase 2 optimizes for converting any AMR drug into a gated evidence packet: strict-MIC census → mechanism/MIC/opacity audit → smoke model → decision artifact.
2. **Label-quality cartography before modeling.** EP1 failed because the cohort wasn't decidable; Phase 2 anchor = "map which BV-BRC AMR tasks are actually modelable under strict phenotype criteria."
3. **Architecture triage, not architecture repair.** Instead of fixing NT immediately, classify drug/mechanism regimes where whole-genome pooling is acceptable, useless, or needs per-locus windows.
4. **Mechanism panel as adversarial oracle baseline.** Treat AMRFinder/curated mechanisms not as the product but as the adversarial oracle every learned model must beat or complement.
5. **Evidence-to-publication runway.** Anchor on generating one credible external artifact, first sprint decides whether that's toolkit paper, architecture note, or portfolio writeup.

## Cross-Domain Analogs

1. **Small clinical ML projects** often pivot from classifier-building to dataset audit tooling after discovering label leakage or phenotype ambiguity. Publishable contribution becomes "how to know when not to train."
2. **Security ML** often starts with a detector, then becomes a benchmark/evaluation harness because the model is less reusable than the adversarial testbed.
3. **Astronomy/bioimage solo projects** often move from "build a model" to "curate a reproducible target-selection and quality-control pipeline" when data quality dominates model choice.
4. **Compiler/performance projects** often pivot from optimization algorithm to diagnostic profiler after Phase 1 shows the real bottleneck is workload classification.

## Sub-Problem Decomposition (top-axis ranking per goal variant)

- **Publish a paper:** D2 > B2 > D1. Audit-framework paper closest; B2 supplies replication; D1 waits for independent cohorts.
- **Ship open-source tool:** A2 > label-quality cartography > A3. Productize the audit gate, not the frozen-NT classifier.
- **Deepen cross-drug architectural claim:** B2 > B1 > B3. 4th mechanism class = cleanest falsifier; tet windows test repair after tet audit; cipro cleanup secondary.
- **Pivot to new organism:** C1 only after A2-style audit abstraction. Port the audit workflow first, not the model.
- **Maximize portfolio impact:** A2 + B2 hybrid > D2 > A1. "Tested AMR cohort-audit + architecture-triage system" stronger than "one classifier attempt."
- **Chase interesting questions:** B2 > B1 > C2. Science curiosity is mechanism-shape vs architecture, not raw classifier shipping.

## Candidate /idea-anchor Sentences

**Candidate 1 (preferred):**

> Build an audit-first AMR evidence-packet system for E. coli that determines, per drug and mechanism class, whether a cohort is label-decidable and whether frozen NT whole-genome pooling is appropriate, using strict-MIC, mechanism, opacity, and smoke-model gates before any scale-up.

**Candidate 2:**

> Validate and operationalize the Phase 1 mechanism-class architectural finding by running low-cost, audit-gated cross-drug experiments that separate cohort quality, curated mechanism signal, and NT pooling failure modes before committing to larger cohorts or new architectures.

**Why Candidate 1 preferred:** absorbs A2 (audit-as-product) + B2 (4th-mechanism falsifier) + Missing Anchors #1 (evidence-packet factory) + #2 (label-quality cartography) without forcing a single output artifact. Gives the BV-BRC strict-MIC 3-drug feasibility census (already named as Phase 2 entry in `project_state/dna-decode-2026-05-11.md`) room to be the first sprint without letting the census become the anchor. Honors the epistemic-infrastructure read of user intent; matches solo + GTX 860M + no-deadline constraints.

## Honest Read on User Intent

**Most likely:** user wants a research-grade, portfolio-visible scientific artifact, not merely a deployed classifier. Evidence: enforces preconditions, decision ledgers, synthesis artifacts, adversarial reviews, publication framing; cares about original biological interpretability goal enough not to accept a black-box or tautological classifier as success.

**Second-order:** drawn to modeling, but project's actual taste is **epistemic infrastructure** — knowing when evidence is valid. Strongest Phase 1 artifacts were the SUSPEND gate, statistical sanity checks, plan-drift detection, cross-drug falsification framing, NOT AUROC.

## Open Tradeoffs (for /idea-anchor session)

**Tradeoff A — classifier-tier vs audit-tier anchor:** original anchor was prediction platform; audit-tier reframe acknowledges Phase 1 didn't ship a model. Is the deliverable a predictive system or an evidence-validity framework?

**Tradeoff B — replication threshold for "audit framework is real":** minimum = cipro + one more drug (~2 weeks solo). Strong = cipro + cef + tet + 4th class with BV-BRC-wide cohorts (months).

**Tradeoff C — strictly E. coli vs portability as success criterion:** stable substrate vs framework generality test.

## Open Questions (must answer before /idea-anchor fires)

1. Is the desired external artifact a paper/blog, GitHub tool, or portfolio case study?
2. Should Phase 2 stay strictly E. coli, or is portability a required success criterion?
3. Minimum replication for calling the audit framework real: cipro + one more drug, or cipro + cef + tet + 4th class?

## Assumptions

- Phase 2 anchor should optimize for solo feasibility, low compute spend, preserving interpretability/science arc.
- A Phase 2 "product" can be an evidence framework or toolkit, not necessarily a predictive API.
- Next /idea-anchor should constrain action but avoid committing to Databricks, N=150, or new organism before the feasibility census.

## Next Step

Fresh session firing `/idea-anchor` with one of the two candidate sentences (or a merge), followed by `/project-init`. No Phase 2 experiment fires from this brainstorm output. Open Questions 1-3 answered as part of /idea-anchor input.
