# v1+ Horizon Framing — Draft for /idea-anchor + /project-init

> Drafts the "AI DNA decoder tool at maturity" framing as input for a future `/idea-anchor` + `/project-init` cycle. NOT auto-invoking those skills; producing the candidate anchor sentences + scoping doc the user can review + run the skills against when ready.

**Status:** DRAFT 2026-05-24.
**Anchors on:** north-star clarification 2026-05-18 ("AI DNA decoder tool, not papers"); v0 closeout 2026-05-24 ("cached-strain cipro predictor"); v0.1 ingestion contract 2026-05-24 (two paths in flight).
**Pending Decisions row:** "Pre-Phase-2: /idea-anchor + /project-init proper entry" (overdue since 2026-05-17).

---

## Why this doc exists

The 2026-05-18 north-star reset said "AI DNA decoder tool, not papers" + locked the v0 spec. But v0 → v0.1 → v0.N is now incremental drift without an articulated destination. Without a v1+ horizon framing:

1. v0.1 decisions feel ad-hoc (Path G vs Path C with no v1+ anchor saying which one is more load-bearing).
2. Future versions inherit whatever architecture v0.1 ships, regardless of whether it scales to the mature tool.
3. The "decoder" label is being used in two contradictory ways: Codex relocked it to mean "cached-strain predictor" for v0, while the original intuition + closeout-handoff's v0.1 recommendation imply "ingests novel genome + emits resistance call."

This doc proposes candidate anchor sentences for the mature tool + scoping decisions that would resolve those ambiguities. Once locked, downstream version planning becomes a backward-induction from the v1+ definition rather than forward drift.

---

## Candidate /idea-anchor sentences (3 framings, user picks)

### Framing 1 — "Honest AMR resistance predictor" (narrow, ship-friendly)

> **An open-source CLI tool that takes a user-supplied E. coli genome assembly (FASTA + optional GFF3) and emits a per-drug resistance prediction with audit-grade provenance: which gene-level regions drove the call, how trustworthy the underlying training labels are, and explicit per-strain attribution-scope confidence.**

Scope:
- Multi-drug (cipro / cef / tet / gent at minimum).
- E. coli only at v1; multi-organism deferred to v2+.
- CLI-first; no UI, no REST.
- Honest output discipline (audit verdict + scope-limit propagation) is the load-bearing product feature.

Trade-offs:
- **Pro:** closes the original "decoder" intuition (user input → decoded output) without overreaching.
- **Pro:** clearly distinct from existing AMR tools (AMRFinder, CARD-RGI, ResFinder) — those emit gene calls; we emit phenotype predictions WITH honest-output framing.
- **Con:** does NOT generalize to "DNA decoder for arbitrary trait." Locks the project into AMR for the foreseeable future.

### Framing 2 — "Bacterial G2P platform" (medium, what was originally pitched)

> **A research platform that maps DNA sequences to phenotypic predictions for bacterial organisms, starting with E. coli AMR + expanding to other phenotypes (growth rate, virulence factors, plasmid stability) and other organisms (Klebsiella, Pseudomonas, etc.) over a 12-24 month horizon. Each prediction carries audit-grade provenance.**

Scope:
- AMR is one of many phenotypes; not the terminal product.
- Multi-organism within bacteria.
- Same "honest output" discipline as Framing 1.
- Still CLI-first at v1; web UI possible at v2+.

Trade-offs:
- **Pro:** the original Phase 1 plan's mid-term + long-term horizons match this framing.
- **Con:** broad enough that "shipping v1" requires a clearer definition of which non-AMR phenotype lands first.
- **Con:** "research platform" pulls toward papers/benchmarks again; conflicts with the 2026-05-18 north-star clarification.

### Framing 3 — "Multimodal genotype-phenotype decoder" (broad, original goal)

> **A multimodal AI system that maps DNA sequences to a wide range of phenotypic outputs — clinical labels (AMR, virulence), continuous measurements (MIC, growth rate), AND eventually paired image phenotypes (colony morphology, microscopy) — for both prokaryotic and eukaryotic organisms over a 3-5 year horizon.**

Scope:
- Multi-organism (bacteria + fungi + eventually eukaryotes).
- Multi-modal output (categorical + continuous + image).
- Research-scale: requires GPU cluster + paired-dataset acquisition.

Trade-offs:
- **Pro:** matches the original Phase 1 plan's long-term vision verbatim.
- **Con:** unbounded scope — research-program territory, not project-shape work (per the 2026-05-11 `/project-init` verdict).
- **Con:** the 2026-05-18 framing reset EXPLICITLY rejected the multimodal framing for paper-driven framing reasons. Reopening that decision needs an explicit acknowledgment.

---

## Recommended framing (subject to user override)

**Framing 1 — Honest AMR resistance predictor.**

Rationale:
- Aligns tightest with v0 (cached-strain cipro predictor) + the v0.1 candidate paths (Path G = genome-input + multi-drug; Path C = cef expansion).
- Honest-output discipline is what differentiates the tool from existing AMR-gene callers — that's the unique product moat.
- Solo-researcher cost is bounded; multi-organism / multi-modal expansion can become Framing 2 / 3 if Framing 1 hits diminishing returns.
- North star ("AI DNA decoder tool, not papers") is satisfied by Framing 1 at v1; Framing 2 reintroduces research-platform ambiguity.

Trade-off explicit: if user wants Framing 2 or 3, the v0.1 path priority changes — Path G (genome-input) is correct for all three, but Path C (cef-cached) is over-specific for Framing 1 (cef is one of many drugs we want), under-specific for Framing 2/3 (cef alone doesn't open multi-organism). The right Path C scoping IS "ship cef the same way cipro shipped, then do tet/gent the same way, then move to other organisms."

---

## v1 success criteria (Framing 1)

If user picks Framing 1, v1 ships when:

1. **Functional:** `pipeline.py predict --genome-fasta X.fna --annotations Y.gff3 --drug <D>` runs end-to-end on a user-supplied novel E. coli genome for any of {cipro, cef, tet, gent}, emitting v0-schema JSON + markdown.
2. **Predictive:** primary CV AUROC ≥ 0.70 on each of cipro/cef (where leakage-safe cohorts exist); tet + gent ship with `attribution_scope_confidence=INDETERMINATE` if cohort feasibility blocks training.
3. **Interpretable:** `top_k_attribution` recovers at least one known-mechanism gene for ≥ 50% of true-positive R predictions on cipro + cef; tet + gent inherit v0's "exploratory attribution" framing.
4. **Honest:** every prediction carries an `attribution_scope_confidence` field; every JSON output is accompanied by a scope-limit doc when below HIGH confidence.
5. **Reproducible:** one-command quickstart in README produces an end-to-end prediction on a public example genome.
6. **Operational:** the cross-machine handoff overhead is bounded — `scripts/cross_machine_sync_check.py` reports IN-SYNC after each session; no silent v0-style divergence.

---

## v1 non-goals (Framing 1)

- Streamlit / web UI / REST API.
- Multi-organism (Klebsiella, Pseudomonas, etc.).
- Continuous MIC head.
- Image / multimodal output.
- AlphaFold-style structural-context modeling.
- Federated / multi-laboratory deployment.
- External benchmark paper / arXiv publication.

These are explicitly v2+ candidates, NOT v1.

---

## v2+ horizon candidates (Framing 1 expansion path)

Ordered roughly by leverage:

1. **2nd organism (Klebsiella).** Reuses the entire v1 architecture; only needs a new cohort + retrain. Tests cross-organism transfer.
2. **Continuous MIC head** instead of binary R/S. Provides finer-grained output; needs MIC regression model + calibrated prediction intervals.
3. **Per-gene NT-window architecture** (replaces mean-pool with attention). The mechanism-class-bounded architectural finding (`wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`) implicates mean-pool dilution on distributed mobile-element resistance; this is the architectural fix.
4. **Federated / multi-laboratory cohort acquisition.** Brings cohort sizes from N≈150 to N≈10K. Mostly a data + ops play.
5. **Multimodal output (paired DNA + image).** v2+/v3+ candidate; gated on dataset existence.

None of these are committed; they're horizon options for after v1.

---

## Open Questions for the user (must answer before /idea-anchor + /project-init fire)

1. Which framing? (1 / 2 / 3 / other)
2. v1 timeline: 3 months / 6 months / 12 months?
3. Compute substrate at v1: stay on Precision 7780 (RTX 3500 Ada) OR formalize Databricks burst as the canonical training environment?
4. Honest-output discipline: keep as v1 hard gate, OR ship without it for parity with existing AMR tools?
5. v1 ship vehicle: GitHub release OR PyPI package OR Docker image OR all three?
6. Open-source license: pick MIT / Apache 2.0 / GPL / other?
7. Co-authorship attribution: solo (you) OR add Codex / Claude / others to README?

---

## Recommended /idea-anchor + /project-init invocation (when ready)

```
/idea-anchor An open-source CLI tool that takes a user-supplied E. coli genome assembly (FASTA + optional GFF3) and emits a per-drug resistance prediction with audit-grade provenance: which gene-level regions drove the call, how trustworthy the underlying training labels are, and explicit per-strain attribution-scope confidence.
```

Then:

```
/project-init
```

The `/project-init` skill will run its 3-sub-gate empirical-concerns / project-vs-research-program / refinement-candidates check + produce a refined goal hierarchy + Bellman frame for the v1 horizon.

---

## What this plan deliberately does NOT cover

- **Auto-invocation of /idea-anchor or /project-init.** Those are user-driven skills; this is INPUT.
- **v0.1 path picks.** Path G vs Path C lives in `plans/v0.1_Ingestion_Contract_Plan.md`; this doc is about v1+, not v0.1.
- **Specific v2+ implementation plans.** Each v2+ candidate would need its own technical plan when it's selected.
- **Co-funding / sponsorship / commercialization.** Solo project, no current external interest; revisit if the tool gains usage.
- **Comparison benchmarks vs existing tools (AMRFinder, RGI, ResFinder).** Useful for v1 validation; not a planning concern.

---

## Bottom line

Pick **Framing 1** (recommended) → run `/idea-anchor` with the candidate sentence → `/project-init` for the v1 horizon → the project has a destination v0.1 can be planned backward FROM, rather than v0.1 being planned forward INTO nothing.
