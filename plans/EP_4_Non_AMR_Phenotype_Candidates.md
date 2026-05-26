# EP-4 Non-AMR Phenotype Candidates — Scoping Memo

> First conceptual jump out of AMR per `plans/Trait_Decoding_Roadmap.md` Phase 4. Surveys candidate bacterial non-AMR phenotypes + ranks for "smallest credible first slice." Stays scoped to E. coli for now (one axis at a time per the cef-cached-first discipline; new organism = Phase 3, not 4).

**Status:** DRAFT 2026-05-26.
**Anchors on:** `plans/Trait_Decoding_Roadmap.md` Phase 4 row; the project's north-star ("DNA input → phenotype + trait identification at gene level"; corrected per 2026-05-26 user clarification).
**Purpose:** identify which non-AMR phenotype is the cheapest credible first step. Final pick locks via `/idea-anchor + /project-init` when Phase 4 actually fires.

---

## Why scoping now

Phase 4 is the first conceptual jump out of AMR. Without a scoping pass, "non-AMR phenotype" stays unbounded and gets deferred forever. With it: when Phase 2 multi-drug AMR ships + the user is ready to widen, this memo names the 1-2 strong candidates + their dataset availability.

This is NOT a `/idea-anchor + /project-init` cycle — Phase 4 isn't firing yet. It's a pre-stage memo so when the cycle fires, the candidate evidence already exists.

---

## Candidate non-AMR bacterial phenotypes (E. coli)

Selection criteria for each candidate:
1. **Paired DNA + phenotype dataset availability** — public + downloadable + ≥50 strains.
2. **Mechanism-class shape** — concentrated signal (likely passes mean-pool architecture per EP1+EP2 2026-05-17 finding) vs distributed signal (would need EP-1.5 architectural fix first).
3. **Phenotype-measurement objectivity** — binary OR continuous OR categorical; how noisy is the measurement.
4. **Clinical / scientific relevance** — does ANYONE care about this prediction.
5. **Cohort size potential** — limited by data; smaller is OK for v0.2-style smoke.
6. **Architectural compatibility with current v0/v0.1 surface** — reuses cached-strain CLI or genome-input CLI; no rewrite.

### Candidate 1 — Growth rate / fitness on a defined medium

- **Definition:** doubling time OR maximum OD₆₀₀ on standard LB / minimal media; continuous regression target.
- **Mechanism shape:** distributed — growth rate is shaped by genome-wide metabolic + regulatory architecture, NOT a small mechanism catalog. Likely FAILS mean-pool architecture per EP1+EP2 pattern (analogous to tet's distributed failure mode).
- **Paired-data sources:**
  - **KEIO collection** (~3,985 single-gene knockouts; growth rate measured) — high signal but knockout-style, not natural variation.
  - **CGSC + E. coli strain collections** — variable growth measurements; cohort assembly hard.
  - **OD₆₀₀ from BV-BRC genomes** — NOT systematically measured; effectively unavailable.
- **Cohort size potential:** KEIO is 3,985 knockouts — large; but knockout-style means feature signal is "gene present / absent," NOT "natural variation." Different problem shape.
- **Clinical relevance:** LOW — growth-rate prediction from genome is academic; no buyer-pain-point.
- **Architecture compatibility:** v0 schema needs a `phenotype: growth_rate` field + regression head. v0.2-class refactor.
- **Verdict:** ❌ Skip. Mechanism is too distributed; data is KEIO-knockout-shaped not natural-variation-shaped; no buyer.

### Candidate 2 — Virulence factor presence / pathogenicity prediction

- **Definition:** binary classification ("pathogenic" vs "commensal") OR categorical ("EPEC / EHEC / ETEC / UPEC / EAEC / commensal" pathotype).
- **Mechanism shape:** **CONCENTRATED** — pathotypes are defined by acquired-gene clusters (stx for EHEC; afa/hra/papC for UPEC; ETEC enterotoxins; LEE for EPEC). Very similar shape to cef plasmid β-lactamases. Mean-pool architecture should work.
- **Paired-data sources:**
  - **EnteroBase E. coli** — large public WGS + curated pathotype labels (>200K isolates).
  - **NCBI Pathogen Detection** — taxonomy + isolate metadata; serotype-derived pathotype proxies.
  - **CGE VirulenceFinder** — gene-call output; comparable to AMRFinder's interface.
- **Cohort size potential:** N ≥ 500 easy; N ≥ 5K plausible.
- **Clinical relevance:** HIGH — clinical microbiology routinely needs pathotype calls; veterinary; food safety; outbreak investigations.
- **Architecture compatibility:** v0 schema: replace `--drug` with `--phenotype` (categorical: pathotype). Same NT mean-pool + XGBoost classifier (multiclass instead of binary). 
- **Verdict:** ✓ **Strong candidate.** Concentrated mechanism + huge public dataset + clinical relevance + minimal architectural change.

### Candidate 3 — Biofilm formation (binary / weak-moderate-strong)

- **Definition:** crystal-violet biofilm assay quantification; categorical (none / weak / moderate / strong) OR binary (former / non-former).
- **Mechanism shape:** **MIXED** — biofilm involves both concentrated (csgA/csgB curli; cellulose synthesis bcs genes; fimbriae fim) AND distributed (regulatory network — csgD, cpxA, OmpR/EnvZ — affecting many distributed genes). Probably borderline for mean-pool; EP-1.5 architecture decision matters.
- **Paired-data sources:**
  - **Subset of EnteroBase** with biofilm assays performed → small (~hundreds).
  - **Published cohorts** (e.g., UPEC biofilm studies) — typically 50-200 strains each; need to aggregate.
  - **Limited public WGS + biofilm pairing** — feasibility uncertain.
- **Cohort size potential:** N ≥ 100 may need aggregation; not as easy as pathotype.
- **Clinical relevance:** MEDIUM — biofilm is clinically relevant (catheter infections, persistent UTIs) but not a routine clinical call.
- **Architecture compatibility:** same as Candidate 2.
- **Verdict:** ⚠️ Defer. Mixed mechanism + smaller cohorts. Possible later candidate; not the first jump.

### Candidate 4 — Plasmid stability / horizontal gene transfer capability

- **Definition:** count or presence of plasmid replicons (rep typing); OR stability of plasmid maintenance.
- **Mechanism shape:** **CONCENTRATED** — plasmid replicons are well-defined; rep typing via PlasmidFinder gives gene calls.
- **Paired-data sources:**
  - **PlasmidFinder + BV-BRC genome metadata** — already in audit infra orbit.
  - **PLSDB** — plasmid sequence database.
- **Cohort size potential:** Large (>10K).
- **Clinical relevance:** MEDIUM — relevant for AMR-resistance-gene transfer prediction; less so as standalone.
- **Architecture compatibility:** trivially v0 schema with `--phenotype plasmid_rep_type`.
- **Verdict:** ⚠️ Skip for v0.2 jump. Too close to AMR (resistance-gene-transfer adjacent). The "jump out of AMR" framing benefits more from a genuinely different phenotype (pathotype).

### Candidate 5 — Lactose fermentation / sugar utilization

- **Definition:** classic E. coli phenotype: lactose+ vs lactose- (lac operon presence/function). Binary.
- **Mechanism shape:** **CONCENTRATED** — single operon (lacZYA + lacI regulator).
- **Paired-data sources:**
  - **CGSC + ATCC E. coli with lac phenotype** — small but high-quality.
  - **Public WGS from food-safety strains** — labeled some.
- **Cohort size potential:** Small (~50-200 explicitly labeled).
- **Clinical relevance:** LOW (academic / educational).
- **Architecture compatibility:** same as Candidate 2.
- **Verdict:** ⚠️ Skip. Too easy + too academic. Lac+ vs lac- could be solved by gene-presence alone; no AI value-add.

### Candidate 6 — Specific metabolic capability (e.g., tryptophan or methionine biosynthesis competency)

- **Definition:** auxotroph vs prototroph for specific amino acid; binary.
- **Mechanism shape:** **CONCENTRATED** — well-defined biosynthesis operons.
- **Paired-data sources:**
  - **ATCC + CGSC auxotroph collections** — clean labels but small.
- **Cohort size potential:** Small.
- **Clinical relevance:** LOW.
- **Architecture compatibility:** same as Candidate 2.
- **Verdict:** ❌ Skip. Same issue as Candidate 5 (too easy + too academic).

---

## Ranked candidates

| Rank | Candidate | Why |
|---|---|---|
| **1** | **Pathotype prediction (Candidate 2)** | Concentrated mechanism + huge public dataset (EnteroBase) + clinical relevance + minimum architectural change + biggest "jump out of AMR" payoff |
| 2 | Biofilm formation (Candidate 3) | Mixed mechanism; defer until EP-1.5 architecture decision lands |
| 3-5 | Plasmid stability / lactose / metabolic auxotrophy | All skipped — too AMR-adjacent OR too easy OR no buyer |

**Recommendation:** Pathotype prediction is the strong first non-AMR phenotype. When Phase 4 fires:
- Build E. coli pathotype cohort from EnteroBase metadata + WGS (≥500 strains; multiclass: EPEC / EHEC / ETEC / UPEC / EAEC / commensal).
- Train multiclass NT + XGBoost.
- Validate side-by-side vs CGE VirulenceFinder gene-call calls (analogous to AMRFinderPlus comparison for AMR).
- Emit v0 schema with `phenotype: pathotype`.

---

## Architecture / spec implications

For Phase 4 to ship, the v0 spec needs:
- `phenotype` field generalizes `drug` (cipro+cef+tet+gent → pathotype + others)
- `pipeline.py predict --phenotype X` CLI flag added
- `mic_tiers.py` extended → `phenotype_catalogs.py` (renamed; per-phenotype mechanism + threshold catalogs)
- `predict_proba` returns multiclass probabilities for categorical phenotypes

These are v0.2 changes. NOT urgent; just naming the eventual rework.

---

## Open questions for Phase 4 firing time

1. Is multiclass classification the right shape for pathotype (vs binary commensal-vs-pathogen)?
2. Does the v0 audit framework (mechanism × MIC × opacity merge) generalize to virulence?  Likely yes via VirulenceFinder analog of AMRFinderPlus.
3. Should Phase 4 ship with E. coli only, OR include a quick cross-organism check (Salmonella pathotype prediction is similar)?
4. Should EP-1.5 architecture decision finish FIRST (so any future distributed-mechanism phenotype has the right architecture)?

---

## Bottom line

**Pathotype prediction is the strong EP-4 entry point.** Defers nicely behind Phase 2 multi-drug AMR completion + EP-1.5 architecture decision. When ready: `/idea-anchor` on "E. coli pathotype prediction from genome + audit-aware multiclass output" → `/project-init` → execute.

For now: this memo is the pre-stage. NO code; NO `/idea-anchor`. Phase 4 doesn't fire until Phase 2 + EP-1.5 settle.
