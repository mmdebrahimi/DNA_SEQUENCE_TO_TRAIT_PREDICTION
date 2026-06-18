---
title: Rough genome map (evidence-tiered functional annotation)
slug: genome-map
status: design
created: 2026-06-17
updated: 2026-06-17
related_plans: []
related_tests: []
related_decisions: [wiki/north_star_distance_brainstorm_2026-06-17.md, wiki/genome_map_idea_anchor_prompt_2026-06-17.md]
related_modules: [dna_decode/eval/amr_rules.py, dna_decode/organism_rules/, dna_decode/data/calibrated_amr_rules.json]
---

# Rough genome map (evidence-tiered functional annotation)

## Summary
A personal tool that turns a single microbial genome into an HONEST, evidence-tiered map of what its
parts likely do — molecular function (homology/annotation) + curated-determinant phenotype overlays —
where every feature carries an explicit evidence tier and the unknown rate is reported. The achievable,
label-free form of the dna_decode north star ("DNA → what its parts do"); explicitly NOT a learned
genotype→phenotype predictor (that arm is a closed negative on free data).

## Problem
Understanding "what the parts of a genome do" today means running several external annotation tools and
reading their raw output, which (a) BLURS confident curated function with weak homology guesses, (b)
OVER-STATES phenotype (association/annotation metadata read as causal), and (c) HIDES how much is
genuinely unknown. The project also has validated determinant/AMR cells but no per-genome view that
situates them among the rest of a genome's features. There is no single honest, confidence-aware view to
explore an unfamiliar genome or to sanity-check the existing determinant calls in genomic context.

## Desired Outcome
Point the tool at ONE microbial genome and get a per-feature map where: each feature carries an explicit
evidence tier; PHENOTYPE claims appear ONLY where a validated determinant cell fires; and the unknown rate
is a reported headline. This enables the developer to (1) understand an unfamiliar genome at a glance and
(2) QC / sanity-check the existing AMR/determinant calls in context. Usefulness is judged against a UX/eval
gate (defined in `design`), NOT a learned-model metric.

## Repo Fit
- REUSE: the project's honesty/evidence-tier discipline (SUSPEND/PLUMBING/abstain rails, lineage
  disclosure) + the determinant PHENOTYPE cells — frozen AMR surface (`dna_decode/eval/amr_rules.py` +
  `dna_decode/data/calibrated_amr_rules.json`) and the non-frozen organism cells
  (`dna_decode/organism_rules/`: TB; the fungal cell). NO edit to the frozen AMR surface.
- NEW: an annotation/homology layer driven by FREE external tools (Bakta/Prokka structural annotation;
  Pfam/hmmer + eggNOG/COG homology-function) following the existing Docker / D:-cache patterns.
- `viz/` already exists with a deferred genome browser (out of v1 scope).

## Scope
**In v1 (the scoping spike):** an evidence-tier schema; a per-genome JSON-first tiered map (+ a flat
feature table); the phenotype overlay from existing determinant cells; molecular-function / homology tiers
from free tools; a reported unknown rate; a UX/eval gate; a prototype on 2-3 genomes; a GO/NO-GO verdict on
whether to invest in catalog integration.

**Out of v1:** any learned/embedding phenotype model (closed-negative arm); a visual genome browser
(deferred `viz/`); catalog integration at scale (gated on the spike verdict); non-microbial / eukaryotic
genomes; phenotype claims outside the validated-determinant wall; multi-genome comparison (the
decision-support intent was NOT chosen — this is an exploration/QC tool).

## Actors
- **Primary:** the developer (Farshad) — exploring / QC-ing a single microbial genome. Solo hobby
  project; NO external user. "Useful" therefore means "lets me understand + sanity-check a genome I don't
  already know," not "ships to a microbiologist."

## Workflow / Behavior
**v1 scope (probe-tightened, 2026-06-17): the "Bakta honesty report".** Built ONLY on the
already-installed Bakta(db-light) + AMRFinderPlus + the EXISTING GFF parser (`annotations.py::parse_gff3`)
+ the EXISTING Bakta runner (`pathotype_laptop_pipeline.bakta_annotate`). NO hmmer/Pfam/eggNOG (uninstalled,
tens of GB). The map RE-TIERS Bakta's own annotation for honesty + overlays determinants — it does not
produce new annotation. This is a function/QC map with determinant call-outs, NOT a phenotype map (for an
arbitrary genome, no-determinant is the MODAL case).

Actor: the developer, one genome at a time. Happy path:
1. **Input** — one microbial genome: a FASTA (+ an optional precomputed GFF). Single genome; no batch/compare.
2. **Annotate** — run Bakta (db-light) → GFF3 via the existing runner, OR accept a provided GFF. **Bakta
   annotation is unproven on this host (smoke deferred 2026-05-15) → a Bakta smoke is the FIRST build action;
   a wedge yields a BLOCKED artifact, never a fake map.**
3. **Tier each feature** — assign every feature its single HIGHEST-confidence tier (4 v1 tiers; table in
   Rules) from Bakta's product/gene_symbol wording-confidence; retain lower evidence as secondary.
4. **Determinant overlay** — features matching a validated determinant cell get a phenotype annotation
   (drug/property + provenance) from the EXISTING cells; the genome-level R/S verdict (`call_resistance` /
   `tb_amr`) is shown SEPARATELY, tier-tagged. The map NEVER promotes a function/homology hit to phenotype.
5. **Output** — a JSON-first per-feature map + a flat table. Headline metrics: per-tier counts + the
   **DB-labelled unknown rate** (`unknown_under_bakta_db_light`) + the determinant-phenotype feature list.

Alternate / failure paths:
- **No annotation tool / offline** → degrade to the tiers achievable from a provided GFF + the determinant
  cells; everything else is `unknown` (offline-safe degradation, mirroring the BLAST/vf_diff pattern). The
  map still emits, with a degraded-coverage flag — never errors out.
- **A determinant cell ABSTAINs/SUSPENDs** on this genome → the feature keeps its annotation tier; the
  phenotype overlay shows `ABSTAIN`, never a forced call.
- **Conflicting tool annotations** for a feature → record both, take the highest-confidence tier, flag the
  conflict.

## Rules And Decisions
**Evidence tiers (precedence high → low). A feature gets ONE primary tier (highest applicable); lower
evidence retained as secondary.**

**v1 tiers = 4 (pathway DEFERRED — no KEGG in v1). Precedence high → low; one primary tier, lower evidence
retained as secondary.**

| tier | fires when (v1, Bakta-only) | may emit a PHENOTYPE claim? |
|---|---|---|
| `determinant-phenotype` | a validated determinant cell matches this feature (AMRFinder `main.tsv` row + `calibrated_amr_rules.json` class→drug; or `tb_amr`/organism determinant) | **YES — the only tier that may** |
| `curated-molecular-function` | Bakta gives a named `gene_symbol` + a SPECIFIC product (a confident call) | no |
| `homology-only-hypothesis` | Bakta low-confidence wording (`putative`/`probable`/`by similarity`/`domain-containing`/`uncharacterized`) | no |
| `unknown` | `hypothetical protein` / empty product | no |
| ~~`pathway-module`~~ | KEGG/MetaCyc module — **DEFERRED to a later version** (no KEGG; a GO-only add) | n/a |

**Load-bearing rules:**
- **Phenotype wall (R1):** a phenotype/property claim attaches ONLY at `determinant-phenotype`.
- **DB-labelled unknown rate (R2):** the unknown rate is REPORTED but its field name carries the DB/version
  coverage caveat — `unknown_under_bakta_db_light` (db-light has reduced functional coverage, so a bare
  "unknown rate" would mislead as biology when it is partly tooling-coverage). Never hidden/back-filled.
- **Highest-confidence primary (R3):** highest tier primary; lower evidence kept as observable secondary.
- **No frozen edit (R4):** the phenotype tier READS the existing cells; never modifies
  `dna_decode/eval/amr_rules.py` / `calibrated_amr_rules.json`.

**UX / EVAL GATE (the make-or-break — probe-tightened so relabelling alone CANNOT pass).** The prototype
passes iff, on the 2-3 test genomes, BOTH hold:
- **G1 Prevent-wrong-inference:** ≥3 concrete features on ≥1 genome where the tiered map PREVENTS a wrong
  inference a reader would make from the raw Bakta TSV — (a) a raw product taken as fact (`putative`/`by
  similarity`) that the map DEMOTES to `homology-only-hypothesis`; or (b) a determinant the map SURFACES as
  `determinant-phenotype` that raw lists as a plain gene. **At least 1 of the ≥3 must be (a)/(b) — "the
  unknown rate exists" alone does NOT satisfy G1** (that would be the catalog-stacking-busywork failure).
- **G2 No-tier-confusion:** in a spot-check, zero non-`determinant-phenotype` features are mistakable as
  phenotype claims.

## Edge Cases
- **Gene that is BOTH curated-function AND a determinant** → `determinant-phenotype` (highest); function
  kept as secondary.
- **Multi-domain protein, mixed evidence** → per-FEATURE primary tier in v1; per-domain tiering is a later
  refinement (recorded, not built).
- **Plasmid vs chromosome** → the map notes the compartment when known (the project's plasmid-marker
  gotcha — a determinant's compartment matters for interpretation).
- **High unknown fraction (>X%)** → a VALID, honest output (high unknown rate reported), NOT a failure.
  The whole point is preserving unknowns.
- **Genome with no determinant hits** → phenotype tier is empty; the map is still useful (function/homology
  tiers + unknowns). No phenotype is fabricated. (This is the MODAL case for arbitrary genomes.)
- **Bakta annotation wedges / unavailable** (the documented Docker-mount-corruption history) → emit a
  `BAKTA_ANNOTATION_BLOCKED` artifact, never a fake map; the offline path (provided GFF) is the fallback.
- **db-light inflates the unknown/homology fraction** → that is expected; the DB-labelled rate
  (`unknown_under_bakta_db_light`) makes it non-misleading. Full Bakta DB / Pfam-eggNOG would shrink it but
  is out of v1 scope.

## Acceptance Criteria
**A. Map production**
- **AC1** — Given one microbial genome (FASTA, or FASTA+GFF), the tool emits a per-feature map (JSON) + a
  flat feature table; every feature has exactly ONE primary evidence tier from the 5-tier schema.
- **AC2** — The output reports headline metrics: per-tier feature counts, the **unknown rate**
  (unknown-tier features / total), and the list of `determinant-phenotype` features.

**B. Honesty rails (the differentiator)**
- **AC3 (phenotype wall)** — No feature outside `determinant-phenotype` carries any phenotype/property
  claim (the phenotype field is non-empty ONLY on `determinant-phenotype` features).
- **AC4 (unknown visibility)** — The unknown rate is a top-level reported field equal to
  `count(unknown)/total`; unknowns are never relabelled to a lower tier.
- **AC5 (highest-confidence primary)** — A feature with multi-tier evidence takes the highest-precedence
  tier as primary, with lower evidence retained as observable secondary evidence.
- **AC6 (no frozen edit)** — A run leaves `dna_decode/eval/amr_rules.py` + `calibrated_amr_rules.json`
  byte-unchanged; phenotype-tier data is read-only from the existing cells.

**C. Phenotype overlay correctness**
- **AC7** — A feature matching a curated determinant is tagged `determinant-phenotype` with the
  drug(s)/property it contributes to AND the provenance (which rule path produced it); the genome-level
  R/S verdict appears separately, tier-tagged.
- **AC8** — When a determinant cell ABSTAINs/SUSPENDs for this genome, the affected feature shows
  `ABSTAIN` (not a forced call) and the genome-level overlay propagates the abstain.

**D. The UX/eval GATE (spike go/no-go)**
- **AC9 (G1 differentiation)** — On ≥1 of the 2-3 prototype genomes, the spike documents ≥3 concrete
  features where the tiered map is more honest/useful than the raw annotation TSV (a demoted homology
  guess, a surfaced determinant-phenotype, or the reported unknown rate raw output omits). Zero such
  features ⇒ NO-GO.
- **AC10 (G2 no-tier-confusion)** — In a spot-check of the prototype maps, zero non-`determinant-phenotype`
  features are mistakable as phenotype claims.
- **AC11 (verdict)** — The spike emits a GO/NO-GO verdict on catalog integration, backed by the AC9
  evidence + the AC10 spot-check + per-genome unknown rates.

**E. Robustness**
- **AC12 (offline-safe)** — With no external annotation tool available, the tool still emits a map (tiers
  from a provided GFF + the determinant cells; the rest `unknown`) + a degraded-coverage flag; it does not
  error.

## Test Scenarios
- **Happy path** — an E. coli cohort strain (rich AMR): map has `determinant-phenotype` features (e.g.
  QRDR) + function/homology/unknown tiers + a reported unknown rate (AC1, AC2, AC7).
- **Edge — no determinant hits** — a genome with no curated determinants: phenotype tier empty; map still
  emits function/homology/unknown; no fabricated phenotype (AC3 on empty set).
- **Edge — high unknown** — a homology-heavy / hypothetical-protein-rich genome: high unknown rate reported
  as a VALID output, not a failure (AC4).
- **Failure / malformed** — malformed/empty FASTA ⇒ a clear error, not a partial map; no annotation tool ⇒
  the offline degradation path (AC12), not a crash.
- **Precedence / overlap** — a gene that is both a curated function AND a determinant ⇒
  `determinant-phenotype` primary, function secondary (AC5).
- **Precedence — abstain** — a determinant cell ABSTAINs ⇒ feature shows `ABSTAIN`, overlay propagates
  (AC8).
- **Regression-sensitive** — frozen AMR surface byte-unchanged after a run (AC6); the phenotype wall holds
  (no phenotype outside the determinant tier) across all 2-3 prototype genomes (AC3).
- **Gate** — the G1/G2/verdict scenarios run on the 2-3 prototype genomes (AC9, AC10, AC11).

## Technical Notes
_(product-facing only; no implementation steps)_
- New code likely under `dna_decode/genome_map/` + a CLI emitting the JSON map + a flat feature table
  (matches the project's `.json` + `.md`-sidecar artifact convention).
- Phenotype tier READS existing surfaces: `dna_decode/eval/amr_rules.py::call_resistance` /
  `cipro_determinants_from_main` (AMRFinder `main.tsv` per-gene rows) + `calibrated_amr_rules.json`
  (class→drug) for AMR; `dna_decode/data/tb_who_catalogue.load_determinants` + `organism_rules/tb_amr` for
  TB; the fungal cell. NO edit to the frozen AMR surface (guarded by the existing leak-guard pattern).
- External annotation (Bakta/Prokka, Pfam/hmmer, eggNOG) via the existing `tools/docker_runner` +
  D:-cache pattern. Offline-safe degradation required (tiers achievable without the heavy tools).
- `viz/` holds a deferred genome browser — OUT of v1 (JSON-first + table only).
- Compatibility: additive new module; touches no frozen file; reuses the determinant catalogs read-only.
- **v1 reuse (probe):** `scripts/pathotype_laptop_pipeline.py::bakta_annotate` (existing Bakta runner) +
  `dna_decode/data/annotations.py::parse_gff3` (existing GFF→AnnotationTable) — the map assembles existing
  pieces. **NO hmmer/Pfam/eggNOG in v1** (uninstalled; the homology tier comes from Bakta's own wording).

**Spec constraints (product contract for /technical-plan):**
- **Prototype genomes (3):** one E. coli cohort strain (rich AMR overlay), one M. tuberculosis H37Rv-relative
  or C. auris (second organism + existing cells), one homology-heavy / hypothetical-protein-rich bacterium
  (stresses the homology + unknown tiers — the honesty stress test).
- **The JSON map IS the product contract:** per-feature `{primary_tier, secondary_evidence[], phenotype
  (determinant-tier only, with drug + provenance)}` + top-level `{per_tier_counts, unknown_rate, verdict}`.
- **Hard requirements:** offline-safe degradation (AC12); frozen-surface byte-unchanged (AC6);
  phenotype-wall (AC3) is non-negotiable.
- **UNRESOLVED product dependency (confirm in /probe BEFORE /technical-plan):** the free annotation/homology
  tool stack (Bakta/Prokka + Pfam/hmmer + eggNOG/COG) must install + run on this Windows/GTX-860M/Docker
  host. If a tool is infeasible, its tier degrades to `unknown` (the offline path) — but the spike's G1
  differentiation test still needs ≥1 working homology source to be meaningful.
- **Open product decision:** exact N for G1 is set at ≥3 features on ≥1 genome — revisit if too strict/loose
  once the prototype runs. Per-domain tiering is OUT of v1 (per-feature primary tier only).

## Repo grounding

### Captured by: intake @ 2026-06-17
- Files read: dna_decode/eval/amr_rules.py (referenced — frozen AMR surface), dna_decode/organism_rules/ (TB cell built this session), dna_decode/data/calibrated_amr_rules.json, wiki/north_star_distance_brainstorm_2026-06-17.md, wiki/genome_map_idea_anchor_prompt_2026-06-17.md
- Key claims:
  - The phenotype tier reuses the EXISTING determinant cells (frozen `amr_rules.py` + non-frozen `organism_rules/`); the feature must NOT edit the frozen surface.
  - The learned genotype→phenotype arm is a CLOSED NEGATIVE on free data (3 de-confounded tests) — out of scope by decision, not omission.
  - The honesty discipline (tier labels, phenotype-behind-the-validated-wall, reported unknown rate) is the differentiated contribution; raw annotation tools (Bakta/eggNOG) already exist, so the spike must test differentiation-over-raw-Bakta.

### Captured by: design @ 2026-06-17
- Files read: dna_decode/eval/amr_rules.py (call_resistance / cipro_determinants_from_main / qrdr_point_determinants — per-gene AMRFinder main.tsv rows), dna_decode/data/calibrated_amr_rules.json (rules: class→drug map), dna_decode/data/tb_who_catalogue.py (Determinant / load_determinants)
- Key claims:
  - The per-FEATURE phenotype overlay is feasible: AMRFinder `main.tsv` rows are per-gene; `calibrated_amr_rules.json.rules` maps determinant class→drug; `call_resistance` gives the genome-level R/S summary. TB uses `load_determinants` (per-variant). So a feature can be tagged `determinant-phenotype` + which drug, with the R/S verdict as a separate genome-level overlay.
  - Calibrated rules are IN-SAMPLE/opt-in (organism-gated); the default `DRUG_RULE` path remains — the map should surface which path produced a phenotype call (provenance), consistent with the project's honesty rails.

### Captured by: design @ 2026-06-17 (re-run, folding probe findings)
- Files read: features/genome-map/feature.md, dna_decode/data/annotations.py, scripts/pathotype_laptop_pipeline.py, tools/docker_runner.py (consulted the probe @ 2026-06-17 grounding below)
- Key claims:
  - Folded the probe tightenings into design: v1 scope = "Bakta honesty report" (Bakta+AMRFinder+existing parser, NO hmmer/eggNOG); 4 tiers (pathway DEFERRED); R2 unknown rate is DB-LABELLED (`unknown_under_bakta_db_light`); G1 tightened to PREVENT-WRONG-INFERENCE (>=3 features, >=1 not "unknown rate exists"); Bakta-smoke-first; framed as a function/QC map (no-determinant is modal).
  - Reuse confirmed: `pathotype_laptop_pipeline.bakta_annotate` (Bakta runner) + `annotations.parse_gff3` (GFF parser) → the map assembles existing pieces.
  - The AC section is now STALE vs this updated design — a `/feature-design spec genome-map` re-run is required to regenerate AC1-AC12 (5-tier→4-tier, DB-labelled unknown rate, tightened G1) before spec-ready.

### Captured by: probe @ 2026-06-17
- Files read: dna_decode/data/annotations.py, dna_decode/eval/amr_rules.py, dna_decode/organism_rules/tb_amr.py, dna_decode/data/tb_who_catalogue.py, tools/docker_runner.py, wiki/stage2_install_artifact_2026-05-15.md, CLAUDE.md
- Key claims:
  - REUSE: `annotations.py::parse_gff3` (:63) is an existing two-pass Bakta-style GFF parser — the annotation-ingestion layer is already built; the map reuses it.
  - INSTALLED: Bakta (db-LIGHT 4GB) + AMRFinderPlus + BLAST+ via tools/docker_runner. NOT installed: hmmer/Pfam, eggNOG, diamond (the dedicated homology-middle-layer stack).
  - [grounded] Bakta ANNOTATION was never smoke-tested on this host (stage2_install_artifact:59 — deferred, CPU-heavy) — the spike's input is unproven-runnable; a Bakta annotation smoke is the first concrete action.
  - [grounded] db-light makes the headline unknown-rate a TOOLING-coverage metric → must be labelled `unknown_under_bakta_db_light` (DB+version in-field).
  - [grounded] G1 (AC9) is gameable by relabelling existing Bakta/determinant output → tighten to "prevent a concrete wrong inference", >=1 not "unknown rate exists".
  - [grounded/inferred] phenotype tier is determinant-only → no-determinant is the MODAL case for arbitrary genomes → reframe v1 as a function/QC map with determinant call-outs, not a phenotype map.

## Changelog
- 2026-06-17 — record created (cascaded from intake; no template present, built from canonical sections)
- 2026-06-17 — intake completed (Summary/Problem/Desired Outcome/Repo Fit/Scope/Actors populated; primary intent ratified = personal exploration/QC tool)
- 2026-06-17 — design completed (Workflow/Behavior + Rules And Decisions [5-tier schema + phenotype wall + the UX/eval gate G1 differentiation-over-Bakta + G2 no-tier-confusion] + Edge Cases; Technical Notes enriched)
- 2026-06-17 — spec completed (Acceptance Criteria AC1-AC12 grouped by outcome [map / honesty rails / overlay / the gate / robustness] + Test Scenarios [happy/edge/failure/precedence/regression/gate] + spec constraints [3 prototype genomes, JSON-as-contract, free-tool dependency]; status spec-ready)
- 2026-06-17 — probe completed (feasibility GREEN-with-caveats: Bakta/AMRFinder/BLAST + existing parser installed, hmmer/eggNOG not; 3 design tightenings surfaced) [grounding appended]
- 2026-06-17 — design RE-RUN (folded probe: v1 = Bakta-honesty-report, 4 tiers [pathway deferred], DB-labelled unknown rate, G1 = prevent-wrong-inference, Bakta-smoke-first, function/QC framing; status spec-ready -> design; AC now stale -> re-spec needed)
