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
Actor: the developer, one genome at a time. Happy path:
1. **Input** — one microbial genome: a FASTA (+ an optional precomputed GFF/annotation). Single genome; no
   batch, no comparison in v1.
2. **Locate features** — structural annotation (genes/CDS/RNA/operons) via a free external tool
   (Bakta/Prokka) OR accept a provided GFF. (Tool feasibility is a `/probe` concern, not design.)
3. **Tier each feature** — assign every feature its single HIGHEST-confidence evidence tier (table in
   Rules), retaining lower-tier evidence as secondary. Sources: the curated determinant catalogs
   (phenotype tier), curated-DB / homology hits (function/pathway/homology tiers), nothing (unknown).
4. **Phenotype overlay** — features that are validated determinants get a phenotype annotation (the
   drug/property they contribute to) from the EXISTING cells; the genome-level R/S verdict (the existing
   decoder's `call_resistance` / `tb_amr` output) is shown SEPARATELY and clearly tier-tagged. The map
   NEVER promotes a homology/pathway hit to a phenotype claim.
5. **Output** — a JSON-first per-feature map + a flat feature table. Headline metrics: per-tier feature
   counts + the **unknown rate** + the phenotype-tier feature list.

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

| tier | fires when | may emit a PHENOTYPE claim? |
|---|---|---|
| `determinant-phenotype` | a validated determinant cell matches this feature (AMRFinder `main.tsv` row + `calibrated_amr_rules.json` class→drug; or `tb_amr`/organism determinant) | **YES — the only tier that may** |
| `curated-molecular-function` | a curated-DB hit names a gene/enzyme (the annotation tool's confident call) | no |
| `pathway-module` | KEGG/MetaCyc module membership inferred | no |
| `homology-only-hypothesis` | Pfam/eggNOG/orthology domain hit, no curated identity | no |
| `unknown` | no confident annotation | no |

**Load-bearing rules:**
- **Phenotype wall (R1):** a phenotype/property claim attaches ONLY at `determinant-phenotype`. No other
  tier may express a phenotype — they are molecular-function / hypothesis / unknown ONLY.
- **Unknown visibility (R2):** the unknown rate is a REPORTED headline metric; unknowns are never hidden,
  collapsed, or back-filled with a weaker tier's guess.
- **Highest-confidence primary (R3):** a feature with multiple evidences takes the highest tier as primary
  and keeps the rest as secondary evidence (transparency, not suppression).
- **No frozen edit (R4):** the phenotype tier READS the existing determinant cells/catalogs; it never
  modifies `dna_decode/eval/amr_rules.py` or `calibrated_amr_rules.json`.

**UX / EVAL GATE (the make-or-break — the feature's go/no-go bar).** The prototype passes iff, on the 2-3
test genomes, BOTH hold:
- **G1 Differentiation-over-raw-Bakta:** for ≥1 genome there are concrete features where the tiered map
  gives a MORE HONEST/useful read than the raw annotation TSV — e.g. a homology guess that raw-Bakta
  states as fact but the map demotes to `homology-only-hypothesis`; or a determinant the map surfaces as
  `determinant-phenotype` that raw-Bakta lists as a plain gene; or the reported unknown rate that raw
  output never surfaces. If the map adds NOTHING over reading Bakta's TSV, the feature FAILS (it would be
  catalog-stacking busywork).
- **G2 No-tier-confusion:** a reader can distinguish phenotype vs molecular-function vs hypothesis vs
  unknown at a glance; in a spot-check, zero features are mistakable as phenotype that are not
  `determinant-phenotype`.

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
  tiers + unknowns). No phenotype is fabricated.

## Acceptance Criteria
_(deferred to `spec`)_

## Test Scenarios
_(deferred to `spec`)_

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

## Changelog
- 2026-06-17 — record created (cascaded from intake; no template present, built from canonical sections)
- 2026-06-17 — intake completed (Summary/Problem/Desired Outcome/Repo Fit/Scope/Actors populated; primary intent ratified = personal exploration/QC tool)
- 2026-06-17 — design completed (Workflow/Behavior + Rules And Decisions [5-tier schema + phenotype wall + the UX/eval gate G1 differentiation-over-Bakta + G2 no-tier-confusion] + Edge Cases; Technical Notes enriched)
