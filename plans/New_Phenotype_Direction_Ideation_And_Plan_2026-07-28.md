# New phenotype direction — ideation + decomposition + technical plan (2026-07-28)

**User goal (verbatim):** "modify a single genotype in any organism and see if we can predict what the
phenotype change will be (or vice versa)" — a NEW decoder target outside the AMR / variant-effect / organism trio.

**Status: PLAN ONLY — awaiting ratification (Planning STOP). Not executed.**

---

## 1. The reframe (what's actually new)

"Single genotype modification → predicted phenotype change" is **precisely the forward variant-effect cell's
paradigm** at the molecular level (one point mutation → molecular fitness, DMS-validated). So the genuinely-new
thing is NOT the capability — it's the **phenotype AXIS**. The question is: *which phenotype has a single-locus,
large-effect, free-labeled, cross-organism genotype→phenotype map that the project doesn't already decode?*

The project's regime framework is the screen:
- **R1 curated-catalog → deterministic wins** (AMR determinants, PGx, ClinVar, serotype). The project's strength.
- **R2 molecular-property → learned wins if fitness-aligned** (the forward DMS cell).
- **R3 organism-polygenic → CLOSED NEGATIVE** (embeddings 0-for-5; the wall is population-structure-in-the-label).

A viable NEW target must be **single-locus / large-effect** (dodges R3) with **free, independently-measured,
provenance-separable labels** (dodges the label wall), and **not already a cell**.

## 2. Ideation — candidates screened against the gates

| candidate | single-edit→phenotype? | free labels | regime | novel? | verdict |
|---|---|---|---|---|---|
| **Gene essentiality** (KO one gene → lethal/viable) | YES (KO = the edit; dead/alive = phenotype) | **YES, large, multi-organism** (DEG, OGEE, Keio E.coli, yeast deletion, DepMap human CRISPR) | R1 conserved-core catalogue + R2 learned tail | **YES** | ✅ **TOP PICK** |
| **Metabolic capability** (gene present → substrate utilization, e.g. lacZ→lactose) | YES | partial (KEGG/BioCyc curated + BIOLOG/phenotype-microarray measured) | R1 — the AMR "determinant→phenotype" paradigm, for metabolism | YES | ✅ strong runner-up |
| Human LoF constraint (gnomAD pLI/LOEUF) | gene-level | YES (gnomAD) | R1-ish, but a derived pop-gen metric (indirect, mild circularity) | YES | ⚠️ weaker (label is derived, not a measured phenotype) |
| Mendelian pathogenicity (ClinVar) | YES | YES | R1 curated | **NO — already `dna-clinvar`** | ❌ covered |
| Complex organism traits (height, yield, biofilm) | NO (polygenic) | — | **R3** | — | ❌ closed-negative regime |

## 3. Recommendation — the "knockout decoder" (gene essentiality)

**Predict: does knocking out gene X kill/impair organism O? (essential vs non-essential) — and vice versa
(given a lethal-KO phenotype, which genes are essential).** Why it's the best fit to the user's goal:

- **Cleanest instance of "single edit → phenotype change, any organism":** the edit is a gene KO, the phenotype
  is binary (viable/lethal), and essentiality is defined for bacteria, archaea, yeast, AND human cell lines.
- **Free, large, independently-measured labels** (KO screens — NOT derived from a sequence model → not circular;
  different organisms/labs → provenance-separable). Dodges the label wall the AMR track hit.
- **Perfect fit to the project's DUAL architecture:** essential genes are dominated by a **universal conserved
  core** (ribosomal, DNA replication, tRNA synthetases, etc.) → an **R1 curated-catalogue decoder wins on the
  core** (the project's strength), with a **learned R2 complement** (ESM/conservation features) for the
  organism-specific tail. This is the same deterministic-catalogue + learned-complement shape as AMR + forward.
- **High benefit:** essential genes = antibiotic/antifungal/anticancer drug targets; minimal-genome + synthetic-
  biology design; the "which genes can I delete" question.
- **Single-locus → dodges the R3 organism-polygenic wall** that closed the multimodal/organism track.

**Runner-up (metabolic capability)** is the smoothest paradigm-extension of AMR (a "metabolic determinant
catalogue") and is a good v0.2 if essentiality validates. Ratifiable alternative.

**Honest risks (H8, named up front):** (a) essentiality is CONDITION/ORGANISM-dependent (a gene essential in
one context isn't in another) — real label noise, the analogue of AMR MIC-label noise, handled by a per-organism
+ conditional split, not a global label. (b) sequence-alone essentiality prediction is modest; the DETERMINISTIC
conserved-core catalogue is expected to carry most of the signal (which is fine — it's the R1 win). These are
hypotheses to test in Family 4, not asserted results.

## 4. Decomposition — project families + critical path

**Family E1 — essentiality DATA layer** *(critical-path root)*
- Assemble free per-gene essentiality labels for ≥3 phylogenetically-spread organisms (proposal: *E. coli*
  K-12 [Keio/DEG], *S. cerevisiae* [yeast deletion collection], one human cell line [DepMap CRISPR]) + each
  gene's protein sequence + COG/ortholog annotation. Provenance-separable by construction (distinct sources).
- Falsifier/gate: data reachable + a clean binary (or graded) essential label for ≥1000 genes/organism.

**Family E2 — the R1 deterministic "conserved-core" decoder** *(depends E1)*
- Curated catalogue: universal-core essential gene set (COG/eggNOG universal single-copy + known essential
  categories) → predict essential by core-membership / cross-organism conservation. Deterministic, offline.
- Falsifier: on held-in organisms, core-membership AUROC ≥ a stated bar (derive-don't-assert; the null =
  base-rate of essentiality). This is the project-strength win to establish first.

**Family E3 — the R2 learned complement** *(depends E1; parallel to E2)*
- ESM2 / conservation / gene-property features → essential (the organism-specific tail the catalogue misses).
- Falsifier (the regime gate): does the learned model BEAT the deterministic conserved-core catalogue on the
  NON-core genes? (the forward-cell discipline: beat the domain-knowledge baseline, not just the null.)

**Family E4 — cross-organism VALIDATION + trust surface** *(depends E2+E3)*
- The "any organism" test: train/curate on 2 organisms, predict essentiality in a HELD-OUT organism (does the
  decoder transfer?). Clonality/phylo-independence disclosure like the AMR lineage layer. Emit a report card
  (mirror the forward/AMR report cards — per-organism honest tier, no aggregate headline).
- Falsifier: held-out-organism AUROC > base-rate null AND the honest tier is stated per organism.

**Critical path:** E1 → {E2, E3} → E4. E2 (deterministic) is the fast first win; E3/E4 test whether a learned
complement + cross-organism transfer earn their keep.

## 5. Technical plan (the --plan) — v0 "knockout decoder"

**MVP bar (checkable, frozen at run start):**
1. `file-exists` a data manifest with essential/non-essential labels + sequences for ≥3 organisms (E1).
2. `file-exists` a v0 result artifact `wiki/essentiality_decoder_v0_<date>.{md,json}` with, per organism:
   deterministic-core AUROC, learned-complement AUROC, base-rate null, and the held-out-organism transfer AUROC.
3. `test-exit-0` a unit-test suite pinning the catalogue lookup + the AUROC computation (offline, no network).

**Ordered steps:**
1. **E1 data (R3 real-surface first):** fetch OGEE / DEG (bacteria+yeast) + DepMap (human) — CONFIRM
   reachability + license before building (the first gate; if a source is walled, swap it, don't fake it).
   Extract {gene → essential?, protein_seq, COG/ortholog}. Commit a manifest (no licensed bulk data in-repo).
2. **E2 deterministic core decoder:** build the conserved-core essential catalogue (COG/eggNOG universal core
   + curated essential categories) → `predict_essential(gene) → {essential, non-essential, uncertain}` with an
   abstain tier. Offline, deterministic. AUROC vs base-rate null per held-in organism.
3. **E3 learned complement (Kaggle T4, reuse the proven ESM pipeline):** ESM2 embedding / conservation features
   → essential; the decisive test = beat the conserved-core catalogue on NON-core genes (regime gate). $0.
4. **E4 cross-organism transfer + report card:** held-out-organism AUROC; per-organism honest tier; emit the
   report card + a `build_essentiality_report_card.py` roll-up (mirror the forward/AMR pattern).
5. **Wire the CLI cell** (4-place contract: console entry + `cli.py` TRAITS + dispatch + `cell_registry`
   CellContract + the two guard tests) — `dna-decode essentiality` — only after E2 validates.

**Falsifiers / pre-committed verdicts (R2 — derive don't assert):**
- E2 FAIL if conserved-core AUROC ≤ base-rate null → the catalogue has no signal (unexpected; would refute the
  conserved-core premise).
- E3 "learned earns keep" only if it beats E2 on non-core genes (else deterministic-only ships, like AMR).
- E4 "cross-organism decoder" only if held-out AUROC > null; else honest "in-organism only" scope.

**Gates:** free data only (no money); no licensed bulk data committed; frozen AMR + forward surfaces
byte-unchanged; Kaggle T4 = $0. New cell is NON-frozen (like the forward/organism_rules lanes).

**Estimated shape:** E1+E2 = a few sessions (data + deterministic catalogue = the fast R1 win); E3 = one
Kaggle run; E4 = one validation pass + report card. First shippable value at E2 (deterministic core decoder).

---

## 6. Ratification asks (the STOP)

1. **Target:** gene essentiality ("knockout decoder") as primary — or the metabolic-capability runner-up? (I
   recommend essentiality: best fit to "single-edit→phenotype, any organism" + cleanest free data + dual-regime.)
2. **Organism set for v0:** E. coli + yeast + one human cell line (DepMap) — or a different spread?
3. **Scope of v0:** deterministic conserved-core decoder first (E1+E2), then decide on E3/E4 — or plan the full
   arc now?

No code written; this is the plan. On ratification I'll run E1 (data reachability first, R3) and build E2.
