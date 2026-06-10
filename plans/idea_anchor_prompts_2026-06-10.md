# Idea-anchor prompts — the 3 deliberate-project recommendations (2026-06-10)

> Drafted by Soraya (conversational-executor) for the user to run through `/idea-anchor` + ratify. These
> are the "genuinely-new, deliberate-project-shaped" moves surfaced when the 4th-kingdom (viral) antiviral
> decoder shipped (commit e3b5711) and the catalog-pattern plateau was called. Each is a SEED problem
> statement — `/idea-anchor` crystallizes it into an anchored goal before `/probe` → `/feature-design` →
> `/technical-plan`. Ratify (confirm/redirect) the framing; don't author from scratch.
>
> Ranked by VOI: #1 closes a blind spot · #2 tests an architectural question · #3 is the highest-rigor
> de-confounder (the only one that could FALSIFY a shipped decoder). All three honor: no money, laptop +
> Docker + free public data, determinism-first, every claim carries a falsifier.

---

## Anchor 1 — `expression_context` signal class (close the EXPRESSION abstain floor)

**One-line idea:** Layer a deterministic *expression-context* feature onto the presence-based AMR decoder
— read IS-element / promoter insertions immediately upstream of a target gene from the SAME assembly the
decoder already has — to convert correct `ABSTAIN → R` on the carbapenem EXPRESSION floor, and generalize
the pattern to other intrinsic-gene overexpression.

**Relaxes binding assumption #1** ("the decoder predicts binary R/S from PRESENCE of curated determinants"
→ it could also read *regulatory context* from the same bytes).

**Grounded state (do not re-derive):** The 2026-06-10 ISAba1→blaOXA-51 falsifier
(`wiki/expression_frontier_isaba1_falsifier_2026-06-10.md`) proved the signal is REAL and *perfectly
specific* on N=30 Acinetobacter (junction-positive R 1/15, S 0/15 — zero false positives) but *partial*
(recovers 1/5 intrinsic-only-R). The ISAba1 ref was a 570 bp partial; the cohort was a single N=30. The
EXPRESSION floor is the decoder's documented #1 limitation (the dominant carbapenem false-negative cause).
This is the ONLY move on the table that *closes a blind spot* rather than adding reach.

**Success / falsifier:** On an INDEPENDENT Acinetobacter cohort (disjoint from the cached 30, label-balanced),
the ISAba1-upstream-of-OXA-51 feature recovers a meaningful fraction of intrinsic-only-R with **≤1 S
false-positive**, AND beats the naive ABSTAIN baseline (which recovers 0). **KILL** if it adds S
false-positives at scale, or if the recovery rate on the independent cohort collapses to ~0 (signal was a
single-cohort artifact). Wiring into the deployed caller (junction → R override) is gated on this passing —
validate-wrapper discipline.

**Hard constraints / first scoping questions for /idea-anchor:** no money; laptop + Docker (blastn on
cached assemblies). Needs (a) a FULL-LENGTH ISAba1 reference (not the 570 bp partial) — and is a fuller ref
even necessary, given the 4 non-rescued FN strains have ISAba1 *not near* OXA? (b) an independent
Acinetobacter cohort with carbapenem AST + downloadable assemblies (the BV-BRC `assembly_accession`
bottleneck applies). (c) generalization target: does ISAba1/ISAba125-upstream-of-ampC/ADC + IS-upstream-of-
efflux-regulators extend the same feature, or is each a separate curation?

**`/idea-anchor` command:**
```
/idea-anchor expression_context signal class: a deterministic IS-element/promoter-upstream-of-target feature layered on the presence-based AMR decoder, to cross the carbapenem EXPRESSION abstain floor (ABSTAIN->R) using the same assembly; validated on an independent Acinetobacter cohort before wiring. Full context: plans/idea_anchor_prompts_2026-06-10.md Anchor 1.
```

---

## Anchor 2 — pfmdr1 directional catalog (opposite-direction selection in the catalog architecture)

**One-line idea:** Extend the antimalarial decoder with **pfmdr1** markers, encoding the DIRECTIONAL
subtlety the architecture has never faced: N86Y / D1246Y *increase* 4-aminoquinoline (chloroquine /
amodiaquine) resistance but the SAME alleles *increase susceptibility* to lumefantrine / mefloquine — so
the catalog must be per-drug-DIRECTION, not a flat marker set. pfmdr1 copy-number (the mefloquine /
lumefantrine driver) is a CNV, invisible to a point-mutation caller — surface it as a named blind spot, do
not pretend to call it.

**Relaxes binding assumption #5** ("the kingdom jump = single-locus unidirectional target-site mechanisms"
→ can the per-drug catalog represent OPPOSITE-direction selection on one allele?).

**Grounded state:** The antimalarial vertical currently ships artemisinin/K13 + chloroquine/pfcrt
(observed + genome, real-validated). The existing `ANTIMALARIAL_RESISTANCE_MUTATIONS` is `drug → gene →
set-of-R-substitutions` — every catalog so far (bacterial / fungal / protozoan / viral) is UNIDIRECTIONAL
(a marker = resistance, full stop). pfmdr1 N86Y is the canonical case where one allele points opposite ways
for two drugs in the same ACT — the antimalarial_amr.py comment already flags it as "deliberately NOT
catalogued" for exactly this reason. This is an ARCHITECTURE test (can the schema encode direction without
a special-case hack?), not just more markers.

**Success / falsifier:** A schema that represents per-(drug,allele) direction calls real pfmdr1 alleles
correctly: N86Y → amodiaquine-contributing R, AND the SAME N86Y → NOT-R (or susceptibility-flag) for
lumefantrine, AND pfmdr1-CNV surfaced as an `undetectable_mechanism`. **KILL** if the only way to represent
direction is a per-drug special-case branch that doesn't generalize (then pfmdr1 is out of scope and the
unidirectional architecture is confirmed as the boundary), or if real directional-phenotype alleles can't
be sourced free.

**Hard constraints / first scoping questions:** no money; pfmdr1 CDS reference + real alleles fetchable
like pfcrt (NCBI). CNV explicitly OUT of scope (point-caller blind spot — surface, don't solve). First
question for /idea-anchor: is the directional-catalog a generic schema upgrade (every kingdom benefits) or
a pfmdr1-only wart? And: is there a free, machine-readable directional-phenotype label to validate against,
or does this stay an observed-substitution + literature-direction tool (no independent label)?

**`/idea-anchor` command:**
```
/idea-anchor pfmdr1 directional catalog: extend the antimalarial decoder to encode per-(drug,allele) DIRECTION (N86Y raises chloroquine/amodiaquine R but lowers lumefantrine/mefloquine R) — an architecture test of whether the per-drug catalog can represent opposite-direction selection without a special-case hack; pfmdr1 CNV stays a named blind spot. Full context: plans/idea_anchor_prompts_2026-06-10.md Anchor 2.
```

---

## Anchor 3 — independent phenotype-label validation (the validate-wrapper de-confounder)

**One-line idea:** Validate at least one shipped deterministic decoder against an INDEPENDENT,
lab-measured phenotype label (e.g. CDC/GISAID NI-assay IC50 fold-change for influenza NA H275Y; NARMS MIC
for a bacterial drug NOT from BV-BRC) rather than against the curated-DB-derived call — the "beat the naive
tool on INDEPENDENT data" discipline the entire product rests on.

**Relaxes binding assumption #6** ("labels come from NCBI / BV-BRC" → an independent label source could
de-confound or falsify a shipped decoder).

**Grounded state:** Every kingdom's validation so far uses the marker's PRESENCE as both feature and (near-)
label, or BV-BRC MIC — the SAME sampling context. The embedding-niche three-part test + the
validate-wrapper memory both say a rule earns trust ONLY when it beats naive tool use on INDEPENDENT labels;
in-cohort accuracy proves the tool works, not that the layer adds value. This is the highest-RIGOR move and
the only one that could legitimately FALSIFY a shipped decoder (high VOI precisely because it can lose).

**Success / falsifier:** Pull a truly independent, free, machine-readable phenotype set for ONE drug and
measure the deterministic rule's sens/spec against it. **The finding is the value either way:** if the rule
holds on independent labels → trust confirmed; if independent-label accuracy is materially below in-cohort
accuracy → the in-cohort number was inflated by the shared-label confound (a real, publishable-internally
correction). **KILL the anchor** only if NO drug has a free independent label source (then document that
the product is inherently in-cohort-validated and say so honestly in every decoder's caveat).

**Hard constraints / first scoping questions:** no money — the HARD part is FINDING a free, machine-readable,
genuinely-independent phenotype source per drug (many are paywalled or PDF-only). The anchor's FIRST job is
a source census: which shipped drug (cipro / cef / fluconazole / oseltamivir / chloroquine) has an
accessible independent label set (NARMS public MIC, CDC flu NI-assay surveillance, WWARN, EUCAST), and is it
machine-readable? Start with the cheapest-to-source, not the most interesting.

**`/idea-anchor` command:**
```
/idea-anchor independent phenotype-label validation: score at least one shipped deterministic decoder against an INDEPENDENT lab-measured phenotype (CDC/GISAID NI-assay IC50 for flu NA, or NARMS MIC for a non-BV-BRC bacterial drug) instead of the curated-DB-derived label — the validate-wrapper de-confounder; the first job is a free-source census per drug. Full context: plans/idea_anchor_prompts_2026-06-10.md Anchor 3.
```

---

## Pipeline note

Per `planning-pipeline.md`: each is class (d) cross-cutting (#1, touches the deployed caller) or class (e)
research-first (#2, #3). After `/idea-anchor`, the chain is `/probe` → `/feature-design` →
`/technical-plan` → pre-exec `/brainstorm` → `/save-plan` → `/execute-plan`. Each step STOPS — none
auto-invokes the next. Anchor 1 is the recommended first pull (closes a blind spot; substrate already on
disk). Anchor 3 is the highest-rigor but gated on a source census that may dead-end (no free label → honest
documentation outcome).
