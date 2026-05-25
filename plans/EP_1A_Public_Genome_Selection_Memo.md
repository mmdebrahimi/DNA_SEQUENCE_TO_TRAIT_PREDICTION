# EP-1A Public Held-Out Genome Selection — Memo

> Selection criteria + candidate sources for the public E. coli genome that fires the genome-input contract test (`plans/Post_V0_EP_Ladder_Plan.md` EP-1A). One genome required; selection decision lives on the Precision 7780 side (Codex has DB access + GPU to validate).

**Status:** DRAFT 2026-05-25. Codex picks the final genome on Precision 7780 once EP-0 close lands.

---

## Why this memo exists

EP-1A's terminal claim: "E. coli cipro genome-input predictor produces v0 JSON output on a public held-out genome." Without pre-staged selection criteria, Codex would re-invent the picking rationale at execution time. This memo pre-stages the rubric so picking takes minutes, not hours.

## Hard selection criteria (all must hold)

1. **Public + freely redistributable.** No clinical-isolate restrictions, no institutional embargo. NCBI / PATRIC / NARMS / ATCC / GenomeTrakr sources are fine.
2. **Not in the N=147 training cohort.** Cross-reference against `data/processed/stage2_n150_cipro_cohort.parquet` (column: `assembly_accession`). The genome's accession + any related strain registrations must NOT appear.
3. **Has known cipro AST phenotype.** Either a public AST result OR a literature-cited resistance/susceptibility class. The phenotype label is the ground truth for the prediction comparison.
4. **E. coli, not a related species.** v0 model was trained on E. coli; testing on Shigella or Salmonella would conflate signal with organism shift.
5. **Genome + GFF3 annotation both available.** GFF3 is required for attribution; without it, top-K attribution would be locus-tag-only.
6. **Reasonable assembly quality.** Contig count ≤ 200 + N50 ≥ 50 kb. Matches the assembly-quality threshold v0 was trained on.

## Soft preferences (pick higher-priority where possible)

7. **Known mechanism for cipro-R cases.** If a R strain is picked: literature + AMRFinderPlus should agree on the mechanism (QRDR point mutation, qnr gene, etc.). Avoids "ground-truth label uncertain" diagnosis if the prediction misses.
8. **Recent deposition.** Genomes deposited 2023+ are less likely to overlap with the N=147 training cohort (which used pre-2023 BV-BRC snapshot).
9. **Reference-quality strain when possible.** ATCC reference strains have well-characterized phenotypes + are commonly used as positive/negative controls. ATCC 25922 (cipro-S, gyrA wild-type) is the canonical S reference; ATCC 35218 is the canonical β-lactamase-positive reference (relevant for EP-2 cef, but cipro-S in most reports).
10. **Mechanism class diversity** is NOT required for EP-1A (only 1 genome). It's required for EP-1B (see separate memo).

## Candidate sources (Codex picks the actual genome from these)

| Source | Pros | Cons | URL |
|---|---|---|---|
| ATCC reference strains | Well-characterized; canonical references | Limited cipro-R diversity | https://www.atcc.org/ |
| NCBI Pathogen Detection | Large E. coli set with AST metadata + AMRFinderPlus calls; freely queryable | Volume — need to filter | https://www.ncbi.nlm.nih.gov/pathogens/ |
| NARMS (US) | Public AST data; per-isolate genome + phenotype linked | US-only; mostly enteric isolates (food/animal sources) | https://www.cdc.gov/narms/ |
| PATRIC / BV-BRC | Large + annotated; already used for N=147 (CAREFUL with overlap) | High overlap risk with training cohort | https://www.bv-brc.org/ |
| GenomeTrakr (FDA) | Foodborne-pathogen focus; public | Limited clinical cipro-R diversity | https://www.fda.gov/food/whole-genome-sequencing-wgs-program/genometrakr-network |
| ResFinderFG / ResFinder reference panel | Curated reference panel for AMR-gene comparison | May overlap with N=147; older entries | https://cge.food.dtu.dk/services/ResFinder/ |
| EuSCAPE / EuropeanCanadian AMR | European AMR surveillance | Access may require registration | https://www.ecdc.europa.eu/ |

## Suggested first pick (placeholder; Codex validates)

**Recommended: ATCC 25922 OR a recent (2023+) NCBI Pathogen Detection isolate with known cipro AST.**

- **ATCC 25922 (Pros):** canonical cipro-S reference; gyrA wild-type; complete published genome + annotation; almost certainly NOT in N=147 cohort (cohort was BV-BRC clinical isolates, not ATCC type strains); tests the S prediction path.
- **ATCC 25922 (Cons):** only 1 mechanism class (wild-type); doesn't test the R prediction path; doesn't stress the model on novel patterns.

If Codex prefers a R-strain test: pick a recent NCBI Pathogen Detection isolate (2024-2025 deposition) with documented gyrA-S83L OR qnrS1 + cipro-R AST label. Source the AMRFinderPlus prediction from NCBI's per-isolate AMR call; that's the ground-truth proxy.

**Suggested approach:** test BOTH a S strain (ATCC 25922 type) AND an R strain (1 NARMS / NCBI candidate). 2 genomes is still scoped enough for "EP-1A thin slice" + dramatically reduces the "we tested only the easy case" risk.

## Selection rubric (Codex executes)

When picking from the candidate set:

```
For each candidate genome:
  1. Does it satisfy ALL hard criteria 1-6? If no, skip.
  2. How many soft preferences 7-10 does it satisfy?
  3. Score = number of soft preferences satisfied (max 4 per genome).

Pick the highest-scoring genome. Tie-break by: prefer recent deposition;
prefer reference-quality strain; prefer mechanism class not already in cohort.
```

## EP-1A success criterion (recap from `Post_V0_EP_Ladder_Plan.md`)

`scripts/pipeline.py predict --genome-fasta <selected.fna> --annotations <selected.gff3> --drug ciprofloxacin --cache <transient> --model-path <retrained.pkl> --audit-merge-json <merge.json> --output <out.json>` runs end-to-end and emits v0 JSON+MD with prediction + calibrated_probability + confidence_tier + attribution_scope_confidence + top_k_attribution + audit_verdict + provenance. **Same-strain parity test:** when run with the FASTA of an N=147 cohort strain, output matches cached-strain path within ε = 0.01 calibrated probability.

## What this memo deliberately does NOT cover

- The actual accession number(s) for the test genome(s). Codex picks on Precision 7780 with real DB access.
- The transient-cache implementation for novel-genome embeddings (lives in EP-1A code work, not this memo).
- Multi-drug genome selection (cef / tet would need their own panels; EP-2 territory).
- The 10-genome external-benchmark panel (EP-1B; see separate memo).

## Open questions for Codex / user

1. Single S strain (ATCC 25922) for EP-1A, OR 1 S + 1 R? The latter is scoped enough but doubles wall-clock.
2. If R strain is picked: which mechanism class? QRDR (matches v0's strongest signal) OR plasmid (qnr) OR something else?
3. Should the EP-1A test also include a Shigella or Salmonella negative-control (assert the model REFUSES to predict, or that confidence falls below threshold)? Out-of-scope for v0 spec, but cheap to add.
