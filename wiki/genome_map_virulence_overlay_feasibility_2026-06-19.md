# Genome-map virulence-overlay tier — feasibility / probe-grade grounding (2026-06-19)

**Decision input for the planning chain.** The user picked "more capable tool → a genome-map
virulence-overlay tier" (extend the honesty map beyond AMR). Per the discipline that caught the AMR
coordinate-join trap, this needs the planning chain (`/idea-anchor` → `/probe` → `/feature-design` →
`/technical-plan` → pre-exec `/brainstorm` → `/save-plan` → `/execute-plan`) BEFORE building. This memo is
the executor-produced grounding so that chain is fast + the load-bearing risk is already resolved.

## The concept

A 5th genome-map overlay: point at ONE genome → besides the AMR `determinant-phenotype` tier, surface a
`virulence-determinant` tier where a curated VirulenceFinder allele is present, behind the SAME validated-
determinant wall + coordinate-join integrity gate. It is an ANNOTATION overlay ("a curated virulence
determinant is present at this locus"), NOT a learned virulence-phenotype predictor — identical honesty
posture to the AMR-determinant tier.

## Load-bearing feasibility finding: GREEN, with ONE required change

The whole feature hinges on: *can a virulence determinant be COORDINATE-joined to a Bakta feature, or is it
symbol-fallback-only (→ auto-NO-GO by the genome-map's own gate)?* Answer:

- **[grounded] The k-mer detector path is NOT viable.** `dna_decode/pathotype/detect.py` returns k-mer-seed
  `percent_coverage` (a presence proxy), NO coordinates. A tier built on it would be ALL symbol-fallback →
  NO-GO. Rule it out.
- **[grounded] The canonical blastn path IS viable — with one small change.** `vf_runner.run_canonical_vf`
  already runs real `blastn` of the VF allele DB vs the assembly, but its outfmt is `6 qseqid pident length
  qlen` — identity + coverage, **no subject coordinates**. blastn emits subject coords for free: extend the
  outfmt to `6 qseqid sseqid sstart send pident length qlen`. The subject DB is `makeblastdb` over the INPUT
  assembly, so `sseqid` = the original FASTA contig name + `sstart/send` = the locus on it.
- **[grounded] The genome-map coord-join is reusable verbatim.** Each VF coord-hit maps directly onto the
  existing `DeterminantHit` (symbol=allele/gene, cls="VIRULENCE"/cluster, method="blastn", contig=sseqid,
  start=sstart, stop=send). It then feeds the SAME `build_contig_name_map` (AMRFinder/blastn use ORIGINAL
  NCBI contig names; Bakta renames → reconcile by length — IDENTICAL to the AMR path) + `join_hits`
  (protein-id → coord-overlap → symbol-fallback; whole-contig `region` exclusion + smallest-span tie-break).
- **[grounded] D:-FREE on this host.** `blastn` resolves natively (`C:/Users/Farshad/ncbi-blast/bin`) +
  the VF allele DB is committed (`data/virulencefinder_db/virulence_ecoli.fsa`, 23 clusters: STX1/STX2/
  LEE/LT/ST/hemolysin/siderophores/…). A live single-genome run needs only an assembly FASTA. Offline-safe:
  no blastn → the overlay degrades to `unavailable` (mirrors the AMR `--allow-degraded` posture).

## Required changes (the technical-plan skeleton)

1. `vf_runner.run_canonical_vf` — extend outfmt to carry `sseqid sstart send`; return per-hit coords
   (additive; existing per_gene/per_cluster unchanged; offline-degrade unchanged).
2. New `dna_decode/genome_map/virulence_overlay.py` — a VF-blast-output → `DeterminantHit`-shaped adapter +
   `join_hits` reuse; emit join-quality counts; a genome whose VF joins are ALL symbol-fallback → NO-GO
   (same guard).
3. `dna_decode/genome_map/__init__.py` — a `virulence-determinant` tier constant; `build_map` assigns it
   when a high-confidence VF coord-join hits a feature (the virulence wall; symbol-fallback excluded).
4. The genome-level overlay = the pathotype RESOLVER call (`pathotype.resolve_call` → EPEC/EHEC/ETEC/UPEC/
   EAEC/ExPEC/commensal), shown SEPARATELY + tier-tagged — the virulence analog of the AMR genome-level R/S
   call. (Per-feature = "this allele is present"; genome-level = "the pathotype call".)
5. CLI + spike wiring + tests (synthetic VF-blast fixtures; offline-degrade test green without blastn).

## Design questions for `/feature-design` — Soraya's DRAFTED answers (ratify or redirect)

- **Q1 — own tier or fold into `determinant-phenotype`?** DRAFT: its OWN `virulence-determinant` tier (the
  phenotype wall is AMR-resistance-specific; virulence is a different claim class). Risk of folding: blurs
  "resistance" with "virulence."
- **Q2 — per-feature claim wording?** DRAFT: `virulence_determinant_present` + cluster + the pathotype the
  cluster contributes to (provenance), NEVER "this strain is pathogenic" (the wall: presence, not phenotype).
  Mirrors the AMR `DETERMINANT_PRESENT` / `drug_rule_counted` discipline.
- **Q3 — organism scope?** DRAFT: E. coli-only in v1 (the committed DB is `virulence_ecoli.fsa`); a
  non-E.coli genome shows NO virulence overlay (honest absence, like cipro's organism scoping). Surface the
  scope explicitly.
- **Q4 — independence caveat?** DRAFT: carry `vf_runner.NON_INDEPENDENCE_CAVEAT` (the VF caller matches the
  VF DB; this is annotation, not independent validation) into every virulence section. Non-negotiable.
- **Q5 — the genome-level overlay = the pathotype resolver call?** DRAFT: yes — `resolve_call` is the
  virulence analog of `call_resistance`; show it tier-tagged + separate. (AUTHORITY-flavored: it widens the
  tool's claim surface to "pathotype" — surface for ratification, don't bury.)

## Drafted `/idea-anchor` sentence (paste-ready)

> An honest genome-map virulence-determinant overlay: point the existing genome-map at ONE E. coli genome →
> a `virulence-determinant` evidence tier where a curated VirulenceFinder allele is present, behind the same
> coordinate-join integrity gate + validated-determinant wall as the AMR tier (presence, never a learned
> pathogenicity claim), with the pathotype resolver call shown separately as the genome-level overlay;
> offline-degrades without blastn; E. coli-scoped in v1.

## Recommended planning chain (user-only — run in order; this memo grounds `/probe`)

```
/idea-anchor <the drafted sentence above>
/probe genome-map virulence-overlay     # this memo IS the probe evidence; confirms feasibility GREEN
/feature-design intake genome-map-virulence   # ratify Q1–Q5 above
/feature-design design genome-map-virulence
/feature-design spec genome-map-virulence
/technical-plan genome-map-virulence     # skeleton in "Required changes" above
/brainstorm                              # pre-exec; target the coord-join + the wall (the integrity crux)
/save-plan
/execute-plan                            # then Soraya can drive it D:-free (blastn + DB present here)
```

**Net:** the major risk (coordinate-join feasibility) is RESOLVED (GREEN via the blastn-outfmt change); the
remaining work is genuine feature-design (the 5 Q's) + a small, well-scoped build. It is D:-free on this host.
