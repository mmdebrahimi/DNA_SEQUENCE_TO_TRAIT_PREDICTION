# dna_decode — project frontier map (2026-06-19)

A step-back, honest current-state map: what's shipped, what's banked, what's closed-with-a-negative,
and the only genuinely-open moves. Written to reconcile drift (the eukaryotic embedding arm was
mis-described as "pending GPU work" in two places — it is a CLOSED NEGATIVE) and to frame the next
decision. North star (unchanged): **an AI DNA decoder TOOL — "DNA → what its parts do" — not papers.**

## Shipped + banked (terminal honest products)

| product | what it is | status |
|---|---|---|
| **Deterministic AMR decoder** | `call_resistance(organism, drug)` over AMRFinder curated determinants; 10 provenance-disjoint SCORED cells, lineage-disclosed, frozen (`wiki/reproducibility_freeze_2026-06-13.md`) | BANKED |
| **Genome-map v1** ("Bakta honesty report") | point at ONE microbial genome → evidence-tiered per-feature map (4 tiers, phenotype behind a validated-determinant wall, DB-labelled unknown rate) | SHIPPED + hardened (C1/C2/M1, verdict-derived drug labels) 2026-06-19 |
| **Fungal AMR cell** | C. auris ERG11 deterministic target-site decoder | G1 reached (LABEL_LIMITED) |
| **Pathotype resolver (EP-4)** | VirulenceFinder-style ExPEC deterministic resolver | SHIPPED |
| **TMP-SMX overlay** | non-frozen AND-across-gene-families experimental cell | SHIPPED (branded EXPERIMENTAL) |
| **TB AMR cell** | RIF + INH on CRyPTIC via the WHO catalogue | BLOCKED-gated pending data runs (D:-/regeno-gated) |
| **External re-validation arm** | Oxford measured-MIC cohort re-validation of the frozen decoder | SHIPPED (live run = manual) |
| **Negative-results map** | the 8 reusable rejection GATES + the verified failure record | SHIPPED (`wiki/negative_results_map_2026-06-13.md`) |

## Closed with a recorded negative — DO NOT reopen

- **Learned/embedding decoder on free public labels** — 3 de-confounded failures across the kingdom boundary: (1) cipro within-lineage = chance (learned lineage, not mechanism); (2) pathotype (sampling-defined label); (3) **Arabidopsis flowering-time (PlantCaduceus, real GPU, 2026-06-12) — embedding within-group r² negative, learned population structure not the causal signal. H2 FALSIFIED.** Scaling to a bigger/paid GPU does NOT help (signal-vs-structure, not window-budget).
- **Public-label AMR expansion** — bounded by LABELS not models; screened by the 8-gate map. Banked 2026-06-13.
- **MIC-continuous head** — infeasible on censored public MICs.

## The binding constraint (the project's scientific finding)

**Honest public bacterial genotype→phenotype decoding is bounded by LABELS, not models.** Every closed
track tripped one of the 8 rejection gates (circular-label / study==class / sampling-defined / surveillance-
domination / assembly-attrition / MIC-censoring / provenance-not-separable / dedup-collapses-balance). The
deterministic decoder + this boundary statement ARE the contribution.

## The only genuinely-open moves (everything else is closed or banked)

Both require a USER decision — neither is a pure-executor task, and neither should be cold-built:

1. **ACQUISITION of a non-public label source** (wet-lab / clinical / partnership / a genuinely-independent
   un-mined public supplement). This is the ONLY path that clears the label gates *by construction* and
   reopens the learned-decoder arm honestly. → a USER sourcing decision; draft anchor
   `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`.
2. **PROSPECTIVE-LOCK** — pre-register the frozen decoder's predictions on a not-yet-labelled cohort, score
   when labels arrive (leakage-impossible by construction). Buildable D:-free, BUT: it leans *rigor/paper*,
   not *tool capability* — weigh against the "tool not papers" north star before investing. Needs the
   planning chain (a new validation subsystem) + a target prospective data stream (itself a user choice).

### Tool-capability extension candidates (warrant the planning chain, NOT cold-build)

- **Genome-map virulence/pathotype overlay tier** — extend the honesty map beyond AMR using the existing
  `dna_decode/pathotype/` resolver. *Tool*-aligned ("what its parts do" for virulence). **Caveat that needs
  design:** the k-mer detector returns coverage, NOT coordinates → it could only symbol-fallback-join (the
  gene-symbol trap the genome-map brainstorm just caught) → it would NO-GO by the genome-map's own gate
  unless built on the BLAST `vf_runner` coords (a tooling dependency like AMRFinder). → run the full chain
  (`/idea-anchor` → `/probe` → … → `/brainstorm`) before building; do NOT cold-build (the integrity crux is
  exactly what the chain protects).
- **Genome-map richer homology tiers** (hmmer/Pfam/eggNOG) / a visual browser — explicitly deferred in v1;
  user investment call.

## Recommendation

The project is at a **terminal honest state on free public data** — the diminishing-returns plateau is real
(further AMR/genome-map increments re-confirm prior findings rather than producing new signal). The honest
move is to **stop manufacturing increments and make the strategic call**:

- If the goal is a MORE CAPABLE TOOL → the genome-map virulence-overlay tier is the strongest candidate
  (run the planning chain first — the coordinate-join contract is load-bearing).
- If the goal is to UNLOCK the learned arm → it requires label ACQUISITION (a sourcing decision only the
  user can make).
- Prospective-lock is available but rigor-flavored; defer unless paper-grade external validation becomes a goal.

Nothing here is blocked by the missing D: drive except the TB regeno run + a live genome-map reconfirmation.
