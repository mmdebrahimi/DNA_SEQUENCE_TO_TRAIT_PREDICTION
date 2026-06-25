# Non-AMR GREEN-cell triage — Round 2 (2026-06-24)

Applies the project's **2-gate GREEN-cell test** (from `plans/Non_AMR_Phenotype_Pivot_Assessment_2026-06-24.md`)
to a fresh slate of candidate non-AMR traits, to pick the **next** cell to build after the Klebsiella
K-antigen (`dna-ktype`) cell shipped. Round 1 picked K-antigen; this round looks past it.

## The gate (unchanged)
```
trait -> [determinant catalog exists?] --no--> RED (emergent/polygenic; closed 0-for-4 negative; PARK)
                 |yes
                 v
         [FREE, INDEPENDENT, MEASURED, isolate-level label?] --no--> GREEN-CALLER (ship, tier=NO_FREE_SOURCE)
                 |yes                                                  (grows the saturation count; low value)
                 v
         GREEN-VALIDATED cell (ship + score vs the measured label + trust badge)  <-- the gold standard
```
**Circularity rail (load-bearing, from [[feedback_validate_wrapper_vs_underlying_tool]] + the VF-diff lesson):**
gate 2 requires a **wet-lab MEASURED** label (Quellung / Kauffmann-White slide agglutination), NOT another
in-silico tool's prediction. Validating a deterministic caller against SeroBA/SeqSero is faithful-to-tool
(in-distribution), NOT the independent GREEN-VALIDATED tier. Every serotype cell must headline this.

## Triage table

| # | Candidate trait | Gate 1: catalog | Gate 2: free MEASURED isolate-level label | Verdict |
|---|---|---|---|---|
| **1** | **S. pneumoniae capsular serotype** | ✓ cps-locus DBs — SeroBA / PneumoCaT / PneumoKITy / GPS-pipeline (curated, free, ~102/107 serotypes) | ✓✓ **GPS project: 11,810 genomes with phenotypic Quellung serotype** records (free, ENA-public). Known ceiling: GPS-pipeline 89.3% in-silico-vs-Quellung concordance | **GREEN-VALIDATED** — richest free measured label found; sibling of the shipped O:H serotype + K-antigen cells |
| **2** | **Salmonella enterica serovar** | ✓ antigen-formula DBs — SeqSero2 / SISTR (curated, free) | ✓ traditional Kauffmann-White serotyping in public validation cohorts (NARMS / PulseNet / the SeqSero2 + Frontiers eval sets) — **must filter to wet-lab-serotyped, not tool-predicted** | **GREEN-VALIDATED** (with circularity filter) — canonical Salmonella identity trait |
| 3 | H. influenzae capsular serotype (a–f) | ✓ hicap (cap locus, curated, free) | ~ slide-agglutination serotype in some public sets (smaller than pneumo) | GREEN-VALIDATED (secondary; thinner label) |
| 4 | E. coli K-antigen (capsule) | ~ partial (group 2/3 K via kps; less complete than Klebsiella KL) | ✗ sparse free measured K-typing | GREEN-CALLER (NO_FREE_SOURCE) — low value |
| 5 | Expanded virulence / toxin profile | ✓ VirulenceFinder clusters | ✗ presence IS the call; no independent measured phenotype → faithful-to-tool only | GREEN-CALLER (already partly covered by `dna-pathotype` + the virulence overlay) |
| 6 | Catabolic capability / biotype (sugar fermentation, urease, oxidase) | ~ operon-presence (partial) | ✗ biochemical (API-panel) results rarely linked to public genomes at isolate scale | GREEN-CALLER / mostly NO_FREE_SOURCE |
| — | flagella length / cell size / "look" / growth rate | ✗ no determinant catalog | ✗ no free measured isolate-level label | **RED** — closed 0-for-4 embedding negative; do NOT reopen |

## Recommendation (ranked)
1. **S. pneumoniae capsular serotype — BUILD NEXT.** Passes BOTH gates with the cleanest, largest free
   measured label in the whole project (11,810 Quellung-typed genomes vs the Klebsiella cell's 733). It is a
   near-exact reuse of the just-built `dna-ktype` / `dna-serotype` caller pattern (BLAST a curated cps-locus
   DB → serotype call) + the trust-surface + report-card infra. It would be the **2nd GREEN-VALIDATED non-AMR
   cell with a real independent number** (alongside pathotype), not NO_FREE_SOURCE growth. Honest ceiling to
   match/beat: ~89% vs Quellung. This is the genuinely high-value next cell.
2. **Salmonella serovar** — strong #2; canonical identity trait, free traditional-serotyping labels exist, but
   the circularity filter (exclude SeqSero-predicted "serovars") is real label-curation work before scoring.
3. Everything else: GREEN-CALLER (NO_FREE_SOURCE) → defer; building them just grows the saturation count the
   report card already flags.

## Why this is a real result, not motion
The pivot memo warned that adding cells without a free label only grows `NO_FREE_SOURCE`. This triage found
**two traits that pass gate 2 with wet-lab measured labels** (pneumo Quellung, Salmonella traditional) — so the
next cell can be GREEN-VALIDATED (scored, trust-badged), the gold-standard shape, not saturation filler.

## Next step (STOP for user direction — triage precedes build)
- **Build the S. pneumoniae serotype cell** (I'd first pin the exact free Quellung-labeled cohort + the cps
  catalog version, then build the caller mirroring `dna-ktype`, then score vs Quellung for the GREEN-VALIDATED
  number), **or**
- pick Salmonella serovar instead, **or**
- bank the triage (the ktype validation is still accruing; the decoder is already trait-general).

No code / no `/idea-anchor` fired — this is the triage memo (class-e research-only), mirroring the EP-4 + Round-1 discipline.

Sources: [GPS Pipeline (Nature Communications 2025)](https://www.nature.com/articles/s41467-025-64018-5) · [GPS Pipeline (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12460886/) · [SeqSero2 (ASM AEM)](https://journals.asm.org/doi/10.1128/aem.01746-19) · [WGS Salmonella serotyping vs gold-standard (Frontiers 2025)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1685741/full)
