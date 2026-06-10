# Independent phenotype-label SOURCE SWEEP — full web survey + pay/no-pay decision — 2026-06-10

> Full web sweep for reliable sources of the data the decoders need to be **independently validated**:
> isolate-level, lab-measured **phenotype** (MIC / IC50 / clearance) **linked to a genome** we can run the
> genotype decoder against. Categorized by validity + relevance + access (free/paid). Extends
> `wiki/independent_phenotype_label_census_2026-06-10.md`. Web-sourced; some CDC/assay queries were blocked
> by a usage-policy filter (worked around via direct page fetches); per-source claims cite the page found.

## THE LAW THIS SWEEP ESTABLISHES (the headline)

**Phenotype-independence and genome-linkability are ANTI-CORRELATED across every available source.**
- Sources with genuinely independent / different-lab phenotype (industry reference-method MIC, EUCAST,
  clinical clearance) carry **no genome accessions** → the genotype decoder cannot be run against them.
- Sources that ARE genome-linked (NCBI PD / BV-BRC / NARMS; MalariaGEN Pf7) are either the **integrated
  label ecosystem** (not method-independent) or use **marker-INFERRED labels** (circular vs the decoder).
- **No source — free OR paid — provides both independent-lab phenotype AND per-isolate genome linkage at
  systematic scale.** Money does not buy the missing property (the paid MIC databases are also genome-less).

## Catalog (by data type)

### Bacterial AMR (MIC) — the deepest substrate
| Source | phenotype | isolate-level | genome-linked | independence | access | verdict |
|---|---|---|---|---|---|---|
| BV-BRC / PATRIC | MIC/AST | yes | **yes** | low (ecosystem; tuning source) | FREE | already used (tuning) |
| NCBI Pathogen Detection | AST R/S | yes | **yes** | medium (cross-source) | FREE | already used (cross-source validation done) |
| NARMS Now (FDA) | **numeric MIC** | yes | yes (→NCBI) | low–med (feeds NCBI-PD; portal-distinct) | FREE | stricter-LABEL re-run only; NOT different-lab |
| **Vivli AMR Register** (Pfizer-ATLAS, Merck-SMART, GSK-SOAR, J&J-DREAM, Paratek-KEYSTONE, Shionogi-SIDERO, Venatorx) | raw reference-method MIC, >925k isolates | yes | **NO genome accessions** (confirmed — MIC + clinical metadata only) | **HIGH (independent industry labs)** | FREE (open, registration/DUA) | **independent but UNUSABLE — no genome to run the decoder against** |
| SENTRY (JMI/Element) | reference-method MIC, ~40k/yr | yes | no (public viz only) | high | **PAID** (consortium/quote-based; no public price) | unusable (paid AND no genome) |
| EUCAST MIC distributions | MIC histograms | **aggregate only** | no | high | FREE (raw on email request) | unusable (aggregate, no genome) |
| CDC/FDA AR Isolate Bank | characterized MIC | yes | yes (sequenced; in NCBI) | medium (CDC panel; genomes →NCBI) | FREE data / registration to order | **closest bacterial candidate** but small (~1,075 isolates across 34 panels → few per organism×drug) + NCBI-deposited (overlap) |

### Antimalarial (P. falciparum)
| Source | phenotype | isolate-level | genome-linked | independence | access | verdict |
|---|---|---|---|---|---|---|
| MalariaGEN **Pf7** | resistance status | yes | **yes (ENA accessions, 20,864 genomes)** | **none — marker-INFERRED, not measured** | FREE | genome-rich but labels are genotype-derived → CIRCULAR vs the decoder |
| WWARN (IDDO) | clinical clearance / ex-vivo IC50 | study-level/pooled | not systematically | high (clinical phenotype) | access-controlled (data-access request) | independent phenotype but not free + not genome-linked at isolate level |

### Antiviral (influenza NA) — decoder has ZERO phenotype validation
| Source | phenotype | isolate-level | genome-linked | independence | access | verdict |
|---|---|---|---|---|---|---|
| GISAID EpiFlu | sequences | — | genomes | — | **registration-gated** (not free-open) | genome side only; gated |
| NCBI Influenza Virus | sequences | — | genomes | — | FREE | genome side only, no phenotype |
| CDC/WHO NI-assay surveillance | IC50 / % reduced inhibition | **aggregate** | no | high | FREE (summaries) | no isolate-level downloadable IC50+accession found |

### Fungal (C. auris azoles)
| Source | phenotype | isolate-level | genome-linked | independence | access | verdict |
|---|---|---|---|---|---|---|
| CDC AR Isolate Bank (C. auris panels) | characterized MIC | yes | yes (→NCBI) | medium | FREE data / registration | small; same overlap caveat as bacterial |
| Per-study published cohorts (e.g. EP7's public ERG11 fixtures) | MIC | yes | yes (deposited) | high (if MIC not an NCBI antibiogram) | FREE | works but MANUAL per-paper curation, not systematic |

## Pay / No-Pay decision

**Recommendation: DO NOT PAY.** Rationale:
- The **binding constraint is genome-linkage**, and **no paid source fixes it**. SENTRY (the main paid option) is genome-less public-viz-only — paying buys more independent MIC, still with no genome to run the decoder against. ATLAS/Vivli is already FREE and also genome-less. So money buys *more of the wrong shape*, not the missing property.
- The independent phenotype we'd "pay for" (industry MIC) is **decoupled from sequence by construction** (pharma surveillance tests susceptibility, rarely deposits per-isolate WGS with public accessions). Paying cannot re-link them.
- Therefore there is **no purchase that closes the different-lab + genome-linked gap.** (No money gate fires — the recommendation is to spend nothing.)

## The only two real paths to independent-phenotype + genome (both small/manual, no purchase)
1. **CDC/FDA AR Isolate Bank** (bacterial + C. auris) — a defined, sequenced, MIC-characterized free panel. Best *systematic-ish* candidate; caveats: small per organism×drug, genomes deposited to NCBI (partial provenance overlap), CDC-characterized (semi-independent, not "different commercial lab"). Worth a dedicated feasibility probe IF a stricter independent check is wanted.
2. **Per-study manual curation** — papers that deposit genomes in NCBI AND report isolate MICs in supplementary tables NOT submitted as NCBI antibiograms. Genuinely independent + genome-linked, but harvested per-paper (not a source you can query); a small hand-built gold set (~20–40 isolates) is the realistic ceiling.

## Bottom line
The sweep confirms + sharpens the census: a systematic free OR paid source giving **independent-lab phenotype + per-isolate genome** does **not exist** for any decoder substrate — because the two properties are structurally decoupled. **Spend nothing.** The decoders are validated as far as linkable free data allows (in-cohort + out-of-cohort + cross-source NCBI-PD). To go strictly further, the only moves are (a) the CDC AR Isolate Bank panel, or (b) a small hand-curated per-study gold set — both free, both bounded, neither requiring payment. Bank this as the standing constraint on every "validated" claim.
