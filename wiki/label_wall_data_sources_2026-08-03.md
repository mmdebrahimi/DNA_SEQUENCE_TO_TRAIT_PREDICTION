# Breaking the label wall — verified downloadable measured-label sources (2026-08-03)

The FBA / non-metabolic cells are stuck at KNOWLEDGE_BASELINE / partial because no free, genome-keyed,
MEASURED label was fetchable. A deep search (2026-08-03) found several that are — some I can pull directly,
most you download manually (registration / SI / bot-blocked servers) and I ingest. Each is mapped to the
cell it unblocks + a join-feasibility note (the load-bearing detail: the label's gene IDs must match the
model/genome's IDs).

## Ranked shopping list

| # | dataset | unblocks | measured label | fetchable? | JOIN feasibility |
|---|---|---|---|---|---|
| 1 | **Madin 2020 `bacteria-archaea-traits`** (GitHub `output/condensed_species_NCBI.csv`) | **motility** cell (+ some metabolic) | per-species motility (yes/flagella vs no), 14,893 species (~6,400 labelled) | **YES — I fetched it** | species → need a genome/feature-table per species to run the presence cell. **NOISY** (mislabels *Shigella flexneri* as motile — a species-aggregation artifact my cell correctly calls non-motile) |
| 2 | **P. aeruginosa GOLD_84 / GOLD_115** (PLOS Comput. Biol. 2026, SI) + PNAS 2015 (352 essential) | **essentiality (P. aeruginosa)** — currently LABEL_WALLED | Tn-seq gold-standard essential genes, PA#### loci | **manual (SI)** | **CLEAN** — PA#### tags join `iJN1463` model gene IDs directly. The best clean cross-organism essentiality unblock |
| 3 | **PMkbase** (PMC12584640) | **carbon-source growth** — currently RECALL-only (specificity walled) | Biolog PM carbon/N substrate pos+neg for **E. coli / P. putida / S. aureus** (41,664 carbon data points) | **manual (web tool / download)** | gives the measured NEGATIVES the Keio positive-only set lacked → full specificity for the carbon cell |
| 4 | **OGEE v3** (v3.ogee.info) / **DEG 15** (tubic) | **essentiality (many organisms)** | aggregated experimental essential/non-essential, locus-keyed, 91 species (OGEE) | **manual (both timed out from this host — you can download)** | per-organism; join depends on the model's ID scheme (PA#### clean; SAUSA300 needs a crosswalk) |
| 5 | **NTML** (ntml.unmc.edu) + Coe 2019 multi-strain Tn-Seq (PLOS Path SI) | **essentiality (S. aureus)** | 579 non-interrupted ORFs = essential, SAUSA300 loci | **manual (SI / web)** | **needs a crosswalk** — `iYS1720` model uses STM#### (Salmonella-style) IDs, not SAUSA300. Join on gene NAME (ArgD…) or swap to a SAUSA300-keyed S. aureus GEM |
| 6 | **BacDive API** (bacdive.dsmz.de) | motility + substrate + many phenotypes, strain-level | 2.6M data points / 97k strains, incl. motility + 8,151 substrate entries | **manual (free registration → API key)** | strain-level (better than Madin's species-level); still needs genome per strain |
| 7 | **metaTraits** (NAR 2026) / **BactoTraits** / **ProTraits** (protraits.irb.hr) | motility + metabolic, broad | integrated multi-source trait matrices | **manual (download)** | species-level, same join shape as Madin; use for coverage breadth |

## What this means for the label wall

**It is partially breakable — cleanly for one cell, noisily/with-a-crosswalk for others:**

- **Cleanest immediate win — P. aeruginosa essentiality (#2).** GOLD_84 is PA####-keyed = a direct join to the
  already-shipped `iJN1463` model. Download the SI → I ingest → the `LABEL_WALLED` P. aeruginosa essentiality
  cell becomes SCORED with a real Tn-seq gold-standard number. No crosswalk, no noise.
- **Fetchable-now but noisy — motility (#1).** I can pull Madin today and upgrade the motility cell from
  literature-anchors to a measured multi-species cohort — but honestly labelled: Madin has a real error rate
  (Shigella), so the headline must be measured-cohort-with-noise, not clean.
- **Specificity fix — carbon growth (#3).** PMkbase gives the measured NEGATIVES the carbon cell lacked →
  turns RECALL-only into full accuracy/specificity for E. coli + S. aureus.
- **S. aureus essentiality (#5)** needs a `STM → SAUSA300` gene-name crosswalk (or a different S. aureus GEM);
  the label now EXISTS (NTML), the blocker is model-ID keying.

## Honest caveats (so a download doesn't over-promise)

- Compiled trait DBs (Madin/BacDive/metaTraits) are **species-aggregated** → real but noisy (Shigella).
  Strain-level (BacDive) is cleaner but still needs a genome per strain.
- Every one of these still requires a **genome / feature-table per organism** to run the presence-based or
  FBA cells against the label — the label breaks the wall, but the validation is a fetch+assemble per cohort.
- **Yeast is NOT in the bacterial trait DBs** — the yeast essentiality tier stays SGD (already have it).
- Independence: these are MEASURED (Tn-seq / Biolog / phenotype assays), so they're genuine phenotype labels
  — but a decoder scored on the same knowledge base a catalog was built from is in-distribution; the strongest
  claim (independent) needs the label to be provenance-separate from the catalog's construction.

## Immediate executable next steps

1. **Motility measured-cohort validation (I can run now):** fetch Madin labels + genomes for the project's
   organisms (E. coli / Salmonella / Klebsiella / P. aeruginosa / S. aureus) → run `dna-decode motility`
   per genome → score vs Madin (with the Shigella-noise caveat). Upgrades the motility cell's tier.
2. **P. aeruginosa essentiality (you download GOLD_84 SI → I ingest):** the cleanest wall-break; `iJN1463`
   single-gene-deletion vs GOLD_84 → a real Tn-seq-validated cross-organism essentiality number.
3. **Carbon specificity (you download PMkbase E. coli carbon PM → I ingest):** full accuracy for the carbon cell.

Sources: Madin 2020 (Sci Data / GitHub `bacteria-archaea-traits`) · P. aeruginosa gold-standard (PLOS Comput
Biol 2026) + Lee 2015 (PNAS) · PMkbase (PMC12584640) · OGEE v3 · DEG 15 · NTML (UNMC) · BacDive · metaTraits (NAR 2026).
