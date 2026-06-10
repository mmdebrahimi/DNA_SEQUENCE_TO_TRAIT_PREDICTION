# Antiviral NA inhibitor decoder — the 4th kingdom (viral) — 2026-06-10

> Extends the shipped deterministic target-site method one domain further: **influenza A virus**
> neuraminidase-inhibitor (NAI) resistance, headline **H275Y / oseltamivir** (N1 numbering). Same proven
> engine as bacteria (AMRFinder), fungi (ERG11), protozoa (K13/pfcrt) — a hand-curated catalog + a
> BLAST(NA-CDS-vs-genome) → codon-map → check-marker call reusing the gene-generic `observed_substitutions`.
> No money, laptop-only, real public data. `/soraya --advance` increment.

## What shipped

| Piece | Path |
|---|---|
| Catalog (4th-kingdom) | `dna_decode/data/antiviral_amr.py` — oseltamivir / peramivir / zanamivir → NA markers |
| Caller (thin wrapper) | `scripts/flu_na_caller.py` — reuses `observed_substitutions` (DRY) |
| CLI route | `dna_decode/amr/cli.py` `_antiviral_main` + `--na-ref` + drug-choices union |
| Reference (real, committed) | `data/antiviral_ref/N1_NA_NC026434_cds.fna` (RefSeq NC_026434.1, A/California/07/2009 H1N1pdm09, 470 aa, WT His275) |
| Real field fixtures (committed) | `Flu_N1_NA_EU716587_1_H275Y.fna` (real R), `Flu_N1_NA_EU716623_1_WT.fna` (real S) |
| Tests | `tests/test_amr_cli_antiviral.py` — 11 tests (catalog + observed + routing + real-BLAST + real-fixture) |

## Catalog (N1 numbering — consistent with the shipped reference)

- **oseltamivir** → NA:{H275Y, I223R, S247N} — H275Y is the globally dominant, highly-reduced-inhibition marker.
- **peramivir** → NA:{H275Y, I223R} — shares the active-site framework; H275Y cross-resistance.
- **zanamivir** → NA:{Q136K, E119G, I223R} — H275Y barely affects zanamivir (deliberately NOT in this set).

## Validation — REAL field isolates, real BLAST, real biological discriminator

12 real H275Y (R) + 44 real WT (S) N1 NA field isolates (2008–2009 seasonal H1N1, when H275Y was circulating)
were scanned from NCBI; one of each committed as a fixture. End-to-end through the CLI (real makeblastdb+blastn):

| Input (real field isolate) | drug | call | driven by |
|---|---|---|---|
| EU716587.1 (H275Y) | oseltamivir | **R** ✓ | NA:H275Y |
| EU716623.1 (WT His275) | oseltamivir | **S** ✓ | — |
| EU716587.1 (H275Y) | **zanamivir** | **S** ✓ | — (H275Y doesn't reduce zanamivir) |

The oseltamivir-R / zanamivir-S split on the SAME isolate is a genuine **biological** validation of the
per-drug directional catalog — not just plumbing. The codon mapper (real BLAST) places His275→Tyr correctly
across a real, non-engineered field allele.

## What this establishes / honest scope

- **Establishes:** the deterministic target-site method now spans **four kingdoms** (bacterial / fungal /
  protozoan / **viral**) on ONE shared engine + uniform `amr-mechanism-call-v1` record. Influenza NA is
  intronless (non-spliced segment) so the colinear single-HSP codon-mapper applies directly.
- **Does NOT establish:** a phenotype-validated surveillance tool. The call = the recognized marker's
  presence; NAI resistance is really an NI-assay IC50 fold-change (no growth MIC). An **S** call surfaces
  the blind spots it cannot rule out (different-class baloxavir PA/PB2 I38X, adamantane M2 S31N, novel/
  uncatalogued NA substitutions, permissive secondary mutations, mixed-population minor variants).
- **Marker scope:** v0 is the well-established N1 oseltamivir/zanamivir markers. Subtype-specific numbering
  (N2 H274Y ≡ N1 H275Y) and broader marker sets are deliberate follow-ons; the reference is N1, so v0 calls
  are N1-numbering-consistent. A genuinely independent phenotype label (NI-assay IC50 from a public
  surveillance set) would be the next de-confounder — same "validate the wrapper vs the naive tool on
  independent labels" discipline as the other kingdoms.

## Bottom line

The proven engine transfers to the viral kingdom cleanly and is validated on real field isolates with a
real drug-specificity discriminator. The pattern is now confirmed across all four kingdoms it can reach with
a single-locus curated target-site mechanism. Honest read: this is **breadth via a proven pattern** — the
4th application of the same method — not a new architectural finding. The decoder's hard blind spot
(everything gene-presence/target-site can't see) is unchanged.
