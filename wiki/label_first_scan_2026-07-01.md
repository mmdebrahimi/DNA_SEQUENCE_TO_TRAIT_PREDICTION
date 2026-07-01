# Label-first scan — off-pathogen datasets vs the admission gate (2026-07-01)

Executed the recommendation from the off-pathogen brainstorm: scan for a dataset clearing the label-first
admission gate (`wiki/off_pathogen_cell_admission_gate_2026-07-01.md`) BEFORE choosing a trait. The gate's
binding AND is **measured/observed (not self-report) free no-DUA label + a KNOWN deterministic (Mendelian)
rule + open individual-level genotype access**. That AND filters out polygenic quantitative traits (the
embedding arm's closed-negative territory).

## Gate screen (5 items: 1 measured-free-label · 2 stable-IDs · 3 genotype-access · 4 deterministic-rule · 5 provenance/privacy)

| Candidate | 1 | 2 | 3 | 4 | 5 | Verdict |
|---|---|---|---|---|---|---|
| **Horse coat colour** (Dryad `10.5061/dryad.3q111`; VGL/Rieder MC1R-E × ASIP-A rule) | ~ | ✓ | ✓ | **✓✓** | ✓ | **BEST — clears 4/5; item 1 unverified (see below)** |
| Arabidopsis 1001 / AraPheno (FRI→flowering, DOG1→germination) | ✓ | ✓ | ✓ | ~ | ✓ | runner-up — rule is major-effect but flowering is QUANTITATIVE (embedding G2 was a closed negative) |
| C. elegans CaeNDR (wild isolates, single-locus large-effect) | ✓ | ✓ | ✓ | ~ | ✓ | runner-up — single-locus but deterministic variant→trait rules are niche/less-curated |
| Yeast 1011 (Peter 2018, SGD; 36–223 measured traits) | ✓ | ✓ | ✓ | **✗** | ✓ | FAIL #4 — traits are polygenic/CNV-driven (GWAS shape, not a curated Mendelian rule) |
| Dog coat colour (Wisdom Panel 65k / MyDogDNA 12k) | **✗** | ✗ | ✗ | ✓✓ | ✓ | FAIL #1/#3 — clean rule, but data is COMMERCIAL-gated (no open individual dataset) |

## The winner: horse coat colour

- **Deterministic rule (item 4) is textbook Mendelian epistasis** — MC1R extension E/e is epistatic to ASIP
  agouti A/a: `e/e → chestnut` (regardless of agouti); `E_ aa → black`; `E_ A_ → bay`. Molecularly pinned:
  chestnut = MC1R C901T; black = ASIP 11-bp exon-2 deletion (Rieder et al. 2001). This is our decoder's exact
  shape and only needs 2 loci (no WGS).
- **Open data (items 2/3)** — Dryad has one downloadable table (`JOH2015154R1_RawData_Updated.xlsx`, N=215
  Tennessee Walking Horses) with per-horse MC1R + ASIP genotypes; Rieder 2001 (120 horses, 11 observed
  colours) + the Noma-horse study (CIE1976 L*a*b* QUANTITATIVELY measured colour) are additional public
  supp-table sources.
- **Non-human, no privacy issue (item 5).** ✓ A step up from the human self-report pilots.

## The one unverified gate item (must check before building)

**Item 1 — is coat colour INDEPENDENTLY OBSERVED, or genotype-DERIVED?** The Dryad TWH table is from a
BEHAVIOR study; its colour may be inferred from the E/A genotype → that would be **circular** (can't validate
a rule against its own output — the circular-label gate). The Dryad file could not be fetched headlessly
(403/401 — needs a browser session), so this is UNVERIFIED. **Rieder 2001 + the Noma CIE-Lab study are the
safer sources** — there the colour was observed/measured independently (that's how the rule was discovered).

## Honest meta-caveat

Coat-colour rules (like IrisPlex/ABO) are ALREADY-DEPLOYED published tools (VGL sells these tests). So a
horse cell is STILL "validate a deployed rule" — but on a MEASURED/OBSERVED (not self-report) label +
non-human. That **upgrades** the off-pathogen evidence tier from PILOT/DEMO (self-report) to
**measured-label validated** — it clears the gate, which the eye/ABO cells do not. It is not a new KIND of
contribution, but it is the honest next rung.

## Recommendation

The gate is **cleared in principle** by horse coat colour (4/5; item 1 pending an independent-colour source).
Next move (a fresh `--until-mvp`, user-greenlit): fetch an independent-colour source (Rieder 2001 / Noma
supp table, or the Dryad xlsx via browser), VERIFY colour is observed-not-derived, build the deterministic
E/A→{chestnut,bay,black} rule + tests, score. If item 1 fails on every source → the gate is binding and
banking at pilot strength (fork #1) is the honest ceiling. Runner-ups (Arabidopsis FRI, C. elegans) are
fallbacks if horse colour proves circular.

## Sources
Yeast 1011: [Nature 2018](https://www.nature.com/articles/s41586-018-0030-5) + SGD. Horse: [Dryad
10.5061/dryad.3q111](https://datadryad.org/dataset/doi:10.5061/dryad.3q111), [Rieder 2001
PMID 11353392](https://pubmed.ncbi.nlm.nih.gov/11353392/), [VGL horse coat](https://vgl.ucdavis.edu/resources/horse-coat-color).
Dog: [VGL dog coat](https://vgl.ucdavis.edu/resources/dog-coat-color) (commercial data). Arabidopsis:
[1001 Genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC2718507/) + AraPheno. C. elegans: CaeNDR.
