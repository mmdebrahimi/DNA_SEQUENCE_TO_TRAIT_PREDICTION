# Eye-colour decoder — first off-pathogen cell (built) + OpenSNP label-source finding (2026-06-28)

> **UPDATE 2026-06-28 (later same day): VALIDATED ON REAL DATA.** The user ratified downloading the archived
> OpenSNP dump (deliberate, colleague-reviewed). The ~20 GB 2017 mirror was fetched to D:, and the rs12913832
> v0 rule scored **accuracy 0.993 / brown-sens 0.98 / blue-spec 1.0 on N=298** binary-scored users — on the
> DECISIVE (homozygote) subset, with honest abstention on heterozygotes. The "live ethical label source
> needed" disposition below is SUPERSEDED for OpenSNP. Full result: `wiki/eye_colour_opensnp_validation_2026-06-28.md`.
> The ethics framing (founders deleted it; this is the archive) stands and is recorded honestly.

The wide-net scan's flagship (`wiki/wide_net_gene_phenotype_source_scan_2026-06-28.md`): generalize the
deterministic gene→trait decoder from AMR to a human visible trait. **The RULE side is built, sourced, and
tested; the live label source (OpenSNP) turned out to be DEAD + ethically-withdrawn** — an external+ethical
wall, surfaced honestly rather than forced.

## What shipped (durable, autonomous, no fabrication)
- **`dna_decode/data/eye_colour.py`** — the rs12913832 (HERC2) single-locus eye-colour rule. v0 = the single
  strongest blue/brown predictor (~74% of blue-vs-brown variance alone). **Strand-agnostic** (the load-bearing
  sourced fact: dbSNP GG→blue/AA→brown, but 23andMe reports the COMPLEMENT C/T; rs12913832 is A/G =
  non-palindromic, so blue allele={G,C}, brown={A,T} — a memory-inversion was caught + corrected by sourcing).
  The IrisPlex 6-SNP model is the documented v0.1 (rsIDs named; coefficients must be SOURCED from Walsh 2011,
  NOT fabricated → deferred).
- **`scripts/eye_colour_opensnp_validate.py`** — OpenSNP ingestion + scoring harness: parses 23andMe/
  AncestryDNA raw files → rs12913832 call, bins free-text eye colour → blue/brown/other, scores blue-vs-brown
  with the honesty caveats (self-reported label = near-independent non-circular; ancestry-confound; intermediate
  excluded). Reports `NOT_FETCHED` (never a synthetic number) when no data.
- **`tests/test_eye_colour.py`** (8 tests) — pins strand-agnosticism (both A/G and C/T give the IDENTICAL
  call), the DTC parser (23andMe + AncestryDNA shapes), and the label binner. FROZEN AMR surface byte-unchanged.

## The label-source finding (the real boundary)
**OpenSNP was sunset 2025-04-10 and ALL user data deleted.** The founders' stated reason: the dissolution of
23andMe + concern about weaponization of genetic data by authoritarian regimes, with vulnerable populations in
the set (The Register / 404media / Wikipedia, 2025). The live API is gone (opensnp.org → DreamHost
"site not found"; reachability + SSL-intercept confirmed). The ONLY remaining copy is an Internet Archive
mirror (`archive.org/details/opensnp_data_dumps`).

**This is a Care-surface, not just a wall.** Downloading a privacy-sensitive human genetic dataset that its
own creators *deliberately deleted for ethical reasons* is a decision for the USER, not an autonomous executor
action — even though the harness is format-ready for the archived dumps. NOT fetched.

## Ethically-clean paths to actually validate this cell (recommended pivot)
The eye-colour CALLER is reusable on ANY DTC-format file with zero ethics issue:
1. **The user's OWN 23andMe/AncestryDNA file** — drop it in `data/raw/opensnp/` (a 2-line phenotype CSV +
   the raw file) → an instant self-test of the rule. Zero fetch, zero ethics.
2. **Pivot the "decoder generalizes off-pathogen" validation to an ethically-clean LIVE source** (shortlist
   #2/#3): **dog coat colour** (MC1R/ASIP/CBD103 + Darwin's Ark, live, citizen-science, no human-privacy issue)
   or **yeast 1011 Mendelian traits** (CuSO₄/anisomycin resistance — cleanest *lab-measured* labels, live,
   no ethics issue). These prove the same generalization claim without the OpenSNP ethics/availability wall.

## Honest disposition
- The off-pathogen decoder RULE is DONE + validated-as-code (the generalization claim is demonstrated at the
  code level; a real-data score awaits an ethically-clean label source).
- OpenSNP (the flagship's label source) is **dead + ethically-withdrawn** → reclassify the flagship from
  "free cell available today" to "rule ready; live ethical label source needed." The wide-net scan's
  conclusion stands, but the FLAGSHIP's data path is now via dog/yeast, not OpenSNP.
- Next: a fresh `--until-mvp` on an ethically-clean live source (dog coat colour or yeast) — a user pick,
  not auto-started (it's a new data-fetch cell + the OpenSNP ethics finding warrants a deliberate source choice).

## Sources
OpenSNP shutdown: [The Register 2025-04-01](https://www.theregister.com/2025/04/01/opensnp_shutdown/),
[404media](https://www.404media.co/open-source-genetic-database-opensnp-shuts-down-to-protect-users-from-authoritarian-governments/),
[archive mirror](https://archive.org/details/opensnp_data_dumps); rs12913832 strand orientation:
[PLOS One 2020 PMC7485777](https://pmc.ncbi.nlm.nih.gov/articles/PMC7485777/). Catalog:
`wiki/wide_net_gene_phenotype_source_scan_2026-06-28.md`.
