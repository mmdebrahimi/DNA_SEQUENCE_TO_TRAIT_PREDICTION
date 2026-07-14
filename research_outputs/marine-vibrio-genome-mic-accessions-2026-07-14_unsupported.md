# Marine / aquaculture bacteria — genome + measured-AST accessions — unsupported claims (V1 invocation)

> Slug: marine-vibrio-genome-mic-accessions-2026-07-14. Captured 2026-07-14. These rows were rejected during V1 intake.
> Reasons: missing audit-floor locator (no verbatim source quote), low-confidence label.

## Rejected rows

| Row content | Rejection reason | Suggested follow-up |
|---|---|---|
| BioProject PRJNA622672 — foodborne *Vibrio* blaCTX-M-14 clade | Low confidence: cited as an **unpublished** study inside another paper; genome + measured-AST pairing unverified | Find the primary deposit for PRJNA622672 on NCBI and confirm whether it carries any measured AST before use |
| Pathogenwatch hosts a public *V. cholerae* genome collection with AMR annotation | Audit-floor fail: no verbatim source quote (the collections page is JS-rendered and did not return content this run); named from domain knowledge only | Open pathogen.watch in a browser, confirm the V. cholerae collection URL + whether it exposes measured vs predicted AST. **Note: Pathogenwatch AMR is largely genome-PREDICTED — likely G1-circular for the decoder; useful as a genome pointer, not a label source** |
| NCBI Pathogen Detection includes *Vibrio* taxgroups with `AST_phenotypes` | Audit-floor fail: no verbatim source quote (JS-rendered browser); named from domain knowledge only | Open ncbi.nlm.nih.gov/pathogens, filter to a Vibrio taxgroup, and check how many isolates have a populated `AST_phenotypes` (measured) field. This is the SAME source the decoder already uses for E. coli/Klebsiella — the exact data shape to show Hamid |

## Summary

- Total rejected: 3
- Reason breakdown:
  - Missing audit-floor locator (no verbatim quote): 2
  - Low confidence: 1
  - Mapping-floor failure: 0
  - Hard-reject banned phrase: 0

**Honest-gap note (not a rejected row, but the biggest gap):** Aeromonas salmonicida/hydrophila, Piscirickettsia salmonis, and other fish/shrimp-farm aquaculture-pathogen genome+measured-AST accessions could NOT be retrieved this run — every resistance-phrased search query tripped a content-filter false-positive. These datasets very likely exist; a follow-up run on an unblocked channel (or a direct Google Scholar / ENA search) is the way to pin them.
