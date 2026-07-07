<!-- memo-schema: 0.4 -->
# CYP2D6-CYP2D7 PSV analysis for hybrid IDENTITY (*13/*36/*68) — supported memo (2026-07-07)

## Research Context
**Problem anchor (verbatim):** `fresh research arc (PSV curation + read-level pileup)  hybrid identity (which of *13/*36/*68 via CYP2D6-vs-CYP2D7 PSV analysis)`

dna_decode ships a CYP2D6 hybrid-PRESENCE detector (elevated CYP2D7 read depth; sens 0.62 / spec 1.0 on 38
1000G CRAMs) but does NOT resolve hybrid IDENTITY (which of *13/*36/*68). This memo curates the field's
method (Cyrius/DRAGEN/PharmVar) for that identity step — what a v0.4 "Cyrius-class" build would need.

**Validation:** 16/16 rows pass the audit floor (value + source + section + stable URL + access date + verbatim
quote each). 0 banned phrases, 0 cite-token noise. 11 high / 5 medium (authoritative Illumina technical docs)
/ 0 low. 0 rows rejected.

## Supported audit table

| Claim / quantity | Value | Units | Source | Year | Section | URL | Quote (≤25 words) | Conf |
|---|---|---|---|---|---|---|---|---|
| CYP2D6/CYP2D7 differentiating bases (PSVs) Cyrius uses to identify hybrids | 117 | bases | Cyrius (Chen et al., Pharmacogenomics J) | 2021 | Methods | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "identified 117 reliable bases that differ between CYP2D6 and CYP2D7" | high |
| Hybrid-detection principle | qualitative | — | Cyrius (Chen et al.) | 2021 | Methods Fig.2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "Hybrids are identified when the CN of CYP2D6 changes within the gene." | high |
| Cyrius concordance vs truth (untrained) | 96.5 | % | Cyrius (Chen et al.) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "five discordant calls in the 144 truth samples, showing a concordance of 96.5%" | high |
| Cyrius trained concordance | 99.3 | % | Cyrius (Chen et al.) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "reaching a 'trained' concordance of 99.3% (143 out of 144 samples)" | high |
| Aldy CYP2D6 concordance | 86.8 | % | Cyrius (Chen et al.) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "Aldy had a concordance of 86.8%" | high |
| Stargazer CYP2D6 concordance | 84 | % | Cyrius (Chen et al.) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "Stargazer had a concordance of 84%" | high |
| Cyrius validation truthset | 144 | samples | Cyrius (Chen et al.) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "compared the CYP2D6 calls ... against 144 truthset samples" | high |
| Benchmark hybrids: NA12878=*68, NA24631=*36, HG01161=*13 | qualitative | — | Cyrius (Chen et al.) | 2021 | Figs.2-4/Table S1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7997805/ | "NA24631 ... *36 ... HG01161 ... *13 ... NA12878 ... *68" | high |
| DRAGEN star-defining small variants | 118 | variants | DRAGEN CYP2D6 Caller | 2023 | Small Variant Calling | https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/CYP2D6_Caller_fDG.htm | "118 small variants that define various star alleles are detected" | medium |
| DRAGEN variants in unique high-MAPQ CYP2D6 regions | 96 | variants | DRAGEN CYP2D6 Caller | 2023 | Small Variant Calling | https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/CYP2D6_Caller_fDG.htm | "96 of these variants are in unique (nonhomologous) regions of CYP2D6 with high mapping quality" | medium |
| Hybrid breakpoints: *36 exon9 / *68 intron1 / *13 exon9,intron4,intron1 | qualitative | — | DRAGEN CYP2D6 Caller | 2023 | Structural Variant Calling | https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/CYP2D6_Caller_fDG.htm | "*36 ... exon 9 ... *68 ... intron 1 ... *13 ... exon 9, intron 4, intron 1" | medium |
| CYP2D6*36 activity | negligible | — | Gaedigk et al. (hybrid tandems) | 2019 | Discussion | https://pmc.ncbi.nlm.nih.gov/articles/PMC6556886/ | "activity encoded by CYP2D6*36 is negligible in vitro and in vivo" | high |
| T-insertion in exon 1 marks nonfunctional CYP2D7-2D6 hybrids | qualitative | — | Gaedigk et al. | 2019 | Results/Discussion | https://pmc.ncbi.nlm.nih.gov/articles/PMC6556886/ | "carried a T-insertion in exon 1 abolishing activity" | high |
| Hybrid-tandem frequency among dup-positive alleles (Caucasian) | 12 | % | Gaedigk et al. | 2019 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC6556886/ | "four out of 33 (12%) of the duplication-positive alleles were hybrid tandems" | high |
| Hybrid-tandem overall allele frequency (Caucasian) | 0.25 | % (upper) | Gaedigk et al. | 2019 | Conclusion | https://pmc.ncbi.nlm.nih.gov/articles/PMC6556886/ | "Hybrid tandem alleles occur infrequently (<0.25%) in Caucasians" | high |
| *36 embedded exon-9 gene-conversion can yield NO CYP2D6-CN change (depth blind spot) | qualitative | — | DRAGEN CYP2D6 Caller | 2023 | Structural Variant Calling | https://support-docs.illumina.com/SW/DRAGEN_v38/Content/SW/DRAGEN/CYP2D6_Caller_fDG.htm | "exon 9 can participate in CYP2D7 gene conversion resulting in an embedded CYP2D7 sequence instead of a true hybrid" | medium |

## Decisions for Human Confirmation

| # | Decision | Candidate use | Verification needed |
|---|---|---|---|
| 1 | Adopt the Cyrius PSV method (117 differentiating bases, per-site CYP2D6-CN profile) as the identity engine | Curate the 117 PSV coordinates from the Cyrius GitHub `CYP2D6.json` config (NOT the paper) + compute per-site CYP2D6-CN from CRAM pileup; a CN-change location → breakpoint → identity | Confirm the Cyrius config license permits reuse; verify the 117 coords lift to GRCh38 |
| 2 | Use the breakpoint-location discriminator (*36 exon9 / *68 intron1 / *13 exon9|intron4|intron1) to assign identity | Map the WITHIN-gene CN-change position to the breakpoint table → the specific hybrid allele | The DRAGEN breakpoint list is a vendor doc; cross-check against PharmVar structural-variation definitions before pinning |
| 3 | Validate against the named benchmark anchors (NA12878=*68, NA24631=*36, HG01161=*13) | These are fetchable GIAB/1000G CRAMs (the same remote-CRAM tooling dna_decode already proved) → a real-surface identity test | Confirm each anchor's CRAM is reachable + the truth is the current PharmVar assignment |
| 4 | Set the accuracy bar at the field ceiling (Cyrius 96.5% untrained / 99.3% trained; beats Aldy 86.8%) | The v0.4 identity caller should target ≥90% hybrid-identity concordance to be worth shipping over "presence-only" | This is a full-diplotype concordance; the hybrid-only subset bar should be derived, not assumed |
| 5 | Accept the *36 gene-conversion blind spot as a documented residual | A pure embedded exon-9 conversion may produce no CN change → undetectable by depth alone (needs read-level PSV base-counting) | Confirm whether read-level PSV allele-fraction (not just depth) is required, raising the build cost |

## Honest gaps
- **The exact 117 PSV genomic coordinates are NOT in the literature** — they are in the Cyrius GitHub `CYP2D6.json` config; PSV *curation* (the topic's first half) is a repo-config task, not a fetch.
- **No single source gives a per-allele PSV→identity lookup table**; identity = (breakpoint) × (fusion direction) × (CN profile), i.e. an algorithm — which is exactly why this is Cyrius-class, not a single-SNP tag (confirms the v0.3 scoping).
- Aldy's method internals not fetched (redirect + budget); its 86.8% concordance captured via the Cyrius comparison.
