# Yeast decoder — mechanistic attribution capstone (2026-07-02): ATTRIBUTION INCONCLUSIVE

> **RESOLVED 2026-07-02 (`wiki/yeast_cnv_attribution_result_2026-07-02.md`):** the open thread below is now
> CLOSED. Attribution was inconclusive *with gene presence/absence* because the canonical mechanisms are
> COPY-NUMBER (CUP1) — the wrong feature type. Given the RIGHT feature (the free `genesMatrix_CopyNumber`),
> a bounded confirmatory test (NOT the full SNP+LMM build) CONFIRMS canonical-gene attribution, de-confounded:
> **CUP1 copy → copper clade-centered ρ +0.73** (perm_p 0.005); **ENA5 copy → sodium ρ +0.25** (perm_p 0.005),
> mechanism-SPECIFIC (ionic Na/Li yes, non-ionic sorbitol null +0.03). So the failure below was a
> feature/mechanism MISMATCH, not a capability gap — the 2nd proof (after DepMap) that the feature type must
> match the mechanism type.

The capstone tested whether the de-confounded within-clade yeast growth signal is attributable to the KNOWN
canonical resistance genes (arsenite→ACR3/ARR, copper→CUP1, benomyl→TUB2). **It is NOT cleanly attributable.**
The win (real, permutation-clean, K-robust de-confounded signal) STANDS; its mechanistic identity is
UNRESOLVED and partly consistent with finer accessory-element structure. This tempers the earlier "MECHANISM"
label to "de-confounded predictive signal, canonical-gene attribution not confirmed."

## What ran
- Recomputed the **within-clade-residualized** top genes per condition (residualize both growth and gene
  presence by clade mean, then correlate) — the genes that actually drive the DE-CONFOUNDED signal (not the
  global-univariate top, which mixes between-clade structure).
- Extracted their pangenome ORF sequences (index-mapped `X<N>.name` ↔ fasta `<N>-name`, verified) and BLASTed
  vs the S288C ACR3-bearing chr XVI + CUP1-bearing chr VIII.

## What the within-clade top genes actually are (the honest catalog)
| condition | top within-clade genes (|r_within| up to) | canonical gene present? |
|---|---|---|
| arsenite | snap/augustus-masked ORFs, YDL247W.A, Q0160 (mito), EC1118 accessory ORFs (0.24) | **No** (ACR3/ARR not in top; one ORF maps chr XVI ~946kb but ACR3 locus unconfirmed) |
| copper | **2μm plasmid genes R0010W–R0040C**, Ty5 element, YJL218W/YJL222W, YER187W (0.31) | **No** (CUP1 absent; copper R is CUP1 COPY-NUMBER, invisible to presence/absence) |
| benomyl | putative iron/AA transporters, augustus-masked ORFs, YNR074C, YOL163W (0.43) | **No** (TUB2 not in top) |

## Why attribution failed (three real reasons, all named)
1. **Wrong feature type for the canonical mechanisms.** Copper resistance = CUP1 tandem **copy number** (CUP1
   is present in ~all strains → zero presence/absence signal). Arsenite = the ARR cluster, largely core/
   systematic-named. Gene PRESENCE/ABSENCE structurally cannot see dosage or SNP-level variation at core loci.
2. **The signal is carried by accessory-element carriage.** Copper's top within-clade genes are 2μm plasmid +
   Ty transposon genes — these co-vary with growth WITHIN clades but are markers of finer accessory-element
   sub-structure, **not** a copper-resistance mechanism. So part of the "within-clade" signal is a subtler
   structure confound that K=18/30 clade de-confounding did not fully remove.
3. **Reference-locus fetching was unreliable** (UniProt P39980 resolved to SIT1 not ACR3; the 61-aa CUP1 did
   not tblastn-locate) — so the one suggestive hit (an arsenite ORF at chr XVI ~946 kb) could NOT be confirmed
   as ACR3/ARR. Reported as suggestive-unconfirmed, not a win.

## Revised honest verdict on the yeast result
- **Stands:** the within-clade r² (0.15–0.37 on 11/35 conditions) is real, permutation-null-clean, and robust
  to K=30 — categorically unlike the embedding failures (Arabidopsis within-group −0.13). Gene content carries
  genuine within-clade predictive signal.
- **Tempered:** it is NOT attributed to the textbook resistance genes; label it **"de-confounded predictive
  signal (mechanistic identity unresolved)"**, NOT "MECHANISM confirmed." Some signal is likely finer
  accessory-element structure (2μm/Ty), not causal resistance loci.

## What would actually resolve attribution (the real next steps — a different feature set)
1. **SNP + copy-number features at the canonical loci** (the full `1011Matrix.gvcf.gz` + CUP1 copy number) —
   presence/absence is the wrong resolution for CUP1/ARR.
2. **Kinship mixed-model** (continuous relatedness matrix from the SNP distance) instead of discrete clades —
   removes the finer accessory-element sub-structure the 2μm/Ty genes track.
3. **Reliable canonical-locus references** (fetch ACR3 YPR201W / CUP1 YHR053C / TUB2 YFL037W CDS correctly)
   to test locus-level hits.

Frozen AMR surface byte-unchanged; this is a read-only analysis. The capstone's value is HONEST SCOPING: it
stops the project from overclaiming a mechanistic win on gene-presence features.
