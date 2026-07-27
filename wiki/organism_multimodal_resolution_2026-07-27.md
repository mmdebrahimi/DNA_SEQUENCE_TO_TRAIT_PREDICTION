# Organism-multimodal — resolution of the "acquire dbGaP/UKB anyway" fork (2026-07-27)

**Directive:** "figure out a way to resolve #3 — pursue organism-multimodal anyway; needs a dbGaP/UKB paired
dataset acquisition = the same authority/money fork as unlocking R3."

**Resolution (decision-grade): the dbGaP/UKB acquisition is the WRONG lever. It is neither NEEDED (a free,
open, individual-level paired substrate exists) nor SUFFICIENT (the 2025-2026 SOTA shows more scale/access
does not move the result). The organism-multimodal question is answerable for $0 — and both the free-repro
answer and the field's newest result confirm the R3 closed negative. The authority/money fork I previously
flagged was a FRAMING ARTIFACT.**

## The R2 pre-bar framing check that flipped it

My prior memos (Family A falsifier, Family C synthesis) asserted: "the decisive organism-phenotype test needs
dbGaP/UKB controlled-access." That is an ASSERTED wall. A cheap real-surface scan falsifies it on two axes:

### Axis 1 — the substrate is FREE, not controlled-access
**GEUVADIS** (Lappalainen 2013, *Nature*) is a fully OPEN paired individual-level genotype→molecular-phenotype
dataset: 465 individuals × LCL RNA-seq × 1000-Genomes genotypes, across 5 populations (CEU/FIN/GBR/TSI/YRI).
Access is public: RNA-seq at ENA **ERP001942**, processed expression at ArrayExpress **E-GEUV-1 / E-GEUV-3**,
genotypes from 1000G Phase 3. This is exactly the paired substrate the "needs dbGaP/UKB" framing claimed was
walled — and it is on the open internet, laptop/Kaggle-feasible, $0. (openSNP is a second free option for
visible traits; the project already uses it.)

### Axis 2 — scale/access is not the binding constraint (N does not fix the failure)
The failure the multimodal ambition runs into is well-mapped by the field's own SOTA on THIS exact question:
- **Nat Genet 2023** ("Personal transcriptome variation is poorly explained by current genomic DL models"):
  Enformer/Basenji2/ExPecto/Xpresso predict expression across genes but FAIL cross-individual — often the
  wrong direction of effect for cis-regulatory variants. DL ties/loses to linear.
- **Genome Biology 2025 + Variformer 2026** (the newest attempt to fix it): fine-tuning a DL model on
  PERSONAL genomes DOES correct the cross-individual shortcoming — but only up to "similar performance and
  limitations of state-of-the-art LINEAR models," and it "does not learn a regulatory grammar that
  generalizes to unseen loci." **Even personal-genome fine-tuning only MATCHES elastic-net/PrediXcan.**

The ceiling is the **linear cis-eQTL model**, not the sample size. UKB's 500k-vs-GEUVADIS's 465 buys N, and N
does not convert a tie-with-linear into a beat-linear (the Arabidopsis R3 negative had n=1003 within-ancestry
and still failed). So acquiring controlled data buys a re-confirmation at high cost + an authority
application, NOT a resolution.

## What "resolving organism-multimodal" actually yields

The multimodal north-star ("DNA-encoder + 2nd modality → phenotype, does the DNA arm add value") inherits the
result above at the organism level:
- If modality-2 = measured expression, the DNA arm adds only the cis-eQTL-predictable part — redundant with
  the measured modality (ties linear).
- If modality-2 = image/organism trait, "DNA adds signal beyond modality-2" IS the R3 organism-phenotype
  question — population-structure-confounded, the project's thrice-confirmed 0-for-5 de-confounded negative.

**The DNA arm demonstrably adds NON-redundant signal in exactly ONE regime: molecular / variant-effect** (the
Family A finding — fine-mapped eQTL sign auROC 0.80; the forward cell). That is Family C, already shipped. At
the organism level the DNA arm ties linear, free-demonstrably.

## The honest options (the fork the user actually has)

1. **Bank organism-multimodal as resolved (recommended).** The acquisition is the wrong lever; the answer is
   the R3 closed negative, now confirmed by the July-2026 SOTA on the exact task, on free-reproducible data.
   No money, no controlled-access application, no build needed. The multimodal north-star is realized in the
   molecular regime (Family C) and closed in the organism regime — for free, not for want of data.

2. **$0 project-owned confirmation on GEUVADIS (reversible, free, if a project-owned number is wanted).** A
   laptop/Kaggle-feasible test: DNA-encoder arm vs elastic-net cis-eQTL vs both, de-confounded WITHIN
   population (the 5 GEUVADIS pops), reporting each within-group vs its own null (the population-structure
   discipline). Expected outcome: the DNA arm ties elastic-net (per the SOTA) — i.e. CONFIRMATORY of the
   negative, not resolving-positive. Worth it only if the project wants to own the number rather than cite
   Variformer + Nat Genet 2023. This is an authorize-able reversible build, NOT auto-run (confirming a
   literature-established negative is motion, not signal — the user decides if the owned number is worth the
   hours).

3. **Acquire dbGaP/UKB — NOT recommended.** It neither unblocks nor resolves; it re-confirms the free-data
   negative at the cost of money + an authority application. Reserve controlled-access acquisition for a
   question where free data genuinely cannot answer it — this is not one.

## H8 falsification rail

Key published claim: **"acquiring dbGaP/UKB does not resolve organism-multimodal — it re-confirms the R3
negative at high cost."** This SURVIVED a cheap falsification attempt: the scan actively searched for a
DNA-encoder arm that BEATS linear cross-individual given more scale/controlled data, and found the OPPOSITE —
the newest fine-tuning SOTA (Variformer 2026) ties linear and does not generalize. NOT "proven" (no owned
reproduction was run); literature-grounded on the field's own SOTA for the exact task, plus the project's
three de-confounded R3 negatives.

## Sources
- GEUVADIS: Lappalainen et al., Nature 501:506 (2013); ENA ERP001942; ArrayExpress E-GEUV-1/E-GEUV-3.
- Nat Genet 2023: "Personal transcriptome variation is poorly explained by current genomic deep learning
  models" (s41588-023-01574-w).
- Genome Biology 2025 / Variformer 2026: "Deep-learning prediction of gene expression from personal genomes"
  (s13059-025-03926-7) + the bioRxiv benchmarking (2023.03.16.532969).
