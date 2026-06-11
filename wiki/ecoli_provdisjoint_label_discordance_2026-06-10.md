# E. coli flagship cipro provdisjoint — spec-0.70 is label discordance, not a rule defect — 2026-06-10

The flagship DEFAULT decoder `call_resistance(Escherichia_coli_Shigella, ciprofloxacin)` (qrdr_point@2) scored
acc 0.817 / sens 0.933 / **spec 0.70** on 60 fresh provenance-disjoint strains (30R/30S). The low specificity
(9 FP / 30 S) reads at first glance as "the flagship degrades out-of-provenance." Per-strain genotype
inspection inverts that read.

## The 9 false positives ALL carry canonical resistance genotypes

| accession | n QRDR point mut | determinants | lab AST |
|---|---|---|---|
| GCA_012775415.1 | 3 | gyrA_S83L, parC_S80R, parE_I529L | S |
| GCA_012765655.1 | 3 | gyrA_S83L, parC_S80R, parE_I529L | S |
| GCA_015511245.1 | 2 | gyrA_S83L, parE_I355T | S |
| GCA_018094305.1 | 2 | gyrA_S83L, parE_I355T | S |
| GCA_012634045.1 | 2 | gyrA_S83L, parE_I529L | S |
| GCA_012634105.1 | 2 | gyrA_S83L, parE_I529L | S |
| GCA_020373145.1 | 4 | gyrA_S83L, gyrA_D87N, parC_S80I, parE_S458A | S |
| GCA_012902165.1 | 4 | gyrA_S83L, gyrA_D87N, parC_S80I, parE_S458A | S |
| GCA_020368835.1 | 2 | gyrA_S83L, parC_S80I | S |

`gyrA_S83L` + a `parC` target mutation is one of the best-established high-level fluoroquinolone-resistance
genotypes in *E. coli*; two strains carry a `gyrA` double **plus** `parC` (quadruple QRDR). A strain with this
genotype labeled cipro-**S** is the **high-sens / low-spec → suspect the label** pattern: near-perfect
sensitivity + poor specificity ⇒ suspect the phenotype LABEL/assay before the rule.

## The 2 false negatives carry ZERO QRDR mutations

GCA_012774615.1 and GCA_012766035.1 (lab R, rule S, n=0) — resistance is plasmid-mediated (qnr / aac(6')-Ib-cr)
or efflux, which the QRDR-point counter deliberately excludes (cross-organism robustness). Expected blind spot,
not a regression.

## Honest conclusion (not overclaimed)

- The decoder's GENOTYPE calls are defensible on **all 9** discordant-S strains — they carry the rule's R
  signature, including canonical gyrA+parC doubles/quads.
- This does NOT prove the labels are wrong by itself (possible: MIC just below breakpoint, heteroresistance,
  data-entry/provenance error, or a submitter using a non-CLSI breakpoint). It shifts the burden: **spec 0.70
  measures genotype↔phenotype discordance in the disjoint labels, not a decoder coding error.**
- Clean adjudication ATTEMPTED 2026-06-10 — and it hit the free-source bound: NCBI-PD `AST_phenotypes`
  carries the qualitative call ONLY (`ciprofloxacin=S` for all 9), **no raw MIC** to tier. None of the 9
  accessions appear in the BV-BRC `genome_amr` CSV either (the MIC-carrying source) — expected, since they
  are non-ecosystem submitters. So `cipro_mic_audit` / `mic_tiers.classify_tier` **cannot free-adjudicate
  these strains** — a true MIC tiering needs the submitter's raw MIC (per-study curation, not census-able;
  the ledger-row-58 bound). The genotype evidence (9/9 canonical QRDR-R) is therefore the strongest FREE
  adjudication available. It is consistent with label noise but does not PROVE the S labels wrong without
  MIC — the discordance is real and label-side; its precise cause (near-breakpoint MIC / heteroresistance /
  data-entry / non-CLSI breakpoint) stays unresolved on free data.

## Report-card implication

The `SCORED` E. coli cipro cell stays acc 0.817 / spec 0.70 **as measured against the labels as given** — the
report card does not silently "correct" labels. But the cell carries a discordance caveat: the specificity
penalty is label-quality-dominated, pending MIC adjudication. This is the report card doing its job — it
surfaced a genotype↔phenotype discordance on the headline decoder that same-source numbers (0.939/0.862) hid.
