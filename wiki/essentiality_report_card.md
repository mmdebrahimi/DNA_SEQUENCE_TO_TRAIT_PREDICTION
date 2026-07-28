# Essentiality decoder report card (standing trust surface)

Single-gene KO -> essential/non-essential, via the conserved-core R1 decoder. Per-organism honest
tier; **no aggregate headline**. E. coli validated by composition (labels walled); human = the
cross-organism TRANSFER AUROC on the citable BAGEL CEG2/NEG reference.

| organism | cell | tier | metric | validation |
|---|---|---|---|---|
| Escherichia coli K-12 | conserved-core v0.1 | AUROC_SCORED | AUROC 0.6952 vs null 0.5 (Goodall-TraDIS gold-standard, n=3783, 351 ess/3432 non, base rate 0.0928); sens 0.3732 spec 0.984 prec 0.7043 | real per-gene AUROC vs the Goodall 2018 mBio Table S1 gold-standard (CC-BY); high-precision moderate-recall -- catches the universal core, misses the E. coli-specific essential tail (the E3 learned-complement target) |
| Homo sapiens | cross-organism transfer (E4) | TRANSFER_SCORED | AUROC 0.5805 vs null 0.50 (BAGEL CEG2 n=681 / NEG n=899); sens 0.1571 spec 0.9978 | universal core (ribosome/tRNA-synth/translation/polymerase) transfers cross-kingdom at high precision; human-specific core (proteasome 0/53, spliceosome 0/49) MISSED -> per-organism catalogue extension is the follow-on |

## Honest scope
- The conserved-core decoder is the R1 PRIOR: high-precision, conservative-recall; captures the
  UNIVERSAL essential core, misses lineage-specific core (the R2/per-organism-catalogue target).
- E. coli per-gene AUROC + a learned E3 complement are gated on gold-standard labels (see
  `wiki/essentiality_label_wall_2026-07-28.md`); human labels (BAGEL CEG2/NEG) ARE available.
- Regenerate: `scripts/build_essentiality_report_card.py` (needs D: gene_info + BAGEL sets).