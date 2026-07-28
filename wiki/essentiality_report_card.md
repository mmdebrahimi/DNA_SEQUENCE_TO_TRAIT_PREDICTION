# Essentiality decoder report card (standing trust surface)

Single-gene KO -> essential/non-essential, via the conserved-core R1 decoder. Per-organism honest
tier; **no aggregate headline**. E. coli validated by composition (labels walled); human = the
cross-organism TRANSFER AUROC on the citable BAGEL CEG2/NEG reference.

| organism | cell | tier | metric | validation |
|---|---|---|---|---|
| Escherichia coli K-12 | conserved-core v0 | COMPOSITION_VALIDATED | 208/4318 predicted essential (known essentialome ~300) | size + composition match the known essentialome (translation/envelope/replication); per-gene AUROC pending gold-standard labels (walled) |
| Homo sapiens | cross-organism transfer (E4) | TRANSFER_SCORED | AUROC 0.5805 vs null 0.50 (BAGEL CEG2 n=681 / NEG n=899); sens 0.1571 spec 0.9978 | universal core (ribosome/tRNA-synth/translation/polymerase) transfers cross-kingdom at high precision; human-specific core (proteasome 0/53, spliceosome 0/49) MISSED -> per-organism catalogue extension is the follow-on |

## Honest scope
- The conserved-core decoder is the R1 PRIOR: high-precision, conservative-recall; captures the
  UNIVERSAL essential core, misses lineage-specific core (the R2/per-organism-catalogue target).
- E. coli per-gene AUROC + a learned E3 complement are gated on gold-standard labels (see
  `wiki/essentiality_label_wall_2026-07-28.md`); human labels (BAGEL CEG2/NEG) ARE available.
- Regenerate: `scripts/build_essentiality_report_card.py` (needs D: gene_info + BAGEL sets).