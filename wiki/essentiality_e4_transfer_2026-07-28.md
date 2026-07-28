# Essentiality E4 — cross-organism transfer (E. coli -> human, 2026-07-28)

**The "any organism" test:** apply the E. coli-tuned conserved-core decoder UNCHANGED to human genes and
measure whether it predicts human essentiality. Uses the citable **BAGEL** human reference (Hart lab):
CEGv2 core-essential vs NEGv1 non-essential — freely fetchable (GitHub), unlike DepMap (portal bot-walled).

## Result

| | value |
|---|---|
| **Transfer AUROC** (E.coli decoder -> human CEG2/NEG) | **0.580** (null 0.50) |
| sens / spec at essential threshold | **0.157 / 0.998** |
| reference | BAGEL CEGv2 (681 essential) vs NEGv1 (899 non-essential), mapped by Entrez ID |

## Finding — the universal core transfers; the lineage-specific core does not

- **Above chance (0.58 > 0.50): the decoder DOES transfer cross-kingdom** — the universal essential core
  (ribosome, tRNA-synthetases, translation factors EEF1A1/EEF2/EIF2, RNA/DNA polymerase) is shared E. coli↔
  human and the decoder recovers it in human (sample: AARS1, DARS1, EEF1A1, EEF2, EIF2S1).
- **HIGH-PRECISION, LOW-RECALL (spec 0.998, sens 0.157):** when the bacterial-tuned decoder fires on a human
  gene it is almost always a true essential (the universal core), but it captures only ~16% of human
  essentiality — because the **human-specific essential core is absent from the E. coli catalogue**:
  proteasome (0/53 caught), spliceosome/splicing (0/49 caught). These are essential in human, not in the
  bacterial function catalogue.
- **-> the E4 conclusion the plan predicted:** the UNIVERSAL core is transferable cross-organism; the
  LINEAGE-SPECIFIC core needs a per-organism catalogue extension (the R2/E3 target). "Any organism" works for
  the shared core out-of-the-box; lineage-specific essentiality is per-organism.

## Bonus: human labels are NOW available (BAGEL)

BAGEL CEG2/NEG is a real, citable human essential/non-essential label set — so the HUMAN arm of v0.1 (real
AUROC) and E3 (train the learned complement on human) is UNBLOCKED (unlike E. coli, whose gold-standard
labels remain walled -- `wiki/essentiality_label_wall_2026-07-28.md`). A human-tuned conserved-core catalogue
(add proteasome/spliceosome/etc.) would raise the transfer recall substantially — the natural v0.1 move.

## Honesty (H8)
- AUROC 0.58 is MODEST-but-real transfer, NOT a strong predictor — it reflects the universal-core overlap
  only; the low recall is the honest signal that lineage-specific essentiality doesn't transfer.
- The decoder was applied to human WITHOUT retuning (a fair transfer test); human-specific tuning is the
  named follow-on, not a claimed result.

## Reproduce
`scripts/build_essentiality_report_card.py` (transfer AUROC + report card; needs D: gene_info + BAGEL sets).
