# Essentiality decoder report card (standing trust surface)

Single-gene KO -> essential/non-essential, via the conserved-core R1 decoder. Per-organism honest
tier; **no aggregate headline**. E. coli validated by composition (labels walled); human = the
cross-organism TRANSFER AUROC on the citable BAGEL CEG2/NEG reference.

| organism | cell | tier | metric | validation |
|---|---|---|---|---|
| Escherichia coli K-12 | conserved-core v0.1 | AUROC_SCORED | AUROC 0.6952 vs null 0.5 (Goodall-TraDIS gold-standard, n=3783, 351 ess/3432 non, base rate 0.0928); sens 0.3732 spec 0.984 prec 0.7043 | real per-gene AUROC vs the Goodall 2018 mBio Table S1 gold-standard (CC-BY); high-precision moderate-recall -- catches the universal core, misses the E. coli-specific essential tail (the E3 learned-complement target) |
| Homo sapiens | cross-organism transfer (E4) | TRANSFER_SCORED | AUROC 0.5805 vs null 0.50 (BAGEL CEG2 n=681 / NEG n=899); sens 0.1571 spec 0.9978 | universal core (ribosome/tRNA-synth/translation/polymerase) transfers cross-kingdom at high precision; human-specific core (proteasome 0/53, spliceosome 0/49) MISSED -> per-organism catalogue extension is the follow-on |
| Escherichia coli | transfer ladder rung (depth 8 from E. coli) | COVERAGE_SCORED | coverage_lift 0.3125 (null p95 0.0298 / MAX 0.0518); coverage(ess) 0.3789, coverage(non) 0.0664; phrasing-adjusted 0.3213; AUROC 0.6952 is SECONDARY and NOT cross-rung comparable | decoder applied UNCHANGED; coverage_lift is class-conditioned so it is base-rate robust and cross-rung comparable, which AUROC is not. Label technology transposon_insertion; class sourcing one_file_two_columns |
| Pseudomonas aeruginosa | transfer ladder rung (depth 5 from E. coli) | WALL_SCHEMA_UNVERIFIED | none -- a walled rung carries no coverage number by construction | route tried: W0 record says unverified: not in D:\dna_decode_cache\essentiality; remote is NOT RECORDED -- reachability is a HYPOTHESIS until fetched |
| Staphylococcus aureus | transfer ladder rung (depth 2 from E. coli) | WALL_SCHEMA_UNVERIFIED | none -- a walled rung carries no coverage number by construction | route tried: W0 record says unverified: not in D:\dna_decode_cache\essentiality; remote is NOT RECORDED -- reachability is a HYPOTHESIS until fetched |
| Saccharomyces cerevisiae | transfer ladder rung (depth 1 from E. coli) | WALL_SCHEMA_UNVERIFIED | none -- a walled rung carries no coverage number by construction | route tried: W0 record says unverified: not in D:\dna_decode_cache\essentiality; remote is http://sgd-archive.yeastgenome.org/curation/literature/phenotype_data.tab |
| Homo sapiens | transfer ladder rung (depth 1 from E. coli) | COVERAGE_SCORED | coverage_lift 0.1633 (null p95 0.0214 / MAX 0.0446); coverage(ess) 0.1689, coverage(non) 0.0056; phrasing-adjusted 0.1787; AUROC 0.5805 is SECONDARY and NOT cross-rung comparable | decoder applied UNCHANGED; coverage_lift is class-conditioned so it is base-rate robust and cross-rung comparable, which AUROC is not. Label technology crispr_ko; class sourcing two_files |

## Cross-organism transfer ladder

Primary metric is **`coverage_lift`** = coverage(essential) - coverage(non-essential),
which is class-conditioned and therefore base-rate robust and cross-rung comparable.
AUROC is a SECONDARY and is **not** compared across rungs (different sampling frames:
BAGEL is two curated extremes at base rate 0.431 vs a genome-wide 0.0928).

Verdict: **`INDETERMINATE_INSUFFICIENT_RUNGS`** - 2 scored rung(s) < min_rungs=3
Scored 2 / walled 3. Frozen bar {"plateau_tol": 0.06, "cliff_drop": 0.12, "min_rungs": 3} - a DESCRIPTIVE-CONSISTENCY LOCK, not a
frozen-before-the-numbers endpoint test (it is derived from the two rungs that
motivated the question). Full detail + the five honest limits:
`wiki/essentiality_transfer_ladder_2026-10-06.md`.

| organism | rung | tier | coverage_lift |
|---|---|---|---|
| Escherichia coli | depth 8 | COVERAGE_SCORED | 0.3125 |
| Pseudomonas aeruginosa | depth 5 | WALL_SCHEMA_UNVERIFIED | - |
| Staphylococcus aureus | depth 2 | WALL_SCHEMA_UNVERIFIED | - |
| Saccharomyces cerevisiae | depth 1 | WALL_SCHEMA_UNVERIFIED | - |
| Homo sapiens | depth 1 | COVERAGE_SCORED | 0.1633 |

## Honest scope
- The conserved-core decoder is the R1 PRIOR: high-precision, conservative-recall; captures the
  UNIVERSAL essential core, misses lineage-specific core (the R2/per-organism-catalogue target).
- E. coli per-gene AUROC + a learned E3 complement are gated on gold-standard labels (see
  `wiki/essentiality_label_wall_2026-07-28.md`); human labels (BAGEL CEG2/NEG) ARE available.
- Regenerate: `scripts/build_essentiality_report_card.py` (needs D: gene_info + BAGEL sets).