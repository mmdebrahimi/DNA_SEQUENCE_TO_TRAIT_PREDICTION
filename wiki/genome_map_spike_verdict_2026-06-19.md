# Genome-map v1 spike — GO/NO-GO verdict (2026-06-19)

**Tiering verdict (Bakta honesty re-tiering): GO**
**Overlay-integrity verdict (determinant->feature join): GO**

- G1 satisfied on >=1 genome; G2 clean; no all-symbol-fallback genome
- determinant-overlay integrity demonstrated on 2 AMR-bearing genome(s)

Honesty contract: phenotype claims appear ONLY behind a high-confidence determinant join (symbol-fallback excluded); the unknown rate is DB-labelled `unknown_under_bakta_db_light` (db-light = reduced functional coverage, not biological unknown).

## GCA_002180195.1 — E. coli ST131 (cipro-R)
- organism (-O): `Escherichia`
- total features: 5185
- per-tier counts: determinant-phenotype=23, curated-molecular-function=4444, homology-only-hypothesis=400, unknown=318
- `unknown_under_bakta_db_light`: 0.061
- join quality: n_main_rows=32 high_confidence=32 symbol_fallback=0 unjoined=0
- all_joins_symbol_fallback: False
- determinant-phenotype features: 23
- genome-level R/S calls (separate from features): ciprofloxacin=R, ceftriaxone=R, tetracycline=R, gentamicin=R
- **G1** prevent-wrong-inference: 423 features (demote=400, surface=23) -> g1_pass=True
    - DEMOTE: raw `Aldehyde oxidase/xanthine dehydrogenase a/b hammerhead domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `putative succinyl-diaminopimelate desuccinylase` -> homology-only-hypothesis (low-confidence wording ('putative'))
    - DEMOTE: raw `Small putative membrane protein` -> homology-only-hypothesis (low-confidence wording ('putative'))
    - DEMOTE: raw `PPC domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Putative pre-16S rRNA nuclease` -> homology-only-hypothesis (low-confidence wording ('putative'))
- **G2** no-tier-confusion: pass=True (violations=0)
- **overlay-integrity**: GO (main_rows=32, high_conf_joins=32, surfaced=23)
- per-genome tiering verdict: GO

## GCA_000417105.1 — K. pneumoniae (carbapenem-R)
- organism (-O): `Klebsiella_pneumoniae`
- total features: 5427
- per-tier counts: determinant-phenotype=11, curated-molecular-function=4329, homology-only-hypothesis=618, unknown=469
- `unknown_under_bakta_db_light`: 0.086
- join quality: n_main_rows=12 high_confidence=11 symbol_fallback=1 unjoined=0
- all_joins_symbol_fallback: False
- determinant-phenotype features: 11
- genome-level R/S calls (separate from features): meropenem=R, ciprofloxacin=R, ceftriaxone=R, gentamicin=S, tetracycline=S
- **G1** prevent-wrong-inference: 629 features (demote=618, surface=11) -> g1_pass=True
    - DEMOTE: raw `CSD domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Major facilitator superfamily (MFS) profile domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Putative 2-Ketogluconate kinase` -> homology-only-hypothesis (low-confidence wording ('putative'))
    - DEMOTE: raw `HTH lacI-type domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Autotransporter outer membrane beta-barrel domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
- **G2** no-tier-confusion: pass=True (violations=0)
- **overlay-integrity**: GO (main_rows=12, high_conf_joins=11, surfaced=11)
- per-genome tiering verdict: GO

## GCF_000171775.1 — Gemmata obscuriglobus (homology stress test)
- organism (-O): `None`
- total features: 9479
- per-tier counts: determinant-phenotype=0, curated-molecular-function=2032, homology-only-hypothesis=1009, unknown=6438
- `unknown_under_bakta_db_light`: 0.679
- join quality: n_main_rows=0 high_confidence=0 symbol_fallback=0 unjoined=0
- all_joins_symbol_fallback: False
- determinant-phenotype features: 0
- genome-level R/S calls (separate from features): (none)
- **G1** prevent-wrong-inference: 1009 features (demote=1009, surface=0) -> g1_pass=True
    - DEMOTE: raw `N-acetyltransferase domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Integrase catalytic domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `DUF3185 domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `HTH marR-type domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
    - DEMOTE: raw `Transposase IS4-like domain-containing protein` -> homology-only-hypothesis (low-confidence wording ('domain-containing protein'))
- **G2** no-tier-confusion: pass=True (violations=0)
- **overlay-integrity**: n/a (determinant-free) (main_rows=0, high_conf_joins=0, surfaced=0)
- per-genome tiering verdict: GO
