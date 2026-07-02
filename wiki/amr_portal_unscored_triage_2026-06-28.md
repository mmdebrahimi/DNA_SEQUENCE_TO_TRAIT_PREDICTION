# EBI AMR Portal — unscored powered-cell triage board (2026-06-28)

The corrected frontier (`wiki/frontier_reassessment_2026-06-28.md`): the AMR Portal has **79 provenance-
disjoint POWERED cells** (≥10 R AND ≥10 S, measured AST, accession-disjoint vs CRyPTIC + our tuning
cohorts). **27 are scored** (`wiki/amr_portal_independent_report_card.md`, 23 SCORED_INDEPENDENT + 4
underpowered-at-strict + Campylobacter added 2026-06-28). That leaves **52 unscored powered cells.** This is
the ranked board for what to do with them — NOT a "score them all" sweep.

**Load-bearing guardrail (the whole reason this is a triage, not a loop):** do NOT spray the generic,
E. coli-derived `DRUG_RULE` across an organism just because the AMR Portal has measured labels for it. Each
organism×drug needs an ENDORSED rule (deployed default validated for that organism, or a curated
organism-specific registry config) + its own independent powering check. Applying the E. coli rule to, say,
Acinetobacter or Enterococcus is the intrinsic-gene over-call trap (`feedback_intrinsic_genes_break_broad_amr_class_rules`).

> **Tier 1 DONE (2026-07-01):** TMP-SMX x5 scored on the AMR Portal -- **4/5 SCORED** (E. coli acc 0.926, Salmonella 0.963, Shigella sonnei 0.874, Shigella flexneri 0.961; strata-reproduced), **Klebsiella INDETERMINATE** (strata gate correctly fired -- sul-AND-dfr pattern doesn't hold there). Namespace-separate artifact wiki/amr_portal_tmpsmx_experimental_2026-07-01.json + wired into the AMR-Portal report card as a branded EXPERIMENTAL section (deployed-surface counts byte-unchanged). Frozen AMR surface untouched.

## Tier 1 — has a REAL rule, just not scored on the AMR Portal yet (the quick wins)

| cells | rule status | action | VOI |
|---|---|---|---|
| **TMP-SMX ×5** — E. coli (2625R/9277S), Klebsiella (2827R/1384S), Salmonella (1667R/24955S), Shigella sonnei (796R/343S), Shigella flexneri (138R/69S) | the **EXPERIMENTAL overlay** exists (`dna_decode/data/experimental_drug_rules.py`, the `sul AND dfr` rule, validated on Oxford+Sci234) — NON-frozen, branded `EXPERIMENTAL_SCORED` | **wire the experimental TMP-SMX scorer to the AMR Portal genotype parquet** (its own namespace-separate `external_validation_*` artifact, NOT the frozen card — the shared-key silent-overwrite trap). Highest-VOI code move on the board. | **HIGH** |

## Tier 2 — ABSTAINS-BY-DESIGN but powered (no action — honest non-cells)

| cells | why |
|---|---|
| Pseudomonas aeruginosa meropenem (718R/503S), Acinetobacter baumannii meropenem (761R/368S) | the deployed rule is `EXPRESSION_FLOOR`/ABSTAIN (intrinsic OXA/oprD-loss carbapenem resistance is not gene-presence-decodable). Powered labels exist, but the honest decoder REFUSES — scoring would force a wrong call. Leave as ABSTAINS_BY_DESIGN; document, don't score. |

## Tier 3 — org has a rule for ONE drug only → the other drugs are guardrail-blocked (new-rule work)

| cells | why blocked |
|---|---|
| Campylobacter jejuni/coli **tetracycline** (2600R/3290S, 923R/844S) + **gentamicin** (80R/5711S, 92R/1671S) | Campylobacter has ONLY the endorsed `Campylobacter\|ciprofloxacin` registry rule (qrdr_point). The generic DRUG_RULE for tet (`tetO`/`tetM` ribosomal-protection — different from the E. coli efflux `tetA` set) + gent is NOT validated for Campylobacter. Needs a curated Campylobacter tet/gent rule first. |

## Tier 4 — NEW-RULE cells (the bulk: ~42 cells across 12 genera) — per-cell research, NOT a sweep

Each genus needs a curated, organism-specific determinant rule + an independent powering check before it can
join the card. Ranked by cell-count (breadth of the labelled surface available):

| genus | cells | drugs with powered labels | note |
|---|---|---|---|
| Enterobacter | 10 | cipro, gent, mero, tet, TMP-SMX | Enterobacterale; may partly reuse E. coli rules BUT has intrinsic AmpC (cef over-call risk) — verify per drug |
| Staphylococcus | 5 | cipro, gent, rif, tet, TMP-SMX | Gram-positive; needs S. aureus-specific determinants (the oxacillin/mecA + cefoxitin-surrogate lesson lives here) |
| Acinetobacter | 5 | cef, cipro, gent, tet, TMP-SMX | intrinsic-gene heavy (OXA-51) — the canonical over-call organism; curated/acquired-only rules essential |
| Streptococcus pneumoniae | 4 | cef, mero, tet, TMP-SMX | Gram-positive; β-lactam resistance is PBP-mutation (not gene-presence) — hard for a determinant caller |
| Neisseria gonorrhoeae | 4 | cef, cipro, tet | mosaic penA / gyrA — partly point-mutation (decodable) |
| Enterococcus | 3 | cipro, gent, tet | intrinsic aminoglycoside resistance — high-level-gent needs specific markers |
| Proteus | 3 | cef, gent, TMP-SMX | Enterobacterale; intrinsic tet/colistin resistance |
| Serratia | 3 | cef, cipro, mero | intrinsic AmpC (SME) |
| Pseudomonas (non-mero) | 2 | cipro, gent | efflux/porin heavy |
| Clostridioides difficile | 1 | cef | anaerobe |
| Haemophilus influenzae | 1 | TMP-SMX | |
| (Mycobacterium cipro/inh/rif | 3 | — | **already scored** via the TB WHO-catalogue pipeline, NOT the AMRFinder scorer — exclude from this board) |

## Recommended sequence
1. **Tier 1 (TMP-SMX ×5)** — the one HIGH-VOI quick win: a real (experimental) rule + powered labels; wire
   it to the AMR Portal parquet as a namespace-separate experimental artifact. The single best next code move.
2. **Tier 2** — document the 2 abstainers on the card as ABSTAINS_BY_DESIGN-powered (no scoring).
3. **Tier 3 + 4** — a curated-rule research backlog, one organism×drug at a time, each gated by the
   intrinsic-gene guardrail + an independent powering check. NOT an autonomous sweep; each is a deliberate
   per-cell effort (and Gram-positive / PBP-mutation / intrinsic-AmpC organisms may simply not be
   gene-presence-decodable — an honest ABSTAIN is a valid outcome).

## Source
`wiki/amr_portal_feasibility.json` (79 powered cells), `wiki/amr_portal_independent_scores.json` (27 scored),
`scripts/amr_portal_score_independent.py` (the scorer + scope contract), `dna_decode/data/experimental_drug_rules.py`
(the TMP-SMX overlay), `dna_decode/data/calibrated_amr_rules.json` (the registry).
