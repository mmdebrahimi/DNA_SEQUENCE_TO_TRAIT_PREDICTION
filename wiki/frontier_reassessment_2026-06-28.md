# Frontier reassessment — the "exhausted / send-the-emails" verdict was STALE (2026-06-28)

A `/brainstorm` adversarial review of the 2026-06-28 "autonomous frontier is exhausted; the only high-VOI
lever is user-authority acquisition (send the TB author emails)" conclusion **overturned it**. The June-28
session anchored on the 2026-06-19 ledger frame + the 2026-06-13 reproducibility freeze — both of which
**predate** the work that changes everything: the **EBI AMR Portal** validation arm (2026-06-23). This is the
project's own recurring failure mode (stale-state / wrong-denominator / under-search), caught before banking.

## What the stale verdict missed (all verified in-repo)
- **The label wall is already broken — for FREE.** `wiki/amr_portal_feasibility_result_2026-06-23.md`: the EBI
  AMR Portal (1.71 M measured-AST rows, free FTP) yields **74 provenance-disjoint POWERED cells** (≥10 R/S),
  measured wet-lab AST, accession-disjoint vs CRyPTIC + our tuning cohorts.
- **The independent TB number ALREADY EXISTS, free, no DUA, no author contact.**
  `wiki/tb_independent_number_2026-06-23.md`: RIF acc 0.937 / INH acc 0.914, N=2,845, WHO rule unchanged.
  → **The TB author emails (`wiki/tb_goldset_author_emails_2026-06-22.md`) are MOOT. Do NOT send them.**
- **21→23 bacterial cells already SCORED_INDEPENDENT** on a standing card
  (`wiki/amr_portal_independent_report_card.md`): E. coli / Salmonella / Klebsiella / Shigella / now
  Campylobacter, measured AST, 0.83–0.995 acc.
- **Track A's "no free isolate-level measured-phenotype DB beyond the Stanford family"
  (`wiki/free_label_source_scan_2026-06-28.md`) is WRONG** — the EBI AMR Portal is exactly that, and was
  already in the repo + exploited 5 days before the scan. The scan missed it (under-search).

## The corrected verdict
**The frontier is NOT exhausted. The live frontier is the EBI AMR Portal: 74 powered cells, ~23 scored.**
The binding constraint thesis ("labels not models, only acquisition unblocks independence") was TRUE at the
2026-06-13 freeze but was SUPERSEDED on 2026-06-23. Independence is now FREE and code-closable for the cells
that have a deployed rule.

## What was actually done this session (recs 1-2, post-correction)
- **Campylobacter ciprofloxacin scored** on the AMR Portal (rec 1): C. jejuni 1533R/4358S acc 0.981,
  C. coli 437R/1331S acc 0.995 — both SCORED_INDEPENDENT, calibrated qrdr_point rule. Card 21→23 / 27.
- **Scope contract documented + routing data-driven** (rec 2): the scorer's `CELLS` now has an explicit
  scope-contract comment (5 full-drug organisms cartesian + Campylobacter CIPRO-ONLY — the generic DRUG_RULE
  is NOT sprayed across organisms without an endorsed rule); the card's calibrated-routing label is now
  data-driven from the CALIBRATED registry keys (was a hardcoded 2-cell set that mislabeled Campylobacter).

## The forward frontier (the ~53 unscored powered cells) — triage
See `wiki/amr_portal_unscored_triage_2026-06-28.md`. The bulk are NEW organism×drug cells with NO deployed
rule (Staph, Neisseria, Strep pneumo, Enterococcus, Enterobacter, Serratia, Proteus, C. difficile,
H. influenzae, …). **Guardrail (load-bearing): do NOT spray the generic E. coli-derived DRUG_RULE across
organisms just because the AMR Portal has labels** — that is the intrinsic-gene over-call trap
(`feedback_intrinsic_genes_break_broad_amr_class_rules`). Each new cell needs a curated, organism-specific
rule + an independent powering check; that is a per-cell research effort, not a card-expansion sweep.

## Honesty residuals (named, not fixed)
- **influenza NA + HCV/HBV** stay `NEEDS-VERIFICATION` (a content filter blocked viral-resistance web search
  this session). The shipped surface still renders them `NO_FREE_PHENOTYPE_SOURCE`; that label is
  conservative-but-unconfirmed. Do NOT treat it as a verified terminal until a filter-free primary-source
  check runs. (Not patched here: `shipped_decoder_surface.py` is in the prospective-lock + .gitattributes
  integrity-pinned set; changing it is a deliberate freeze amendment, not a doc fix.)
- The TB author-emails doc + the acquisition-strategy memo + the Track-A scan are top-bannered as
  SUPERSEDED/CORRECTED (this commit), not deleted (audit trail).

## Sources
`wiki/amr_portal_feasibility_result_2026-06-23.md`, `wiki/tb_independent_number_2026-06-23.md`,
`wiki/amr_portal_independent_report_card.md`, `wiki/cross_kingdom_validation_summary.md`,
`scripts/amr_portal_score_independent.py`, `scripts/build_amr_portal_report_card.py`.
