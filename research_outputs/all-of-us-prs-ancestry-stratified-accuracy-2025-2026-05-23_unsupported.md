# All of Us research program polygenic risk score ancestry-stratified accuracy (2025) — unsupported claims (V1 invocation)

> Slug: all-of-us-prs-ancestry-stratified-accuracy-2025-2026-05-23. Captured 2026-05-23. These rows were rejected during V1 intake.
> Reasons fall into: missing locator, mapping-floor failure, hard-reject banned phrase, low-confidence label.

## Rejected rows

| Row content | Rejection reason | Suggested follow-up |
|---|---|---|
| Atrial Fibrillation PRS — European AUROC (AoU validation) = 0.646 (Jurgens S et al., PMC12622184, Figure 2) | **Audit floor failure (locator #4 verbatim quoted excerpt).** Quote field reads "(same paper Table S3 ancestry stratification)" — meta-reference, not a verbatim excerpt containing the value 0.646. Mechanical audit-floor check requires the quote to contain the numeric value. | Re-extract the actual verbatim AUROC sentence from Figure 2 / Table S3 of the Research Square preprint (e.g., a sentence stating "AUROC = 0.646 in European-ancestry individuals"). Then resubmit with a clean quote field. |
| Atrial Fibrillation PRS — Asian ancestry AUROC = 0.637 (Jurgens S et al., PMC12622184, Figure 2) | **Audit floor failure (locator #4 verbatim quoted excerpt).** Quote field reads "(same paper)" — not verbatim, does not contain 0.637. (Note: the paired OR/SD row for Asian ancestry passed intake because its quote read "Asian (OR/SD = 1.76; AUROC = 0.637)" — that quote contains the AUROC value 0.637 also, so the row could be re-validated by reusing that quote string for this row.) | Reuse the verbatim quote "Asian (OR/SD = 1.76; AUROC = 0.637)" from the paired row as the locator-quote for this row, with rationale explicitly stating the quote yields AUROC 0.637. |
| Atrial Fibrillation PRS — Admixed American AUROC = 0.595 (Jurgens S et al., PMC12622184, Figure 2) | **Audit floor failure (locator #4 verbatim quoted excerpt).** Quote field reads "(same paper)" — not verbatim. | Reuse the verbatim quote "Admixed American (1.45; 0.595)" from the paired OR/SD row; rationale must state that the quote yields AUROC 0.595. |
| Atrial Fibrillation PRS — African ancestry AUROC = 0.573 (Jurgens S et al., PMC12622184, Figure 2) | **Audit floor failure (locator #4 verbatim quoted excerpt).** Quote field reads "(same paper)" — not verbatim. | Reuse the verbatim quote "African ancestry groups (1.39; 0.573)" from the paired OR/SD row; rationale must state that the quote yields AUROC 0.573. |

## Summary

- Total rejected: 4
- Reason breakdown:
  - Missing audit-floor locator(s): 4 (all four are locator #4 — non-verbatim "(same paper)" quote-field shorthand on the AUROC pair-rows of the AF PRS ancestry stratification)
  - Mapping-floor failure: 0
  - Hard-reject banned phrase: 0
  - Low confidence: 0

## Lesson for future /research orchestrator runs

When extracting paired metrics (OR/SD + AUROC) from the same source-text sentence, do NOT shorthand the second row's quote-field as "(same paper)" — the audit-floor check is row-local and treats each row independently. Either reuse the full verbatim sentence on both rows, or split the source-quote across the two rows such that each row's quote-field contains the value claimed in that row. This is the same locator-shorthand anti-pattern that the V1 intake checklist's audit-floor was designed to catch.
