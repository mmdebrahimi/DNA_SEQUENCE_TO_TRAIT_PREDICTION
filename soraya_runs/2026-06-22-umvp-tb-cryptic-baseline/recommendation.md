# recommendation — tb-cryptic-baseline
Baseline DONE + honest. Three named residuals (all surfaced in the result memo + CLAUDE.md):
- INDEPENDENT TB number (deliverable-b): hand-curate a post-2023 gold set (gated; not executor-closable).
- INDEL normalization (code-closable): add delins/LoF/indel matching to the parquet adapter to lift the
  lineage-collapsed lower bound (esp. INH katG_LoF). A real but larger build; do it only if the tighter
  number is wanted.
- CALLABILITY (ABSTAIN vs S-by-absence): the ~1.6 TB regeno run is now needed ONLY for this, not the
  baseline number.
Recommend banking the baseline here; the higher-value next move is the independent gold set (a user/data
decision), not more local indel-plumbing.
