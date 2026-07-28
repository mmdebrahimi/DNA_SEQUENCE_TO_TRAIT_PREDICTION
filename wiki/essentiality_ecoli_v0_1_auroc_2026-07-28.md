# Essentiality v0.1 — the real E. coli AUROC (2026-07-28)

**The gated number, delivered.** The conserved-core decoder now has a hard per-gene accuracy score against a
gold-standard label, replacing v0's composition-only validation. Label: **Goodall 2018 mBio Table S1** (TraDIS
genome-wide essentiality classification, CC-BY open access; user-downloaded, the DB routes were bot-walled).

## Result (n=3783 E. coli genes: 351 essential / 3432 non-essential; base rate 9.3%)

| metric | value |
|---|---|
| **AUROC** | **0.695** (null 0.5) |
| sensitivity @ threshold | 0.373 |
| specificity @ threshold | 0.984 |
| precision @ threshold | 0.704 |

## Reading it

- **AUROC 0.695 ≫ 0.5 null** — the conserved-core function catalogue has genuine predictive signal for E. coli
  essentiality (a label-independent decoder, validated now against real knockout data). This is the honest
  "beats base-rate" number: with essentials only 9.3% of genes, a naive "always non-essential" guesser is 91%
  *accurate* but 0.5 AUROC; the decoder's 0.695 is real signal above that floor.
- **High precision (0.70) + high specificity (0.98), moderate recall (0.37):** when the decoder predicts
  essential it is right ~70% of the time, but it catches only ~37% of essential genes — it captures the
  **universal core** (ribosome/replication/transcription/division) and misses the **E. coli-specific essential
  tail**. That tail is exactly the **E3 learned-complement** target (now trainable: this same Goodall label is
  the training set that was previously walled).

## What this unblocks
- **E3 is now trainable for E. coli** — the Goodall Table S1 is a genome-wide essential/non-essential training
  label. The E1 label wall (row 583) is resolved via the user-downloaded open-access supplement.
- The report card's E. coli row is upgraded COMPOSITION_VALIDATED → **AUROC_SCORED**.

## Honesty
- "Non-essential" = TraDIS-non-essential (the screen assessed all genes; "Unclear" rows excluded). AUROC is vs
  this single screen's calls — a strong gold standard, though essential-gene lists differ ~10% across screens
  (Keio/PEC/TraDIS); the 248-consensus would be a slightly cleaner (smaller) positive set.
- Data on D:/dna_decode_cache/essentiality/. Reproduce via `scripts/build_essentiality_report_card.py`.
