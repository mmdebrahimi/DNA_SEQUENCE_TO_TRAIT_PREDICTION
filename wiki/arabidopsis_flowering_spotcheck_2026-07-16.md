# Flowering cell vs measured AraPheno — the attempt to SCORE it, and why it can't be (yet) (2026-07-16)

Goal: turn the flowering cell from *faithful-to-literature* into a **SCORED** cell (the project's gold
standard) against free measured data. **Outcome: cannot be scored on free data this session** — for two
independent reasons, both worth recording precisely. `scripts/flowering_arapheno_spotcheck.py`.

## What worked

**The measured phenotype is free and live:** AraPheno DTF1 (days to visible flowering buds, 1001 Genomes) —
`HTTP 200`, **934 accessions** with values, CSV via the free REST API. The join harness is **built and runs**.

## Wall 1 — the allele mapping EXISTS but every free route is walled

A real score needs per-accession FRI/FLC allele calls for a powered cohort. That mapping **exists**:

| source | scope | route | result |
|---|---|---|---|
| Zhang 2020 *Plant J* Table S1 | **1016** accessions × FRI allele | Wiley | **HTTP 402 Payment Required** — a **money gate; not paid** |
| " (preprint) | same | bioRxiv supplementary | figures only (`media-1.pdf` = Suppl. Figs 1–2); **no Table S1** |
| Shindo 2005 *Plant Physiol* | 176 accessions: FRI haplotype + flowering time + FLC level | PMC | **CAPTCHA-blocked** |
| Werner 2005 *Genetics* | 145 accessions genotyped for both common FRI lesions | PMC / OUP | CAPTCHA / abstract-only |
| a CSV / repo mirror | — | targeted search | **does not exist** |

**Wall class: EXTERNAL (publisher access), not code.** The harness is done; it needs **one spreadsheet**.

## Wall 2 — the joinable set is DEGENERATE (the more interesting one)

Falling back to the accessions whose FRI/FLC status the cell's own catalog carries, **5 joined** — and the
null-baseline gate killed the result:

| accession | FRI | FLC | predicted | measured (d) | observed |
|---|---|---|---|---:|---|
| Col-0 | LoF | strong | early | 28.00 | early |
| Van-0 | functional | nonsense | early | 27.00 | early |
| Bur-0 | functional | null-behaving | early | 54.25 | early |
| **Sf-2** | **functional** | **strong** | **late** | **40.00** | early |
| Wil-2 | LoF | strong | early | 36.00 | early |

Cohort median (the **derived**, not asserted, early/late split) = **61.25 days**. **All five joined accessions
fall below it** → the observed labels carry **only one class**.

> **Cell 4/5. Constant-"early" null 5/5.** A binary call cannot be scored on a set with no label variation —
> the null attains the majority rate *by construction*, so **"4/5 = 80%" would be a meaningless number**.
> **Verdict: `DEGENERATE_NO_LABEL_VARIATION`.**

This is the `feedback_threshold_vs_null_baseline_sanity_check` lesson firing on my own work: the agreement
count looked publishable until the null was computed. The Sf-2 "miss" is likewise uninterpretable on a
degenerate set (and DTF1's growth/vernalization regime is a live confound for a winter-annual's day-count).

## Honest status

**The flowering cell remains faithful-to-literature, NOT scored.** Nothing about its correctness changed —
its four literature anchors still reproduce (incl. the Da(1)-12 anchor a naive FRI-only rule mis-calls). What
changed is that the *validation* claim is now precisely bounded rather than aspirational.

**Exactly what would close it (a named, one-step external dependency):** any ONE of —
- Zhang 2020 *Plant J* **Table S1** (Wiley supporting information; institutional access or manual download), or
- Shindo 2005 supplementary (PMC, needs a human past the CAPTCHA), or
- Werner 2005 supplementary.

Drop that spreadsheet in and the harness scores immediately: it already parses AraPheno, derives the
threshold, joins by accession, and gates on the null. **This converts an external wall into a code wall in
one file.**

Run: `uv run python scripts/flowering_arapheno_spotcheck.py` (fetches free AraPheno; `--dtf1 <cached.csv>`
for offline). Artifact: `wiki/arabidopsis_flowering_spotcheck_2026-07-16.json`. Frozen decoder surface
byte-unchanged.
