# Horesh 2021 File F1 — Per-Record Label-Provenance Audit

**Date:** 2026-05-29
**Action:** Short-term Action 1-new (EP-4 pathotype CLI ledger)
**Gates:** H1 (label-circularity) status on the Tier-2 substrate; H2 per-class floor feasibility; H3 fold-construction feasibility
**Source file:** `F1_genome_metadata.csv` (4.83 MB, 10,146 records) — Figshare DOI 10.6084/m9.figshare.13270073, file id 25552514. Local copy: `data/external/horesh2021_F1_genome_metadata.csv`
**Method paper:** Horesh et al. 2021, *Microbial Genomics* (PMC8208696)

---

## TL;DR verdict

| Question | Verdict |
|---|---|
| **H1 on full 10,146-record collection** | **FAIL** — only 20.5% of records carry v0-independent label provenance; 52.2% are explicitly gene-rule-derived (`(predicted)`), tripping the H1 FAIL sub-criterion (gene-rule-derived > 40%). |
| **H1 on the non-predicted subset (N=2,077)** | **PASS by construction** — these labels are isolation-source-rule (ExPEC) or dedicated-diarrheagenic-study curated (von Mentzer / Hazen / Ingle), independent of v0 marker rules. |
| **Horesh Tier-2 role** | **CONFIRMED for 3 of 6 v0 target classes** (ExPEC/UPEC, EPEC, ETEC) — clear H2 floors AND H3 fold-diversity. **INSUFFICIENT** for EHEC/STEC-LEE (N=47 + lineage leakage), EAEC (N=2), commensal (N=2). |
| **Net** | Horesh is a usable Tier-2 substrate **only when filtered to non-`(predicted)` labels**, and only covers the extraintestinal + EPEC + ETEC classes. EAEC / commensal / EHEC-STEC must come from Tier-1 (DECA / von Mentzer) or another source. **Do NOT train on the full 10,146 — the gene-predicted majority is circular with the v0 resolver.** |

---

## 1. The provenance flag: `(predicted)` suffix

File F1's `Pathotype` column encodes per-record provenance directly via a `(predicted)` suffix. This is the per-record metadata the H1 audit needed (the methods summary alone hid it — cf. per-record-label-provenance-audit lesson).

Split of all 10,146 records:

| Provenance class | N | % | v0 independence |
|---|---|---|---|
| **`(predicted)`** — gene-rule-derived (ariba + VirulenceFinder DB markers) | 5,295 | 52.2% | **CIRCULAR** with v0 |
| **Clean label** — no `(predicted)` suffix | 2,077 | 20.5% | **INDEPENDENT** (see §2) |
| **`Not determined`** | 2,774 | 27.3% | no label |

Confirmation that `(predicted)` == gene-only: of the 5,295 predicted records, **4,379 (83%) have `Isolation = Unknown`** — i.e. no isolation-source signal at all, the label rests entirely on marker-gene calls.

## 2. Why the clean labels are independent

The clean (non-`predicted`) 2,077 break into two independence mechanisms:

**(a) ExPEC by isolation-source rule (N=1,574).** Horesh method: *"if the source of isolation was blood or urine the assignment was ExPEC."* 1,567 of the 1,574 clean ExPEC are blood/urine isolates. This label is derived from clinical isolation site, **not** from virulence genes → independent of the v0 marker resolver.

**(b) Diarrheagenic by dedicated source-study (N=501 clean diarrheagenic).** Source-publication breakdown:

| v0 class | von Mentzer (ETEC study) | Hazen (DECA) | Ingle (diarrheagenic) | generic surveillance (PHE) | other study |
|---|---|---|---|---|---|
| ETEC | 181 | — | — | — | 2 |
| EPEC | — | 126 | 141 | — | 2 |
| EHEC | — | 1 | — | — | 5 |
| NON-O157EHEC/EPEC | — | 10 | — | — | — |
| STEC | — | — | — | 31 | — |
| EAEC | — | 1 | — | — | 1 |

ETEC, EPEC, EHEC, NON-O157 clean labels come overwhelmingly from **dedicated diarrheagenic / pathotype studies** where the pathotype was established by the original authors via clinical / epidemiological / serotyping means — literature-curated, independent of v0's marker rules.

**Caveat (weaker independence):** the 31 clean **STEC** labels come from PHE generic surveillance, not a dedicated study. PHE STEC surveillance uses PCR/serotyping (partially independent), but conservatively these should be flagged lower-confidence-independent.

## 3. H1 fraction math

- Independent / all records = 2,077 / 10,146 = **20.5%** (FAIL vs ≥70% target)
- Gene-rule-derived = 5,295 / 10,146 = **52.2%** (FAIL — H1 fails if labels explicitly gene-rule-derived > 40%)
- Provenance-missing (`Not determined`) = 27.3% (< 40%, so the "missing" sub-criterion alone is OK; the gene-derived criterion dominates)

→ **H1 FAILS on the full collection.** The non-predicted subset is independent by construction, so the operational reading is: *Horesh is Tier-2 only after filtering to clean labels.*

## 4. H2 per-class floor feasibility (on independent subset, N=2,077)

| v0 class | Independent N | H2 ship floor | Verdict |
|---|---|---|---|
| ExPEC / UPEC-compatible | 1,574 | ≥50 | ✓ **PASS** (massive) |
| EPEC | 269 | ≥50 | ✓ **PASS** |
| ETEC | 183 | ≥50 | ✓ **PASS** |
| EHEC / STEC-LEE | 47 (6 EHEC + 31 STEC + 10 NON-O157) | ≥50 | ✗ **marginal FAIL** (+ STEC independence weak) |
| EAEC | 2 | ≥50 | ✗ **HARD FAIL** |
| commensal / low-marker | 2 | ≥75 | ✗ **HARD FAIL** |

## 5. H3 fold-construction feasibility (ST concentration per class)

| v0 class | N | distinct ST | top-ST share | top-3 ST share | H3 read |
|---|---|---|---|---|---|
| ExPEC/UPEC | 1,574 | 231 | ST73 = 17% | 43% | ✓ foldable |
| EPEC | 269 | 84 | ST328 = 8% | 23% | ✓ **best** — well-distributed |
| ETEC | 183 | 82 | ST443 = 9% | 20% | ✓ well-distributed |
| EHEC/STEC-LEE | 47 | 6 | **ST21 = 60%** | 94% | ✗ severe lineage leakage (ST21/O157 shortcut — the exact H3 confound) |
| EAEC | 2 | 2 | 50% | 100% | ✗ unfoldable |
| commensal | 2 | 2 | 50% | 100% | ✗ unfoldable |

ExPEC, EPEC, ETEC support ≥5 clade-balanced folds without lineage leakage. EHEC/STEC additionally fails H3 (ST21 = 60% would let an O157/ST21 serotype shortcut explain class success). EAEC + commensal cannot be folded at N=2.

## 6. Implications for the ledger

- **U6-new RESOLVED:** Horesh independent-label fraction = 20.5% full / 100% on filtered subset.
- **H1 status:** keep `under-investigation` overall, but record the per-substrate finding — **FAIL on full collection, PASS on the non-predicted subset for ExPEC/EPEC/ETEC only.**
- **Substrate-tier strategy:** Horesh is **Tier-2 confirmed** for ExPEC + EPEC + ETEC (independent, H2-clearing, H3-foldable). It is **NOT** a source for EAEC, commensal, or (cleanly) EHEC/STEC — those stay dependent on Tier-1 (DECA / von Mentzer) + the still-open Whittam contact + NCBI host_disease (commensal proxy).
- **Action ordering nudge:** EAEC + commensal are now the binding constraint on H2 across the whole project, not EPEC/ETEC. The Whittam DECA contact (Short-term Action 2) and a commensal source matter more than previously weighted.

---

## Provenance & reproducibility

- File F1 enumerated via Figshare REST API (`GET /v2/articles/13270073`) — WebFetch on the Figshare HTML 403s; the API does not.
- All counts reproducible: `data/external/horesh2021_F1_genome_metadata.csv`; analysis is the three pandas cross-tabs in this session (Pathotype × `(predicted)`, clean-diarrheagenic × Source, per-class ST concentration).
- Independence interpretation rests on the Horesh method text ("marker virulence genes … refined by the source of isolation: if … blood or urine … ExPEC") + the `(predicted)` suffix semantics + the 83%-Unknown-isolation confirmation for predicted records.
- Residual uncertainty: the 31 PHE-surveillance STEC clean labels — independence not fully established (flagged lower-confidence). Does not change any verdict (EHEC/STEC fails on N + H3 regardless).
