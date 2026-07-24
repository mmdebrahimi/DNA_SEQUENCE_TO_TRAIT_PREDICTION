<!-- memo-schema: 0.4 -->
# [CORRECTED 2026-07-24] Cross-lab RBP test set — the DATA-BLOCKED verdict was WRONG

> **CORRECTION:** this memo originally concluded the cross-lab RBP number was DATA-BLOCKED. That was an
> under-search error — it MISSED the **phageReceptor** database (Zhang et al., Bioinformatics 2020; Peng
> lab, Hunan U), surfaced by a user-relayed DeepSeek pointer. phageReceptor has 37 E. coli phages with
> OMP-class receptors (24 independent of LBNL/classic). The cross-lab number WAS computed:
> **0.364 (4/11 called)** — see wiki/phage_rbp_crosslab_result_2026-07-24.md. The Kaneko-vocab-mismatch
> finding below still stands, but it was not the whole picture. Lesson: an absence-of-evidence verdict
> needs a wider search than I ran.

---

# Cross-lab independent RBP receptor test set for E. coli phages (research memo, 2026-07-24)

> Topic (verbatim): "cross-lab independent RBP number". Source: Claude Code `/research` orchestrator (run by hand). Slug: cross-lab-phage-receptor-rbp-dataset-2026-07-24.
> Web search via WebSearch + WebFetch, single-pass. The underlying question: is there a SECOND, free, public E. coli phage dataset with per-phage MEASURED host receptors + RBP annotation + genomes, from a lab OTHER than LBNL/Arkin-Mutalik and BASEL/Maffei, usable to score the phage receptor cell's RBP caller cross-lab (its 0.975 is within-LBNL leave-one-out)?

## Research Context (problem anchor)

The phage receptor cell ships an RBP-level caller validated at 0.975 leave-one-out on the LBNL Phage Datasheets (Moriniere et al. 2026). That is a WITHIN-study method validation. A genuinely INDEPENDENT (cross-lab) RBP number needs a second public E. coli dataset with (a) measured per-phage receptors, (b) in the SAME receptor-class vocabulary the caller uses (the outer-membrane proteins OmpA/OmpC/OmpF/BtuB/FhuA/Tsx/LamB/LptD/FadL + LPS_core/NfrA), (c) genomes, ideally (d) RBP annotation.

## Candidate second sources (audit table)

| Dataset | Phages (n) | Receptor vocabulary | Genomes deposited | RBP annotated | Independent of LBNL/BASEL? | Usable as cross-lab test of covered OMP classes? | Source | Confidence |
|---|---|---|---|---|---|---|---|---|
| Kaneko et al. 2025 "From phenotype to receptor" (Waseda/Japan) | 13 (+3 validation) | **LPS R-core biosynthesis genes** (WaaV/WaaW/WaaT/WaaY/WaaG) + **NfrB** + flagella (FihD) + inner-membrane (YhaH/TolA) | Yes — DDBJ LC739530–LC739542, DRX534190/534192 | tail-fiber sequences extracted for phylogeny (not deposited separately) | **YES** (different lab) | **NO — ~0 class overlap** (paper states NfrB is "distinct from OmpA/OmpC/OmpF/BtuB/FhuA/Tsx/LamB/LptD/FadL"; receptors are LPS-biosynthesis genes + flagella, a different vocabulary; only the T5→FhuA reference overlaps) | J Virol 2025 jvi.01061-25 (PMC12548472) | high |
| LBNL Phage Datasheets (Moriniere et al. 2026) | 260 | OMP + LPS + ECA (matches caller) | .gbk in repo; GenBank "submitted" | YES (Table_S1 col 17) | **NO — this IS the set the caller is built + LOO-validated on** | n/a (the incumbent) | github.com/mjohnson11/PhageDataSheets | high |
| BASEL collection (Maffei 2021 / Humolli 2025) | 106 | OMP + LPS (matches catalogue) | GenBank MZ501046– / PRJNA1207239 | genus-level, not per-phage RBP CDS | **NO — the catalogue is curated FROM it; BASEL phages are also INSIDE the LBNL set** | n/a (the catalogue basis) | PLOS Biology 3001424 / 3003063 | high |
| EP75 / EP335 (E. coli O157 phages) | 2 | O-antigen tailspike (serotype-specific) | yes | RBP structurally characterized | yes | NO — n=2, O-antigen not OMP-class | PMC8217332 | high |
| Klebsiella phage RBP sets (e.g. Beamud/Pilar-Ferrer; pbio.3003515) | tens | capsule (K-type) via RBP | yes | yes | yes | NO — different HOST (host-range generalization test, not same-organism E. coli cross-lab) | PLOS Biology 3003515 | medium |

## Verdict (decision-grade)

**The cross-lab independent RBP number, in the caller's covered receptor-class vocabulary, is DATA-BLOCKED.**

- The two LARGE public E. coli measured-receptor datasets in a compatible vocabulary are **LBNL** and **BASEL** — and they are NOT independent of each other or of the cell: the BASEL phages are literally *inside* the LBNL set, and the catalogue is curated from BASEL. So they cannot cross-validate each other.
- The one clearly **cross-lab** E. coli measured-receptor + genome + tail-fiber dataset found — **Kaneko et al. 2025 (Waseda)** — defines receptors in a **different vocabulary** (LPS R-core biosynthesis genes WaaV/W/T/Y + NfrB + flagella), which the paper explicitly notes is *distinct from* the OMP classes the caller predicts. Overlap with the caller's classes ≈ 0 scoreable phages, so it cannot produce a cross-lab number for the covered classes (only the trivial T5→FhuA reference overlaps).
- Klebsiella RBP datasets exist but test a **different host** (a host-range-generalization question, not an E. coli cross-lab number).

**So a cross-lab independent RBP number is not achievable from current public data without either:** (a) a NEW lab publishing an E. coli phage set with measured **outer-membrane-protein** receptors + genomes at scale (does not exist publicly today beyond LBNL/BASEL), or (b) re-deriving Kaneko's LPS-biosynthesis/flagella receptors into the OMP vocabulary — not possible, it is different receptor biology.

## Decisions for Human Confirmation

| # | Candidate use | Verification needed | Source |
|---|---|---|---|
| 1 | Accept the within-LBNL LOO (0.975) as the terminal RBP validation; mark the cross-lab RBP number DATA-BLOCKED (no compatible-vocabulary second source) | none — this is a documented negative | this memo |
| 2 | OPTIONAL: run a Kaneko-based cross-lab check on the NfrB/NfrA overlap only (ΦWec189/191/193/196/270/272 → NfrB ≈ the caller's NfrA class) — a tiny (~6-phage), single-class cross-lab spot-check | fetch the 6 DDBJ genomes + extract tail fibers + map NfrB→NfrA; underpowered (1 class) but genuinely cross-lab | Kaneko 2025 (DDBJ LC739530–) |
| 3 | OPTIONAL: a Klebsiella host-range-generalization test (different question) using the Klebsiella RBP sets | separate build; different host + capsule receptors | pbio.3003515 |

## Honest gaps

- No public E. coli phage dataset with measured **OMP-class** receptors + genomes independent of LBNL/BASEL was found. The field's two large measured-receptor E. coli resources are LBNL + BASEL (non-independent of each other).
- Kaneko RBP sequences are not deposited as separate accessions (only whole-genome DDBJ); extracting them needs the same .gbk/annotation route used for LBNL.
