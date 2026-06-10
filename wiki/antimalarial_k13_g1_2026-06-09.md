# Antimalarial K13 — G1 real-data validation — 2026-06-09

> Gate G1 for the protozoan vertical: validate the Pfkelch13 caller (`scripts/pf_kelch13_caller`,
> `dna-amr --drug artemisinin`) on REAL GenBank *P. falciparum* K13 sequences with documented variants —
> not just the synthetic planted-C580Y (G0) or the real-3D7-reference numbering check (G0-completion).
> Mirrors the fungal EP-7 G1. Real 3D7 reference: `data/antimalarial_ref/Pf3D7_K13_cds.fna`
> (NCBI XM_001350122.1, 726aa, WT Cys@580). Caller = BLAST(K13-CDS-ref vs allele) → gap-aware codon-map →
> WHO-validated marker catalog (`dna_decode/data/antimalarial_amr.py`).

## VERDICT: G1 MET — caller decodes real K13 sequences correctly (genotype-fidelity)

| GenBank allele | documented | caller call | determinant | correct? |
|---|---|---|---:|---|
| OQ064537.1 (voucher "C580Y"/Lazai50, propeller) | C580Y *(voucher label)* | **R** | **K13:P574L** | ✅ (see below) |
| PV890521.1 (isolate pfk13A-46, partial CDS) | WT propeller | **S** | — | ✅ |

Plus 3 more real pfk13A field isolates (PV890522/523/524) all → **S** (WT propeller) in the exploratory run.

## The "suspect-the-label" finding (G1's real catch)

OQ064537.1's GenBank **voucher/isolate name is "C580Y"**, but the deposited amplicon is **colinear with the
3D7 reference at 99.88% — exactly ONE SNP, at codon 574 (ref Pro → Leu), i.e. P574L**. The reference has
Cys at 580 and the amplicon **matches the reference there** (no C580Y). So:

- The caller correctly decoded the **actual** single real variant (P574L) at the correct 3D7-numbered codon
  — validating the caller's numbering on REAL data (not just the synthetic/planted check).
- **P574L is itself a WHO-validated artemisinin partial-resistance marker**, so the **phenotype call (R) is
  correct** regardless of the title.
- The "C580Y" in the record name is a **strain/voucher label, not this amplicon's variant** — a textbook
  metadata-label-vs-sequence discrepancy. The genotype (the sequence) is the trustworthy signal; the caller
  reads it correctly. (Same discipline as the mecA/oxacillin + C. auris F126L label cautions.)

## What G1 establishes vs what it doesn't

- **Establishes:** the K13 caller + 3D7 reference numbering + WHO marker catalog correctly identify real
  WHO-validated propeller markers in real GenBank sequences (P574L detected at the right codon), and call
  WT field isolates S. The protozoan kingdom jump decodes real data. Genotype fidelity confirmed.
- **Does NOT establish (honest scope):** these records carry no artemisinin *clearance phenotype* (RSA /
  half-life), so this validates GENOTYPE (marker detection), not phenotype PREDICTION on a clinically-typed
  cohort — exactly the label limitation the fungal G1 also hit. N=5 real alleles (1 marker-positive, 4 WT),
  propeller-domain partial sequences. A clinically-phenotyped WGS cohort (e.g. Pf3k/MalariaGEN with
  documented ART-R status) is the future G1-completion. C580Y specifically is validated by NUMBERING (the
  G0-completion test on the real reference), not yet by a real C580Y-carrying sequence (this voucher
  mislabel means a genuine C580Y amplicon still needs sourcing).

## Fixtures + tests
- `data/antimalarial_ref/Pf_K13_OQ064537_P574L.fna` (real, P574L) + `Pf_K13_PV890521_WT.fna` (real, WT).
- `tests/test_pf_kelch13.py::test_caller_on_real_genbank_k13_alleles` (skip-if-no-BLAST/fixtures). 8 tests green.
