# Soraya --until-mvp — result (new gear: antimalarial kelch13)

- Run: 2026-06-09-ruv8-antimalarial-kelch13 · verdict: mvp-reached (+ stretch G0-completion)
- Gear chosen by user (over EXPRESSION-frontier / antiviral): antimalarial vertical.

## MVP bar — all met
- C1 file-exists dna_decode/data/antimalarial_amr.py — MET
- C2 file-exists scripts/pf_kelch13_caller.py — MET
- C3 test-exit-0 pytest tests/test_pf_kelch13.py — MET (6 tests; full suite 955 passed, 0 regressions)
- STRETCH (G0-completion): real 3D7 K13 reference committed + numbering validated on the REAL reference — MET

## What shipped — the 3rd kingdom (protozoan)
P. falciparum artemisinin partial resistance via WHO-validated Pfkelch13 markers. Reuses the fungal caller's
gene-generic BLAST+codon-map machinery (DRY); only new surface = catalog + thin wrapper + real reference.
bacteria (AMRFinder) -> fungi (ERG11) -> protozoa (K13).

## Next (follow-on, not this MVP)
1. Wire into dna-amr --drug artemisinin (CLI route, like the fungal productization).
2. Real C580Y isolate cohort validation (G1-style).
3. pfcrt/pfmdr1 (chloroquine / partner-drug) as additional antimalarial loci.
