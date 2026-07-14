# TB full 12-drug CRyPTIC-parquet panel — the concentration rule, measured across a kingdom's drug classes (2026-07-14)

Completes the TB second-line build: the shipped RIF/INH cell + the flagship moxifloxacin/bedaquiline pair
(`wiki/tb_secondline_cryptic_parquet_baseline_2026-07-14.md`) are now joined by the **8 remaining panel
drugs**, all scored in a **single parquet pass** (union every drug's determinant+barcode positions, stream
`VARIANTS.parquet` once, score each drug from the shared calls — `scripts/score_tb_cryptic_parquet.py
--drugs all-remaining`, commit `776f234`; 8×15min → 1×~18min). WHO mutation catalogue v2 (grade-1/2)
applied UNCHANGED; scorer (`tb_amr.score_drug`) drug-generic.

## The full panel (10 cells scored; 12-drug `DRUG_CODE` wired)

| Drug | Class / target | n | raw sens/spec | **lineage-collapsed sens/spec** | R/S/disc lineages |
|---|---|---:|---:|---:|---:|
| **levofloxacin** | fluoroquinolone — gyrA/gyrB QRDR | 7773 | 0.816/0.984 | **0.333 / 0.997** | 45/299/45 |
| **ethambutol** | embB (M306 hotspot) | 6614 | 0.855/0.953 | **0.293 / 0.991** | 75/213/45 |
| **moxifloxacin** | fluoroquinolone — gyrA/gyrB QRDR | 6784 | 0.841/0.961 | **0.279 / 0.988** | 43/259/40 |
| **amikacin** | aminoglycoside — rrs a1401g | 8972 | 0.797/0.996 | **0.263 / 1.0** | 38/366/36 |
| **kanamycin** | aminoglycoside — rrs/eis | 9332 | 0.771/0.983 | **0.229 / 0.997** | 48/376/41 |
| ethionamide | ethA LoF + inhA | 8251 | 0.638/0.983 | 0.083 / 1.0 | 48/301/51 |
| linezolid | rrl / rplC | 7140 | 0.303/1.0 | 0.0 / 1.0 | **5**/317/22 |
| clofazimine | Rv0678 LoF (BDQ-cross) | 7762 | 0.009/0.998 | 0.0 / 1.0 | **6**/359/24 |
| delamanid | ddn/fbiA-C/fgd1 LoF | 8094 | 0.059/1.0 | 0.0 / 1.0 | **5**/386/21 |
| bedaquiline | Rv0678/atpE/pepQ LoF | 8535 | 0.183/0.991 | 0.0 / 0.938 | **5**/406/17 |

(RIF raw 0.916→lineage 0.41, INH 0.889→0.349 from the shipped cell — same shape.) Every cell `status =
WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`.

## The finding, sharpened — the concentration rule holds, along THREE axes not one

The project's recurring **concentrated-mechanism-works / distributed-mechanism-fails** boundary (E. coli
cipro-QRDR passes, tet distributed-mobile fails; HIV ESM antagonistic-selection blindness) is now measured
across an entire pathogen's drug panel. The 10 cells separate cleanly into three regimes:

**1. Concentrated target-site → WORKS (5/5).** Every drug whose resistance concentrates at defined
point-mutation hotspots — fluoroquinolone gyrA/gyrB QRDR (LEV, MXF), aminoglycoside rrs a1401g (AMI, KAN),
embB M306 (EMB) — lands at lineage-collapsed sens **0.23–0.33** with robust spec **≥0.99**, the SAME band
as the validated first-line RIF/INH cell. This is the fluoroquinolone/aminoglycoside analog of the E. coli
cipro cell transferring into M. tuberculosis with an honest in-distribution number.

**2. Diffuse loss-of-function → FAILS.** Resistance reached by *diverse* LoF across a gene (ethA for ETH;
Rv0678 for CFZ/BDQ; ddn/fbiA-C for DLM) is structurally invisible to a point-mutation catalogue.
**Ethionamide is the cleanest evidence:** 48 R-lineages (NOT rare), a 2107-row determinant set, and 729
exact indel targets — and it *still* collapses to lineage-sens **0.083**. The catalogue's size can't rescue
a mechanism whose realizations aren't enumerable. (CFZ/DLM/BDQ are diffuse *and* rare — see axis 3 — so ETH
is the load-bearing diffuse-mechanism cell.)

**3. Rare-signal under-powering → a DISTINCT failure mode (LZD).** Linezolid is *not* diffuse — rrl/rplC
carry real hotspots, and its **raw** sens is 0.303 (the isolate-level signal IS there). But with only **5
R-lineages** in CRyPTIC, lineage-collapse zeroes it. That is a *rarity/power* limit, not a mechanism limit —
and it must be reported separately from the diffuse-LoF failures it superficially resembles (CFZ/DLM/BDQ
share the ≤6-R-lineage rarity, confounding their diffuseness). **Reusable rule: a lineage-collapsed 0.0 has
two causes — un-cataloguable mechanism vs too-few-lineages — and the raw-sens + R-lineage-count disambiguate
them.** ETH (48 lineages, raw 0.638, lineage 0.083) is mechanism; LZD (5 lineages, raw 0.303, lineage 0.0)
is power.

## Indel-matching marginal effect (the exact WHO-indel → CRyPTIC-string realization)

Indel determinants matter only for the LoF drugs: ETH +17 true-positive S→R flips (15 FP), BDQ +9 TP (61 FP),
CFZ +1 TP (13 FP). For the concentrated drugs there are 0 indel determinants (point mutations) — the caveat
is moot. Even with indels ON, the diffuse drugs stay low-sens: indel matching lifts the floor but does not
close the un-cataloguable gap.

## Honesty tier (identical to the shipped cell — do not overclaim)

1. **In-distribution, NOT independent.** The WHO catalogue was built partly from CRyPTIC → a CRyPTIC-scored
   number is a KNOWLEDGE baseline, not the independent validation the AMR-Portal arm gives RIF/INH.
2. **Raw is clonality-inflated → lineage-collapsed is the headline.** TB is heavily clonal; raw counts one
   vote per isolate. Lineage-collapsed (barcode-on-VCF, one vote per lineage) is the honest number.
3. **SNV+exact-indel determinant match; complex delins not normalized → a LOWER BOUND.** The diffuse-LoF
   sens numbers are lower bounds; the concentrated-drug numbers are point-mutation-complete.

## What this banks

- **10 real TB AMR cells** (2 first-line shipped-and-independently-validated + 8 second-line in-distribution
  baselines), spanning fluoroquinolone / aminoglycoside / ethambutol / thioamide / oxazolidinone /
  Rv0678-LoF / nitroimidazole mechanism classes — the broadest single-organism decoder coverage in the repo.
- **The concentration rule is now a measured, three-axis screening tool**, not a two-way heuristic: before
  building a determinant-catalogue cell, classify the mechanism as concentrated-hotspot (expect the RIF/INH
  band) vs diffuse-LoF (expect low-sens) vs rare-in-substrate (expect under-powered, raw>lineage), and read
  raw-sens + R-lineage-count to tell diffuse from rare.
- **The independent (non-CRyPTIC) arm remains the next real bar** for any of these — the AMR-Portal
  provenance-disjoint runner already validated RIF/INH; the same substrate could lift a concentrated
  second-line cell (LEV/EMB/AMI) to a genuine out-of-distribution number.

Frozen decoder surface (`amr_rules` / `calibrated_amr_rules` / `mic_tiers` / `shipped_decoder_surface` /
`cohort_manifest`) byte-unchanged (`verify_lock OK`); `tb_amr` + the parquet adapter are NON-frozen.
