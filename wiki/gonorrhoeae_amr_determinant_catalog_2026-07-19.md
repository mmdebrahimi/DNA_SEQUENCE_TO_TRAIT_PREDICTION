# N. gonorrhoeae AMR determinant catalog — sourced, for the new gonococcal cell

**Date:** 2026-07-19 (cell BUILT 2026-07-20) · **Status:** CATALOG SOURCED + **CELL BUILT + real-symbol-
validated** (`dna_decode/organism_rules/neisseria_amr.py`, commits 115b6cb + b1e017c); the 94-genome
measured-MIC SCORING is the remaining step (Kaggle-native AMRFinder). · **Purpose:** the determinant→
phenotype catalog + rules to score the AR Bank's 94 gonococcal isolates. · **Grounding:** Pathogenwatch
community catalog + AMRFinderPlus + Eyre 2017 + WHO 2016 reference panel. Frozen decoder surface untouched.

## Cell build status (2026-07-20)

The gonococcal cell was NOT a new file — `organism_rules/neisseria_amr.py` already existed (cipro gyrA
QRDR + tet tet(M)). It was EXTENDED (non-duplicating) to the full AR Bank panel: `call_ng_azithromycin`
(23S A2045G/C2597T → R, mtrR accessory), `call_ng_ceftriaxone`/`call_ng_cefixime` (any curated penA point
→ R, ponA/porB/mtrR accessory), `call_ng_penicillin` (blaTEM → R, chromosomal accessory), `call_ng_gentamicin`
(**ABSTAIN** — no determinant), + a `call_ng_amr(drug, symbols)` dispatcher. **R3 real-surface validated:**
ran AMRFinder `-O Neisseria_gonorrhoeae` on a real gono genome (AR#0165, GCA_042036815.1) — this CAUGHT a
bug (a hard-coded penA codon set missed real mosaic positions 504/510) now fixed to match any curated penA
point; the full real symbol vector is pinned in `test_ng_real_amrfinder_symbols_AR0165` (12 tests green).
**Remaining:** wire an AR-Bank-gono scorer (route to `call_ng_amr` instead of the frozen `call_resistance`)
+ run the 94 genomes via Kaggle-native AMRFinder → measured-MIC one-sided validation (spec≥0.85 falsifier
per drug, like the cipro/tet rules).

## The two-layer answer (mirrors the E. coli/Klebsiella cells)

1. **Detection layer — AMRFinderPlus `-O Neisseria_gonorrhoeae`** (VERIFIED present in the pinned
   `ncbi/amr:4.2.7-2026-03-24.1`: `amrfinder -l` lists `Neisseria_gonorrhoeae`). Calls the gonococcal
   point mutations (penA, gyrA, parC, 23S rRNA, mtrR, ponA, porB, rpsJ, folP) + acquired genes
   (blaTEM, tetM) — the same `main.tsv` shape the E. coli/Klebsiella cells already consume.
2. **Interpretation catalog — Pathogenwatch `485.toml`** (taxid 485; `github.com/pathogenwatch-oss/
   amr-libraries/blob/main/485.toml`, ~62 KB, CGPS-curated, EuroGASP + Grad-lab/PHE/PHAC provenance).
   The authoritative, open, machine-readable **determinant → S/I/R** library; validated on Pathogenwatch's
   12,000-genome collection with **specificity >99% for azithromycin, ciprofloxacin, ceftriaxone** and
   **sensitivity ~99% for penicillin, tetracycline** (Genome Medicine 2021).

**VERIFY-IN-BATCH CATCH (honesty):** a web result asserted `666.toml` was "the current library" — it is
NOT the gonococcal file. `666.toml` extends `gram_neg_esbl`, uses **E. coli gyrA numbering (S83/D87)**,
and carries blaCMY/tetA acquired genes with **no** penA/gono-gyrA(S91)/23S — it's an Enterobacterales
library. The gonococcal file is **`485.toml`** (NCBI taxid 485), confirmed by its content (gyrA S91F,
23S A2045G, penA mosaic). Do not wire 666.toml.

## Drugs — catalog vs the AR Bank panel

- **Pathogenwatch 485.toml predicts 8:** AZM, CRO (ceftriaxone), CFM (cefixime), CIP, PEN, SMX
  (sulfamethoxazole), SPT (spectinomycin), TCY.
- **AR Bank gono panel measures 7** (confirmed from a real cached isolate page, AR Bank #1280
  `SAMN35332250`): **Azithromycin, Cefixime, Ceftriaxone, Ciprofloxacin, Gentamicin, Penicillin,
  Tetracycline.** 52+ gono isolates already enumerable from the cached panel pages (of the 94 total).
- **Scorable overlap = 6** (azithromycin, cefixime, ceftriaxone, ciprofloxacin, penicillin, tetracycline).
- **⚠ GAP — gentamicin:** the AR Bank measures it, but **there is no validated gonococcal genetic
  determinant for gentamicin resistance** (a newer/alternative agent; no established marker in Pathogenwatch
  or the literature). The gono cell must **ABSTAIN** on gentamicin — do NOT invent a rule. (Same discipline
  as the fungal/antiviral abstainers.)

## Per-drug determinant rules (verbatim-sourced from 485.toml)

**Ciprofloxacin (CIP) — the clean single-marker drug** (Eyre 2017: gyrA_S91F 32/32 = 100%):
- `gyrA`: **S91F → R**, D95N → R, D95A → R, D95G → I
- `parC`: D86N → R, S87R → R; S87I / S87N / S88P / E91K → I
- Practical rule: gyrA S91F (or D95N/A) ⇒ R; wild-type gyrA-91 ⇒ S (the validated single-marker call).

**Azithromycin (AZM):**
- `23S rRNA`: **A2045G → R**, **C2597T → R** (gonococcal nomenclature; = E. coli A2059G/C2611T).
  ⚠ short-read assemblies collapse the 4× 23S copies into one → detection is consensus-base-dependent
  (a mixed-copy genotype may under-call). Note this as a known sensitivity caveat.
- `mtrR`: promoter mosaics 1–3 + `a-57del` ⇒ R (SUPPRESSED if `mtrC` disrupted); `a-57del` alone ⇒ I (PEN/TCY).

**Ceftriaxone (CRO) + Cefixime (CFM) — multilocus, mosaic-driven:**
- `penA` mosaic combos: `I312M+V316T+G545S` ⇒ R (CFM); `I312M+V316P+T483S+G545S` ⇒ R (CFM/CRO);
  single A501P / T483S / V316P ⇒ I–R (CRO).
- `ponA` L421P ⇒ I (PEN); R only combined with porB + mtrR.
- `porB1b`: G120K + A121N/D ⇒ I (PEN alone), R with mtrR_promoter `a-57del`.

**Penicillin (PEN):**
- `blaTEM` (whole-gene presence) ⇒ R (effect SUPPRESSED if disrupted).
- `penA` combos (e.g. `ins346D+L421P+G542S`) + `porB1b` A121D ⇒ R; mtrR_promoter `a-57del` ⇒ I alone.

**Tetracycline (TCY):**
- `tetM` (gene presence) ⇒ R (SUPPRESSED if disrupted).
- `rpsJ` V57M ⇒ I alone; R when combined with mtrR (A39T / G45D / disrupted) or promoter `a-56c`.

## The catalog schema (what the build must implement)

485.toml is **materially more complex than the frozen E. coli `DRUG_RULE`** (which is simple
`threshold + subclass_any` count/OR). Each mechanism is:
```toml
{phenotypes = [{effect = "RESISTANT|INTERMEDIATE", profile = ["AZM", ...],
   modifiers = [{gene = "mtrC", variants = ["disrupted"], effect = "SUPPRESSES"}]}],
 members = [{gene = "penA", variants = ["G542S"]}, {gene = "porB1b", variants = ["G120K","A121D"]}]}
```
Combination logic: **ALL `members` (genes + their variants) must be present simultaneously** for the
phenotype; a `modifier` with `effect="SUPPRESSES"` cancels it if the named variant is detected. So the
build needs a small **combination interpreter** over the AMRFinder determinant set, NOT a re-derivation
from AMRFinder's `Class`/`Subclass` alone (which can't express the multi-gene AND + suppressor logic).

## Recommended build architecture (for the next run)

1. **Vendor `485.toml`** (pin a tagged release + sha256, like the WHO TB catalogue pin) into a NON-frozen
   `dna_decode/organism_rules/gono_amr.py` + `data/gono_catalogue/`. Do NOT touch the frozen E. coli surface.
2. **Detection:** reuse the shipped AMRFinder path with `-O Neisseria_gonorrhoeae` (the scorer is already
   organism-parameterized; add the gono determinant-parsing).
3. **Interpret:** a small TOML combination-rule engine (members-all-present + SUPPRESS modifiers) →
   per-drug S/I/R. Map I→(R or S) explicitly per the AR Bank/CLSI intent (document the choice).
4. **Score** the 6 overlapping drugs against AR Bank measured MIC via the external-cohort arm (hardened
   preflight → provenance-disjoint → one-sided powering as before); **ABSTAIN on gentamicin**.
5. **Honesty:** Pathogenwatch's catalog is a curated genotype→SIR library (a *knowledge* baseline, like the
   WHO TB catalogue on CRyPTIC) — scoring AR Bank (independent CDC MICs, provenance-disjoint) is a real
   external test of the catalog, but it is NOT methodology-independent. Label it
   `PATHOGENWATCH_CATALOGUE_ON_AR_BANK` accordingly.

## Compute note

94 gonococcal genomes >> the ~50-genome local-Docker threshold (the Klebsiella run wedged repeatedly at
63). Per the captured lesson, **run AMRFinder on the Kaggle-native-bioconda path** (`ncbi-amrfinderplus`,
no Docker) → pull `main.tsv` → apply the gono catalog locally. Don't babysit local Docker for this one.

## Sources

- Pathogenwatch gono catalog (machine-readable): `github.com/pathogenwatch-oss/amr-libraries/blob/main/485.toml`
  (+ `FORMAT.md`); GitLab mirror `gitlab.com/cgps/pathogenwatch/amr-libraries`.
- Pathogenwatch gono resource + validation (spec >99% AZM/CIP/CRO): Genome Medicine 2021,
  https://link.springer.com/article/10.1186/s13073-021-00858-2 (PMC8054416).
- WGS→AST validation (gyrA_S91F 32/32; single-marker CIP): https://pubmed.ncbi.nlm.nih.gov/36735316/ +
  the 481-isolate global collection https://www.omicsdi.org/dataset/biostudies-literature/S-EPMC10319951.
- NG-STAR 7-locus scheme (penA, mtrR, porB, ponA, gyrA, parC, 23S) — PubMLST `pubmlst.org/neisseria`.
- Determinant meta-analysis (penA/mtrR/gyrA/parC/tetM): https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2024.1414330/full
- AMRFinderPlus organism support: verified locally (`amrfinder -l` in ncbi/amr:4.2.7 lists Neisseria_gonorrhoeae).
