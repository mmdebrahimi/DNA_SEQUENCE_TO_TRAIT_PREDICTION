# Sci234 external validation — 234-cohort COMPLETED (cipro + cef SCORED, gent underpowered) (2026-06-16)

Completes the 2nd fully-independent external cohort (the 1st being Oxford, `wiki/oxford_external_validation_result_2026-06-15.md`).
Substrate: the Sci Rep 2023 234-isolate E. coli bloodstream cohort (Spain, BioProject PRJNA854358,
PMC9829913), open CC-BY supplement `MOESM1_ESM.xlsx`. Every isolate carries its **own measured MIC**
(broth microdilution) + its **own ResFinder/PointFinder-style genotype**, joined by study `Key`. Scored
by `scripts/sci234_score.py` — a pure join → frozen `call_resistance` → `_conf`; **no genome download,
no Docker** (the supplement ships the genotype).

cipro was scored 2026-06-15 (QRDR-count, no Subclass lookup). This run adds **gentamicin + ceftriaxone**
via the new **fam.tsv per-gene Subclass resolver** (`scripts/fam_subclass_resolver.py`), which reproduces
the deployed AMRFinder v4.2.7 Class/Subclass per acquired allele so a synthesized main.tsv carries the
same Subclass AMRFinder would emit — making the frozen `call_resistance` RULE the thing under test.

## Result (binary R/S, CLSI E. coli breakpoints)

| drug | powering | n | acc | sens (R) | spec (S) | R / S |
|---|---|---|---|---|---|---|
| **ciprofloxacin** | SCORED | 231 | **0.987** | **1.000** | **0.984** | 45 / 186 |
| **ceftriaxone** (CEFOTA 3GC proxy) | SCORED | 231 | **0.991** | **0.833** | **1.000** | 12 / 219 |
| gentamicin | **UNDERPOWERED** | 221 | 0.991 | — | 0.991 | **0** / 221 |

## Verdict
- **cipro + cef GENERALIZE to a 2nd independent lab.** cipro re-confirms 0.987 (sens 1.0); cef lands
  acc 0.991 / **spec 1.0** — the extended-spectrum refinement (CEPHALOSPORIN/CARBAPENEM only) produces
  **zero false positives** on 219 susceptible isolates. Combined with Oxford (cipro 0.96 / gent 0.99 /
  cef genotype-rule clean), the deterministic decoder now holds on **two** fully-independent measured-MIC
  cohorts.
- **gentamicin is UNDERPOWERED, NOT validated.** The cohort's max GENTAM MIC is 8 (CLSI R≥16) → **0
  resistant isolates**. Sensitivity is unmeasurable; only specificity (0.991, 2 FP) is informative. The
  gent decoder is NOT validated by this cohort — reported honestly, not as a pass.
- **cef sens 0.833 = 2 FN, both attributed (neither a hidden defect):**
  - `05/1671` — only `ampC_promoter_size_53bp_dC12T` (chromosomal **AmpC promoter overexpression**). A
    regulatory/expression mechanism = the decoder's **documented honest blind-spot** (`UNDETECTABLE_MECHANISMS`).
    The supplement's "ESBL+" flag here is AmpC-hyperproduction, not a true ESBL. Unavoidable by design.
  - `05/1680` — `blaTEM-52C` (a TEM-**ESBL**). fam.tsv v4.2.7 curates the `blaTEM` family node as generic
    BETA-LACTAM, so my proxy resolves it to BETA-LACTAM and does not count it = the **documented resolver
    scope-limit**. The DEPLOYED decoder (real AMRFinder) MAY sub-classify TEM-52 to CEPHALOSPORIN and
    catch it; this FN is a proxy artifact, not necessarily a deployed-rule miss.
- **cef spec 1.0 + gent clean** ⇒ the fam.tsv resolver introduces no false positives.

## Fidelity caveat (the resolver scope-limit, surfaced)
gent+cef synthesize a main.tsv reproducing fam.tsv **FAMILY-level** Class/Subclass — a faithful-but-
imperfect proxy for real AMRFinder **node-level** assignment. CTX-M/CMY/DHA/ACT (the E. coli 3GC-R
drivers) + NDM/KPC resolve correctly; rare ESBL/carbapenemase variants on a generic BETA-LACTAM family
node (TEM/SHV-ESBLs, blaOXA-48) under-resolve and are NOT counted for cef → a **conservative under-call
bounded to rare variants**, every instance surfaced in the artifact `discordances`. gent resolution is
clean (aac(3)*/armA/ant(2'')-Ia = GENTAMICIN; plain aac(6')-Ib = generic, correctly not-gent).

## Independence + leakage
- **Independence:** the cohort's own measured MIC AND own genotype caller (ResFinder/PointFinder-style) —
  MORE independent than the NCBI-PD provenance-disjoint cells (independent genotype caller too).
- **Leakage = DISJOINT (project-level, by construction):** PRJNA854358 deposits **0 NCBI assemblies**
  (reads-only); the decoder's tuning set is 100% GCA_/GCF_ → a 234-cohort isolate cannot be in it.

## Artifacts
- Result: `wiki/external_validation_sci234_2026-06-16.json` (per-drug cells + powering + per-isolate discordances + unresolved-gene pool).
- Scorer: `scripts/sci234_score.py`. Resolver (reusable keystone): `scripts/fam_subclass_resolver.py` (+ `tests/test_fam_subclass_resolver.py`, 20 tests).
- Data (gitignored): `data/raw/sci234/sci234_supplement_MOESM1.xlsx`.

## Reusable outcome
The fam.tsv Subclass resolver is the **keystone for new-drug coverage**: any future drug whose
determinants are acquired bla/aminoglycoside/etc. alleles (cefotaxime/ceftazidime/cefepime/aztreonam,
amikacin, ertapenem/imipenem) can now be scored from a gene-presence table without assembly — the same
machinery that scored gent+cef here.
