# EP-PGx — independent functional-evidence + co-segregation layer (decompose + --plan, 2026-06-25)

From `/hypothesise` (8 grey-zone hypotheses) → the subset useful for what we're actually doing (the PGx
decoder, currently CYP2C19/CYP2C9/VKORC1 validated vs GeT-RM). **Planning artifact — NOT executed.** Per the
planning-STOP discipline, this is decompose + a technical plan; build is a separate greenlight.

## The thinking — which hypotheses are useful, and why
The PGx cells report the diplotype-CALLING as independently validated (GeT-RM 72/72, 73/73) but the
metabolizer PHENOTYPE as **"faithful-to-CPIC"** — and the per-allele **function assignments ARE CPIC's**.
That function-assignment step is the one **circular** link (the caller and the "truth" both come from CPIC).
The hypotheses that attack exactly that, with free data:

- **H3 (ΔΔG / structure → independent functional label)** + **H8 (eQTL measured-expression label)** —
  together an **independent functional-evidence layer**: a NON-CPIC signal for each allele's function.
  Missense alleles → variant-effect predictors (REVEL/AlphaMissense/CADD via dbNSFP); regulatory alleles
  (promoter/expression: CYP2C19 *17, VKORC1 -1639) → GTEx liver cis-eQTL direction (a genuinely MEASURED
  independent signal). Breaks the circularity gate — the project's prized discipline.
- **H7 (pedigree co-segregation)** — a cheap free internal QC: the 1000G 602 trios let us check Mendelian
  consistency of the diplotype calls (a child's alleles must descend from the parents'). Strengthens the
  existing GeT-RM validation at near-zero cost.
- **H1 (population-control negatives)** — banked as a TECHNIQUE for the next S-starved cell (PGx isn't
  S-starved — GeT-RM gives both classes), not built now.
- H2/H4/H5/H6 — marginal / high-noise / already-have / animal-catalog-only (the dog wall). Not in scope.

## Grounded data sources (all free, verified 2026-06-25)
| Source | What | Access |
|---|---|---|
| **dbNSFP 4.9a** (academic) | per-variant REVEL / AlphaMissense / CADD / SIFT / PolyPhen, GRCh38, tabix-indexed | dbnsfp.org (Box/Amazon mirror); the academic "a" branch (the "c" branch strips REVEL/CADD/PolyPhen) |
| **GTEx V8** `Liver.v8.signif_variant_gene_pairs.txt.gz` | significant liver cis-eQTLs (variant→gene, effect size/slope) | gtexportal.org/home/datasets (open). IDs: CYP2C19 ENSG00000165841, CYP2C9 ENSG00000138109, VKORC1 ENSG00000167397. eQTL Catalogue API = alt |
| **1000G** `20130606_g1k_3202_samples_ped_population.txt` | 602 trios (family/paternal/maternal/sex/pop) | 1000genomes EBI FTP / S3 (open) |

## Decompose — three units
- **Unit A (PRIMARY — the circularity-break): `dna_decode/pgx/functional_evidence.py` + a report.**
  For each catalogued PGx allele, attach the independent signal and a concordance verdict vs CPIC's function:
  - missense (CYP2C9 *2 R144C / *3 I359L / *5 / *8 / *9 / *11; CYP2C19 *2 is splice, *3 is stop — handled
    by category, see below) → dbNSFP score at the GRCh38 coord (REVEL/AlphaMissense; "damaging" ⇒ supports
    decreased/no-function).
  - regulatory/expression (CYP2C19 *17 rs12248560 ↑; VKORC1 rs9923231 ↓) → GTEx liver cis-eQTL presence +
    SLOPE direction vs the expected ("increased function" ⇒ +slope on the ALT; "sensitive/↓expr" ⇒ -slope).
  - splice/stop (CYP2C19 *2 splice, *3 stop; CYP2C9 *6 frameshift) → LoF-by-consequence (no predictor needed;
    consequence class is itself the independent mechanism).
  - Output `wiki/pgx_functional_evidence_2026-06-25.{md,json}`: per-allele {CPIC function, independent signal,
    AGREE / DISAGREE / NO-SIGNAL}. **Honest tier:** dbNSFP predictors are ML (independent of CPIC curation,
    NOT ground truth → a cross-check); GTEx eQTL is a real MEASURED signal for the regulatory alleles. This is
    an evidence-ANNOTATION layer (small-N, per-allele), not a concordance-% claim.
- **Unit B (FAST adjacent — co-segregation QC): `scripts/pgx_trio_concordance.py`.**
  Parse the 1000G ped → the trios whose 3 members are in the cohort VCFs → run the PGx caller on each →
  assert Mendelian consistency (child allele set ⊆ {one allele from each parent}); count violations.
  Output `wiki/pgx_trio_mendelian_2026-06-25.{md,json}`. A clean trio-consistency rate strengthens the
  CALLING claim independently of GeT-RM. (Uses the CYP2C19 + CYP2C9 region VCFs already on disk.)
- **Unit C (BANKED — not built): population-control negative-class helper** for the next S-starved cell
  (`dna_decode/eval/`-side; H1). Documented technique, deferred.

## Plan (Unit A primary; Unit B fast win) — ordered steps
1. **Data fetch (one-time, gitignored):** dbNSFP 4.9a (or just the PGx-region slices via tabix to keep it
   small — the genes are a few kb each), GTEx Liver signif-pairs (filter to the 3 Ensembl IDs → a tiny TSV,
   committable), the 1000G ped (tiny, committable). Care-check disk before the dbNSFP download (it is large;
   prefer tabix-slicing the gene regions, not the full ~30 GB).
2. **Unit A catalog wiring:** add an `independent_evidence` block per allele to the catalogs (or a sidecar
   `functional_evidence.py` keyed by rsID) — variant consequence class (missense/splice/stop/regulatory) +
   the looked-up dbNSFP score / GTEx eQTL slope + the AGREE/DISAGREE verdict vs `ACTIVITY_VALUE`/function.
3. **Unit A report builder + tests:** `build_functional_evidence()` → the report card; unit tests pin the
   verdict logic on committed mini-fixtures (a damaging-REVEL missense ⇒ AGREE with no-function; a +slope
   eQTL ⇒ AGREE with increased-function for *17). Offline-safe (skip live dbNSFP/GTEx fetch when absent).
4. **Unit B trio QC:** `scripts/pgx_trio_concordance.py` + 6 tests (ped parse, trio extraction, Mendelian
   rule, a planted violation). Run on the on-disk CYP2C19 + CYP2C9 region VCFs → the trio-consistency number.
5. **Verify-in-batch:** run both → inspect the per-allele evidence table + the trio number; confirm the
   *17→+eQTL and VKORC1→-eQTL directions actually resolve in GTEx (if GTEx liver lacks the eQTL, report
   NO-SIGNAL honestly — do not force agreement).
6. **Docs + ship:** CHANGELOG + report cards; FROZEN AMR surface untouched; commit. (No PyPI bump required —
   this is a validation/evidence layer, not a new CLI surface; optional minor bump if a `--functional-evidence`
   flag is added.)

## Success bar / falsifier
- **Bar (Unit A):** every catalogued PGx allele gets an independent-evidence verdict (AGREE/DISAGREE/NO-SIGNAL)
  vs its CPIC function, with ≥1 regulatory allele confirmed by a real GTEx liver eQTL in the expected direction
  (the measured-signal proof-point). **Bar (Unit B):** trio-Mendelian-consistency computed on the on-disk
  cohorts with the violation count reported.
- **Falsifier / honest red flag:** if the independent signals DISAGREE with CPIC on a core allele
  (e.g. REVEL benign for *3-no-function, or no *17 eQTL), that is a REAL finding to surface — it does not
  invalidate the cell (CPIC is the clinical standard) but it flags where "faithful-to-CPIC" rests on softer
  evidence. Do NOT tune the thresholds to manufacture agreement (that would re-introduce circularity).

## Why this is the right next move
- It attacks the **deepest unvalidated link** in the PGx cells (the function assignment) with **free,
  independent, non-CPIC** signals — the circularity-break the project repeatedly prizes (and the reason
  several prior tracks were closed). It needs NO new wet-lab data and NO money.
- Unit B is a near-free QC that strengthens the calling claim using trios already in the 1000G data on disk.
- Both are small, bounded, offline-testable, and leave the FROZEN surfaces untouched.

## Status: PLANNED (decompose + technical plan). Build (Units A+B) is a separate greenlight.
