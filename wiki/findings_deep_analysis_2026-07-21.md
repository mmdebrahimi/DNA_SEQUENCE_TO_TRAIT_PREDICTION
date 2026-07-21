# dna_decode — Deep synthesis of all findings (2026-07-21)

**Scope:** a cross-cutting analysis of every finding across the platform — the validated cells (bacterial /
viral / fungal / molecular), the closed negatives, the infrastructure, and the strategic state. Synthesized
from 4 parallel reads over the report card, lineage metrics, the negative-results map, CLAUDE.md, and this
session's AR-Bank artifacts. **Honesty discipline preserved: there is deliberately NO aggregate "X% works"
headline** — every cell keeps its own tier. Numbers are cross-checked across sources.

---

## 0. The one-paragraph honest state

dna_decode is a **genotype→phenotype decoder that works — cell by cell, in the regimes where a curated
determinant catalog is the right tool, and only where a free wet-lab-measured label exists to prove it.**
The validated product is a **deterministic determinant-scan engine + a supervised learned complement**, NOT
a foundation-model embedding predictor (that bet is a closed negative across the kingdom boundary). The
single engine transfers across **bacteria → fungi → viruses → protozoa**. Its strongest evidence is HIV
(25 cells against a free, independent, isolate-level wet-lab label) and cross-kingdom fungal C. auris. Its
binding constraint — reconfirmed at every wall — is **LABELS, not models or compute.** The platform's real
moat is not any model; it is a **mature, self-skeptical validation stack** that makes every claim inspectable,
leakage-clean, and clonality-disclosed, so the project can extend honestly wherever a clean free label opens.

---

## 1. Unified scorecard — the validated cells, tiered by rigour

Ordered by **independence tier** (the honest axis), not by headline number. Independence tiers:
**IND-EXT** = free independent isolate-level wet-lab label · **IND-of-interp** = independent of the
competitor tool's *interpretation* but in-distribution of its knowledge base · **AR-Bank** = CDC-measured
BMD-MIC, provenance-disjoint but not methodology-independent · **IN-DIST** = catalog-derived (lower rigour) ·
**MARKER** = marker-presence validated, no phenotype label.

### A. The flagship — HIV-1 (IND-of-interp; the wall that fell)

The central "labels-not-models" wall fell the moment a free, independent, isolate-level wet-lab
genotype↔phenotype label existed: **Stanford HIVDB PhenoSense fold-change** (non-circular vs HIVDB's own
Sierra rules-engine). **25 cells, 5 classes, 4 genes.** First-generation drugs are excellent:

| drug | class | AUC | drug | class | AUC |
|---|---|---|---|---|---|
| nevirapine | NNRTI | **0.985** | atazanavir | PI | 0.957 |
| lamivudine | NRTI | **0.975** | nelfinavir | PI | 0.948 |
| abacavir | NRTI | **0.974** | indinavir | PI | 0.929 |
| efavirenz | NNRTI | **0.962** | lopinavir | PI | 0.933 |
| didanosine | NRTI | 0.919 | raltegravir | INSTI | 0.905 |

The **NRTI→PI→INSTI mutant-specific deconfounding arc is COMPLETE** (v0.1 improves 8/8 PI +0.056, 5/5 INSTI
+0.087, lifting specificity without losing sensitivity). Within-subtype de-confounding cleared all 4 classes
→ **the catalogs decode MECHANISM, not subtype population structure.** The weak spots are honest and
predicted, not hidden: 2nd-gen DTG 0.745 / BIC 0.846 / doravirine 0.555 are lower *exactly* as the class-level
over-call model predicts; cabotegravir AUC 1.0 is UNDERPOWERED (n_S=4). **Caveat that must ride every HIV
claim:** independent of interpretation, but in-distribution of the HIVDB knowledge base — not a provenance-
disjoint external cohort.

### B. Cross-kingdom fungal — C. auris (IND-EXT via CDC AR-Bank measured MIC)

Genuine cross-kingdom independent validation — a **2nd mechanism AND 2nd gene**, not just a 2nd azole:
- **ERG11 / fluconazole (azole):** POWERED; binary 1.00/0.714; **9/9 perfect on the mechanism-attributable
  subset.** Surfaced a clade-IV lineage confound (haplotype K177R/N335S/E343D identical in 1R+2S → zero
  discriminative signal) — the *fungal analogue of the bacterial QRDR-vs-lineage trap*, handled identically
  (confidence tier, R preserved, weak calls made visible).
- **FKS1 / micafungin (echinocandin):** POWERED + **SCORED_ENDORSED** (sens 0.60 / spec 1.00 / acc 0.846).
  3/5 R caught via canonical S639F/P/Y; the uncatalogued-variant disclosure separated a *real catalog gap*
  (F635C — a substitution at a catalogued hotspot where the catalog lists only F635del; would lift sens
  0.60→0.80, deliberately NOT auto-added on n=1) from an *isolate-specific quirk* (W691L, non-hotspot).

### C. New species this session — N. gonorrhoeae (AR-Bank; the headline win)

The first whole-new-species validation; the Kaggle AMRFinder factory transferred cleanly to a new
Gram-negative:
- **ciprofloxacin** (gyrA/parC QRDR): **sens 1.00 / spec 1.00 / acc 1.00** (11R/9S) — perfect out of the box.
- **cefixime** (penA mosaic-34): **sens 0.917 / spec 1.00 / acc 0.95** (12R/8S) after a **v0.1 fix that lifted
  spec 0.0→1.0** — v0 fired on the shared A510V/F504L markers (which the reduced-susceptibility S isolates
  also carry); v0.1 requires ≥3 of the 4 mosaic-penA-34 core markers {I312M, V316T, N512Y, G545S} (11/11 R,
  0/8 S). One disclosed FN (a non-mosaic penA_D346DD high-MIC path). Mirrors the ceftriaxone v0.1 narrowing.

### D. The frozen 10 SCORED cells — NCBI-PD provenance-disjoint (the banked baseline)

The terminal honest product of the public-label AMR track (frozen 2026-06-13, commit b3761c8). **Read the
lineage-collapsed column, not the raw** — raw sens/spec is clonality-inflated ~2–2.5×.

| organism | drug | raw sens/spec | lineage sens/spec [Wilson CI] | honest grade |
|---|---|---|---|---|
| Campylobacter | cipro | 1.00/1.00 | 1.0/1.0 (15R/14S eff) | **cleanest — endorsed both raw+lineage** |
| E. coli | ceftriaxone | 0.967/0.967 | 1.0/1.0 (11R/17S) | endorsed |
| E. coli | tetracycline | 0.933/0.933 | 0.882/1.0 (17R/19S) | endorsed (robust both sides) |
| E. coli | cipro | 0.933/**0.70** | **0.5**/0.8 (**4R** eff) | SCORED, weak — R collapses to 4 lineages |
| E. coli | gentamicin | 0.90/1.00 | 0.6/1.0 (**5R** eff) | SCORED, R clonal |
| Klebsiella | ceftriaxone | 1.00/0.90 | 1.0/0.95 (16R/21S) | endorsed (ESBL) |
| Klebsiella | tetracycline | 0.80/0.967 | 0.842/0.963 (19R/27S) | endorsed |
| Klebsiella | gentamicin | 0.933/0.933 | 1.0/0.857 (11R/7S) | endorsed (S-side small) |
| Klebsiella | cipro | 0.967/0.967 | **0.5**/1.0 (**2R** eff) | SCORED — raw is clone inflation (2 eff lineages) |
| Klebsiella | meropenem | **0.467**/0.90 | 1.0/0.952 (6R eff) | **raw NOT_ENDORSED** — 16 FN, determinant-invisible carbapenem-R |

Plus INDEPENDENTLY_VALIDATED (EBI AMR-Portal, metadata-only amendment): Salmonella cipro (acc 0.959, N=24,972)
+ Klebsiella cipro (spec 0.994, N=4,385).

### E. M. tuberculosis (NON-FROZEN; WHO catalogue v2 rule)

The first genuinely-INDEPENDENT (out-of-CRyPTIC-build) TB number:

| drug | in-dist (CRyPTIC) lineage | **independent (AMR-Portal) lineage [CI]** |
|---|---|---|
| RIF (rpoB) | 0.41 | **sens 0.444 [.246–.663] / spec 0.979 [.889–.996]** |
| INH (katG+inhA) | 0.349 | **sens 0.321 [.179–.507] / spec 0.972 [.858–.995]** |

**Key finding: the independent lineage number ≈ the in-distribution lineage number → the WHO-catalogue rule
does NOT degrade out-of-distribution at the lineage level; spec is high + robust (~0.97) in both.** The raw
0.92/0.88 was clonality inflation (2,845 isolates → ~67 barcode-lineages).

### F. Molecular — forward / inverse DMS (the one regime where LEARNED wins)

The fitness-aligned regime where a learned model beats a catalog:
- **forward ESM2-650M:** ProteinGym median Spearman **0.484–0.490** (217 assays); TEM-1 β-lactamase 0.7315.
  **ESM2 PEAKS at 650M — 3B (0.467) / 15B (0.438) REGRESS.**
- **forward hybrid (rank-average):** ESM2+GEMME+ProSST **+0.056 vs ESM2-650M, win 90.5%** — the headroom is
  **MODALITY (seq⊕evolution⊕structure), not parameters.**
- **inverse `propose_edits`:** validated as a **RANKER** (blaTEM +53.0% on measured DMS) — ranks, not doses
  (calibrators can't transfer); Regime B only, never clinical resistance.

### G. Lower-rigour / honestly-tiered
- **SARS-CoV-2 Mpro** — IN-DIST (catalog + fold both from CoV-RDB) AND underpowered (37R/5S, spec 0.0
  uninformative). Not the independent win HIV is. (Mutant-level is load-bearing: Omicron P132H correctly S.)
- **Influenza NA** — MARKER-only (osel-R/zana-S drug-specificity on real isolates, no IC50 phenotype).
  Breadth-via-proven-pattern, explicitly not a new finding.

---

## 2. The meta-pattern — three G2P regimes (classify BEFORE building)

The single most useful organizing principle the project produced. **Classify a candidate into one of three
regimes before spending any labor:**

1. **Curated-catalog regime → DETERMINISTIC RULES WIN.** A known set of causal loci + a wet-lab label. All
   the SCORED cells above live here. The determinant-scan engine transfers across kingdoms unchanged.
2. **Organism-polygenic regime → NEITHER WORKS.** Distributed/polygenic phenotype, no compact catalog.
   Learned is 0-for-5, embeddings 0-for-4 (they learn lineage/population structure, not mechanism).
3. **Molecular-property regime → LEARNED WINS, but ONLY when the phenotype is fitness-aligned.** ESM2 on
   ProteinGym DMS works. **Antagonistically-selected drug resistance INVERTS it** — ESM2 0.454 *below chance*
   on HIV resistance; mutation-count beats ESM on 10/11 drugs; resistance is reached via chemically
   conservative single-nt substitutions that every likelihood/exchangeability scorer (even BLOSUM62, which
   never saw HIV) mis-calls as benign. The blindness is a property of the **phenotype**, not model capacity.

---

## 3. The two dominant caveats that govern how to read EVERYTHING

### Caveat 1 — Clonal inflation (the pervasive one)
Every SCORED bacterial R-class is clonally dominated (≤17 effective lineages; several at 2–6). **Raw sens/spec
is one vote per ISOLATE, not per lineage → inflated ~2–2.5×.** Klebsiella cipro R = 60 isolates ≈ 3 lineages;
TB raw 0.92 → lineage 0.44. This is *disclosed, never demoted* — the lineage-collapsed column with its Wilson
CI is the real evidence, and effective-N is often tiny, so read the CI. **Any headline sens/spec quoted
without its lineage-collapsed companion is misleading by ~2×.**

### Caveat 2 — Labels, not models (the binding one)
Reconfirmed at every wall. Models are good enough; compute is not the wall (a negative de-confounded metric
is a signal-vs-structure problem, NOT a window-budget one — do NOT scale on a bigger GPU). **The wall is a
free + independent + isolate-level + wet-lab-measured + provenance-separable label, and that quadrant is
narrow.** Every closed expansion is label- or regime-blocked, not model-blocked.

### The ONE nuance that must NOT be over-generalized (user-corrected twice)
**"AI failed / 0-for-N" is scoped to ZERO-SHOT / frozen embeddings ONLY.** The shipped architecture is a
**hybrid: deterministic catalogue + SUPERVISED learned complement.** A supervised model *rescues* the
catalog's blind spot (0.81 leave-study-out vs ESM zero-shot 0.449) — but ONLY on **convergent-evolution**
pathogens (HIV); on **clonal** organisms (TB) it collapses to 0.51 ≈ chance. Carry the scope in the sentence
every time: zero-shot vs supervised; convergent vs clonal.

---

## 4. The closed negatives (do-not-reopen) + the 8 gates

**Do-NOT-reopen (recorded so labor isn't re-spent):** foundation-model embeddings (0-for-4, learn lineage
not mechanism); Arabidopsis flowering-time embedding (H2 falsified, the designed best-case); self-supervised
catalog (circular); learned variant scorer to fill the catalog blind spot (phenotype-property ceiling); ESM2
3B/15B scaling (peaks at 650M); proximity/one-step DRM forecaster (base-rate ceiling). **Label-blocked (needs
a non-public source):** pathotype (G1/G2/G3), MIC-continuous (G1/G6/G8), Salmonella grid (G4).

**The 8 reusable rejection GATES — screen every candidate dataset BEFORE building:**

| Gate | Trips when |
|---|---|
| **G1 Circular label** | phenotype is produced by a genomic tool the decoder would compete against |
| **G2 Study==class** | label confounded with source study/submitter |
| **G3 Sampling-defined** | label IS the sampling context, not a measurement |
| **G4 Surveillance domination** | excluding the surveillance ecosystem collapses the positive pool <20/class |
| **G5 Assembly attrition** | label-bearing records lack downloadable assemblies |
| **G6 Phenotype censoring** | quantitative label interval-censored at the breakpoint |
| **G7 Provenance not separable** | metadata too thin for a leakage-clean provenance-disjoint split |
| **G8 Dedup collapses balance** | clonality correction drops a class below ~3 effective lineages |

---

## 5. The infrastructure moat (the platform's real asset)

Not any model — a **mature anti-self-deception stack** that lets the project extend cell-by-cell without
fooling itself:
- **Leakage registry** (`cohort_manifest.py`) — exact-identity exclusion across all cohorts, fail-closed.
- **Clonality/lineage disclosure** (`clonality.py`) — greedy-representative Mash clustering (chaining-
  resistant), mixed-label clones excluded as DISCORDANT, Wilson CI + effective-N on every weighted point.
- **`canonical_cell_key`** — the single join key; the reason the external arm got a separate namespace (a
  shared-key overwrite of the frozen E. coli cell was caught in pre-save brainstorm).
- **Prospective-lock** (`prospective_lock.py`) — sha256-pins the frozen surface, `lock_date` cutoff makes any
  post-lock isolate leakage-free by construction; tamper-evident. First sweep = a *genuine accruing zero*
  (ingestion lag, not decoder failure).
- **External-cohort arm** — namespace-separate, BioSample-level leakage check, powering gate, operator-aware
  MIC censoring.
- **This session's addition:** the **empty-assembly integrity gate** (`MIN_DETERMINANTS_FOR_VALID_ASSEMBLY=3`)
  — an empty assembly is INDETERMINATE, never scored S-by-absence. It caught a *false SCORED_ENDORSED*
  (vancomycin sens had cratered to 0.2 from empty-assembly R isolates counted as S). A new instance of the
  recurring "verify the numbers mean what they claim" discipline.

---

## 6. Strategic frontier — what's genuinely open vs gated

**The public-free-label AMR track is BANKED (frozen).** The platform is NOT frozen — it extends sideways into
new species/kingdoms via free reference-lab banks (CDC AR-Bank), which is forward-path #1 on its free tier,
not a reopening of closed negatives.

**Genuinely open, free + reversible, executor-runnable (ranked):**
1. **NARMS overlap join** — highest-VOI/cheapest. Access confirmed reachable (FDA MIC file + ENA/NCBI WGS,
   not paywall). The whole question is the net-new count vs the 744-accession footprint (NARMS flows into
   NCBI-PD → G4 risk). A scoped metadata-join, not a guess.
2. **CDC AR-Bank non-Enterobacterales panels** — reuse the existing scraper (Pseudomonas / Acinetobacter /
   Aspergillus / Candida); nearly-free given current infra.
3. **MaveDB de-dup + AMR-enzyme DMS** — expands the molecular forward/inverse cell (ProteinGym is a subset).
4. **GISP/Euro-GASP gonococcus** — deepen the just-launched N. gonorrhoeae cell from public ENA/PubMLST.
5. **Prospective-lock re-sweep** — periodic; accrues as post-2026-06-13 isolates appear.

**Authority/money-gated (user decision, NOT executor):** any physical AR-Bank isolate request (biosafety
attestation); WWARN malaria IPD (DAC-gated — the one true acquisition); pharma programs (ATLAS is the wrong
shape — measured MIC but no WGS).

**Genuinely outside-the-box (needs a new cell design):** BASEL phage-susceptibility — a wet-lab EOP plaque
phenotype that dodges G1 AND G3 by construction. The one candidate that escapes the AMR label quadrant
entirely; a different modeling shape (phage×host matrix).

---

## 7. Tensions & things to watch (the honest loose ends)

1. **The "SCORED_ENDORSED" label carries very different weight across cells.** Campylobacter cipro (15 eff
   lineages, endorsed raw+lineage) and Klebsiella cipro (2 eff R-lineages, raw is clone inflation) both read
   "SCORED" — the lineage-effective-N is the disambiguator and must be read every time. Consider whether the
   report card should surface effective-N more prominently than the raw headline.
2. **Klebsiella meropenem is the honest scar** — raw sens 0.467 (16 determinant-invisible carbapenem-R FN).
   The lineage-rescue to 1.0 rests on 6 lineages. This is the hard FN ceiling of gene-presence decoding
   (expression/porin-loss resistance) and should not be smoothed over.
3. **AR-Bank cells are one-sided** (resistance-enriched bank → powers one class per drug); "sens 1.0" on a
   0-S-scored drug (Enterococcus doxycycline) is not a two-sided validation.
4. **The new-species cefixime v0.1 rule is derived + validated on the SAME cohort** (like ceftriaxone v0.1).
   It's literature-grounded (mosaic penA-34) but a held-out gonococcal cohort (GISP) would upgrade it from
   "cohort-fit + literature-aligned" to "externally confirmed."
5. **Registry-count vs report-card-count divergence:** the frozen AMR report card = 27 cells; the
   cross-kingdom certification-capstone registry = ~84 cells (AMR+viral+PGx+typing). These are different
   scopes — don't conflate. (The capstone file is being actively edited by a parallel session; its exact
   count is live.)
6. **The whole track is provenance-disjoint stress-testing, NOT methodology-independent** — most labels are
   CLSI BMD and the caller is the same AMRFinder `-O` + frozen `call_resistance`. A truly method-independent
   validation (a different phenotyping method + a different caller) has not been done and would be the next
   rigour tier.

---

## 8. Bottom line

The platform has quietly become **the honest thing it set out to be: a working DNA decoder tool** — a single
determinant-scan engine, validated cell-by-cell across four kingdoms, wrapped in a validation stack rigorous
enough to catch its own false victories (three this session alone: the empty-assembly gate, the cefixime
over-call, the FKS1 blind-spot-vs-quirk separation). Its strongest evidence (HIV, 25 cells against a free
independent label) proves the thesis that the wall was never the model — it was the label. Everything
downstream of a clean measured-phenotype⋈genotype set is already built. **The frontier is finding/fetching
more of the narrow "free wet-lab measured, provenance-separable" label quadrant — and the two cheapest doors
(NARMS overlap, CDC AR-Bank non-Enterobacterales) are open and free right now.**
