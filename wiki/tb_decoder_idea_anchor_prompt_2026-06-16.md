# `/idea-anchor` prompt — extend the deterministic decoder to M. tuberculosis (CRyPTIC) (2026-06-16)

**Why this exists:** the acquisition fork's feasibility probe came back GREEN
(`wiki/cryptic_feasibility_probe_result_2026-06-16.md`): CRyPTIC is a free, all-gates-passing TB
measured-MIC substrate (12,287 isolates × 13 drugs, reference UKMYC broth-microdilution MIC + per-isolate
VCF), all 13 drugs powered, and the POINT-mutation determinant rule transferred in a PoC (RIF rpoB-RRDR
sens/spec 1.0 on 15R/15S). The next step is a real build: the first TB decoder cell. `/idea-anchor` is
user-only → this file holds the paste-ready command + Soraya's drafted framing for you to RATIFY.

---

## Paste-ready command

```
/idea-anchor Extend the deterministic AMR decoder to a new organism, M. tuberculosis, using the free CRyPTIC measured-MIC compendium.

PROBLEM: The deterministic E. coli/Klebsiella/S. aureus/C. auris decoder counts CURATED resistance
determinants -> R/S and is strongest on POINT-mutation mechanisms (cipro QRDR 0.925, Klebsiella 1.0). The
public-label AMR track on existing organisms is saturated; the binding constraint is LABELS. The CRyPTIC
consortium provides a FREE, deep, sampling-independent TB substrate that no prior decoder organism matched:
12,287 M. tuberculosis isolates, 13 drugs, reference broth-microdilution MIC + binary phenotype + quality
flag + per-isolate VCF (variant calls vs H37Rv NC_000962.3) — all on the EBI FTP, already downloaded to
data/raw/cryptic/. A feasibility probe confirmed: all 13 drugs clear >=20/class at HIGH quality (rifampicin
3448R/5507S, isoniazid 4467R/5051S), the VCF genotype is fetchable + parseable with NO assembly/Docker, and
a PoC rpoB-RRDR rule scored RIF at sens/spec 1.0 on a 15R/15S sample.

GOAL: Build the FIRST TB decoder cell — a curated TB determinant rule (from the WHO M. tuberculosis mutation
catalogue) scored on a CLONALITY-CORRECTED CRyPTIC cohort, starting with rifampicin (rpoB) then isoniazid
(katG + inhA promoter) — extending the deterministic decoder to a new organism in its strongest mechanism
class, with an HONEST validation that accounts for the catalogue-vs-cohort circularity below.

THE LOAD-BEARING RISK (must be designed around, not ignored): the WHO TB mutation catalogue was built PARTLY
FROM CRyPTIC DATA. So scoring a WHO-catalogue rule on CRyPTIC is NOT a clean independent validation — it is
the KNOWLEDGE-BASELINE measuring how well the curated catalogue recovers its own training phenotypes (the
project's recurring "validate the wrapper vs the underlying tool on INDEPENDENT data" + "circular label"
lessons). The cell is honest ONLY if this is named and mitigated: either (a) validate the rule on an
INDEPENDENT non-CRyPTIC TB cohort, or (b) hold out CRyPTIC isolates/mutations not used to build the
catalogue, or (c) frame the CRyPTIC number explicitly as the in-distribution knowledge-baseline (like
AMRFinder-on-its-own-DB) with the independent number deferred.

CONSTRAINTS: solo hobby project; free data only (no money); the deterministic decoder + reproducibility
freeze (2026-06-13) must stay byte-unchanged unless a deliberate frozen-surface decision is made (prefer an
organism-routed extension / overlay first, like the TMP-SMX experimental cell). North star = a DECODER TOOL.
HONESTY: this expands the DETERMINISTIC decoder to a new organism; it does NOT reopen the learned-embedding
niche (TB AMR has a curated catalog, so embeddings would lose to the determinant rule here too, as on E.
coli). Do NOT propose an embedding/learned approach for this cell.

DELIVERABLE the downstream chain should produce (NOT this idea-anchor): a scored first TB cell (RIF, then
INH) on a clonality-corrected CRyPTIC cohort, with the catalogue-circularity caveat explicit and the
determinant rule living in a non-frozen organism-routed extension until a promotion decision.
```

---

## Soraya's drafted anchor framing (ratify or redirect)

**Formal rephrase.** Stand up the first M. tuberculosis decoder cell by porting the curated-determinant
rule to TB (WHO mutation catalogue → VCF-position/codon match), scoring rifampicin then isoniazid on a
clonality-corrected CRyPTIC cohort, with an honest treatment of the WHO-catalogue-built-from-CRyPTIC
circularity, and without touching the frozen E. coli surface.

**Fundamental clarifications (the ~3 the skill will likely ask — drafted answers):**
1. **Determinant source: WHO catalogue vs tb-profiler vs AMRFinder-TB?** → *Draft: the WHO M. tuberculosis
   mutation catalogue v2 (2023)* — the authoritative, expert-graded TB resistance-mutation list (use the
   "Associated with resistance" grades 1+2). It's the TB analogue of the curated-AMRFinder approach. Risk
   of alternatives: tb-profiler bundles its own variant caller (wrapper-trap — we'd be testing their caller,
   not our rule); AMRFinder-TB coverage is thinner. **Carries a hidden authority decision** — the WHO
   catalogue's CRyPTIC provenance is the circularity risk above; ratifying "WHO catalogue" also means
   ratifying the mitigation (independent cohort / held-out / baseline-framing). *Yours.*
2. **Clonality metric for TB (no assemblies — only VCF + reads)?** → *Draft: SNP-distance between the
   per-isolate VCFs* (the standard TB lineage/transmission metric) feeding the existing
   `clonality.py` greedy-representative clustering — NOT Mash (which needs assemblies; SNP-distance from the
   aligned-to-H37Rv VCFs is the natural, cheaper TB clonality signal). *Technical — my call to draft.*
3. **First-cell scope + acceptance bar?** → *Draft: rifampicin first (rpoB, single-gene, 3448R — cleanest),
   then isoniazid (katG + inhA promoter — tests the 2-locus path).* Acceptance bar = lineage-weighted
   sens/spec with a Wilson CI on the clonality-corrected cohort, reported alongside the in-distribution
   baseline number + the circularity caveat. **The bar value is an authority decision.** *Yours.*

**Current assumptions (flagged for test):**
- The WHO catalogue is obtainable machine-readable (it is — WHO publishes an Excel/CSV; needs a fetch).
- VCF positions align to catalogue coordinates (both vs H37Rv NC_000962.3 — verify, don't assume).
- TB resistance is determinant-visible. *True for RIF/INH/fluoroquinolones; WEAKER for newer drugs
  (bedaquiline/clofazimine — Rv0678 efflux + novel mechanisms) = the same expression/efflux blind-spot the
  E. coli decoder has. Scope the first cells to the determinant-visible drugs.*
- "Provenance-disjoint" — TRUE for leakage (no TB tuning set), but the real gate-2 risk is the
  **catalogue-built-from-CRyPTIC circularity**, not accession overlap. Don't conflate the two.

**Blunt opinion.** Strong move: free data, the decoder's best mechanism class, power that dwarfs E. coli,
and the rule already transferred in the probe. The one thing that can make it dishonest is the
**circularity** — a "0.99 on CRyPTIC" headline would be the WHO catalogue recovering its own training
phenotypes, not an independent validation. Name it, and either find an independent TB cohort or frame the
CRyPTIC number as the knowledge-baseline. Also: keep it in a non-frozen organism-routed extension (the
TMP-SMX pattern) — do not edit the frozen E. coli surface for a new organism. And it's worth restating
plainly: this grows the *deterministic* tool's breadth; it is not the learned-decoder experiment.

**Recommended next step.** `/probe` — this is code-touching (a new organism route + a VCF-position/codon
determinant rule, distinct from the AMRFinder-Class plumbing + the clonality-on-VCFs path) AND needs the
circularity-mitigation cohort design grounded before a plan. First probe actions: how the decoder routes
organisms today, where a TB rule should live (frozen vs organism-routed extension), the VCF-vs-AMRFinder
determinant plumbing, and whether an independent (non-CRyPTIC) TB measured-MIC cohort exists for the
honest validation.

```
/probe tb-decoder-cryptic
```
