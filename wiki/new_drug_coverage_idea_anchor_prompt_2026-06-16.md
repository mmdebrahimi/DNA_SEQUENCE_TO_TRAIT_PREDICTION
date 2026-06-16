# Next-move `/idea-anchor` prompt — expand decoder drug coverage from in-hand measured MIC (2026-06-16)

**Why this exists:** Distinct from the 2026-06-13 next-epoch anchor (which targets ACQUIRING a *new label
source* — labels are the constraint there). This move is the opposite: the **labels already exist on
disk**. Oxford (8 measured-MIC drugs) + Sci234 (19) carry measured MIC + complete genotype for ~2-3× more
drugs than the decoder covers. The shipped fam.tsv Subclass resolver (`scripts/fam_subclass_resolver.py`,
2026-06-16) just proved out by scoring gent+cef on Sci234 from a gene-presence table with no assembly. So
NEW drug coverage = **rule-definition + validation on data we already have**, not acquisition. This file
holds the drafted `/idea-anchor` (user-only skill → paste-ready) + Soraya's drafted framing for the user
to RATIFY (draft-then-ratify), since I can't self-invoke it.

---

## Paste-ready command

```
/idea-anchor Expand the deterministic AMR decoder's drug coverage using measured MIC already in hand.

PROBLEM: The frozen v0.5.0 deterministic decoder covers 6 drugs (ciprofloxacin, ceftriaxone, tetracycline,
gentamicin, meropenem, oxacillin). Two fully-independent measured-MIC cohorts already on disk carry MIC +
complete genotype for many MORE drugs: Oxford (ampicillin, cephalexin, cefuroxime, cefepime, levofloxacin,
co-amoxiclav, trimethoprim-sulfa, piptaz) + Sci234 (amikacin, amox/amox-clav, aztreonam, cefepime,
cefotaxime, ceftazidime, ceftazidime-avibactam, colistin, ertapenem, imipenem, piperacillin, pip-tazo,
temocillin, tobramycin, trimethoprim-sulfa). The labels are NOT the constraint here. The question is which
new drugs are AMRFinder-rule-TRACTABLE (acquired-gene / POINT-driven, in the curated determinant set the
decoder reads) vs decoder-BLIND (efflux / porin / expression / enzyme-inhibitor combos), and which have
enough resistant isolates IN THESE COHORTS to actually validate.

GOAL: Add the tractable, adequately-powered subset of new drugs as honestly-validated decoder cells —
each with a CLSI/EUCAST breakpoint, an AMRFinder class/subclass filter, a determinant-count rule, and a
binary R/S validation against the in-hand Oxford/Sci234 MIC — while abstaining explicitly on the
intractable ones, all WITHOUT regressing the 6 existing frozen cells.

TRACTABILITY READ (to be confirmed/ranked downstream, not decided here):
  - CLEAN (acquired-gene, curated, no intrinsic trap): trimethoprim-sulfamethoxazole (sul/dfr — needs NO
    fam.tsv resolver, the cleanest first win, MIC in BOTH cohorts); cefotaxime/ceftazidime/cefepime/
    aztreonam (ESBL blaCTX-M/CMY — the validated ceftriaxone rule PATTERN, via the resolver); amikacin
    (aac(6')/rmt — aminoglycoside-subclass refinement, via the resolver).
  - CARBAPENEM-EXTEND: ertapenem/imipenem (carbapenemase, the meropenem rule extends) — but carbapenem-R
    is RARE in E. coli, so likely born UNDERPOWERED in these cohorts (cf. the gent 0-R outcome).
  - CAVEATED / PARTIAL: colistin (mcr acquired BUT mostly chromosomal mgrB/pmrB → FN-prone, half-blind);
    co-amoxiclav / pip-tazo (β-lactamase + inhibitor interaction → NOT a clean gene count).
  - ABSTAIN-BY-DESIGN: ampicillin/amoxicillin (intrinsic blaEC makes ~all E. coli R → uninformative, the
    'intrinsic genes break broad class-rules' trap); ceftazidime-avibactam / temocillin (novel/complex).

DELIVERABLE the downstream chain should produce (NOT this idea-anchor): a RANKED candidate-drug table —
drug × in-hand-R-count(Oxford+Sci234) × AMRFinder-tractability × proposed rule × powering verdict
(SCORED-feasible / UNDERPOWERED / ABSTAIN) — plus the frozen-surface decision (below) made explicitly.

KEY ARCHITECTURE DECISION (surface it, don't drift): adding a drug touches the REPRODUCIBILITY-FROZEN
surface — DRUG_RULE in dna_decode/eval/amr_rules.py + DRUG_BREAKPOINTS/DRUG_AMRFINDER_CLASSES in
dna_decode/data/mic_tiers.py (frozen 2026-06-13). Two options: (a) extend the frozen files ADDITIVELY +
re-validate the 6 existing cells don't regress (the cef-fix discipline), or (b) build a SEPARATE
non-frozen drug-rule overlay/catalog the scorer reads alongside the frozen one (keeps the freeze intact).

CONSTRAINTS: solo hobby project; free/in-hand data only (no money, no new acquisition); north star is a
DECODER TOOL — a new drug cell only counts if it decodes that drug HONESTLY (a SCORED cell on an
uninformative drug like ampicillin, or an underpowered one, is worse than no cell). Failure-tolerant:
naming a drug ABSTAIN/UNDERPOWERED is a valid, expected outcome. Do NOT pad the cell count.
```

---

## Soraya's drafted anchor framing (ratify or redirect)

**Formal rephrase.** Use the just-shipped fam.tsv Subclass resolver + the two in-hand measured-MIC cohorts
to add the *tractable, adequately-powered* subset of new antibiotics as honestly-validated decoder cells,
abstaining explicitly on the intractable/underpowered ones, without regressing the 6 frozen cells.

**Fundamental clarifications (the ~3 the skill will likely ask — drafted answers):**
1. **Frozen-surface: extend additively vs separate overlay?** → *Draft: separate non-frozen overlay
   catalog.* Keeps the 2026-06-13 reproducibility freeze byte-intact (no re-validation of 6 cells needed
   per new drug), and the report card already unions a shipped-surface registry. Risk of the alternative
   (additive edit): every new drug forces a full no-regression re-run of the frozen cells. *Authority
   decision — yours.*
2. **Scope: clean tier only, or push into caveated (colistin)?** → *Draft: clean tier first
   (trimethoprim-sulfa → β-lactam siblings → amikacin), defer colistin/inhibitor-combos.* The clean tier
   is high-confidence; the caveated tier risks shipping half-blind cells that need their own honesty
   framing. *Authority — yours.*
3. **Powering bar for a SCORED new cell?** → *Draft: ≥10 R AND ≥10 S in the pooled in-hand MIC (mirrors
   `independent_cohort_validate` MIN_PER_CLASS); below that → UNDERPOWERED, reported not hidden* (the gent
   0-R outcome is the precedent). *Authority — yours.*

**Current assumptions (flagged for test):**
- The in-hand cohorts have enough R per new drug to validate. **Likely FALSE for several** — gent already
  hit 0 R; carbapenems are rare in E. coli. *The R-count census per drug is the FIRST action, before any
  rule is written* (threshold-vs-null-baseline discipline).
- fam.tsv family-level resolution is faithful for new β-lactams. *True for CTX-M/CMY (the drivers); the
  surfaced scope-limit on TEM/SHV/OXA-ESBL variants applies (cef sens 0.833 already showed it).*
- AMRFinder curates determinants for these drugs. *True for sul/dfr/bla/aac; verify the trimethoprim
  folate-pathway Class naming (TRIMETHOPRIM/SULFONAMIDE) in the deployed fam.tsv.*
- Adding a drug is "additive" and can't regress existing cells. *Must be VERIFIED — DRUG_AMRFINDER_CLASSES
  and the resolver are shared; a new class filter could perturb a shared-token match.*

**Blunt opinion.** This is the cheapest high-value move on the board: labels are already on disk, the
resolver just proved out, and several drugs are one-rule-each. **trimethoprim-sulfa is the obvious first
win** (no fam.tsv, no intrinsic trap, MIC in both cohorts, big acquired sul/dfr signal). BUT roughly half
the candidate list is a trap — **ampicillin is intrinsically blaEC-R in ~all E. coli** (the decoder would
just always say R: an uninformative cell), inhibitor combos aren't gene-count-decodable, colistin is
chromosomally half-blind. The dominant *real* risk is **powering**: gent showed these cohorts can have 0
R for a drug, so several "new drugs" may be born UNDERPOWERED. Do the R-count census FIRST; let it kill
the unpowered/uninformative candidates before a single rule is written. And make the frozen-surface
decision deliberately — retrofitting the reproducibility freeze for a marginal drug is a real cost.

**Recommended next step.** `/probe` — this is code-touching (frozen `amr_rules.py` / `mic_tiers.py`) and
needs the in-hand-MIC R-count census + the extend-vs-overlay grounding before any plan. The first
concrete probe action is the **per-drug R/S count over Oxford + Sci234** — it decides which drugs are even
worth a rule.

```
/probe new-drug-coverage-from-in-hand-mic
```
