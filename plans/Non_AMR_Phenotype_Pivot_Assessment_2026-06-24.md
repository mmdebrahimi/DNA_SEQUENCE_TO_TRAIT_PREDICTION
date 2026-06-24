# Non-AMR phenotype pivot — state assessment + plan (2026-06-24)

> User question (verbatim intent): "make a plan for other phenotypes we can link to gene/DNA -- flagella length,
> enzymes they produce, their look, length, what they eat, other qualities. We've mostly concentrated on drug
> resistance. Are we in a good place to pivot? What do you think?"

This memo answers it grounded in the project's own hard-won prior art (the embedding 0-for-4 negative + the
EP-4 scoping memo + the "labels not models" wall). **TL;DR: you've already partly pivoted; the rest splits
sharply into a GREEN path (do it) and a RED path (a closed negative — don't).**

---

## 1. State assessment — where we actually are

**The deterministic decoder is NOT an AMR tool. It's a general `genotype -> catalogued-trait` engine** that
happens to have AMR as its deepest cell. It already decodes **8 traits**, most of them NON-AMR:

| Already shipped | Trait class | This IS one of the user's examples |
|---|---|---|
| `dna-pathotype` | virulence pathotype (EPEC/EHEC/ETEC/UPEC/EAEC) | "other qualities" — the EP-4 #1 pick, done |
| `dna-serotype` | O:H antigen (incl. **`fliC` = flagellar antigen type**) | "their look" / flagellar identity |
| `dna-mlst` | sequence type (clonal identity) | "other qualities" |
| `dna-plasmid` | plasmid replicon type | mobile-element capability |
| `dna-resfinder` / `dna-pointfinder` / `dna-disinfinder` | acquired genes / point mutations / biocide tolerance | "enzymes they produce" (β-lactamases etc.) |
| `dna-amr` | antibiotic R/S (bacteria + TB + fungi + 3 viruses) | the deep cell |

Plus the validation infra (provenance-disjoint scoring, trust-surface honesty tiers, report cards) is built
and reusable. **So the "pivot" is smaller than it feels — the engine is already trait-general.**

**The one thing that was tried as a LEARNED predictor and FAILED is the relevant warning:** the
foundation-model-embedding bet is a **closed 0-for-4 negative** (cipro within-lineage, pathotype-via-embedding,
Arabidopsis flowering-time x2) -- in every de-confounded test the embedding learned **population structure, not
the causal trait signal**. See `wiki/negative_results_map_2026-06-13.md` +
`wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md`.

---

## 2. The decisive split — every trait the user named falls into one of two classes

The project's own gate (the embedding-niche three-part test + determinant-catalog tractability) sorts ANY
proposed phenotype cleanly:

### GREEN — determinant / catalog-tractable (the deterministic engine transfers directly)
A KNOWN gene/allele drives the trait; you detect presence/identity, not a quantitative emergent value.

| User's example | GREEN form | Why it works |
|---|---|---|
| "enzymes they produce" | a SPECIFIC enzyme's GENE present? (β-lactamase, urease, catalase, specific metabolic enzyme) | gene-presence = the AMR/resfinder pattern exactly |
| "what they eat" | a SPECIFIC catabolic capability (lac operon, specific sugar-utilization operon, auxotrophy) | single operon presence -> deterministic call |
| "their look" (antigenic/structural identity) | O:H **serotype** (DONE), K capsule TYPE, fimbriae/flagella TYPE, pili presence | curated antigen/structure catalogs exist (serotype already does this) |
| "other qualities" | pathotype (DONE), biotype, toxin profile, virulence-factor set | curated marker catalogs (VirulenceFinder etc.) |

GREEN traits are FREE, on-brand, low-risk, and reuse the trust-surface + honesty framework. **This is where a
real pivot lives.**

### RED — emergent / quantitative / polygenic (the closed-negative territory — DO NOT reopen)
Shaped by genome-wide architecture, NOT a catalog; usually NO free isolate-level MEASURED label; and exactly
where the learned-embedding bet died 0-for-4.

| User's example | Why it's RED |
|---|---|
| **flagella LENGTH** | a quantitative morphological measurement; no curated "length gene"; no public per-isolate flagella-length DB; polygenic + growth-condition-dependent |
| cell **size / length** | same -- emergent morphology, no catalog, no free measured label |
| "their **look**" (morphology broadly) | emergent; would need imaging + a learned model = the multimodal track, which is parked and label-blocked |
| **growth rate** / full metabolism ("what they eat" broadly) | EP-4 Candidate 1 already verdict-ed SKIP: distributed mechanism, KEIO-knockout-shaped not natural-variation, no buyer, no free measured label |

RED traits fail the SAME gate everything else fails: **no free, independent, MEASURED, isolate-level label** +
no determinant catalog. Trying them re-hits the label wall AND the embedding-learns-structure-not-signal
problem. The project has a RECORDED negative here -- re-litigating is the decision-avoidance pattern.

---

## 3. The binding constraint is unchanged: LABELS, not models

For ANY new phenotype, the question is NOT "can DNA encode it" (often yes) -- it's:

1. Is there a curated **determinant catalog** (gene/allele -> trait)?  -> if yes, GREEN deterministic cell.
2. Is there a **FREE, independent, MEASURED, isolate-level label** to VALIDATE the cell?  -> if no, you can
   BUILD the caller but you can't honestly SCORE it (the fungal cells are exactly this: shipped, but
   `NO_FREE_PHENOTYPE_SOURCE`).
3. Is it already shipped?

Most of the user's vivid examples (flagella length, look, size) fail #1 AND #2. The tractable ones
(enzymes-as-genes, antigen types, catabolic operons) pass #1; whether they pass #2 decides if the cell gets a
real validation number or ships as a curated-but-unvalidated caller.

---

## 4. The plan — a phenotype triage + a GREEN slate (if you want a new cell)

### 4a. Triage gate (apply to ANY proposed trait before building)
```
trait -> [determinant catalog exists?] --no--> RED (emergent/learned -> closed negative; PARK unless a free
                                       |          measured label appears)
                                       yes
                                       v
         [free independent measured isolate-level label?] --no--> GREEN-caller (ship, tier=NO_FREE_SOURCE,
                                       |                            like the fungal cells)
                                       yes
                                       v
         GREEN-validated cell (ship + score + trust badge)  <-- the gold standard (AMR/pathotype shape)
```

### 4b. GREEN slate worth considering (catalog-tractable, on-brand) -- ranked by validate-ability
1. **Capsule / K-antigen typing** (Klebsiella K-locus; E. coli) -- curated DB (Kaptive-style); the natural
   sibling of the shipped O:H serotype; high clinical relevance. *Validation label availability = the gating
   question.*
2. **Expanded virulence / toxin profile** beyond the current pathotype resolver (more VirulenceFinder
   clusters surfaced as discrete trait calls).
3. **Specific catabolic capability** (a named, clinically-used biochemical, e.g. a specific sugar fermentation
   or a urease/oxidase-type marker) where a gene catalog + a labeled cohort coexist.

NONE of these is urgent; all are smaller than they sound (the engine + trust-surface already exist).

### 4c. What NOT to do
- Do NOT build a learned predictor for flagella length / cell size / morphology / growth rate from DNA -- that
  is the closed 0-for-4 negative; it will re-hit the label wall + learn structure not signal.
- Do NOT reopen the multimodal "DNA -> image / look" track -- parked + label-blocked; the north star
  explicitly is NOT a stepping stone to "DNA -> animal image".

---

## 5. My honest opinion (the user asked "what do you think?")

- **Are we in a good place to pivot?** For GREEN (catalog-tractable) traits: **YES, excellent** -- the engine,
  trust-surface, and honesty framework all transfer, and you've already shipped ~7 non-AMR cells. For RED
  (emergent/quantitative -- flagella length, look, size, growth): **NO, and it's a recorded closed negative**;
  pivoting there re-fights a battle the project already lost on principle (labels + structure-vs-signal).
- **The deeper read:** the most VALUABLE next move is probably NOT a new phenotype cell at all. The validation
  surface is rich, the tool is honest, and the standing open decision is the **productization fork** (packaging
  gate vs editable-only) + possible banking. A new GREEN cell is a fine, low-risk increment IF a specific trait
  has a FREE validation label -- but adding cells without a validation label just grows the
  `NO_FREE_PHENOTYPE_SOURCE` count, which the report card already shows is the project's saturation point.
- **Bottom line:** the honest "pivot" is to recognize the decoder is ALREADY trait-general, add a GREEN cell
  only when a free-labeled tractable trait appears, and NOT chase emergent morphology/metabolism from DNA.

---

## 6. Recommended next step (STOP for user direction)
Pick one:
- **(i) Build a GREEN cell** -- name a catalog-tractable trait + I'll first run the LABEL gate (is there a free
  measured isolate-level validation set?) before any build. Best first candidate: **Klebsiella K-antigen
  capsule typing** (serotype's sibling).
- **(ii) Bank the phenotype question** -- the decoder is already trait-general; resolve the productization fork
  / new free-label cell instead.
- **(iii) Deeper scan** -- I run `/research` for free isolate-level labeled datasets for a specific GREEN trait
  you care about, and report feasibility before committing.

NO code / no `/idea-anchor` fired here -- this is the scoping memo (mirrors the EP-4 discipline).
