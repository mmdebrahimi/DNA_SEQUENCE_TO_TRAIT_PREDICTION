# The GLM: first-principles design

**Date:** 2026-10-06
**Status:** DESIGN — written before the SOTA scan, deliberately, so prior art checks the design rather
than anchoring it.
**Mandate (user, 2026-10-06, emphatic):** *build the GLM.* Do not benchmark it against the deterministic
catalog — that comparison is a category error and is now recorded as a standing lesson
(`memory/feedback_glm_is_the_product_not_the_catalog.md`). The catalog was a BRIDGE built so there would
be any number at all; the GLM is the product.

---

## 0. What the thing has to do

From the north star, verbatim: **name a trait, get the genome edits.** Four steps:

1. identify the organism — SHIPPED (`dna-identify`)
2. read the genotype → which sections JOINTLY make a phenotype
3. **act in bacteria** — *"make a minor modification on a bacteria to produce a protein we want it to produce"*
4. climb to animals

Step 3 is the load-bearing one and it is stated as a **generative** request. It has a precise form:

> **input:** a host genome + a desired product (a protein, an enzyme activity, a metabolite)
> **output:** an actual DNA edit — sequence, and where to put it — with a calibrated expectation of what
> it will do and an honest interval on that.

Nothing in the repo does this today. What the repo has is a **scoring stack** and two **narrow symbolic
proposers**. The gap is specific and nameable, which is the point of this document.

---

## 1. What already exists (derived from the code, 2026-10-06, not from prose)

| Capability | Module | Status | Validation it carries |
|---|---|---|---|
| protein variant → effect | `forward/variant_effect.py` | SHIPPED | DMS, Spearman **0.35–0.76** across two β-lactamases; modality hybrid wins 84–90% of ProteinGym proteins |
| genome coordinate → effect | `forward/genome_edit.py` | SHIPPED | 0.7611 over 1,715 variants |
| **effect → candidate edits** | `forward/inverse.py` | SHIPPED | falsifier PASS (blaTEM +53.0%); **ranks, does not dose** |
| **score → calibrated magnitude** | `forward/dosage.py` | BUILT | split-conformal + `interval_narrowing` honesty check |
| gene KO → organism growth | `fba/model.py` | SHIPPED | Keio accuracy **0.954**; per-condition MCC **0.70–0.74** |
| **product → KO design set** | `fba/design.py` | SHIPPED | reproduces the literature anaerobic succinate design `PFL+LDH_D+ALCD2x` @ 9.26 mmol/gDW/h |
| **missing reaction → repair** | `fba/gapfill.py` | SHIPPED | measured: sucrose growth **0.000 → 1.7798 /h** by adding one donor reaction (`FFSD`) |
| missense → cell trait | `fba/compose.py` | SHIPPED | chains the two above, honest about the heuristic join |
| hard biophysical law | `constraints/` | SHIPPED | frame / codon reachability / mass balance / clade genetic code; **71.6%** of *M. genitalium* CDS are false-nonsense under the wrong table |
| genome annotation + tiering | `genome_map/` | SHIPPED | 5 evidence tiers, honest unknown rate |

**Read that table again as an architecture, not an inventory.** It is a scoring committee with two
validated oracles at *different levels of biology* (molecular and organism), plus a hard-constraint layer.
That is exactly the substrate a generative design loop needs — and it is the reason the label wall does
not block this work (§3).

---

## 2. The four missing pieces

### M1 — There is no generative model over DNA *sequence*
Every proposer above operates on a **symbolic** edit space: a mutation string (`S83L`), a reaction id
(`PFL`), a gene id. **Nothing emits nucleotides.** `models/foundation.py` wraps NT / DNABERT-2 / GENA-LM,
which are **encoder-only** — they score, they do not generate — and Evo was recorded as not callable on
this host's GPU. So the literal GLM does not exist yet. This is M1 and it is the centre of the build.

### M2 — There is no UP-regulation edit class
`design.py` can only knock **out**. `inverse.py` can only substitute **within** a CDS. To make a host
*produce more of* something you need the expression levers — promoter, RBS/translation initiation, codon
usage, copy number, locus context — and **no proposer in the repo can express an edit of that kind.**
This is the gap that sits directly on top of the user's step 3.

### M3 — There is no bridge from a model-level edit to a buildable construct
`gapfill` says *"add reaction FFSD"*. It does not say *"here is the cassette, this is its sequence, insert
it here."* That translation — **reaction → gene → codon-optimized CDS + promoter + RBS + integration
site** — is the single highest-value missing link, because it is what converts the repo's existing
organism-level design capability into an actual answer to "give me the genome edits".

### M4 — Trait → target resolution for non-metabolic traits
"Produce more succinate" resolves to an exchange reaction (`resolve_target` does this). "Pink hair"
does not resolve to anything the current stack can score. This is step 4 and it is deliberately **last** —
not because it is uninteresting but because M1–M3 are its prerequisites.

---

## 3. Why this is buildable now — the label wall does not bound it

The repo's negative-results corpus is real, well-measured, and **scoped to one thing: discriminative
zero-shot scoring of variants in natural populations.** Specifically it says — correctly — that you cannot
learn genotype→phenotype supervised from natural cohorts without ancestry confounding it, and that
zero-shot likelihood does not rank resistance variants.

**Generation scored by a validated simulator is a different problem.** The generator needs an *oracle*,
not a *label set*. This stack already has two oracles with real validation attached (`forward`: DMS;
`fba`: Keio 0.954 / MCC 0.70–0.74), plus a hard-constraint layer that is correct **by construction**
(mass balance and reading frame are not statistical claims). So:

> **The wall blocks supervised g→p on natural populations. It does not block design-and-score.**
> Treating one wall as if it bounded both is what stalled this for months.

The precedent class is protein design: RFdiffusion/ProteinMPNN do not learn from phenotype labels, they
generate and let a structure oracle filter. The analogue here is stronger, because an organism-level
oracle (FBA) exists and is validated per-condition.

### The named risk: oracle hacking
A generator optimizing against FBA **will** find stoichiometric loopholes; one optimizing against ESM2
will find high-likelihood nonsense. This is the central failure mode and the design answers it structurally:

1. **Hard constraints run FIRST and are non-negotiable** — frame, codon reachability from the host's own
   genetic code, mass balance. These cannot be gamed because they are not learned.
2. **Committee, not a single scalar.** An edit must satisfy the molecular oracle *and* the organism oracle
   *and* the constraints. Loopholes rarely survive two independent levels of biology.
3. **Trust region.** Cap edit distance from wild type. Far-from-distribution proposals are where learned
   oracles are least reliable, and `design.py` already encodes the matching instinct (fewer knockouts wins
   ties, because each is real bench work).
4. **A held-out oracle that never scores during search** — used only to audit the selected set.
5. **Report the committee's disagreement** rather than averaging it away. Disagreement is the honest signal
   that an edit is in oracle-hacking territory.

---

## 4. The architecture

```
                        TRAIT SPEC
         (a molecule · a protein · an enzyme activity · a host)
                             │
   ┌─────────────────────────┴──────────────────────────────┐
   │ L0  TRAIT → TARGET RESOLVER                            │
   │  molecule → exchange reaction      [HAVE: resolve_target]│
   │  protein  → gene / CDS / pathway   [PARTIAL]            │
   │  trait    → ???                    [M4, deferred]       │
   └─────────────────────────┬──────────────────────────────┘
                             │
   ┌─────────────────────────┴──────────────────────────────┐
   │ L1  PROPOSER — the generative core, FOUR edit classes   │
   │  (a) coding substitution      HAVE  inverse.py          │
   │  (b) knockout / down-tune     HAVE  design.py           │
   │  (c) EXPRESSION up-tune       M2    ← step 3 lives here │
   │  (d) heterologous insertion   M3    gapfill → construct │
   │        ── all four emit DNA, not symbols ──  M1          │
   └─────────────────────────┬──────────────────────────────┘
                             │  candidate edits (sequence + locus)
   ┌─────────────────────────┴──────────────────────────────┐
   │ L2  ORACLE COMMITTEE                                    │
   │  HARD    constraints/   frame · codon reach · balance   │
   │  MOLEC   forward/       effect ρ 0.35–0.76  + dosage    │
   │  ORGAN   fba/           growth · flux · coupling        │
   │  EXPR    [M2 oracle]    predicted expression level      │
   │  → committee verdict + DISAGREEMENT, never a mean       │
   └─────────────────────────┬──────────────────────────────┘
                             │
   ┌─────────────────────────┴──────────────────────────────┐
   │ L3  SELECT → k proposals + assay plan + honest interval │
   │  propose-k / assay-k / keep-best (inverse's measured    │
   │  workflow: top-1 is ~4× worse than best-of-5)           │
   └────────────────────────────────────────────────────────┘
```

**Why this shape and not an end-to-end model.** An end-to-end trait→edit network is the thing the label
wall genuinely blocks: it would need (edit, measured phenotype) pairs at scale, which is exactly what does
not exist for free. The committee architecture needs **no such pairs** — it needs a generator (self-
supervised on sequence, which is free) and oracles (already validated). That is the whole argument for this
decomposition, and it is why it can start today rather than after an acquisition.

---

## 5. The first target, stated as a falsifiable build

Of the four edit classes, **(c) expression up-tuning** is where step 3 lives and where nothing exists.
So v0 of the GLM is:

> **Given a host genome and a target gene already in it, emit a ranked set of DNA edits predicted to
> RAISE that gene's product level, each with a calibrated expectation and an interval — then verify the
> predictions against a held-out measured expression set.**

That is a complete, checkable, label-honest unit of work. It needs: a nucleotide generator (M1), an
expression oracle (M2), and a real measured expression dataset for the held-out check. Whether such a
dataset is free and downloadable is the single most important unknown and is being scanned now.

**Explicitly NOT in v0:** the catalog comparison (category error), the climb to animals (M4), and any
wet-lab claim. A proposal is a hypothesis for the bench — the same honesty `design.py` already stamps
into every record.

---

## 6. Open questions this design cannot settle by itself

1. **Which nucleotide generator?** Needs: downloadable weights, a licence that permits this use, and
   inference inside 4 GB VRAM locally or 16 GB on a free T4. Evo2's smaller variants are the obvious
   candidate; unverified until the scan lands.
2. **Is there a free, measured, large bacterial expression dataset** (MPRA-style: regulatory sequence →
   measured expression)? If yes, M2 has a real supervised target and this stops being oracle-only. If no,
   M2 falls back to a biophysical oracle (RBS-Calculator-class) and the honest claim weakens accordingly.
3. **Does anyone already compose a sequence generator with a genome-scale metabolic oracle?** If yes,
   adopt their guardrails. If no, that composition is this project's novel contribution and should be said
   so plainly — once, with evidence, not as a boast.
4. **Trust-region radius** — an empirical question, answered by measurement, not assertion.
