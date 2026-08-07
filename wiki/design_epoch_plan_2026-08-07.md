# The design epoch — from decoding organisms to designing them (2026-08-07, DRAFT for ratification)

User direction: *"the AI/learned route is not a blocked path... a model that predicts missing parts of a
genome, then genotype→phenotype, would let us design new varieties of organisms — E. coli that produce
new proteins/enzymes."*

**The user is right, and the correction is sharper than "AI isn't blocked."** This plan states the
correction with evidence from this repo, revises the proposed critical path, and names the cheap
falsifier that could sink it.

---

## 1. The correction: I over-generalized a scoped negative (3rd occurrence)

Every learned-decoding **failure** in this project was on **natural populations**, where genotype tracks
ancestry:

| failure | substrate | root cause |
|---|---|---|
| Arabidopsis flowering (G2) | 1,003 wild accessions | learned population structure, not the trait |
| cipro within-lineage | clinical isolates | chance within lineage; learned the clade |
| pathotype | clinical isolates | sampling-defined label |
| organism-multimodal | wild/clinical panels | ties the linear cis-eQTL ceiling |
| ESM zero-shot vs AMR catalog | clinical | antagonistic phenotype; ESM **below chance** (0.454) |

Every learned **success** was on **constructed variation**:

| success | substrate | result |
|---|---|---|
| **Bloom-2013 yeast segregant cross** | 1,008 segregants of ONE cross | **12/12 traits decode, r 0.46–0.80 vs null p95 0.03–0.12** |
| **forward cell (DMS)** | designed mutant libraries | ESM2+ProSST+GEMME beats best-single on **84–90%** of proteins |
| inverse cell | constructed candidate edits | blaTEM +53.0% on measured wet-lab DMS |

The dividing line is not "AI vs rules". It is **natural variation (confounded) vs constructed variation
(confound-free by construction)**.

> ### The load-bearing consequence
> **Design IS constructed variation, by definition.** You build the variants, so there is no ancestry for
> a model to secretly track. **The design goal therefore lives entirely inside the regime where this
> project's learned models already win.** The "learned route is closed" verdict was drawn from the
> natural-population regime and **does not transfer to design.**

This is the third time I compressed a scoped negative into "AI failed" — flagged by my own standing note
`feedback_dont_overgeneralize_scoped_negative_zeroshot_vs_supervised`. Recorded so it stops recurring.

---

## 2. Revising the proposed critical path

Proposed chain: **genome infilling → g→p map → organism design.**

Honest revision: **infilling is not on the critical path.** It is a *generator*, useful but not
foundational, and putting it first creates a long serial dependency before any design capability exists.
The real path is three questions in dependency order:

| # | question the designer must answer | status in this repo | gap |
|---|---|---|---|
| **Q1** | Will the **part** work? (protein → stability / activity) | **SHIPPED** — `forward` + `inverse` cells, DMS-validated | none for variant-level design |
| **Q2** | Will the **host express** it? (promoter/RBS/codon → protein level) | **MISSING** | the real gap; free well-powered data exists |
| **Q3** | Can the host be **rewired to yield more**? (edits → flux to product) | **PARTIAL** — `fba` does edit→growth | missing the *design direction* (which edits maximize product) |

Q3's missing half is **growth-coupled strain design** (production envelope / OptKnock-style): given a
target product, find the gene edits that make producing it *necessary* for growth. This is deterministic,
native to cobrapy, needs **no labels**, and is the single largest step from "we decode organisms" to
"we design them."

### Where "predict missing parts of the genome" genuinely belongs

Not nucleotide infilling — the highest-value reading is **filling the missing FUNCTIONS**, and this repo
has already *measured* that gap:

- `genome_map` dark matter: **6.1% / 8.6% / 67.9%** of genes unknown under Bakta db-light (E. coli ST131 /
  K. pneumoniae / Gemmata).
- The FBA cell's weakness is **model gaps**: the sucrose false-negative was a missing transporter, not a
  wrong prediction (`fba_carbon_growth_validation_2026-08-03`); essentiality is weak off E. coli
  (yeast MCC **0.252**; P. aeruginosa recall **0.57**, worst on the metabolic class).

So "missing parts" → **metabolic model gap-filling driven by protein-function prediction**, which has an
*already-instrumented* metric: does essentiality/growth accuracy improve against Keio / SGD / GOLD_115?

---

## 3. Data reality — verified today, not assumed

The label wall that closed the AMR epoch **does not hold here**, because designed variation is measured in
the lab by construction. Verified from this host:

| source | what it gives | verified status |
|---|---|---|
| **Kosuri 2013** (PNAS 110:14024) | **12,563 constructed promoter×RBS combos**, DNA + RNA + protein measured | landing page OK; **SI download 403 (bot-block, not paywall)** → browser download or GEO/SRA or a cleaned CSV in a reuse repo |
| **JBEI ART** `Limonene_data_for_ART.csv` | engineered E. coli, 9 pathway-enzyme proteomics + measured limonene titer | **VERIFIED downloadable, no login — but N = 30 strains** |
| **Evo 2** (Nature 2026, `arcinstitute/evo2`) | 9T bp gLM, 1M ctx, designs bacterial-genome-scale sequence | repo reachable; **7B/40B — GTX 860M (4 GB) cannot run it**; Kaggle T4 16 GB fits ~1B, borderline 7B |

**Powering is split, and it dictates the allocation:** Kosuri is well-powered (12,563) → a learned model
is appropriate. Strain→titer public data is **small-N (30)** → a learned model there would be
underpowered; use it as a **validation set** for deterministic design, not as training data.

---

## 4. Recommended sequencing

**Track A — FBA growth-coupled strain design (deterministic, no label wall, buildable now).**
Add production-envelope + growth-coupled design to the `fba` cell: target metabolite → candidate gene
edits → predicted product flux vs growth. Validate against the N=30 limonene strains (real measured
titers) and against published OptKnock designs. *Highest VOI: biggest capability jump, zero label
dependency, entirely in the open mechanistic regime.*

**Track B — sequence→expression predictor on Kosuri (the learned track that can actually work).**
12,563 confound-free constructs, measured protein levels. **Falsifier:** held-out r² must beat a simple
sequence-feature baseline (GC / codon-adaptation / RBS free-energy). Closes Q2 and is the decisive test
of the learned-design thesis. Blocked only on getting the SI table.

**Track C — model gap-filling as the honest "missing parts."**
pLM function prediction over dark-matter genes → candidate missing reactions → measure the FBA accuracy
delta on the existing gold standards. Metric already wired.

**Evo 2** enters as an optional *generator* of candidate regulatory/coding sequences for Track B to rank —
compute-gated to Kaggle, and explicitly **not a prerequisite**.

---

## 5. The epoch's falsifier (pre-registered)

**Track B's held-out r² vs the simple-feature baseline.** If a learned model cannot beat GC/codon/RBS
features on 12,563 confound-free constructs with directly measured protein output, then the
learned-design thesis is in serious trouble and the epoch should fall back to deterministic Track A only.
Cheap, decisive, and it tests the thesis rather than assuming it.

Secondary honesty rail: a design that FBA says is growth-coupled has **not** been shown to work in a
cell. Every output of Track A is a *hypothesis for the bench*, never a validated strain — the same wall
that separates our decoders from clinical claims.

---

## 6. What needs the user (authority, not technical)

1. **Ratify the epoch + the ordering** (recommend A → B → C).
2. **Kosuri SI table** — one browser download unblocks Track B, or I search for a mirrored/cleaned copy.
3. **Compute** — Evo 2 is optional; if wanted, it is a Kaggle T4 job, not a local one.

Nothing here spends money. Tracks A and C are fully reversible and need nothing from outside.
