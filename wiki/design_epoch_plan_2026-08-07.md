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
  (yeast **iMM904** MCC **0.377** on a label-matched rich medium — the only cross-organism cell that is
  SCORED at all; the 0.252 first reported here was a MEDIUM-MISMATCH artifact, see
  `wiki/fba_label_matched_medium_2026-08-11.md`).

  > **Stale-citation correction (2026-08-11).** This line previously also read *"P. aeruginosa recall
  > **0.57**, worst on the metabolic class."* **Withdrawn.** That number came from a run where
  > `paeruginosa` was mapped to `iJN1463` — a *P. putida* model — the wrong-organism defect fixed in
  > `78ea4ba`. The artifact was corrected to `MODEL_WALLED` in `0d8f9ad` (there is no P. aeruginosa GEM in
  > BiGG, and substituting another organism's model was refused), but this plan was not updated with it.
  > There is **no** P. aeruginosa essentiality number. S. aureus is `LABEL_WALLED` (iYS854 uses
  > `USA300HOU_####` ids; NTML is JE2 `SAUSA300_####` — a crosswalk away). So the honest statement of the
  > cross-organism weakness rests on **yeast alone**, which weakens the evidence for this bullet without
  > changing its direction.

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

## 5. The epoch's falsifier (pre-registered — SHARPENED 2026-08-07 with the paper's own baseline)

Original wording: "beat a simple sequence-feature baseline". **Kosuri 2013 supplies a better, published
one**, so the bar is now sourced rather than invented:

| model (Kosuri 2013, their Fig. 4 / ANOVA) | RNA R² | **protein R²** |
|---|---|---|
| simple multiplicative (promoter strength × RBS strength) | 0.92 | **0.76** |
| ANOVA (promoter + RBS, both affecting both levels) | 0.96 | **0.82** |

Also reported: 80% of RNA and **64% of protein** levels land within twofold of prediction, while **the
worst 5% deviate by 13-fold on average**.

**The bar: a learned sequence model must beat protein R² ≈ 0.82 on HELD-OUT constructs.**

Two honesty conditions on that comparison, both load-bearing:
1. The paper's R² values come from element strengths **fit on the same data**, so they are effectively
   in-sample. The fair test re-fits the element-strength baseline on the training split only and scores
   both models on the same held-out constructs.
2. Split by **element**, not at random — held-out *promoters* and *RBSs*, not just held-out combinations.
   Predicting an unseen combination of seen parts is a much easier problem than predicting a new part,
   and only the latter is what a designer actually needs.

The paper's own conclusion makes this a genuinely decisive test rather than a formality: in 2013 the
authors judged prediction good enough that they recommended *screening libraries instead of predicting*.
If a modern sequence model clears 0.82 on held-out elements, that is a real advance and Track B is live.
If it does not, that is an honest negative and the epoch falls back to deterministic Track A.

### Getting the data (verified 2026-08-07 — needs ONE browser download)

Automated fetch is blocked, and it is a **bot-block, not a paywall** (the article is Free access).
**DIAGNOSIS IS DEFINITIVE — do not re-try these; four independent routes are exhausted:**

| route | result |
|---|---|
| PNAS `suppl_file` + `downloadSupplement`, correct `.xls` names, proper referer | first **HTTP 403**; later **HTTP 200 returning 59 211 bytes of HTML** — identical for `sd01`/`sd03`, confirmed to contain `cloudflare` / `captcha` / `challenge` / `verify` markers. A **Cloudflare CAPTCHA**, not a transient error |
| PMC mirror (`PMC3752251`, real paths parsed from the article HTML) | **JavaScript proof-of-work** "Preparing to download…" interstitial |
| GitHub repository search (4 query forms) | **0 repos** |
| reuse-paper vendored copy (arXiv 2506.10271, which uses this dataset) | **no repo link** on the abstract page |

A browser passes the Cloudflare and PMC gates invisibly; nothing scriptable from this host does.

- Article: <https://www.pnas.org/doi/10.1073/pnas.1301301110> → Supporting Information
- **`sd03.xls` (15.8 MB) is the one that matters** — *"All values used in intermediate and final
  calculations are enumerated in Dataset S3"* (per-construct DNA/RNA/protein). `sd01`/`sd02` are the
  per-promoter and per-RBS strength summaries.
- PMC direct paths (same files, if the PNAS page is awkward):
  `https://pmc.ncbi.nlm.nih.gov/articles/instance/3752251/bin/1301301110_sd0{1,2,3}.xls`

Drop any of them anywhere on disk and point me at the path; ingestion is a short build from there.

Secondary honesty rail: a design that FBA says is growth-coupled has **not** been shown to work in a
cell. Every output of Track A is a *hypothesis for the bench*, never a validated strain — the same wall
that separates our decoders from clinical claims.

---

## 6. What needs the user (authority, not technical)

1. **Ratify the epoch + the ordering** (recommend A → B → C).
2. **Kosuri SI table** — one browser download unblocks Track B, or I search for a mirrored/cleaned copy.
3. **Compute** — Evo 2 is optional; if wanted, it is a Kaggle T4 job, not a local one.

Nothing here spends money. Tracks A and C are fully reversible and need nothing from outside.
