# The hybrid decoder: CALL / DOUBT / EVIDENCE

**Status:** DRAFT for ratification, 2026-08-31. Decomposition + plan only — nothing implemented.

Assembled from the measured record, not from recollection. Every number below is traceable to a committed
artifact; scope figures come from `scripts/project_status.py` at read time.

---

## Part 1 — Assessment

### What is actually deployed

`dna-decode` v0.13.1 on PyPI. **48 console entry points · 44 CLI traits · 110 registered cells** across 7
tracks (typing 34, viral 29, amr 25, pgx 14, finder 6, hla 1, mendelian 1). Evidence tiers: **28
`INDEPENDENT_MEASURED`**, 24 `NEAR_INDEPENDENT`, 33 `KNOWLEDGE_BASELINE`, 13 `FAITHFUL_TO_TOOL`, 11
`NO_FREE_SOURCE`, 1 `NOT_CENSUSED`.

### Five findings that constrain any design

**1. The deployed core is right, and it is not a model.** Curated deterministic rules beat learned scoring
wherever a catalog exists: **0.926–0.962** vs ESM2's **0.454 — below chance**. Independently confirmed in
the literature (rule-based ResFinder superior on divergent genomes; *mBio* 2020: *"the prediction model
matters far less than the dataset"*).

**2. The catalog's failure mode is COMPLETENESS, not accuracy — and it has now been found twice, the same
way.**

| gap | mechanism | cost | how it was found |
|---|---|---|---|
| gentamicin `rmt` | rule matches AMRFinder `Subclass == GENTAMICIN`; `rmt*` files under generic `AMINOGLYCOSIDE` | sens **0.523** vs 0.893; 24/31 FN carry `rmt` | an independent label set |
| HIV NNRTI | catalog scoped to the 8 Stanford major positions | **53** resistant isolates missed; drivers at 179/98/221/227/108 | an independent label set |

Neither is a model failure. Both were invisible until independent labels arrived. **This is the real
product risk and it is systematic.**

**3. The learned layer's honest role is SELF-AWARENESS, not prediction.** The position-novelty flag
recovers **60.4%** of the HIV blind spot (lift 4.69) at zero cost — and its own artifact insists it is
*"a 'catalog call may be incomplete' self-awareness flag, NOT a resistance predictor."* That framing is
the finding. Every learned *predictor* attempt in the natural-population regime is a closed negative
(0-for-5, de-confounded).

**4. Where learned methods do work is narrow and known.** Constructed variation, molecular endpoint,
**orthogonal modalities — not scale**. `ESM2+GEMME+ProSST` beats ESM2 on **90.5%** of proteins paired,
while 650M > 3B > 15B. Extending this to organism-level natural populations is the arm that keeps dying.

**5. Our differentiator is the evaluation discipline, and the literature says so by omission.** DART-Eval,
Kedzierska (*Genome Biology* 2025) and Ahlmann-Eltze (*Nature Methods* 2025) are, in effect, three papers
about the absence of exactly what this repo does by default.

### The design conclusion

**The hybrid is NOT "catalog + ML predictor".** That is the framing that has failed five times and that
today's framing sweep scored at **zero survivors**. The measured shape is different:

> A deterministic **CALL**, a cheap **DOUBT** layer that qualifies the call without competing with it,
> and an **EVIDENCE** layer that makes both auditable.

---

## Part 2 — The system

| layer | what it does | status | key property |
|---|---|---|---|
| **L1 CALL** | curated deterministic rules → R / S / INDETERMINATE | **shipped**, 110 cells | high precision, honest abstention |
| **L2 DOUBT** | *"this call may be incomplete, and here is why"* | **mostly missing** — one prototype, unwired | **never emits a competing call** |
| **L3 EVIDENCE** | de-confounding, nulls, denominators, leakage, provenance | **built**, under-exposed | what the field lacks |
| **L4 LEARNED (narrow)** | forward/inverse, orthogonal-modality hybrid | **shipped**, bounded | molecular + constructed only |

**The load-bearing rule for L2:** a doubt signal may qualify a call, request abstention, or explain
itself. It may **never** overrule L1 or emit a resistance prediction of its own. That constraint is what
keeps it out of the failed-predictor regime, and it is exactly how the position-novelty artifact already
frames itself.

**Why L2 is the innovation:** every paper in the field ships an L1 competitor or an L4 predictor. Nobody
ships *"my catalog might be wrong here."* It is cheap (the prototype needs no model, no GPU, no
structures), it targets the measured product risk, and it composes with the deterministic core rather
than fighting it.

---

## Part 3 — Decomposition

| family | scope | blocked_by | authority? |
|---|---|---|---|
| **F-A · doubt-layer** | generalize + surface the completeness signals | — | no |
| **F-B · catalog-curation** | close the measured gaps; define the procedure | F-A | **YES** — edits a shipped catalog |
| **F-C · evidence-surface** | expose L3 as a first-class product surface | — | no |
| **F-D · learned-narrow** | hold L4 to its regime; no extension | — | no (a *restraint*) |
| **F-E · acquisition** | the label wall (PEAR; condition-resolved expression) | — | **YES** — external/money |

### Requirements flow-down

- **F-A → F-B.** Curation must be *measured*, not asserted: without a doubt-layer baseline you cannot say
  what curation recovered. F-B is blocked on F-A by construction.
- **F-A** consumes L1's existing `AbstentionVocab` and the existing flag — **no new evidence-tier concept
  needed** (executed-verified today).
- **F-C** is independent; it exposes machinery that already exists.
- **F-D** is a restraint family: its deliverable is a boundary that stays enforced, not a build.
- **F-E** is external. Nothing downstream may be planned as if it will land.

**Critical path: F-A → F-B.** Everything else is parallel or held.

---

## Part 4 — Plan for F-A (the critical path)

**Checkable bar** — draft, for ratification:

1. `test-exit-0` — a doubt-layer module exists with tests, and its signals are pure functions over
   already-computed determinant calls (no model, no network, no structures).
2. `test-exit-0` — the HIV decoder record carries a `doubt` block; a guard test asserts the block can
   never contain a resistance call.
3. `file-exists` — an artifact reports doubt-layer sensitivity **per cell**, against the flag's measured
   0.604 on the EFV blind spot as the incumbent baseline.

**Ordered steps** (steps 1–4 carry no authority; step 5 is the fork):

1. **Generalize the completeness signal.** The `rmt` and HIV gaps are one shape: *a determinant family
   present in the data but unrepresentable by the rule*. Write the screen that detects it from cached
   determinant calls — for AMR, determinants whose class/subclass the rule cannot match; for target-site
   cells, the existing position-novelty logic. **This is the object F4-same-class-as-rmt showed does not
   exist.**
2. **Measure it per cell.** Run across the 10 SCORED AMR cells + the HIV cells. Report sensitivity and
   false-positive rate on labelled data, against the 0.604 incumbent. Report per-cell, never pooled.
3. **Wire it into the record**, additively, as a `doubt` block behind a hard guard test: the block may
   carry a reason and a tier, never a call. `hiv_amr.py` is **not** in the prospective-lock manifest
   (executed-verified), so this does not touch the frozen surface.
4. **Register the cells** so a doubt signal is visible on the trust surface without changing any cell's
   evidence tier — the same augment-only discipline as the lineage, prospective and source-concentration
   layers.
5. **Then, and only then, F-B** — with a measured baseline for what curation must beat.

**What this does NOT do:** it does not predict resistance, does not touch the frozen AMR surface, does not
require a GPU or a structure, and does not change a single existing call.

---

## Part 5 — What is deliberately NOT planned

- **A better scorer for the blind spot.** Today's framing sweep: incumbent framing, **0 survivors, 2
  executed kills** — including my own ΔΔG-first recommendation, disproved by a zero-tool flag that
  already recovers 60.4%.
- **Extending L4 to organism-level natural populations.** Closed negative, 0-for-5, de-confounded, and
  independently confirmed at 24,000 genomes (*PLOS Biology* 2025, where **more data does not rescue it**).
- **Pretraining anything.** Scale is measured dead in this regime; the field is crowded; `gLM2-650M` is a
  download.
- **Unpinning the frozen AMR surface.** The gentamicin `rmt` fix invalidates the prospective lock and the
  reproducibility freeze. It needs its own validation and a **new** lock — a user authority call, not a
  step.

## Part 6 — Open authority calls (yours, not mine)

1. **F-B — edit a shipped catalog?** The HIV drivers are named and countable and `hiv_amr.py` is outside
   the lock, so this is *cheap*; whether curated biology enters a shipped catalog is still a scope call.
2. **The gentamicin v2 lock** — measured at **+0.369 sensitivity** at zero measured specificity cost, but
   deploying it invalidates the lock and the freeze.
3. **F-E acquisition** — PEAR (~23,000 *E. coli* strains, single-copy `blaCTX-M-14` variants, measured
   growth) is the highest-value target found; constructed variation at scale on an AMR target.
4. **Whether L2 becomes a headline product claim** or stays an internal diagnostic. It changes what the
   tool *is*.

## Note on process

`decompose` would normally self-invoke `/project-init` per family. I did **not** — creating five ledgers
for an unratified decomposition is the sprawl the population cap exists to prevent, and the standing
Planning-STOP rule says wait for direction. Say the word and the ledgers get seeded.
