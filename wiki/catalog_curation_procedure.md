# What is allowed to enter a shipped catalog

Written 2026-09-01 from **two executed instances**, one shipped and one declined — not from principle.
The generalization was deliberately deferred until both existed, for the same reason the rejection-gate
screen was: a procedure written before any real case encodes opinion.

| instance | outcome | what settled it |
|---|---|---|
| Gentamicin `rmt*` / `npmA` rescue | **SHIPPED** under a v2 lock | +0.369 sensitivity on 131 disjoint isolates, at no measured specificity cost |
| HIV NNRTI catalog additions | **DECLINED** | every variant recovered less of the blind spot than the free doubt layer already does |

---

## The four conditions

An entry may enter a shipped catalog when **all four** hold. They are ordered by how often they are the
one that fails.

### 1. Every entry is sourced to a named external authority

Not to an OLS coefficient, not to recall. The authority is named in the module, and the entry is
traceable to it. **This is the fabrication gate**, and it is the reason the 40 unrecorded colour loci
remain uncurated: writing curated biological facts into a shipped catalog from anything other than a
source is fabrication with a citation-shaped hole.

*Enforced, weakly:* `tests/test_catalog_provenance.py` fails a catalog module that names no authority.

### 2. The recovery is MEASURED against the doubt-layer baseline, not asserted

A curation competes with the L2 doubt layer, which is free, needs no model, and touches no shipped
surface. If a proposed addition recovers less of the measured gap than the incumbent flag, **the flag
wins and the catalog is left alone**.

This is why the curation family was declared `blocked_by` the doubt-layer family: the baseline had to
exist before curation could be judged. It is also what decided the NNRTI case — every variant tested
recovered 0.000–0.500 of the blind spot against the position-novelty flag's 0.604.

### 3. The lock and freeze consequence is stated explicitly, and authorized

Two different consequences, never conflated:

- `hiv_amr.py` is **not** in the prospective-lock manifest → curating it is lock-safe.
- `amr_rules.py` **is** → changing it invalidates the prospective lock *and* the reproducibility freeze,
  and requires an unfrozen revision, its own validation, and a **new** lock.

The second is a user-authority call. Conflating them would smuggle a lock-invalidating change past a
lock-safe one.

### 4. A derivation is reviewed as biology, not accepted as a metric

Two failures from the NNRTI derivation, both of which a good balanced-accuracy would have hidden:

- **Read the derived entry list.** The first run admitted `L234L`, `K238K`, `M230M`, `R72R` — WT == MUT,
  not substitutions at all but mixture markers correlated with treatment experience. Excluding them
  changed the result *against* my favour, so the headline would have been wrong in the flattering
  direction.
- **An automated rule that deletes a canonical entry is measuring the wrong thing.** The best-scoring
  variant dropped `Y181C/A/I` from the EFV catalog, because Y181C co-occurs with K103N which absorbs the
  coefficient. Co-occurrence absorption is a property of the fit, not of the biology.

---

## What a test can and cannot enforce

Measured, not assumed. **Per-entry citation is not representable.** The shipped catalogs are bare
collections — `NNRTI_RT_MAJOR_DRMS: set[str] = {"L100I", ...}` — with no field to hang a source on.
Adding one means restructuring the frozen surface: large, risky, and buying nothing that was measured.

So the split is:

| condition | status |
|---|---|
| 1. named authority — **per module** | **enforced by test** |
| 1. named source — **per entry** | **review discipline** (not representable without restructuring) |
| 2. measured against the doubt layer | review discipline (the measurement is per-case) |
| 3. lock consequence stated | review discipline; the lock itself is tamper-evident by hash |
| 4. derivation reviewed as biology | **irreducibly** review discipline |

**Three of four conditions are review discipline.** Saying so is the point: a green suite that checked
one weak condition would read as more assurance than it is. The enforced test catches an *unsourced
catalog*; it cannot catch a fabricated entry inside a correctly-cited one.

## The standing default

**A catalog is not edited because an edit is available.** Both instances began as "close the measured
gap"; one earned it and one did not, and the difference was measurement against a free incumbent — not
whether the fix was technically possible. It always was.
