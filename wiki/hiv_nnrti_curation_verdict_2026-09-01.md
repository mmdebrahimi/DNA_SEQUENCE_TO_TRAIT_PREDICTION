# Curating the NNRTI catalog does NOT earn its place — measured, three ways

**Verdict: do not curate `hiv_amr.py`.** The blind spot is real (53 resistant isolates carry no
catalogued DRM), the fix is technically available, and it **loses** to the free zero-tool doubt layer
already shipped. The position-novelty flag stays the right instrument.

This is F-B, and it is the payoff of F-B having been **blocked_by F-A**: curation had to be *measured
against a baseline*, not asserted. The baseline won.

---

## What was tried

The other three HIV classes each shipped a deconfounded mutant catalog (`hiv_{nrti,pi,insti}_mutant_catalog.py`).
NNRTI was the one that never got one, so the natural move was the same method: multivariate OLS on
log10 fold-change, keep mutants whose **independent** coefficient clears a threshold after controlling
for co-occurrence, ≥5 carriers, 5-fold cross-validated held-out.

Two adaptations were forced, and both matter:

**Candidates span the whole RT.** The shared helper restricts to a class's catalogued positions. Reusing
it would let the fit re-derive only inside the 8 positions the catalog already has — and could never
reach the blind spot, whose drivers sit elsewhere. Candidates therefore come from the dataset's own
`CompMutList`, which carries full `<WT><pos><MUT>` strings across the whole RT and needs no
consensus-WT inference of ours.

**NNRTI is the mirror image of the other three.** Their v0 was *position-based* and over-called, so
deconfounding won by lifting **specificity**. NNRTI's v0 is already mutant-level; its failure is **false
negatives**. The question was whether sensitivity could be bought without paying for it.

**And NNRTI is scored ABSOLUTELY, not delta-honestly** — the one place it is stronger than the PI/INSTI
arc. Real Stanford `DRMcv.R` clinical cutoffs exist for EFV/NVP/ETR/RPV. Doravirine postdates that
script → reported `CUTOFF_UNAVAILABLE`, never guessed.

## Three measurements, one answer

**1. Inherited threshold (1.5×) — catastrophic over-call.**

| drug | v0 sens/spec | v0.1 sens/spec |
|---|---|---|
| efavirenz | 0.947 / **0.904** | 0.989 / **0.579** |
| nevirapine | 0.906 / **0.991** | 0.985 / **0.454** |
| rilpivirine | 0.822 / 0.544 | 0.975 / **0.150** |

75 additions. Balanced accuracy falls hard. The 1.5× threshold was tuned for *position-restricted*
candidates; spanning the whole RT multiplies the comparisons and it stops holding.

**2. Swept threshold (3×) — a hold, bought at an unacceptable price.**

5 additions, EFV balacc 0.932 vs v0 0.926 (**+0.006**). But it **drops `Y181C`, `Y181A`, `Y181I`** from
the EFV catalog. Y181C is among the most firmly established NNRTI resistance mutations there is; its
independent EFV coefficient is modest because it co-occurs heavily with K103N, which absorbs the effect.
**Any automated rule that removes Y181C from a shipped catalog is measuring something other than
resistance biology.** Blind-spot recovery 0.415 — below the incumbent.

**3. Additive-only, minimal, biologically anchored — still negative.** The most favourable framing
available: keep every shipped entry, add only deconfounding-survivors at *already-catalogued* positions
(`K101E`, `K101H`, `K103S`) plus `V179D`, a named blind-spot driver that survived deconfounding
independently. `G335D` excluded (connection domain, no NNRTI association).

| drug | Δ balanced accuracy | blind-spot recovered |
|---|---|---|
| efavirenz | **−0.004** | 0.453 |
| nevirapine | +0.005 | 0.323 |
| etravirine | +0.007 | 0.500 |
| rilpivirine | **−0.021** | 0.286 |

Mean **≈ −0.003**. Two drugs lose. And these are **in-sample optimistic** — the four entries were chosen
knowing this data — so the true out-of-sample effect is no better than shown.

**Every blind-spot recovery figure across all three variants (0.000–0.500) sits below the position-novelty
flag's measured 0.604.** On the final artifact, etravirine recovers **zero**.

## A defect found in my own script, and what it cost

The first derivation admitted `L234L`, `K238K`, `M230M`, `R72R` as candidates — **WT and MUT are the same
letter**. These are not mutations; `CompMutList` uses them to encode a mixture containing wild-type or an
ambiguity. The OLS happily assigned them coefficients, because they are markers of *sequencing/mixture
status*, which correlates with treatment experience — a pure confound, and four of them landed in the
first derived catalog.

Excluding self-to-self entries changed EFV's 3× result from 8 additions to 5 and lifted balanced accuracy
0.929 → 0.932. **The headline would have been wrong in the flattering direction had I not read the entry
list.** Pinned by test.

## Why this is the right outcome, not a failure

The blind spot is genuinely there and the doubt layer already addresses it — **without touching the
shipped catalog, without a model, and at a higher recovery rate**. L2's whole design claim is that a
signal which *qualifies* a call can beat one that *changes* it; here that is now measured rather than
argued.

**This does not close curation in general.** What is closed is *data-derived* curation from this dataset
at any threshold tested. A **literature-anchored** addition — each entry sourced per-mutation to a named
authority rather than to an OLS coefficient — is a different proposal with different provenance, and
`V179D` emerging independently from the data is corroboration for it. That would be a curation project
with a fabrication hazard to manage, not a decoder change; it is not recommended on these numbers, and it
is not foreclosed either.

## Honest limits

In-distribution against HIVDB-PhenoSense, not provenance-disjoint. Censored folds are kept at their
numeric bound (inherited v0 convention). The additive-only numbers are in-sample. `G335D` was excluded on
a biological-association judgement, not a measurement. Doravirine is unscored for lack of a sourced
cutoff.

## Artifacts

`scripts/hiv_nnrti_mutant_catalog.py` → `wiki/hiv_nnrti_mutant_catalog_2026-09-01.json`.
Label: Stanford HIVDB PhenoSense fold-change (independent wet-lab, **not** Sierra). Cite Rhee 2003.
`hiv_amr.py` is **unmodified** — no catalog entry was added or removed.
