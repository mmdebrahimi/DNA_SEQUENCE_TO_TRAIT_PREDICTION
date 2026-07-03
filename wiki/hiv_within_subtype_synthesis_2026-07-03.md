# HIV within-subtype de-confounding — the rail now cleared across ALL classes (2026-07-03)

**Finding: the deterministic HIV resistance catalogs are MECHANISM, not subtype structure — verified
within subtype B for all four target classes (NRTI, NNRTI, PI, INSTI).** This closes the within-subtype
gap the class validators explicitly flagged, and it is the HIV analogue of the project's central
within-lineage rail (overall AUROC conflates population structure with mechanism; the honest test is *within*
a single group). Script: `scripts/hiv_within_subtype.py`.

## Why this matters (the de-confounding question)

The NNRTI/PI/INSTI cells were validated on a subtype-MIXED set (~96% subtype B + a non-B tail). A high
mixed AUC could, in principle, reflect the catalog tracking subtype structure rather than the resistance
mechanism. The clean test: restrict to a SINGLE subtype (B, the well-powered arm) and re-compute the same
cutoff-free AUC. If it holds within B, the catalog is decoding mechanism; if the pooled number is materially
higher than within-B, it was riding structure. NRTI cleared this on 2026-06-21; this run closes the other
three.

## Result — all four classes HOLD within subtype B

| Class | call mode | median within-B AUC | pooled − within-B (subtype-structure contribution) | verdict |
|---|---|---|---|---|
| NRTI (2026-06-21) | position-based v0 | balacc 0.63–0.80 within B | ≈ 0 (B ≈ non-B) | HOLDS (prior) |
| **NNRTI** | mutant-level | **0.795** | **−0.0005** | **HOLDS_WITHIN_SUBTYPE** |
| **PI** | position-based v0 | **0.921** | **−0.006** | **HOLDS_WITHIN_SUBTYPE** |
| **INSTI** | position-based v0 | **0.898** | **+0.043** | **HOLDS_WITHIN_SUBTYPE** |

(Cutoff-free AUC = P(fold of a called-R isolate > fold of a called-S isolate), PhenoSense fold label —
independent wet-lab IC50, NOT HIVDB's own Sierra interpretation. `Method=PhenoSense AND Type=Clinical`.
Per-drug tables in `wiki/hiv_{nnrti,pi,insti}_within_subtype_2026-07-03.md`.)

**The pooled − within-B delta is ≈ 0 for every class** (NNRTI −0.0005, PI −0.006, INSTI +0.043; all far
below the 0.10 subtype-inflation threshold). The class-mixed numbers were NOT subtype-inflated — restricting
to a single subtype barely moves the AUC. For NNRTI and PI the within-B AUC is even *slightly higher* than
pooled (the small, noisier non-B tail drags the pooled number down, not up).

## Honest scope (unchanged rails)

- **Within-B is the powered arm; non-B is under-powered** (the free HIVDB data is ~96% B: non-B pooled
  n ≈ 39–154 per class). A within-B result is a de-confounding test (signal survives inside one subtype);
  it does NOT prove non-B generalisation at scale — that needs a deeper non-B cohort (an external/label
  wall, not a code wall).
- **PI/INSTI are POSITION-BASED v0** (deliberate over-call at major positions) → their AUC has a built-in
  ceiling below a mutant-specific catalog. The de-confounding readout is the within-vs-pooled DELTA (near
  zero), not the absolute level. The PI/INSTI mutant-specific v0.1 catalogs (2026-06-23 / 2026-06-27) lift
  the level; this check confirms the level that exists is mechanism.
- **CAI (capsid / lenacapavir)** is deferred: single-drug, resistance-enriched (n≈140), too thin for a
  subtype split.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); this is a READ-only validation
  over the frozen `dna_decode.data.hiv_amr` dispatch. HIV's own trust surface is namespace-separate from the
  bacterial NCBI-PD card.

## What it adds to the independent-validation frontier

The HIV cell was already the project's first free-independent-label win (PhenoSense fold ≠ Sierra). This
strengthens it: the win is now shown to be de-confounded (within-subtype) across the entire HIV target-class
surface, not just NRTI. It is a genuinely new number on the one track where free independent labels exist —
not a re-confirmation.

## Reproduce

```bash
# .Full HIVDB datasets carry the Subtype column (regular filtered sets do not):
#   curl -L -o data/raw/hiv/<CLASS>_DataSet.Full.txt \
#     https://hivdb.stanford.edu/download/GenoPhenoDatasets/<CLASS>_DataSet.Full.txt   # NNRTI / PI / INI
uv run python scripts/hiv_within_subtype.py --class all
uv run pytest tests/test_hiv_within_subtype.py -q    # 8 offline tests (no network)
```
