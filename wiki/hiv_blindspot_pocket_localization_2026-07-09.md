# ΔΔG go/no-go pre-check: is the EFV blind spot pocket-mediated? — VERDICT: GO (qualified)

**Date:** 2026-07-09 · **Zero-tool, D-free** (no FoldX/Rosetta/structures/embeddings) ·
`scripts/hiv_blindspot_pocket_localization.py` · data: `data/raw/hiv/NNRTI_DataSet.txt`

## The question

Before spending a FoldX/Rosetta install on a ΔΔG pilot, decide cheaply: are the EFV **catalog-blind-spot
resistant** isolates resistant via **NNRTI-pocket (binding-site)** mutations a binding-ΔΔG scorer *could*
see — or via non-pocket mechanisms it's structurally blind to? A prior brainstorm predicted the blind
spot would be **adversarial** to binding scorers, because it is *defined* by the absence of any major DRM
at the 8 primary pocket positions {100,101,103,106,181,188,190,230}.

## Result — the brainstorm's prediction is FALSIFIED

The blind-spot R are **strongly enriched for known-functional NNRTI-pocket mutations** vs blind-spot S:

| metric | R (n=53) | S (n=1058) | R/S |
|---|---|---|---|
| ≥1 known **in-window** secondary-pocket NNRTI mut (V108I/V179x/P225H/F227x/L234I/P236L) | **41.5%** | 5.0% | **8.3×** |
| ≥1 known **functional** NNRTI mut (+ ext 98/138/221/238/318) | **60.4%** | 14.4% | **4.2×** |
| ≥1 **any** window-position mut (naive — burden control, ~11% chance/mut) | 92.5% | 67.1% | 1.38× |

**Burden-adjusted functional enrichment = 4.2 / 1.38 = 3.05×.** Blind-spot R carry more mutations
(median 15/isolate), but even after dividing out that burden effect, they remain **~3× enriched for
functional pocket mutations specifically** — this is not a burden artifact.

The specific drivers (count in R): **V179D ×12, A98G ×10, H221Y ×7, F227C ×5, V108I ×4, V179E ×3** — all
genuine NNRTI accessory/secondary mutations that reduce drug binding. So **uncatalogued binding-site DRMs
DO exist on the blind spot**, and a binding-ΔΔG scorer has real targets there.

**VERDICT: GO** — a FoldX/Rosetta pilot is warranted (the blind spot is *not* the adversarial test set the
brainstorm feared).

## Why "GO" is qualified — three honest caveats

1. **Enrichment ≠ ddG will succeed.** These are *secondary/accessory* mutations (V179D, A98G, V108I),
   subtler than the majors — their binding effects are smaller, and cheap-ddG RMSE (1–1.9 kcal/mol =
   5–24× fold-change error) may not resolve them above the >1 kcal/mol hotspot bar. GO means "the targets
   exist," not "the tool will clear the bar." The pilot must still prove ddG *recovers* these.
2. **Coverage ceiling ~60%.** ~40% of blind-spot R carry **no** known-functional pocket mutation → still
   ddG-blind (non-pocket mechanism or a truly novel mutation outside the curated set). A successful ddG
   pilot addresses at most the pocket-mediated ~60%.
3. **Passenger risk.** A single weak accessory (A98G, V108I) may co-occur with the real non-pocket cause
   rather than drive resistance itself; the 3× burden-adjusted enrichment argues they contribute, but
   single-weak-accessory isolates are ambiguous.

Two of the enriched drivers (A98G at 98, H221Y at 221) fall just *outside* the user's literal windows
(<100, <225); the strict in-window-only reading is still **41.5% / 8.3×** — GO holds either way.

## Recommended next move

GO to the ddG pilot **when tools + D: return**, but frame it correctly: first prove ddG **recovers the
cataloged majors** (K103N/Y181C) as >1 kcal/mol hotspots (validity gate), THEN test whether it flags the
enriched *secondary* pocket mutations (V179D/F227C/V108I) on the blind-spot R. If cheap ddG can't resolve
the secondary mutations (likely, given their subtlety vs its RMSE), that's the real wall — not the test
set. The strongest single validation target is **V179D** (12 blind-spot R carry it).

## Provenance

- `wiki/hiv_blindspot_pocket_localization_2026-07-09.json` — machine-readable.
- Curated secondary/functional NNRTI mutation set: Stanford HIVDB accessory NNRTI positions.
- This reverses the pessimistic fatal-risk call in the prior ΔΔG brainstorm; the cheap pre-check earned
  its keep by overturning an assumption with data before any tool install.
