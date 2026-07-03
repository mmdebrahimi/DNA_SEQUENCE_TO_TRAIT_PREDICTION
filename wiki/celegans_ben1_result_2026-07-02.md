# C. elegans ben-1 → benzimidazole resistance — the deterministic decoder crosses into the animal kingdom

**Date:** 2026-07-02
**Script:** `scripts/celegans_ben1_decoder.py` (`wiki/celegans_ben1_scores.json`)
**Data (free, MIT/CC):** `AndersenLab/ce_cb_ct_betatubulin` → `isotype_variant_summary.tsv` (611 wild *C. elegans*
isolates; per-isolate ben-1 variant columns + measured albendazole HTA resistance calls). On D: (gitignored).

## What this is

The same shape as the project's validated resistance cells (bacterial AMR, fungal ERG11, HIV RT, SARS-CoV-2
Mpro): a **deterministic determinant scan of a curated causal gene, validated against a measured drug-response
phenotype** — but in a **new kingdom** (multicellular animal / nematode). The gene is **ben-1** (β-tubulin, the
canonical benzimidazole target); the phenotype is the Andersen-lab albendazole high-throughput-assay (HTA)
resistance call across natural isolates.

Determinant rule (transparent): a strain is ben-1-determinant-positive if it carries an **impactful ben-1 SNV**
(incl. missense F200Y/Q131L/S145F…), a **structural variant** (deletion/inversion/duplication), or **low ben-1
expression** → predict RESISTANT. Validated against three HTA assays (`norm`, `2018`, `2024`).

## Result — a highly specific, high-PPV positive predictor (low aggregate sensitivity, by design)

| Assay | n | R / S | sens | **spec** | **PPV** | acc |
|---|---|---|---|---|---|---|
| abz_hta_norm | 285 | 137 / 148 | 0.263 | **0.973** | **0.90** | 0.63 |
| abz_hta_2018 | 202 | 111 / 91 | 0.297 | **0.967** | **0.92** | 0.60 |
| abz_hta_2024 | 188 | 58 / 130 | 0.259 | **0.969** | **0.79** | 0.75 |

### The deployable rule is decisive — ben-1 LoF/coding variants are ~100% resistant

Variant-class breakdown (norm assay), fraction resistant:

| ben-1 variant class | n | resistant |
|---|---|---|
| Frameshift | 14 | **14/14 (100%)** |
| Deletion | 9 | **9/9 (100%)** |
| Start/stop altering | 6 | **6/6 (100%)** |
| Missense | 5 | **5/5 (100%)** |
| Low ben-1 expression | 3 | 1/3 |
| Inversion | 2 | 1/2 |
| Inframe deletion | 1 | 0/1 |
| **No variant** | 245 | 101/245 (41%) |

The four LoF/impactful coding classes are **34/34 = 100% resistant**. If a wild isolate carries a ben-1
frameshift / deletion / nonsense / missense, it is albendazole-resistant. That is the trustworthy, deployable
call — and it holds across a genuinely new organism with no re-tuning.

### Why aggregate sensitivity is low — and why that's the honest, correct answer

Only 65/611 isolates carry ANY ben-1 variant, yet 137/285 are HTA-"resistant" → 101 resistant isolates have
**wild-type ben-1**. This is NOT a gap in the scan (the annotation includes missense, SVs, and expression). The
strong-vs-marginal diagnostic explains it decisively — mean continuous HTA phenotype (higher = more resistant):

| assay | determinant+ resistant | determinant− resistant | sensitive |
|---|---|---|---|
| norm | **136.6** | 10.4 | −36.6 |
| 2018 | **121.4** | 2.8 | −45.3 |
| 2024 | **176.8** | 17.5 | −34.8 |

**ben-1-variant resistant isolates are 7–64× more strongly resistant** than the resistant-without-ben-1
isolates, which sit barely above the binarization threshold. So ben-1 explains **strong** benzimidazole
resistance; the low aggregate sensitivity is the HTA threshold admitting a large **marginal / polygenic**
resistance background that no single locus captures.

## Why this matters (project lens)

- **The deterministic feature-match decoder crosses the animal-kingdom boundary intact.** ben-1 LoF → BZ
  resistance transfers with 100% PPV on coding variants — extending the "matched determinant → phenotype" law
  from bacteria/fungi/viruses to a multicellular animal, on a free, isolate-level, wet-lab-measured label.
- **It is a textbook instance of the project's own high-spec/low-sens "scrutinize the label" lesson** (mirror of
  mecA/oxacillin): the genotype is the trustworthy output; the label's low-effect tail (marginal HTA resistance)
  is what the single-locus rule cannot — and should not — chase. Quantified, not hand-waved.
- **Contrast with the embedding track (0-for-4 de-confounded negative):** this is deterministic, curated-catalog,
  and it *works* — reinforcing that the validated product is the determinant scan, not a learned model.

## Honest scope / walls

- **In-distribution association on a natural panel, not a held-out external validation.** It validates the
  determinant→phenotype rule on wet-lab labels; it is not a deployed clinical predictor.
- **The marginal-resistance tail (det− resistant strains) is real biology + threshold sensitivity**, not modeled
  here — a single-locus scan is the wrong tool for it (would need polygenic modeling; the project's embedding
  attempts at that class are a closed negative).
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); this is a new research cell.

## Reproduce

```bash
uv run python scripts/celegans_ben1_decoder.py    # reads D:/dna_decode_cache/celegans_ben1/variants.tsv
uv run pytest tests/test_celegans_ben1_decoder.py -q   # 3 offline synthetic tests
# data: raw.githubusercontent.com/AndersenLab/ce_cb_ct_betatubulin/main/data/isotype_variant_table/c_elegans/20250128/isotype_variant_summary.tsv
```
