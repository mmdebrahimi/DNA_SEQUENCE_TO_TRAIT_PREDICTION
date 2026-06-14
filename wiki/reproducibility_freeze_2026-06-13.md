# Reproducibility freeze — deterministic AMR decoder (2026-06-13)

**Status: BANK-AND-SHIP.** This is the terminal honest product of the dna_decode project's AMR track. Every
expansion beyond it is closed (see `wiki/negative_results_map_2026-06-13.md`): embeddings 0-for-4,
MIC-continuous infeasible, pathotype label-blocked, provdisjoint AMR grid saturated. This freeze makes the
validated state **inspectable, rerunnable, and challengeable** so the tool is bankable as-is and serves as
a clean base for any future prospective-lock validation.

Frozen at commit **b3761c8** (origin/main, 2026-06-13).

## What is frozen (the load-bearing units of release)

| Unit | Path | What it is |
|---|---|---|
| **Decoder rules** | `dna_decode/data/calibrated_amr_rules.json` | the deployed per-(organism,drug) call_resistance config (counters, thresholds, intrinsic-family exclusions, EXPRESSION_FLOOR abstentions) |
| **Deployed-claim surface** | `dna_decode/data/shipped_decoder_surface.py` | the authoritative grid of what the tool claims to decode (coverage-tested vs the CLI drug catalogs) |
| **Per-drug catalogs** | `dna_decode/data/mic_tiers.py` | CLSI/EUCAST breakpoints, mechanism→loci catalogs, AMRFinder class filters |
| **Trust surface** | `wiki/decoder_validation_report_card.{md,json}` | the standing report card: 10 SCORED + lineage-disclosure table + honest tier strings |
| **Scored-cell artifacts** | `wiki/provenance_disjoint_validation_*.json` (10) | per-cell acc/sens/spec + confusion + leakage-control provenance |
| **Lineage metrics** | `wiki/provdisjoint_lineage_metrics.json` | clonality-corrected cluster-weighted sens/spec + Wilson CI per cell |
| **Accession lists** | `data/raw/*_provdisjoint_*/selected.tsv` (10 cohorts, 580 accessions) | the exact FRESH leakage-clean isolates each cell was scored on |
| **Leakage registry** | `dna_decode/eval/cohort_manifest.py` | data-driven exact-identity exclusion across all cohorts (fail-closed on incomplete) |
| **Env lockfile** | `uv.lock` (+ `pyproject.toml`) | exact dependency pins for `uv sync` |

## The 10 SCORED cells (the validated grid)

campylobacter·cipro · escherichia_coli_shigella·{cipro, ceftriaxone, gentamicin, tetracycline} ·
klebsiella·{cipro, ceftriaxone, gentamicin, meropenem, tetracycline}

Honest tier: isolate-level provenance-disjoint stress test; R classes clonally dominated; lineage-effective
N + cluster-weighted metrics (with Wilson CI) disclosed. NOT lineage-independent external clinical
validation. Read the report card cell-by-cell — there is deliberately no aggregate "X% validated" headline.

## One-command reproduction path

```bash
# 1. Environment (exact pins)
uv sync

# 2. Rerun the full validated test suite (the rerunnability guarantee)
uv run pytest tests/ -q --ignore=tests/test_models_foundation.py
#    test_models_foundation.py is excluded: it loads torch and fails to COLLECT on a memory-constrained
#    host (OSError WinError 1455 paging-file) — a host limit, NOT a decoder regression. All decoder /
#    validation / clonality / report-card tests run without it.

# 3. Rebuild the trust surface from the on-disk scored artifacts (read-only roll-up; no network/Docker)
uv run python scripts/build_validation_report_card.py

# 4. (Docker host only) re-derive lineage metrics from the frozen cohorts
uv run python scripts/compute_lineage_metrics.py      # needs Docker for Mash

# 5. Decode a strain with the shipped tool
uv run python -m scripts.pipeline predict --drug ciprofloxacin --strain-id <id> --cache <h5> ...
#    or the installed console scripts: dna-decode / dna-amr / dna-pathotype
```

To re-derive a SCORED cell end-to-end from scratch (Docker + network):
```bash
uv run python scripts/provenance_disjoint_validate.py --group <Group> --amrfinder-organism <Org> \
  --drug <drug> --per-class 30 --registry-organism <Group>
# Re-selection is forced onto FRESH accessions by the leakage registry (exact-self cohort identity),
# so a re-run cannot silently reuse a prior cell's isolates.
```

## What this freeze is NOT

- NOT a claim of clinical/external validation — it is a provenance-disjoint stress test (see tier strings).
- NOT lineage-independent — the cipro QRDR cells are clonality-inflated at isolate level (disclosed in the
  report card's lineage table; weighted metrics + CIs tell the honest story).
- NOT the end of the project's vision — it is the terminal point of what *free public labels* can honestly
  support. The two non-foreclosed forward paths (both gated on a user decision, neither an executor task):
  1. **Acquisition** — a non-public wet-lab/clinical label source in hand (clears the label gates by
     construction). Draft anchor at `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`.
  2. **Prospective-lock** — this frozen decoder + a pre-registered protocol scoring later-arriving
     independent isolates as they appear. Needs no new label today; this freeze is its baseline.

## Provenance
- Negative-results map (why expansions were rejected): `wiki/negative_results_map_2026-06-13.md`
- Embedding closure (0-for-4): `wiki/embedding_niche_cross_domain_synthesis_2026-06-12.md`
- Strategy review that recommended this freeze: brainstorm 2026-06-13 (this session).
