# DNA_SEQUENCE_TO_TRAIT_PREDICTION

Genotype-to-phenotype (G2P) inference platform — predicts phenotypic traits from genomic DNA sequences AND identifies which genomic regions are most strongly associated with those predictions. Biologically interpretable, not causal-claim-making.

## Status: pre-Phase-1

No code yet. Planning artifacts only:

- `plans/Ecoli_G2P_Platform_Technical_Plan.md` — Phase 1 technical plan (17 steps in 8 waves; E. coli antibiotic resistance for ciprofloxacin + ceftriaxone + tetracycline)
- `project_state/dna-decode-2026-05-11.md` — project planning ledger (Bellman-inspired decision frame; 11 hypotheses; 4 decisions made; 12 action-log entries)

## Phase 1 scope

| Aspect | Value |
|---|---|
| Organism | E. coli |
| Phenotypes | Ciprofloxacin, ceftriaxone, tetracycline binary resistance |
| Foundation models | Evo (primary), DNABERT-2, Nucleotide Transformer, GENA-LM (leaderboard) |
| Baseline | Frozen embeddings + XGBoost per drug |
| Attribution | In-silico mutagenesis (gene-level + nucleotide-level saturation) |
| CV | Leave-one-Mash-clade-out + clade-only baseline + per-clade reporting |
| Target | AUROC ≥0.80 SLO / ≥0.85 stretch; clade-baseline-gap ≥0.10 on ≥75% of held-out clades |
| Horizon | 3 months Phase 1; 12 months Phase 1+2+3 |
| Compute | Single RTX 4090 with 4-bit Evo quantization (bitsandbytes); A100 fallback |

## Long-term vision

Multimodal genotype-phenotype platform — start with bacterial AMR (Phase 1), expand toward eukaryotes + image-paired phenotype data in later phases. NOT a direct stepping stone to "DNA → animal image" prediction; that would require a parallel multimodal track.

## Setup

```bash
# Install dependencies (after Step 1 bootstrap fires via /execute-plan)
uv sync

# Install Mash CLI (Step 10 phylogeny clustering dependency)
# Linux/WSL: apt install mash
# Windows: download binary from https://github.com/marbl/Mash/releases

# Run real-data pilot gate (Step 0.5 — HARD gate)
uv run python -m scripts.pilot_gate --drugs ciprofloxacin,ceftriaxone,tetracycline --target-per-drug 150

# Run smoke pipeline (synthetic-fixture end-to-end, no real downloads)
uv run python scripts/smoke_pipeline.py
```

## Project workflow

Built using a personal Claude Code skill ladder for project planning:
- `/idea-anchor` → `/project-init` → `/brainstorm` → `/technical-plan` → `/execute-plan`
- Project ledger maintained via `/project-state` skill
- All planning artifacts captured as audit trail

See `project_state/dna-decode-2026-05-11.md` for full decision history.
