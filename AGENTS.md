# AGENTS.md

A README **for coding agents** (Claude Code, Codex, Cursor, OpenHands, …). Setup, smoke checks, runnable
examples, and the scientific guardrails an agent must preserve. The human landing page is `README.md`.

## Project overview

`dna-decode` is a **research-use** Python CLI + library for deterministic, interpretable microbial
**genome → trait** decoding. For each call it reports:

- antibiotic / antiviral / antifungal **resistance R/S** where supported (bacteria, *M. tuberculosis*, *Candida auris*, HIV-1, SARS-CoV-2, influenza),
- the **exact genes or mutations** driving the call,
- typing / pathotype calls (E. coli pathotype, MLST, serotype, Inc-plasmid, K-type, …),
- a **validation tier + caveats + provenance** for every call, and
- **abstention** when the mechanism is outside the validated scope.

Mechanism-feature / rule based — **not an embedding black box, not a clinical diagnostic tool.**

> The authoritative, always-current list of supported traits + their validation tiers is the command
> `dna-decode list` — prefer it over any hand-written table (including this file), which can drift.

## Install

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
# or, from PyPI:  pip install dna-decode
```

This project also supports `uv` (recommended): `uv sync` then prefix commands with `uv run`.

## Smoke checks

```bash
dna-decode list          # the authoritative supported-trait + validation surface
pytest tests/ -q         # ~1770 tests; the deterministic decoder needs no torch/GPU
```

## Minimal runnable examples (all verified, no Docker / no GPU)

```bash
# 1. Zero-dependency: an observed mutation -> R/S + provenance (wheel-only, no genome/DB)
dna-decode amr --drug efavirenz --observed RT:K103N
#   -> CALL: R  driven by RT:K103N  validation: INDEPENDENT_WETLAB (Stanford HIVDB PhenoSense)

# 2. From a committed AMRFinder run (bacterial; no Docker needed for this path)
dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_000492655.1
#   -> CALL + driven-by determinants + validation tier (INDEPENDENT_MEASURED, EBI AMR Portal)

# 3. List everything the tool supports + each cell's honest validation tier
dna-decode list
```

Genome-FASTA input (e.g. `dna-decode amr --drug ceftriaxone --genome-fasta X.fna`) additionally needs
Docker (AMRFinder) for the bacterial path or BLAST+ for the target-site engines — see `docs/quickstart.md`.

## Rules for agents (load-bearing — the project's scientific caution is explicit, not inferred)

- **Never** describe this as a clinical / diagnostic / medical tool. It is research-use only.
- **Never** claim organisms or drugs that `dna-decode list` does not show as supported.
- **Preserve the validation tier + caveats + provenance** in any generated docs or summaries — never
  average distinct tiers into one aggregate score (the project forbids an aggregate headline by design).
- Prefer **deterministic / rule-based** language. Do **not** call this "AI", "predicts phenotype from DNA",
  a "breakthrough", or "better than AMRFinder" (over-broad) — see `README.md` "What not to say".
- Honor **abstention**: when the tool abstains or says NOT_CENSUSED / no-free-source, report that, don't
  invent a call.
- When changing decoder logic: update tests **and** the validation docs; **run `pytest tests/ -q`** before
  proposing a final patch. The frozen AMR surface (`dna_decode/eval/amr_rules.py` +
  `dna_decode/data/calibrated_amr_rules.json`) is sha256-pinned — do not edit it without intent.

## Project map (where things live)

- `dna_decode/amr/` `dna_decode/pathotype/` `dna_decode/mlst/` … — per-trait decoders (CLI: `dna-<trait>`).
- `dna_decode/cli.py` — the unified `dna-decode <trait>` dispatch (TRAITS).
- `dna_decode/data/` — curated catalogs + the frozen rule surface + the Evidence-Contract Registry (`cell_registry.py`).
- `wiki/*report_card*.{md,json}` — per-domain validation surfaces (the provenance the calls cite).
- `tests/` — pytest; `docs/` — quickstart + validation; `CLAUDE.md` — maintainer internal notes (not agent onboarding).
