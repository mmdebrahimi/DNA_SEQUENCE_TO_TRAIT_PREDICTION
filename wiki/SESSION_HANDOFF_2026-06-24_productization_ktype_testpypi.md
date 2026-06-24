# Session handoff — 2026-06-24 (productization + ktype cell + TestPyPI)

Pick up cleanly in a fresh session. Git is **clean + fully pushed** (origin/main `0/0`), HEAD = **`2abaddd`**,
ledger at **row 193**, project ledger `project_state/dna-decode-2026-05-11.md`.

---

## 0. One-paragraph orientation
`dna_decode` is a **deterministic, interpretable genome->phenotype decoder** (NOT a learned/embedding model —
that bet is a closed 0-for-4 negative). It calls antibiotic/antiviral/antifungal **resistance R/S** from DNA
+ names the exact genes/mutations driving each call, across **bacteria, M. tuberculosis, fungi, HIV-1,
SARS-CoV-2**, each cell carrying its own **honest validation tier** (independent-measured / in-distribution /
faithful-to-tool / no-free-source — never conflated). It is now a **pip-installable tool** (the productization
arc completed this session). North star: "an AI DNA decoder **tool**, not papers." Binding constraint (still):
a **free, independent, measured** phenotype label.

## 1. What THIS session shipped (the arc)
1. **Cross-kingdom validation summary** (`wiki/cross_kingdom_validation_summary.md`) — one legible view of all
   validation surfaces, no aggregate headline.
2. **HIV PI v0.1** deconfounded refinement (`scripts/hiv_pi_mutant_catalog.py`) — 8/8 PIs improve, mean +0.056.
3. **Productization moves 1-3 + the packaging gate** (the big one):
   - **Inline trust-surface** (`dna_decode/data/trust_surface.py`) — every `dna-amr` call carries its
     `validation:` tier badge inline (never fabricated, never borrowed across organisms/genus).
   - **Unified `dna-decode profile`** — one genome -> AMR R/S + typing decoders + trust badges in one report.
   - **Dependency split** — torch/transformers/xgboost -> `[ml]` extra; the deterministic CLI is the lean default.
   - **Packaging gate (Path A, user-chosen)** — the wheel now **ships the trust cards as package data**
     (`pyproject` force-include -> `dna_decode/report_cards/`; `trust_surface._card_path` dual-path loader).
     Verified: a fresh-env wheel install serves correct badges from the PACKAGED cards.
   - **Shigella exact-organism-match fix** + honesty-wording cleanup.
4. **First GREEN non-AMR cell BUILT: Klebsiella K-antigen (wzi) capsule typing** (`dna_decode/ktype/`) — the
   9th decoder, `serotype`-sibling. Label gate CLEARED (free measured serological label exists). Caller +
   tests + CLI + namespace-separate `wiki/ktype_report_card.md`.
5. **TestPyPI dry-run SUCCEEDED** — `dna-decode 0.5.0` uploaded + installed-from-index + ran correctly.
   Package metadata made real-PyPI-ready (author **Eman <zoghalefaal@gmail.com>**, URLs, **bumped to 0.5.1**).

## 2. OPEN ITEMS — what to continue (ranked)
| # | Item | State | Next action |
|---|---|---|---|
| **A** | **Revoke the TestPyPI token** | a sandbox token is in the prior transcript | test.pypi.org -> Account -> API tokens -> revoke. **Do this.** |
| **B** | **Real PyPI publish** (0.5.1) | user-authority; package is publish-READY (`twine check` PASSED, metadata done) | user creates a **real** PyPI token, then `cd` project + `uv run --with twine python -m twine upload dist/dna_decode-0.5.1*` (username `__token__`). Runbook: `wiki/pypi_publish_runbook.md`. NEVER auto-publish — explicit user go-ahead + their token. |
| **C** | **Full ktype validation** (the real number) | caller shipped + measured-label CEILING computed (0.745 naive / 0.845 curated, N=733); the genuine wzi-caller-vs-serology number is NOT run | fetch the 731 ENA-run genomes (`WGS_accession` in the Zenodo xlsx) -> targeted wzi read-mapping (reuse `scripts/assemble_sra_cohort.py --method map`, the project's single-locus pattern) -> run `dna-ktype` -> concordance vs `Serological_Ktype`. Multi-hour cohort op -> best in a long-run/unattended window. Stage on **D:**. Feasibility + plan: `plans/Klebsiella_K_antigen_capsule_typing_feasibility_2026-06-24.md`; validator scaffold `scripts/ktype_validate.py`. |
| D | More GREEN cells | K-antigen was the first; the triage gate is in `plans/Non_AMR_Phenotype_Pivot_Assessment_2026-06-24.md` | apply the gate [catalog? -> free measured label?] to any new trait BEFORE building |

## 3. How to continue (commands + gotchas)
- **The venv is LEAN** (default deps only — no torch/transformers/xgboost; they were split to `[ml]`). The
  **deterministic decoder + its tests run fine** on it. To run the FULL test suite or do embedding/foundation
  work: `uv sync --extra ml --extra dev` (restores torch etc. — multi-GB; C: is disk-tight, ~18 GB free).
- Run a decoder: `uv run dna-decode list` · `uv run dna-amr --drug efavirenz --observed RT:K103N` ·
  `uv run dna-decode ktype <assembly.fna> --db-dir data/ktype_db`.
- Tests (deterministic): `uv run --no-sync pytest tests/ -q` (lean venv; ml-track tests need `--extra ml`).
- Verify the wheel ships the cards: `uv run python scripts/verify_wheel_ships_cards.py --fresh-env`.
- Verify the quickstart paths: `uv run python scripts/verify_quickstart.py`.
- ktype DB is **gitignored** (`data/ktype_db/wzi.fasta` + `wzi.txt`); re-fetch from Kleborate raw if missing
  (URLs in `dna_decode/ktype/cli.py` docstring). Also staged at `D:/dna_decode_cache/ktype_build/`.

## 4. Standing discipline (do NOT break)
- **FROZEN AMR surface**: `dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` are
  byte-frozen + sha-pinned (leak-guard `tests/test_tb_leak_guard.py` + prospective-lock manifest). Any change
  = a ratified freeze-amendment re-pinning the shas. Most work is ADDITIVE around it.
- **Honesty tiers**: never relabel in-distribution as independent; each surface has its OWN namespace-separate
  report card; the trust-surface never fabricates a metric or borrows across organisms.
- **Gates** (Soraya): money -> hard pause; publish/send -> user-authority (TestPyPI was explicitly authorized
  + token-provided); destructive-local -> pause. Verify-in-batch every built artifact (it caught a TOML bug
  twice this session).
- **Cross-machine**: commit to `main` is the sync channel; another machine (Precision 7780 / "workhorse")
  occasionally pushes — `git fetch` + rebase before push if divergent (happened cleanly twice this session).

## 5. Key artifacts map (this session)
- Decoders: `dna_decode/{amr,profile,ktype,...}/`; trust surface `dna_decode/data/trust_surface.py`.
- Validation cards: `wiki/{cross_kingdom_validation_summary,ktype_report_card,*_report_card}.md`.
- Plans: `plans/Non_AMR_Phenotype_Pivot_Assessment_2026-06-24.md` (GREEN/RED triage),
  `plans/Klebsiella_K_antigen_capsule_typing_feasibility_2026-06-24.md` (the cell's GO + scale-up plan).
- Publish: `wiki/pypi_publish_runbook.md`; build `pyproject.toml` (0.5.1, author Eman, `[ml]` extra,
  force-include report cards); artifacts in `dist/` (gitignored).
- Ledger: `project_state/dna-decode-2026-05-11.md` rows 180-193.

## 6. Recommended first move next session
Either (B) **publish to real PyPI** (1 user-token step — the tool is ready), or (C) **the full ktype
validation run** (a long-run window). Both are scoped + ready. Or bank — the session delivered a complete,
honest, pip-installable tool + the first validatable non-AMR cell.
