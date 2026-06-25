# Session handoff — 2026-06-24 (real PyPI PUBLISHED + ktype validation running)

Pick up cleanly in a fresh session. Git is **clean + pushed** (origin/main `0/0`), HEAD = **`4f3504f`**,
ledger at **row 195**, project ledger `project_state/dna-decode-2026-05-11.md`.

> Supersedes `wiki/SESSION_HANDOFF_2026-06-24_productization_ktype_testpypi.md` (the TestPyPI-era handoff).
> The two open items it listed (real-PyPI publish + full ktype validation) are now (B) DONE and (C) running.

---

## 0. One-paragraph orientation
`dna_decode` is a **deterministic, interpretable genome->phenotype decoder** (NOT a learned/embedding model —
that bet is a closed 0-for-4 negative). It calls antibiotic/antiviral/antifungal **resistance R/S** from DNA
+ names the exact genes/mutations driving each call, across **bacteria, M. tuberculosis, fungi, HIV-1,
SARS-CoV-2** + typing decoders (serotype, MLST, plasmid, the new Klebsiella **K-antigen/ktype**), each cell
carrying its own **honest validation tier** (independent-measured / in-distribution / faithful-to-tool /
no-free-source — never conflated). It is now **PUBLICLY pip-installable**: `pip install dna-decode`. North
star: "an AI DNA decoder **tool**, not papers." Binding constraint (still): a **free, independent, measured**
phenotype label.

## 1. What THIS session shipped
1. **`dna-decode 0.5.1` PUBLISHED to REAL public PyPI** — https://pypi.org/project/dna-decode/0.5.1/
   - Verified: a fresh-venv `pip install dna-decode` straight from pypi.org installs + runs + serves the
     `INDEPENDENT_WETLAB` trust badge from the **packaged** card (not the repo wiki/).
   - The prior **403** was diagnosed + fixed: the name was FREE on pypi.org (not squatting) -> the 403 was
     the TOKEN (interactive-prompt paste mangled control chars + TestPyPI!=PyPI separate accounts). Fix was
     a real pypi.org token + the **env-var upload method** (`TWINE_USERNAME=__token__` + `TWINE_PASSWORD=...`
     exported, then `twine upload`), which avoids the prompt entirely.
2. **Ledger rows 194-195** record the publish + verification.

## 2. OPEN ITEMS — what to continue (ranked)
| # | Item | State | Next action |
|---|---|---|---|
| **A** | **REVOKE both API tokens** | a TestPyPI token AND a real-PyPI token are in prior transcripts | **Do this.** Real PyPI: https://pypi.org/manage/account/token/ -> revoke `dna-decode` token. TestPyPI: https://test.pypi.org/manage/account/token/ -> revoke. The real-PyPI one is live + can publish publicly under the account. For future publishes make a fresh **project-scoped** token each time. |
| **B** | **Fold the full ktype number into the report card** | the cohort validation run is IN FLIGHT (background; **178/447** at handoff time) | When `D:/dna_decode_cache/ktype_build/COHORT_DONE` exists: read `wiki/ktype_cohort_validation.json` (written by the runner on completion) -> the real wzi-caller-vs-serology concordance (N up to 447) -> fold it into `wiki/ktype_report_card.{md,json}` next to the existing measured-label ceiling (0.745 naive / 0.845 curated). Commit + ledger row. **Honesty caveat baked into the runner:** concordance is naive `KL#==K#` (a lower bound; curated KL<->K equivalence lifts it) + wzi single-gene v0 (~94% of full Kaptive). |
| C | Next version bump if you change anything | 0.5.1 is published + permanent | a published version can NEVER be re-uploaded; bump `[project] version` in `pyproject.toml` before the next `twine upload`. |
| D | More GREEN non-AMR cells | K-antigen was the first; triage gate in `plans/Non_AMR_Phenotype_Pivot_Assessment_2026-06-24.md` | apply the gate [curated catalog exists? -> free measured label exists?] to any new trait BEFORE building. |

## 3. The background ktype run (how to check / resume)
- **Runner:** `scripts/ktype_cohort_validate.py` (per ERR isolate: `fetch_reads` -> `map_erg11(WZI_REF)`
  targeted single-locus read-mapping -> `call_ktype` -> concordance vs `Serological_Ktype`). **Checkpointed**
  to `D:/dna_decode_cache/ktype_build/cohort_results.jsonl` (restartable — skips done accessions); reads are
  deleted per-isolate to bound disk.
- **Self-healing wrapper:** `scripts/ktype_cohort_finisher.sh` loops the runner until 447 scored, then writes
  `D:/dna_decode_cache/ktype_build/COHORT_DONE`. Log at `D:/dna_decode_cache/ktype_build/cohort.log`.
- **Check progress:** `wc -l < D:/dna_decode_cache/ktype_build/cohort_results.jsonl` (count scored).
- **If it died:** just re-run `bash scripts/ktype_cohort_finisher.sh` (or the runner directly) — it resumes
  from the checkpoint. Docker Desktop must be up (targeted-map uses minimap2/samtools/SRA biocontainers).
- **Cohort scope:** only the 447 isolates with an ERR accession are fetchable (283 'TBD' + 3 unavailable are
  skipped — reported, not hidden) out of the N=733 measured-serology set.

## 4. How to continue (commands + gotchas)
- **The venv is LEAN** (default deps only — torch/transformers/xgboost split to `[ml]`). The deterministic
  decoder + its tests run fine on it. Full suite / embedding work: `uv sync --extra ml --extra dev` (multi-GB;
  C: is disk-tight).
- Run a decoder: `uv run dna-decode list` · `uv run dna-amr --drug efavirenz --observed RT:K103N` ·
  `uv run dna-decode ktype <assembly.fna> --db-dir data/ktype_db`.
- Tests (deterministic): `uv run --no-sync pytest tests/ -q`.
- Verify the wheel ships the cards: `uv run python scripts/verify_wheel_ships_cards.py --fresh-env`.
- **Publish runbook:** `wiki/pypi_publish_runbook.md`. The working upload method (avoids the 403):
  ```bash
  cd C:/Users/Farshad/PythonProjects/dna_decode
  export TWINE_USERNAME='__token__'
  export TWINE_PASSWORD='pypi-...'   # a REAL pypi.org token (NOT test.pypi.org)
  uv run --with twine --no-sync python -m twine upload --non-interactive dist/dna_decode-<version>*
  ```
  (Real PyPI = plain `twine upload`, no `--repository`. TestPyPI = `--repository testpypi`. Separate
  accounts/tokens — a TestPyPI token gives 403 on real PyPI.)

## 5. Standing discipline (do NOT break)
- **FROZEN AMR surface**: `dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` are
  byte-frozen + sha-pinned (leak-guard `tests/test_tb_leak_guard.py` + prospective-lock manifest). Any change
  = a ratified freeze-amendment re-pinning the shas. Most work is ADDITIVE around it.
- **Honesty tiers**: never relabel in-distribution as independent; each surface has its OWN namespace-separate
  report card; the trust-surface never fabricates a metric or borrows across organisms/genus.
- **Gates** (Soraya): money -> hard pause; publish/send -> user-authority (the 0.5.1 publish was explicitly
  authorized + token-provided); destructive-local -> pause. Verify-in-batch every built artifact.
- **Cross-machine**: commit to `main` is the sync channel; Precision 7780 ("workhorse") occasionally pushes —
  `git fetch` + rebase before push if divergent.

## 6. Key artifacts map (this session)
- Publish: `dist/dna_decode-0.5.1*` (gitignored); `pyproject.toml` (0.5.1, author Eman, `[ml]` extra,
  force-include of 4 report cards -> `dna_decode/report_cards/`); `wiki/pypi_publish_runbook.md`.
- ktype: `dna_decode/ktype/`; `scripts/ktype_cohort_validate.py` + `scripts/ktype_cohort_finisher.sh`;
  `wiki/ktype_report_card.{md,json}`; `wiki/ktype_cohort_validation.json` (written on cohort completion);
  plan `plans/Klebsiella_K_antigen_capsule_typing_feasibility_2026-06-24.md`.
- Trust surface: `dna_decode/data/trust_surface.py`. Cards: `wiki/*_report_card.{md,json}`.
- Ledger: `project_state/dna-decode-2026-05-11.md` rows 180-195.

## 7. Recommended first move next session
1. **Revoke the two tokens** (item A — 1 min, do it).
2. When `COHORT_DONE` appears, **fold the ktype number into the report card** (item B) + commit.
3. Otherwise bank — the session delivered a **publicly pip-installable, honest, multi-pathogen decoder**.
