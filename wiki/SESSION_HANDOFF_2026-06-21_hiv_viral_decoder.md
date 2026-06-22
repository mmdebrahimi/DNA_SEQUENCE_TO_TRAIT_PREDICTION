# Session handoff — 2026-06-21: HIV viral decoder (Wave B) + genome-map virulence overlay

Canonical fresh-session entry point for everything landed in the long 2026-06-21 session. Repo:
`C:\Users\Farshad\PythonProjects\dna_decode` (origin `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`).
**Everything is committed + pushed to `origin/main` (synced, 0 ahead/behind).** Nothing is owed.

---

## ⇒ START HERE (1-paragraph state)

The session's headline: **the bacteria/virus phenotype→trait tool now has a VALIDATED virus half.** A
deterministic **HIV-1 drug-resistance decoder** was built, validated against a free independent wet-lab
label (Stanford HIVDB PhenoSense), baseline-compared against Stanford's own regression, refined, subtype-
checked, given its own trust surface, and **wired into the `dna-amr` CLI as a first-class tool.** The
project's central finding held cleanly: *the moment a free, independent label existed, the viral cell
validated — labels were the wall, not models.* Earlier in the same session the **genome-map virulence-
determinant overlay** shipped (the 5th tier) and **Wave A** (cross-kingdom decoder unification) was closed
as already-shipped. **The frozen AMR surface (`dna_decode/eval/amr_rules.py` + `calibrated_amr_rules.json`)
is byte-unchanged across all 15 commits.** Full suite: **1537 passed, 0 regressions** (excl.
`tests/test_models_foundation.py` host-torch limit). The project-level goal is, in my assessment,
**DELIVERED — pending the user's explicit sign-off** (Soraya never auto-declares).

## Verify it's sound (fresh session, ~90s)
```
.venv/Scripts/python.exe -m pytest tests/test_hiv_amr.py tests/test_hiv_nnrti_validate.py \
  tests/test_hiv_nnrti_baseline.py tests/test_hiv_nrti_mutant_catalog.py \
  tests/test_hiv_nrti_within_subtype.py tests/test_build_hiv_report_card.py tests/test_amr_cli_hiv.py -q
# -> 32 passed. (use .venv/Scripts/python.exe -m pytest, NOT uv run — uv cache wedges with os error 183.)
git status --short -- dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json   # EMPTY
```

---

## What shipped this session (the 15 commits, oldest→newest)

| Commit | What |
|---|---|
| `878c5e3` | **Genome-map virulence overlay** — 5th `virulence-determinant` tier (the session's first build; see `executed_plans/Genome_Map_Virulence_Overlay_Plan/` + `wiki/wave_a_cross_kingdom_unification_closeout_2026-06-21.md` is NOT this — that's Wave A). Live-verified on E. coli ST131: 27 virulence features, ExPEC_COMPATIBLE. |
| `bd838ec` | Strategy: tool-completion assessment + Wave-A closeout + viral-label research |
| `a6036ff`→`4f07743` | Wave B grounding: HIV cell `/probe`-by-hand (the circularity-trap catch) → HIVDB license **GO** → dataset confirmed + DRM source (Bennett 2009) |
| `47d8a98` `830208c` | **HIV NNRTI cell built + VALIDATED** (first validated viral cell) |
| `2a5c3db` | **Validate-vs-underlying-tool** — Python reimpl of Stanford `DRMcv.R` OLS baseline |
| `be92cd6` `7f17201` | **NRTI class added + validated**; plateau call (banked at 2 classes) |
| `412eb21` | **NRTI v0.1 mutant-specific** (deconfounded, held-out CV) |
| `fe92df5` | **NRTI within-subtype** transfer check |
| `022871b` | **HIV report card** (own trust surface) |
| `091f375` | **HIV wired into the `dna-amr` CLI** (first-class invokable tool) |

---

## The HIV viral decoder — what it is + how to run it

**Module:** `dna_decode/data/hiv_amr.py` — RT major-DRM catalog (NNRTI mutant-specific + NRTI position-based)
+ `call_from_observed_substitutions` (NNRTI) + `call_nrti_from_observed` (NRTI). Same `HIVCall` shape as
the fungal/antiviral cells.

**CLI (the deployed tool):**
```
dna-amr --drug efavirenz  --observed RT:K103N            # -> R  (amr-mechanism-call-v1 record)
dna-amr --drug zidovudine --observed RT:T215Y --json-only
# NNRTI: efavirenz nevirapine etravirine rilpivirine doravirine
# NRTI:  lamivudine abacavir zidovudine stavudine didanosine tenofovir
# --organism relabels Escherichia -> HIV-1. Genome-FASTA mode is DEFERRED to v0.1 (errors, no fake call).
```

**Validation scripts** (re-run with `.venv/Scripts/python.exe scripts/<name>.py`; need the gitignored data):
- `scripts/hiv_nnrti_validate.py` → `wiki/hiv_nnrti_v0_validation_*.{md,json}`
- `scripts/hiv_nnrti_baseline.py` → the OLS baseline (`DRMcv.R` reimpl)
- `scripts/hiv_nrti_validate.py` → NRTI catalog + baseline
- `scripts/hiv_nrti_mutant_catalog.py` → NRTI v0.1 (deconfounded)
- `scripts/hiv_nrti_within_subtype.py` → subtype transfer
- `scripts/build_hiv_report_card.py` → `wiki/hiv_decoder_report_card.{md,json}` (the roll-up)

## The results (key numbers — do NOT re-derive, they're committed)

| Class | Drug | Headline | Note |
|---|---|---|---|
| NNRTI | efavirenz | **AUC 0.962** | mutant-specific, excellent; **ties** the OLS regression (Δbalacc +0.015) |
| NNRTI | nevirapine | **AUC 0.985** | catalog **beats** the OLS regression (Δ −0.023) |
| NNRTI | etravirine/rilpivirine/doravirine | AUC 0.75/0.70/0.56 | honest class-level degradation on 2nd-gen NNRTIs (K103N spares them) |
| NRTI | (position-based v0) | sens ~0.98-1.0, spec 0.41-0.64 | the T215-revertant over-call |
| NRTI | (deconfounded mutant-specific **v0.1**) | +0.06..+0.14 balacc (held-out) for 5/6 | **ddI regresses** -0.201 (low-signal drug → keep position-based for ddI) |
| NRTI | within-subtype | non-B balacc 0.80-0.93 (transfers) | UNDER-POWERED: data ~96% subtype B |

---

## ⚠️ Design decisions + methodological catches a fresh session MUST NOT undo

1. **Circularity-safe label.** Validate against the **PhenoSense fold-change** (independent wet-lab IC50),
   NEVER HIVDB's own **Sierra/GRT-IS interpretation** (rule-vs-rule = circular). Every script enforces this.
2. **The Y188-not-G188 catch.** A web fetch garbled the consensus-B wild-type at RT-188 to "G188"; the
   Stanford dataset page confirms **Y188C/L/H**. The catalog + a regression test pin Y188. Do not "fix" it.
3. **The deconfounding catch (NRTI v0.1).** A naive "carriers' median fold ≥ cutoff" rule is **CONFOUNDED**
   (revertants ride on resistant lineages). The fix is the **multivariate OLS coefficient** rule (independent
   ≥1.5× effect). A test pins that a co-occurring revertant is excluded. Do not revert to median-fold.
4. **Report-card namespace separation.** HIV has its **OWN** report card (`hiv_decoder_report_card`),
   deliberately NOT folded into the bacterial NCBI-PD provenance-disjoint card — in-distribution ≠
   provenance-disjoint; conflating them is the exact trap the Oxford-arm lesson forbids.
5. **`amr/cli.py` is the dispatcher, NOT the frozen surface.** Adding `_hiv_main` is additive (same pattern
   as `_antiviral_main`); the frozen `amr_rules.py` + `calibrated_amr_rules.json` are byte-unchanged.
6. **Provenance discipline.** The catalog is sourced VERBATIM from the Stanford dataset page (NNRTI majors
   mutant-level; NRTI majors positions-only → position-based v0). NRTI v0.1 mutants are data-derived. No
   fabrication-from-memory. Cite **Rhee et al. 2003 Nucleic Acids Res 31:298-303**.

---

## Data situation (gitignored — re-fetch if missing)

The HIVDB datasets live at `data/raw/hiv/` (gitignored under `/data/*`; **NOT committed** — HIVDB is "© All
Rights Reserved" so research-use download is fine but re-hosting isn't; the project's standard external-data
pattern). Re-fetch:
```
curl -L -o data/raw/hiv/NNRTI_DataSet.txt      https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt
curl -L -o data/raw/hiv/NRTI_DataSet.txt       https://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt
curl -L -o data/raw/hiv/NRTI_DataSet.Full.txt  https://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.Full.txt
# Unfiltered (has Subtype + Type): {CLASS}_DataSet.Full.txt. Filtered: {CLASS}_DataSet.txt.
# PI/INSTI/CAI: PI_DataSet.txt / INI_DataSet.txt / CAI_DataSet.txt (+ .Full variants).
# The Stanford DRMcv.R OLS script is at C:\Users\Farshad\Downloads\DRMcv.R (R not installed; reimplemented in Python).
```
Clinical fold cutoffs (from `DRMcv.R`): NNRTI EFV/NVP/ETR/RPV = 3; NRTI 3TC=5, ABC=2, AZT=3, D4T/DDI/TDF=1.5.

---

## The forward-options menu (NONE auto-started — the user's call)

The project-level goal is DELIVERED (bacteria validated + virus validated). Remaining moves, ranked honestly:

1. **More HIV drug classes (PI / INSTI / CAI)** — **MECHANICAL** repetition of the two proven shapes (PI =
   position-based like NRTI; INSTI = mutant-level like NNRTI, Stanford gave its mutant majors 148H/R/K etc.;
   CAI underpowered n=140). Datasets + cutoffs available. LOW marginal info (the plateau). ~1 increment each,
   full machinery reuse.
2. **HIV genome-mode caller** — accept a novel HIV FASTA (not just `--observed`). Needs a committed
   HXB2/consensus-B RT CDS reference + a BLAST caller, mirroring `scripts/flu_na_caller.py`. MEDIUM effort.
3. **A fresh strategic direction** — a new organism / new phenotype / the parked Evo 2 zero-shot probe
   (`wiki/evo2_zeroshot_vep_lead_2026-06-19.md`, GPU-gated). A new `/idea-anchor` cycle.
4. **Sign-off + bank** — accept completion of the bacteria/virus tool goal.

My recommendation: **option 4 (sign-off)**, or option 2 if a genome-input HIV tool is wanted. Option 1 is
the plateau (I'd only do it on an explicit "do all classes").

---

## Working conventions (this project)
- `.venv/Scripts/python.exe -m pytest` (NOT `uv run` — uv cache wedges, `os error 183`). `PYTHONIOENCODING=utf-8`
  for scripts that print non-ASCII (the cp1252 console else errors on the print; the utf-8 file write is fine).
- Commit straight to `main` (= the cross-machine sync channel; user syncs ~weekly). Footer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- **Leave dirty, NOT mine:** `uv.lock`, `wiki/ciprofloxacin_mechanism_audit_2026-06-05.{json,md}`,
  `bash.exe.stackdump`, the stale `plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md`. Stage feature files only.
- Project ledger: `project_state/dna-decode-2026-05-11.md` (Action Log up to **row 146**; rows 132-146 are
  this session). The genome-map virulence overlay's own handoff is
  `wiki/SESSION_HANDOFF_2026-06-21_virulence_overlay_commit.md` (now consumed).

## Parked / not-this-thread
- The learned/embedding decoder is a CLOSED NEGATIVE on free data (4 de-confounded failures). Do NOT reopen
  without a label-clean substrate. (HIV did NOT reopen it — HIV is a DETERMINISTIC catalog, not a learned model.)
- Public-label bacterial AMR expansion is banked (`wiki/negative_results_map_2026-06-13.md`, the 8 gates).
- TB RIF/INH on CRyPTIC: v1b SCORED run still BLOCKED on the ~1.6 TB regeno fetch + a hand-curated gold set.
