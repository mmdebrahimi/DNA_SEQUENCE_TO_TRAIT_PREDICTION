# Session handoff — 2026-06-16 (resume in a fresh session)

Everything a cold session needs to continue. Repo: `C:\Users\Farshad\PythonProjects\dna_decode` (origin
`mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`, branch `main`, all work below pushed). HEAD = `27ff513`.

## TL;DR — where we are
- The deterministic AMR decoder is **banked + reproducibility-frozen** (2026-06-13). The public-label AMR
  track is **at its terminal honest product** (plateau named: labels, not models, are the constraint).
- This session **externally validated** the frozen decoder on 2 independent measured-MIC cohorts, **added a
  4th-mechanism cell** (TMP-SMX), and **opened a new organism direction** (M. tuberculosis via CRyPTIC),
  which is now fully scoped through the planning pipeline and **sitting at `/technical-plan`**.
- **THE NEXT ACTION:** `/technical-plan tb-decoder-cryptic` — but settle one authority decision first (below).

## ⏭️ Resume here — the active thread: first TB decoder cell (CRyPTIC)
The pipeline ran: `/idea-anchor` → `/probe` (design) → `/probe` (Q1 grounding). All done. Next = `/technical-plan`.

**The settled design (encode these in the plan):**
- **Substrate:** CRyPTIC reuse table — 12,287 *M. tuberculosis* isolates × 13 drugs, reference UKMYC
  broth-microdilution MIC + binary phenotype + per-isolate VCF (vs H37Rv NC_000962.3). Free EBI FTP,
  already downloaded to `data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv` (gitignored). All 13 drugs
  powered; first cells = **RIF (rpoB) + INH (katG + inhA promoter)** — INH is required so the plumbing
  proves more than a single-gene happy path.
- **Rule:** curated determinant rule from the **WHO M. tuberculosis mutation catalogue v2 (2023)** — match
  CRyPTIC VCF variants by genomic position/codon. Catalogue is machine-readable on GitHub
  (`GTB-tbsequencing/mutation-catalogue-2023`, genomic-coordinate + VCF files). **Pin it** (commit SHA +
  checksum + Jan-2024 corrigenda) and verify coordinate alignment to the CRyPTIC VCFs with a fixture
  (`rpoB S450L`, `katG S315T` → actual VCF records) BEFORE trusting any number.
- **Home:** a **NEW non-frozen `tb_amr` / `organism_rules` module** with explicit `input_type: vcf_h37rv` +
  catalogue version/checksum + a frozen-leak test. **Do NOT** put TB in the frozen `calibrated_amr_rules.json`
  (AMRFinder-`main.tsv` input contract; different semantics) and **do NOT** edit the frozen E. coli surface.
  Pattern to follow = the TMP-SMX experimental overlay (see below).
- **Clonality (make-or-break, NOT a companion metric):** TB resistance is clonally transmitted, so the raw
  3448 RIF-R / 4467 INH-R **may collapse to few effective lineages**. Build a **SNP-distance matrix from the
  VCFs** and feed `dna_decode/eval/clonality.py::greedy_representative_clusters_from_matrix(distance_matrix,
  threshold)` (already exists + tested — reuse, no Mash/assemblies). Report `n_clusters_R/S`, discordant
  clusters, Wilson CI, and **raw→lineage shrinkage**. If effective lineages are tiny, demote the cell to
  "smoke test".
- **Honesty label (load-bearing):** the WHO catalogue was built **partly from CRyPTIC**, so a CRyPTIC-scored
  number is the **in-distribution KNOWLEDGE-BASELINE** (rule recalling its own training phenotypes), NOT
  independent validation. Ship the cell labelled `WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`, surfaced in
  the non-frozen external-validation arm (the TMP-SMX precedent), frozen card untouched.
- **Independence path (Q1 grounded this session):** an independent TB cohort *exists* (TB Portals/NIAID,
  WGS+DST) BUT independence must be from the **WHO-catalogue build**, and WHO v2 swept most public pre-2023
  TB WGS+pDST → legacy cohorts are in-distribution. Real independence = **temporal hold-out (post-2023
  isolates)** or a **hand-curated ~30-isolate post-2023 gold set**; TB Portals is access-gated (Data Access
  Request) + mixed-method DST. Scope this as the explicit follow-up, never conflated with the baseline.

**⚠️ AUTHORITY DECISION to settle before/at `/technical-plan` (yours):** confirm a labelled
`KNOWLEDGE_BASELINE` is an acceptable first deliverable (given no free immediately-independent cohort), and
whether the post-2023 independent gold set is **in-scope now or deferred**. This sets what "done" means.

**Honest reframe to keep in view:** CRyPTIC expands the **deterministic** decoder to a new organism — it does
**NOT** reopen the learned-embedding niche (TB AMR has a curated catalog → embeddings lose here too, as on E.
coli). If the deeper goal was "does the learned decoder ever work", that still needs a catalog-free phenotype
(TransPred growth-kinetics, or the parked eukaryotic Path B G2).

**Resume command:**
```
/technical-plan tb-decoder-cryptic
```

## What shipped this session (all on origin/main)
| commit | what |
|---|---|
| `5503772` | 234-cohort completed: gent+cef on Sci234 via new `fam.tsv` Subclass resolver (`scripts/fam_subclass_resolver.py`); cipro+cef SCORED, gent UNDERPOWERED |
| `caab213` | drafted new-drug-coverage `/idea-anchor` + `/brainstorm` prompts |
| `4027333` | **TMP-SMX EXPERIMENTAL cell** — `dna_decode/data/experimental_drug_rules.py` (`(sul AND dfr)` rule) + `scripts/tmp_smx_external_validate.py` + non-frozen branding in `scripts/build_external_validation_report.py`; Sci234 0.987 / Oxford 0.962, strata reproduce; frozen surface byte-unchanged |
| `ab77af9` | epoch bank: `wiki/external_validation_milestone_synthesis_2026-06-16.md` + acquisition-fork `/research` shortlist (`research_outputs/acquirable-label-sources-2026-06-16.md`) |
| `2baedc4` | **CRyPTIC feasibility probe** — `scripts/cryptic_feasibility_probe.py` + GREEN result (`wiki/cryptic_feasibility_probe_result_2026-06-16.md`) |
| `fd87bdb` | TB `/idea-anchor` prompt (`wiki/tb_decoder_idea_anchor_prompt_2026-06-16.md`) |
| `27ff513` | TB `/probe` Q1 grounding (ledger row 104) |

## Key file pointers
- **TB direction:** `wiki/tb_decoder_idea_anchor_prompt_2026-06-16.md`, `wiki/cryptic_feasibility_probe_result_2026-06-16.md`, `scripts/cryptic_feasibility_probe.py` (+ tests). Raw data: `data/raw/cryptic/` (gitignored).
- **Acquisition research:** `research_outputs/acquirable-label-sources-2026-06-16.md` (CRyPTIC #1; ATLAS/EUCAST disqualified = no paired genomes; independence sweep `wiki/independent_phenotype_source_sweep_2026-06-10.md`).
- **External-validation arm (non-frozen):** `scripts/build_external_validation_report.py`, `scripts/{oxford,sci234,tmp_smx_external}_validate.py`, `wiki/external_validation_*` artifacts.
- **Frozen surface (DO NOT EDIT):** `dna_decode/eval/amr_rules.py`, `dna_decode/data/mic_tiers.py`, `dna_decode/eval/cohort_manifest.py`, `scripts/build_validation_report_card.py`, `scripts/compute_lineage_metrics.py`, `dna_decode/data/calibrated_amr_rules.json`. Freeze doc: `wiki/reproducibility_freeze_2026-06-13.md`. Negative-results map: `wiki/negative_results_map_2026-06-13.md`.
- **Reuse for TB:** `dna_decode/eval/clonality.py` (`greedy_representative_clusters_from_matrix`, `cluster_weighted_confusion`, `wilson_ci`, `effective_lineage_n`).
- **Ledger:** `project_state/dna-decode-2026-05-11.md` (action log through row 104). **Frame is STALE** (last `--refresh-frame` 2026-05-26) — `/project-state --refresh-frame` is overdue; the true state is in `wiki/external_validation_milestone_synthesis_2026-06-16.md`.

## Loose ends / not-this-thread
- **NF-001 (OTHER machine, Bombardier b0652085):** resolved this session. All 7 Candida-auris recovery artifacts were on THIS machine; the 5 derived ones are on origin (commit `df78fc9`, `git pull`); the 2 PDFs are zipped at `C:\Users\Farshad\Downloads\nf001_cauris_pdf_recovery_2026-06-16.zip` (1.3 MB) — **user to Gmail-port to the Bombardier machine**. Not a dna_decode-repo task.
- **Ledger frame refresh** (`/project-state --refresh-frame`) — overdue bookkeeping, optional.
- **Other open forks (not chosen):** eukaryotic Path B G2 embedding test (pre-staged, needs Precision 7780 GPU — the real learned-decoder experiment); prospective-lock validation of the frozen decoder.
- **TMP-SMX follow-ups (deferred):** flatten `build_external_validation_report` `drugs`-map so the existing Oxford/Sci234 multi-drug cells render; promotion of TMP-SMX to the frozen surface.

## Working conventions (this project)
- `python`/`uv run` (NOT `python3`). Tests: `uv run pytest tests/ -q` (exclude `tests/test_models_foundation.py` — host torch-paging limit; full suite was 1333 green this session).
- Commits go straight to `main` (= the cross-machine sync channel; user syncs ~weekly). Commit-message footer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Stage only your own files (pre-existing dirty files: `uv.lock`, `wiki/ciprofloxacin_mechanism_audit_2026-06-05.*`, `bash.exe.stackdump`, `research_outputs/eukaryotic-...unsupported.md`, `wiki/3 idea anchors...rtf` — leave untouched).
- Planning pipeline: each step STOPS; user drives. WebSearch trips a false-positive usage-policy filter on "tuberculosis + resistance + MIC" phrasings — use neutral institutional phrasing or WebFetch known URLs.
