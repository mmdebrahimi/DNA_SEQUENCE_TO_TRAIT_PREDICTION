> **⚠️ SUPERSEDED 2026-06-17 (later same-day session).** This is an EARLIER snapshot. A subsequent session
> (a) ratified the independent post-2023 gold set **IN-SCOPE NOW** (no longer "deferred" — see line ~20),
> and (b) ran the pre-exec `/brainstorm`, which **overturned three technical choices below**. Do NOT drive
> the build from this doc — it would re-introduce caught bugs. Corrections:
> - **Ratification C (VCF call policy):** there is **NO `GCP` FORMAT key** in CRyPTIC VCFs; the floor is the
>   `MIN_GCP`/`MIN_DP`/`MIN_FRS` FILTERs, so `FILTER==PASS` + `GT` non-reference is the rule (no extra floor). [verified]
> - **Stage 3 aggregation:** drop "representative-isolate dedup"; **reuse `clonality.cluster_weighted_confusion`**
>   (already excludes/counts DISCORDANT clones).
> - **Stage 3 distance:** masked-VCF union-site SNP distance is untrustworthy for TB lineage → use a
>   **deterministic lineage-barcode caller** (Coll/Napier SNPs) as the v1b collapse key; SNP-distance/REGENOTYPED path → v1c, BLOCKED until validated.
> - **Data state advanced:** WHO catalogue **pinned** (commit 0bb39143 + CHECKSUMS), coordinate alignment **VERIFIED**,
>   **300** RIF VCFs staged (the 150R/150S cache is plumbing substrate, NOT a v1b cohort).
> Current source of truth: the regenerated plan (re-run `/technical-plan tb-decoder-cryptic`) + ledger row 105
> + `wiki/tb_coordinate_alignment_probe_2026-06-16.md`. A/B/E ratifications still pending.

# dna_decode — Session Handoff (2026-06-17): TB AMR decoder PLAN ready, awaiting ratification

Pick-up doc for a fresh session. Driven by Soraya (attended). Goal this session: point Soraya at dna_decode
and set up the next bounded mission. Outcome: the TB-decoder technical plan is **written + brainstorm-hardened**,
**not yet executed** — it's gated on 4 domain ratifications.

## ⇒ START HERE
1. `cd C:\Users\Farshad\PythonProjects\dna_decode`
2. Read the plan: **`plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md`** (status: candidate; pre-exec /brainstorm incorporated).
3. **Ratify the 4 open tradeoffs** (table below) — a simple "defaults are fine" unblocks the build.
4. Then: add a `wiki/plans-index.md` entry (the /save-plan step) → drive the build with
   **`/soraya --until-mvp`** against a TB-decoder family ledger whose MVP bar = the plan's v1a + v1b acceptance predicates.
5. **Do NOT start coding before ratification** — the buried substrate-build was the real risk; the plan now surfaces it.

## What this goal is
The first **M. tuberculosis AMR decoder cell** — **RIF (`rpoB` RRDR) + INH (`katG` + `inhA` promoter)** —
deterministic WHO-catalogue determinant rules, scored over the **CRyPTIC reuse cohort** (~12,287 isolates),
**lineage-collapsed** metrics, honestly labeled **`WHO_CATALOGUE_ON_CRYPTIC_KNOWLEDGE_BASELINE`** (the WHO
catalogue was partly built from CRyPTIC → knowledge-baseline, NOT independent validation). Independent
validation (post-2023 / gold set) is the **deferred** next milestone.

This aligns with the prior dna_decode direction (2026-06-16 strategic fork: acquire new label sources; TB is
the first concrete one). The public-label AMR track is DONE + frozen (2026-06-13); binding constraint = labels.

## The staged build (why staged — pre-exec /brainstorm findings)
The naive "write rules + score the local CSV" framing is WRONG: the CRyPTIC CSV is only phenotypes/MICs/IDs
(no variant data); `data/raw/cryptic/vcf_cache/` has only **30** masked VCFs (RIF PoC), not 12,287; and there's
**no** lineage-collapse substrate (`eval/clonality.py` builds its matrix via Mash-on-FASTA through
`eval/phylogeny.py`; TB has only masked VCFs, no FASTAs). So four gated stages:
- **Stage 0** — VCF acquisition + a **genotype/filter/callability** parser (`GT`/FILTER-aware; "any-ALT=present" mis-scores; handle `0/0`-with-ALT, `./.`, multi-allelic, masked, MNV/indel normalization for `inhA`/`katG`). Unit-tested.
- **Stage 1** — WHO catalogue (`GTB-tbsequencing/mutation-catalogue-2023`) **join + pin by commit SHA + per-file checksums** (NO releases exist; files changed Feb/May 2024). Excel-grades ↔ VCF-coords JOIN, sentinel-tested.
- **Stage 2 — v1a "plumbing cell"** — 20–50 sentinel VCFs (known determinants + negatives), end-to-end scoring, **NO metric claims**. Acceptance: `file-exists organism_rules/tb_amr.py` · `test-exit-0` cell tests · plumbing artifact.
- **Stage 3 — v1b** — direct VCF→SNP distances → **representative-isolate dedup** (phenotype-blind; NOT R/S majority-vote — intra-cluster discordance may be real) → lineage-collapsed sens/spec, **`BLOCKED`-gated** (`LINEAGE_COLLAPSE_BLOCKED_NO_DISTANCE_MATRIX` if no provenance-stamped matrix; never fake metrics). Prefer sparse threshold-neighbor clustering over a dense 12k×12k (~1.2 GB) matrix.
- **Leak test** — assert no CRyPTIC phenotype columns read during rule construction + catalogue checksum-pinned + frozen E. coli surface byte-untouched. (Cannot prove biological independence — that's the deferred follow-up.)

## ⚠️ 4 ratification points (domain-science calls — PENDING)
| # | Decision | Drafted default |
|---|---|---|
| A | INH scope | `katG` + `inhA`-promoter grade-1/2 only, coverage-reported (honest partial) |
| B | v1 cohort cut | bounded VCF subset for v1b, scale later (not all 12,287 up front) |
| C | VCF call policy | `FILTER==PASS` + `GT` non-reference + a min-quality floor (ratify FRS/DP/GCP values) |
| D | Canonical VCF + SNP threshold | masked VCF (matches cache); threshold ~5–12 SNPs + sensitivity check |

## Hard boundaries
- **FROZEN — never edit:** `dna_decode/eval/amr_rules.py`, `dna_decode/data/calibrated_amr_rules.json` (freeze 2026-06-13).
- New code only: `dna_decode/organism_rules/tb_amr.py` (+ `tb_vcf.py`, `tb_snp_distance.py`, `data/tb_who_catalogue.py`) + tests. Mirror the TMP-SMX overlay pattern (`dna_decode/data/experimental_drug_rules.py`).
- Deterministic rules only (no learned/embedding model). Frozen E. coli test suite must stay green (0 regressions).

## Data state
- `data/raw/cryptic/CRyPTIC_reuse_table_20240917.csv` — 12,288 rows, phenotypes/MICs/quality + IDs + VCF paths only.
- `data/raw/cryptic/vcf_cache/` — 30 masked VCFs (PoC subset). Full cohort NOT fetched.
- WHO catalogue — NOT fetched locally yet (Stage 1 fetches + pins it).
- Reuse: `dna_decode/eval/clonality.py:greedy_representative_clusters_from_matrix` (pure-logic dedup, reusable).

## Pointers
- Plan: `plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md`
- Prior probe/design: `wiki/cryptic_feasibility_probe_result_2026-06-16.md`, `wiki/tb_decoder_idea_anchor_prompt_2026-06-16.md`, `wiki/SESSION_HANDOFF_2026-06-16.md`
- Stale Bellman ledger (reference only): `project_state/dna-decode-2026-05-11.md` (true state lives in wiki).

## Not-yet-committed
The plan file + this handoff are **uncommitted** in the dna_decode repo (`mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`).
Commit them before/at session start if you want them durable.
