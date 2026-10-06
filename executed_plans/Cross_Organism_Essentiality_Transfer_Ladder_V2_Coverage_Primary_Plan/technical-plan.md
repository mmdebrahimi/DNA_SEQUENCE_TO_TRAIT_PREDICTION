# Cross-Organism Essentiality Transfer Ladder — Technical Plan (v2, coverage-primary)

## Lens status
Dispatch: conversation-driven. **This is a v2 re-plan, not a fresh plan** — v1 was drafted in conversation, adversarially reviewed (deep, 2 rounds), and the review's 3 critical issues plus a follow-up measurement changed the metric, deleted a step, and added a step. v1 was deliberately never saved (the save-once rule: don't persist a plan a review just mutated).
Project root: `C:/Users/Farshad/PythonProjects/dna_decode`.
Optional context: `wiki/` present; `DESIGN.md` / `AUDIT_REPORT.md` / `project-rules.md` ABSENT → no Project Rules Check section.
Repo index: absent (`update-index` CLI not installed) → direct file reading, no degradation.
Tools: `gh` + `codex` present; `sentrux` ABSENT → architecture gate `n/a — sentrux not installed`.
Wiki replay: `wiki/plans-index.md` read. Nearest prior entries `Dna_Identify_Closed_Set_Organism_Router_Plan/` and `Meropenem_FN16_Species_Composition_Audit_Plan/` (both executed). Three of their ratified decisions are adopted rather than re-litigated: **thresholds frozen only after the distribution analysis**, **disclosure is namespace-separate + augment-only verified by diff**, and **reconcile before you attribute**.

## Problem Statement

The project's north star ends at *climb to higher life forms*, and **the entire empirical basis for that step is one measurement**: the E. coli-tuned conserved-core essentiality decoder applied UNCHANGED to human (BAGEL CEGv2 vs NEGv1) at AUROC 0.5805 against a 0.5 null. One point at the maximum available phylogenetic distance cannot separate **(A)** continuous decay with distance from **(B)** a cliff at the eukaryote boundary, and those imply different work: under (A) climbing is a gradient and per-rung catalogue extension buys proportional recall; under (B) there is one nameable qualitative barrier and the extension is a finite curation job.

**What changed since v1, and why this is a different plan rather than an edited one:**

1. **AUROC is the WRONG primary metric, measured.** The decay is not a discrimination failure, it is **silence**: **83.1% of human essentials score exactly 0** (E. coli 62.1%), and in human **11 of 13 core patterns have P(hit | non-essential) = 0.0000** — the patterns never misfire, they never fire. AUROC on a distribution that is 83% tied is mostly a tie-mass statistic. Worse, its two components move in **opposite directions** across the two existing rungs, so AUROC hides both:

   | rung | coverage(essential) | coverage(non-essential) | **coverage_lift** | precision where it fires |
   |---|---|---|---|---|
   | E. coli (n=351/3432) | 0.3789 | 0.0664 | **0.3125** | 0.3684 |
   | human (n=681/899) | 0.1689 | 0.0056 | **0.1633** | 0.9583 |

   **`coverage_lift = coverage(essential) − coverage(non-essential)` is the new primary.** Both terms are conditioned on the class, so it is base-rate robust and therefore *genuinely cross-rung comparable* — which AUROC is not (spectrum bias: BAGEL is two curated extremes at base rate 0.431; E. coli is genome-wide at 0.0928) and sens/spec/precision are not. It has an obvious permutation null (shuffle labels → expected lift 0).

2. **A cheap pre-test was run and it does NOT pre-empt the ladder — but it bounds it.** `wiki/essentiality_missed_vocabulary_2026-10-05.{md,json}` asked whether the human miss is really a catalogue *phrasing* gap, which would make a ladder the wrong instrument. Verdict `INDETERMINATE_BAR_SENSITIVE_VOCABULARY_CANNOT_ADJUDICATE`, with two floors that survive every bar: **mechanism (iii) host-specific 141/566 = 0.2491** vs **mechanism (ii) phrasing 57/566 = 0.1007**. Host-specific absence is ~2.5× the demonstrable phrasing gap, so the ladder measures something real — **but ~10% of the human miss has nothing to do with phylogeny and must be subtracted before any rung is scored, or it is silently attributed to distance or to the domain boundary.**

3. **One v1 step is DONE and is removed from scope.** v1 Step 5 (the `0.911` within-human-CV over-claim in three shipped surfaces) shipped as commit `c8db73b`, with `tests/test_essentiality_overclaim.py` guarding the prose relationship.

4. **One v1 step is DELETED on evidence.** v1 Step 4 proposed a functional-class vocabulary with per-rung synonyms. The review flagged it as a second adjudication system and a post-hoc rescue table; the measurement then showed that vocabulary overlap **cannot adjudicate the middle at all** (three defensible bars, two verdict bands). So the broad class audit is dropped and replaced by a narrow, provenance-gated **phrasing-floor subtraction** using only the 6 words already verified to appear in the decoder's own regexes.

5. **A W0 schema probe is ADDED.** v1 wrote parsers from prose for files not present on disk. That failed twice on me inside one session: a candidate second-annotation source died on inspection (`uniprot_ecoli_p1.tsv` is 501 rows of `Disruption phenotype`, no product column), and my first Goodall parse silently mis-joined to 3936/372 against the committed 351/3432 because row 0 of that xlsx is a title.

**Non-goals.** Retuning the decoder per organism (that destroys the transfer test). Extending `_CORE` with eukaryote-specific classes (the named follow-on; it would end the experiment). Any learned/supervised model. Touching the five sha256-pinned frozen AMR files. Re-scoring any existing cell. Adjudicating the vocabulary middle — measured as not adjudicable.

## Codebase Context

**The decoder under test** — `dna_decode/essentiality/core_decoder.py` (69 lines, stdlib + `re`). `score_gene(gene, product, threshold=2.0)` regex-matches **product text** against a 13-entry weighted universal-core catalogue (`_CORE`), minus 2.0 on a negative-signal match. Organism-agnostic by construction (inputs are a gene symbol and a product string), which is exactly why "apply unchanged" is a coherent experimental condition and why a per-organism synonym layer would break it.

**THE ENABLING FINDING (unchanged from v1, still a hypothesis until fetched).** `dna_decode/fba/essentiality_labels.py` records S. aureus as `LABEL_WALLED` ("iYS854 ids are `USA300HOU_####`; NTML is USA300 JE2 (`SAUSA300_####`) -> crosswalk") and P. aeruginosa as `MODEL_WALLED` ("no P. aeruginosa GEM in BiGG … the gold standard DOES exist and is fetchable (PLOS Comput Biol 2026 `pcbi.1013945.s011`, sheets GOLD_84 / GOLD_115) but is keyed to **PA14** locus tags — so the blocker is the MODEL, not the label"). **Both are MODEL-GENE-JOIN walls.** This decoder joins a label to a genome **annotation** (`locus_tag -> product`), never to a metabolic model, so two rungs the FBA arm files as blocked are reachable here. The review correctly noted the wall does not vanish, it **moves** to "label ids must join the exact annotation table" — which is what Steps 2 and 3 exist to verify and gate.

**Rung ordering** must be DERIVED, not asserted: E. coli and P. aeruginosa are both Gammaproteobacteria while S. aureus is a different phylum, so the distance order is not the order a human would guess from familiarity. Step 1 derives it from committed NCBI Taxonomy lineage strings, each with a `source_url`.

**Reusable, confirmed present (verified this session):**
- `scripts/essentiality_missed_vocabulary.py` — `human_rows()` (BAGEL CEG/NEG → `(symbol, description)`), `ecoli_rows()` (NCBI feature table → `(symbol, product)`), `content_words()`, `classify_word()`, `STOP`, `ROBUST_II`. **All four loaders/helpers the ladder needs already exist and are tested.**
- `scripts/essentiality_e3_learned.py::load_labels` — the **correct** Goodall xlsx parser (header row 1, skips `Unclear` and ambiguous rows). Reuse; do not write a second one.
- `dna_decode/essentiality/cli.py::_iter_feature_table` — NCBI feature-table reader; promote to the shared module.
- `scripts/build_essentiality_report_card.py` — the read-only roll-up (exit 0 always, per-organism tier, **no aggregate headline**); already computes the human rung.
- `scripts/oxford_w0_probe.py` — the W0 pattern verbatim: *"pin the … schema + crosswalk feasibility BEFORE building the ingester against it (the 'build against real data' gate)"*, with PURE summarize helpers split from the fetch so they unit-test offline.
- `dna_decode/eval/regime.py` — `ORGANISM_TRANSFER_LEVELS`, `organism_transfer_is_unmeasured()` (returns all 8 regimes), and `_transfer_gap_conditions`' string "no regime in this map carries one today".
- `dna_decode/data/ecoff_catalog.py` — the provenance-gate precedent (`source_url` + `verbatim_quote` or raise).

**On disk, verified:** `D:/dna_decode_cache/essentiality/` holds `CEGv2.txt`, `NEGv1.txt`, `Homo_sapiens.gene_info.gz`, `goodall_TableS1_essential.xlsx`, `ecoli_k12_feature_table.txt.gz`. **Absent:** SGD `phenotype_data.tab`, any NTML file, the PLOS supplementary.

**Recorded negative that bounds reachability:** `wiki/essentiality_label_wall_2026-07-28.md` documents SEVEN dead routes for essentiality labels and says "do NOT re-hunt without a new lead". NTML and the PLOS supplementary are not in that table, so they are new leads — but a recorded "is fetchable" is a hypothesis until fetched.

### Reusable-Code Survey

Candidates found and adopted (searched: `dna_decode/essentiality/`, `dna_decode/fba/`, `dna_decode/eval/`, `scripts/` [grep `essential|transfer|w0_probe|ladder`], `wiki/` [`ls | grep essential`], `tests/` [`ls | grep essential` → 8 files], `D:/dna_decode_cache/essentiality/`; `graphify-out/GRAPH_REPORT.md` absent; no `src/`, `lib/`, `utils/`, `common/`, `helpers/` dirs exist in this repo):

1. `essentiality_missed_vocabulary.{human_rows,ecoli_rows,content_words,classify_word,STOP,ROBUST_II}` — the loaders + the phrasing-floor vocabulary. Adopted by import; promoted into the package in Step 4 so a `scripts/` module is not a library dependency.
2. `essentiality_e3_learned.load_labels` — the correct Goodall parser. Adopted; Step 5 must not reimplement it.
3. `essentiality/cli.py::_iter_feature_table` — feature-table reader, promoted to `annotation_join` (Step 3).
4. `build_essentiality_report_card.human_transfer_auroc` — its Mann-Whitney AUROC becomes the SECONDARY metric, computed once in a shared place so the card and the ladder cannot drift.
5. `oxford_w0_probe.py` — the structural template for Step 2 (pure summarizers split from fetch).
6. `ecoff_catalog.py` — the provenance-gate shape reused for lineage strings (Step 1) and the phrasing-floor vocabulary (Step 4).
7. `compute_lineage_metrics.reconcile_raw_metrics` / `ReconcileMismatch` — the reconcile-before-attribute gate reused in Step 6.

## Pre-Change Baseline

- **Test suite:** `uv run pytest tests/ -q --collect-only` → **5777 collected**, 0 collection errors. Last full run: **5773 passed / 12 skipped / 0 failed** (2026-10-05).
- **The two existing rungs, on the NEW primary metric** (computed this session, reproducing the committed 351/3432 join exactly): E. coli `coverage_lift` **0.3125** (cov_ess 0.3789, cov_non 0.0664); human **0.1633** (cov_ess 0.1689, cov_non 0.0056). Spread **0.1492**.
- **On the old secondary metric:** E. coli AUROC 0.6952, human 0.5805 (spread 0.1147).
- **Phrasing floor already measured:** mechanism (ii) 57/566 = 0.1007; mechanism (iii) 141/566 = 0.2491 (`wiki/essentiality_missed_vocabulary_2026-10-05.json`).
- **Regime map:** `organism_transfer_is_unmeasured()` returns **8 of 8** keys.
- **Report card:** `wiki/essentiality_report_card.json` carries **2 organism rows**, schema `essentiality-report-card-v1`, no aggregate headline.
- **Frozen five** byte-unchanged and `verify_lock` reports `ok=True` against the 2026-08-31 lock. This plan touches none of them.

## Verification Signal

1. **A ladder artifact** at `wiki/essentiality_transfer_ladder_<date>.{md,json}` carrying, per rung: organism, NCBI taxid, derived shared-rank depth from E. coli, label source + URL + technology, `n_essential` / `n_nonessential` / base rate, **join rate AND joined-vs-unjoined composition**, **`coverage_lift` with its permutation null (≥200 draws, mean + p95)**, `coverage(essential)`, `coverage(non-essential)`, `precision_where_fires`, the **phrasing-floor-subtracted `coverage_lift_adjusted`**, and AUROC as a clearly-labelled secondary.
2. **Every declared rung accounted for** — a rung that cannot be scored emits a NAMED wall (`WALL_SCHEMA_UNVERIFIED` / `WALL_LABEL_UNREACHABLE` / `WALL_ANNOTATION_UNREACHABLE` / `WALL_JOIN_RATE_BELOW_FLOOR`) with the route tried. A test asserts the artifact's rung set equals the declared rung set.
3. **The two existing rungs reproduce**: the runner recomputes E. coli and human and matches `coverage_lift` 0.3125 / 0.1633 to 4 dp **and** the legacy AUROC 0.6952 / 0.5805 to 3 dp, **raising** on mismatch before anything downstream is trusted.
4. **A verdict from a bar frozen before any new rung lands**, naming one of: `SMOOTH_DECAY_WITH_DISTANCE`, `DOMAIN_CLIFF`, `NO_DECAY_DISTANCE_INSENSITIVE`, `INDETERMINATE_INSUFFICIENT_RUNGS`, `INDETERMINATE_NULL_NOT_CLEARED`, `INDETERMINATE_BAR_SENSITIVE`.
5. **The decoder is byte-unchanged** — a test asserts `sha256(core_decoder.py)`; a retuned decoder voids the whole design.
6. **No W0-unverified parser runs** — a test asserts every label source consumed by the runner has a `w0_verified` schema record, and that the runner refuses otherwise.
7. **Report-card extension is augment-only, verified by diff** — existing 2 rows keep every field value.
8. **Suite green:** ≥ 5777 + new tests, **0 failed**; frozen five byte-unchanged; `verify_lock ok=True`.

## Implementation Steps

### Step 1: Freeze the design with coverage_lift as primary and thresholds as parameters
Files: dna_decode/eval/transfer_ladder.py, tests/test_transfer_ladder.py
Depends on: none

**What changes:**
- `PREREGISTERED` dict naming the question, the primary vs secondary metric, all six verdict names, and the falsification condition — mirroring `transfer_benchmark_falsifier.PREREGISTERED`.
- **`coverage_lift(scores_essential, scores_nonessential) -> float`** = frac non-zero among essentials minus frac non-zero among non-essentials. Pure. Docstring states WHY it replaces AUROC: 83.1% of human essentials are tied at 0 so AUROC is a tie-mass statistic, and the two components move in opposite directions (0.3789→0.1689 vs 0.3684→0.9583) so AUROC hides both. Both terms are class-conditioned ⇒ base-rate robust ⇒ cross-rung comparable, which AUROC is not.
- `RUNGS` as data: organism key, NCBI taxid, label-source id, label **technology** (`transposon_insertion` / `deletion_collection` / `crispr_ko`), and a `lineage` string with `source_url`. Construction RAISES without both, mirroring `ecoff_catalog`.
- `shared_rank_depth(lineage_a, lineage_b) -> int` — pure; DERIVES the ladder ordering instead of asserting it (E. coli and P. aeruginosa are both Gammaproteobacteria while S. aureus is a different phylum, so the order is not the intuitive one).
- `PRIMARY_TEST = "bacterial_technology_matched_subladder"` — the 3 bacterial rungs are all transposon-insertion screens, so distance varies with label technology held constant. The full 5-rung ladder is SECONDARY and labelled technology-confounded AND sampling-frame-confounded.
- `classify_ladder(rungs, *, cliff_drop, plateau_tol, min_rungs)` — pure; every threshold a PARAMETER, deliberately unfrozen here. Step 6 freezes them.

**Test strategy:**
- Unit: `coverage_lift` on synthetic score vectors incl. all-tied-at-zero (lift 0), perfect separation (lift 1.0), and the real two-rung values reproducing 0.3125 / 0.1633.
- Unit: `coverage_lift` is INVARIANT to class-size reweighting (duplicate the non-essential set; lift unchanged) — this is the base-rate-robustness claim, pinned rather than asserted.
- Unit: `shared_rank_depth` on synthetic lineages; assert the derived order of the committed lineages puts P. aeruginosa NEARER E. coli than S. aureus (the counter-intuitive case).
- Unit: a `RUNGS` entry without `lineage` + `source_url` RAISES.
- Unit: `PREREGISTERED` names every verdict `classify_ladder` can return.

### Step 2: W0 schema probe — pin every label file's real schema before any parser exists
Files: scripts/essentiality_ladder_w0_probe.py, tests/test_essentiality_ladder_w0_probe.py
Depends on: none

**What changes:**
- New probe modelled structurally on `scripts/oxford_w0_probe.py`: PURE summarize helpers split from fetch/read so they unit-test offline.
- Per declared label source, records to `wiki/essentiality_ladder_w0_probe_<date>.json`: resolved URL + fetch outcome, file kind, **sheet names and the actual header ROW INDEX** (not assumed 0), the key column, key-space shape with representative ids, row counts, label polarity/columns, duplicate-key count, and the strain/reference the keys belong to (PA14 vs PAO1; USA300 JE2/FPR3757 vs TCH1516/HOU; S288C).
- Emits `w0_verified: true|false` + a reason per source. **This record is the gate Step 5's parsers and Step 6's runner read.**
- Two failures from this session are encoded as explicit probe checks, because both already happened to me: a candidate source that is the wrong shape entirely (`uniprot_ecoli_p1.tsv` — 501 rows, `Disruption phenotype`, no product column), and a header that is NOT row 0 (the Goodall xlsx row 0 is a title; reading it as the header silently mis-joins to 3936/372 against the true 351/3432).

**Test strategy:**
- Unit (offline): each pure summarizer against small synthetic fixtures — an xlsx-shaped table whose header is row 1, a tab file, and a header-on-row-0 table; assert the detected header index is right in all three.
- Unit: a source whose fetch fails yields `w0_verified: false` + a reason and does NOT abort the other sources.
- Unit: the Goodall-shaped fixture with a title row is detected as header-row-1; a test asserts that reading row 0 as header would have produced a DIFFERENT label count (the mis-join reproduced, so the check is non-vacuous).
- Live (manual, network): run against the real sources; its artifact is Step 5's input.

### Step 3: Annotation join, coverage math, and a gated, composition-reported join rate
Files: dna_decode/essentiality/annotation_join.py, dna_decode/essentiality/cli.py, tests/test_essentiality_annotation_join.py
Depends on: none

**What changes:**
- `load_annotation(path, kind) -> dict[str, str]` for `ncbi_feature_table` / `sgd_features` / `ncbi_gene_info`. The feature-table reader is PROMOTED from `essentiality/cli.py::_iter_feature_table` (single definition; the CLI imports it).
- `join_labels(annotation, essential_ids, nonessential_ids) -> JoinResult` returning scored rows plus `join_rate`, `n_unjoined`, sample unjoined ids, **and `joined_vs_unjoined_composition`** — the review's point that a 70% join can be enriched for conserved genes in a way the rate cannot show.
- `MIN_JOIN_RATE = 0.70`, an **asserted** floor flagged as asserted. Below it the rung emits `WALL_JOIN_RATE_BELOW_FLOOR` and **no coverage key at all**. Rationale: the repo's most expensive identifier trap gave 0% vocab overlap and a confident AUROC 0.000.
- Per-key-space identifier normalisation is explicit (`SAUSA300_` zero-padding, locus-tag case, ORF casing) — never a generic fuzzy fallback.

**Test strategy:**
- Unit: `join_rate` at 1.0, 0.70 (boundary — assert inclusive as documented) and 0.05.
- Unit: a 0.05-rate rung returns the wall and carries **no** coverage key (absence, not a None beside a number).
- Unit: `joined_vs_unjoined_composition` detects a deliberately conserved-enriched joined subset at an acceptable join rate — the bias the rate alone cannot see.
- Unit: the promoted reader yields byte-identical output to the current CLI helper on the real E. coli feature table (skip if D: absent).
- Unit: `essentiality/cli.py` imports the shared reader rather than defining its own.

### Step 4: Phrasing-floor subtraction — promoted, provenance-gated, narrow by evidence
Files: dna_decode/essentiality/phrasing_floor.py, scripts/essentiality_missed_vocabulary.py, tests/test_essentiality_phrasing_floor.py
Depends on: none

**What changes:**
- Promotes `human_rows` / `ecoli_rows` / `content_words` / `STOP` / `ROBUST_II` out of `scripts/essentiality_missed_vocabulary.py` into the package, so the ladder does not import a `scripts/` module. The script then imports from the package (single definition).
- `phrasing_floor_genes(missed_rows) -> set[str]` using ONLY `ROBUST_II` (`polymerase`, `helicase`, `replication`, `division`, `topoisomerase`, `primase`) — the 6 words measured to survive every bar.
- `coverage_lift_adjusted` = `coverage_lift` recomputed with phrasing-floor genes excluded from the essential set, so the ladder's decay is not inflated by a catalogue bug. Reported BESIDE the raw lift, never replacing it.
- **Provenance gate:** every word in `ROBUST_II` must appear in the decoder's own `_CORE` patterns, enforced at import (raise), not documented. This is why the broad per-rung synonym vocabulary from v1 is NOT here: vocabulary overlap was measured unable to adjudicate the middle (three bars, two bands), so only the floor ships.

**Test strategy:**
- Unit: every `ROBUST_II` word appears in a `_CORE` pattern; a planted non-pattern word RAISES at import (non-vacuous).
- Unit: the floor reproduces 57 genes on the real human missed set (skip if D: absent) — the published number, pinned.
- Unit: `coverage_lift_adjusted` ≤ `coverage_lift` is NOT asserted blindly — a test constructs a case where excluding floor genes RAISES the lift and asserts the function reports it honestly rather than clamping.
- Unit: assert `scripts/essentiality_missed_vocabulary.py` imports from the package (guards re-divergence) and that its existing 5 guards still pass.

### Step 5: Label parsers, each gated on its W0 schema record
Files: dna_decode/essentiality/label_sources.py, tests/test_essentiality_label_sources.py
Depends on: Step 2

**What changes:**
- Decoder-side label registry, namespace-separate from `fba/essentiality_labels.py` (model-gene-keyed). Docstring states the enabling finding and names which repo file records each FBA wall.
- Entries for all 5 rungs with `url`, `key_space`, `technology`, `recorded_in`.
- Parsers: **import and reuse** `essentiality_e3_learned.load_labels` for Goodall and `fba.essentiality_labels.parse_sgd_essential` for yeast (identity asserted by test — no copies). New pure parsers only for NTML, the PLOS GOLD sheets, and BAGEL, each written **against the W0 record's recorded header index and key column**, never against prose.
- `parse_or_refuse(source)` RAISES `SchemaUnverified` when the W0 record is missing or `w0_verified: false`, and `LabelParseError` on a zero-id parse (an empty label set yields a clean-looking coverage of 0 on no data — a plumbing failure wearing a result's costume).

**Test strategy:**
- Unit: each parser against a synthetic fixture in its W0-recorded shape; assert `n_parsed` and that empty/malformed RAISES.
- Unit: yeast parser **is** `fba.essentiality_labels.parse_sgd_essential` and Goodall **is** `essentiality_e3_learned.load_labels` (identity, not equality).
- Unit: `parse_or_refuse` RAISES when the W0 record is absent, and when present-but-unverified. This is the gate that makes Verification Signal 6 real.
- Unit: registry key set equals `transfer_ladder.RUNGS` key set — a rung cannot exist without a declared source.

### Step 6: The ladder runner — reconcile first, freeze thresholds, report every rung
Files: scripts/essentiality_transfer_ladder.py, tests/test_essentiality_transfer_ladder_runner.py
Depends on: Step 1, Step 3, Step 4, Step 5

**What changes:**
- Per rung: resolve label (fetch into `D:/dna_decode_cache/essentiality/`, **never C:**) → `parse_or_refuse` → load annotation → join (+composition) → score with **unchanged** `score_gene` → `coverage_lift` + permutation null (≥200 draws, mean + p95) → `coverage_lift_adjusted` → AUROC as labelled secondary.
- **Reconcile-before-attribute gate, first and blocking:** recompute E. coli and human and assert `coverage_lift` 0.3125 / 0.1633 to 4 dp AND AUROC 0.6952 / 0.5805 to 3 dp, **raising** on mismatch (reuses `compute_lineage_metrics.reconcile_raw_metrics` / `ReconcileMismatch` discipline). Nothing downstream is trusted until it passes — this is the gate that caught my Goodall mis-join.
- **Every rung reported** with a named wall + route tried. Exit 3 `INDETERMINATE_INSUFFICIENT_RUNGS` below `min_rungs`.
- **Freeze `cliff_drop` / `plateau_tol` HERE, from the two EXISTING rungs only** (spread 0.1492) plus their own permutation bands, committed **before any new rung is fetched**, with the derivation recorded in the artifact. **Honestly labelled:** this is a descriptive-consistency lock, NOT a frozen-before-the-numbers test of the endpoint hypothesis — those two numbers are the motivating observation, so a bar derived from them cannot also test them.
- **Single-rung smoke before the full loop**, unbuffered, so a broken fetch/join trips in seconds rather than after five fetches.

**Test strategy:**
- Unit (offline): whole per-rung pipeline on synthetic fixtures; assert coverage_lift, null, join composition, adjusted lift and the secondary AUROC all land in the record.
- Unit: the reconcile gate RAISES on a deliberately perturbed committed value (proves it is live, not decorative).
- Unit: a rung whose `parse_or_refuse` raises `SchemaUnverified` yields `WALL_SCHEMA_UNVERIFIED` and does not abort other rungs; 2 scorable rungs → exit 3.
- Unit: `sha256(core_decoder.py)` matches the plan-start value.
- Unit: the frozen thresholds appear in the artifact together with the evidence used to set them.
- Live (manual, D: + network): full run; output feeds Step 7.

### Step 7: Verdict, memo, and honest limits baked into the template
Files: wiki/essentiality_transfer_ladder_memo_template.md, scripts/essentiality_transfer_ladder.py, tests/test_essentiality_transfer_ladder_verdict.py
Depends on: Step 6

**What changes:**
- `--emit-memo` renders the per-rung table (coverage-primary), the PRIMARY bacterial technology-matched verdict, and the SECONDARY full-ladder verdict labelled technology- AND sampling-frame-confounded.
- **Honest limits written into the TEMPLATE so they cannot be omitted:** (a) label technology is confounded with distance across the full ladder and is irreducible with these sources — which is why the bacterial sub-ladder is primary; (b) **cross-rung AUROC is a DIFFERENT SAMPLING FRAME** (BAGEL two curated extremes, base rate 0.431, vs E. coli genome-wide 0.0928) so the memo refuses cross-rung AUROC and sens/spec comparison and leads with `coverage_lift`; (c) **~10% of the human miss is a catalogue phrasing gap, not phylogeny** — the adjusted lift ships beside the raw one; (d) a rung is one organism, so a between-rung difference is not separable from that organism's annotation quality; (e) `MIN_JOIN_RATE` and the frozen thresholds are asserted bars, and the bar is a consistency lock not an endpoint test.
- The verdict stays the mechanical output of the frozen bar; any post-hoc reading ships beside it under a separate key (the `tb_implied_boundary` `quality_gate` pattern).
- **`INDETERMINATE_BAR_SENSITIVE` is reachable**: if the raw and phrasing-adjusted lifts fall in different verdict bands, the memo says so and refuses a single verdict — the lesson from the vocabulary probe, encoded.

**Test strategy:**
- Unit: rendering for each of the six verdicts; assert all five honest-limit headings appear regardless of verdict.
- Unit: the memo contains **no** cross-rung AUROC or sens/spec comparison (string-level refusal check on rendered output).
- Unit: raw and adjusted lifts in different bands → `INDETERMINATE_BAR_SENSITIVE`, proven with a constructed case.
- Unit: `verdict` and any post-hoc reading occupy distinct JSON keys.

### Step 8: Report-card extension, augment-only
Files: scripts/build_essentiality_report_card.py, tests/test_essentiality_report_card_ladder.py
Depends on: Step 7

**What changes:**
- `load_ladder()` + ladder rows appended to `organisms`, each with its own honest tier (`COVERAGE_SCORED` / `WALL_<reason>`), plus a separate `## Cross-organism transfer ladder` section leading with `coverage_lift`.
- Existing invariants preserved explicitly: no aggregate headline, exit 0 always, per-organism tiers only.
- `human_transfer_auroc`'s AUROC computation replaced by a call to the shared scorer from Step 6 so card and ladder cannot report two different numbers for one rung.

**Test strategy:**
- **Augment-only verified by diff:** snapshot before/after; the 2 pre-existing rows retain every field value, only new rows/keys appear. Non-vacuity: simulate a merge bug overwriting the human row's `metric` and assert the test fails.
- Unit: a wall rung renders its named wall, carries no metric, card still exits 0.
- Regression: the pre-existing transfer-test assertions (`schema`, no `aggregate`/`headline`, transfer above null, spec > 0.9) still pass **with that test file unmodified** — it exercises the scorer this step swaps, so editing it would hide the regression it guards.

### Step 9: Register the measured transfer as a regime row and derive the stale string
Files: dna_decode/eval/regime.py, scripts/regime_map.py, tests/test_regime_boundary.py
Depends on: Step 7

**What changes:**
- Add a `Regime` row for the essentiality conserved-core transfer cell with `split_unit=frozenset({"organism"})` and `organism_transfer="held_out_organism"`, evidence quoting the measured ladder, artifact pointing at the Step 7 memo. `organism_transfer_is_unmeasured()` then SHRINKS — the event its own docstring names as the one worth noticing.
- `_transfer_gap_conditions` emits "no regime in this map carries one today", which becomes FALSE once a row carries that evidence; rewritten to be DERIVED from the live table so it cannot go stale again (the `resolve_active_lock` lesson).
- The new row's verdict is set FROM the measurement, not chosen: `NO_DECAY_DISTANCE_INSENSITIVE` and `DOMAIN_CLIFF` imply different verdict strings and both branches are written before the number lands.

**Test strategy:**
- Unit: `organism_transfer_is_unmeasured()` no longer contains the new key and the count dropped by exactly one (not that it is empty — the other 8 are genuinely unmeasured).
- Unit: the gap-condition string is DERIVED — a table where a row carries `held_out_organism` does not emit the "no regime carries one today" wording, and the converse.
- Unit: the three historical-compression guards in `test_regime_boundary.py` still fail loudly.
- Unit: `_transfer_rank` ordering unchanged; a proposal claiming `held_out_clade` against the new row still returns the gap condition.

## Execution Preview

| wave | steps | parallelism |
|---|---|---|
| 0 | 1, 2, 3, 4 | 4 |
| 1 | 5 | 1 |
| 2 | 6 | 1 |
| 3 | 7 | 1 |
| 4 | 8, 9 | 2 |

- **Total waves:** 5. **Max parallelism:** 4. **Critical path:** Step 2 → Step 5 → Step 6 → Step 7 → Step 8|9 (5 deep).
- **Zero intra-wave file overlap.** Wave 0: Steps 1/2/4 create new modules + own tests; only Step 3 touches a pre-existing file (`essentiality/cli.py`). Wave 4: Steps 8 and 9 share no file. v1's Step-5/Step-3 collision on `essentiality/cli.py` is gone because that step shipped already.
- **The critical path is longer than v1's (5 vs 4) on purpose:** the W0 probe is a serialising gate. That is the point — v1 wrote parsers from prose and the cost of being wrong is a confident number from a silent mis-join.
- Steps 2 and 6 are the only steps needing network + D:; everything else is offline-testable.

## Risk Flags

- **Label reachability remains a HYPOTHESIS.** The FBA registry records the PLOS gold standard as "fetchable" and NTML as existing, but `wiki/essentiality_label_wall_2026-07-28.md` documents seven dead routes for essentiality labels generally. Step 2 verifies rather than assumes. **The plan's value does not require all five rungs** — three scorable rungs including one bacterial pair still discriminate (A) from (B) at reduced power.
- **If only E. coli + human score, the plan delivers no new discrimination.** Pre-committed verdict is `INDETERMINATE_INSUFFICIENT_RUNGS` plus a named-wall table, not a write-up.
- **The full ladder confounds distance with label technology** and cannot be de-confounded with these sources. Mitigated by making the technology-matched bacterial sub-ladder PRIMARY — but that sub-ladder has at most 3 rungs, so its power is low. Disclosed, not solved.
- **n=1 organism per rung.** A between-rung difference is not separable from that organism's annotation quality or curation depth. `coverage_lift` removes the base-rate confound, not this one.
- **The frozen bar is a consistency lock, not an endpoint test** (Step 6) — derived from the two numbers that motivated the question, so it cannot also test them. Stated in the memo rather than quietly implied.
- **`MIN_JOIN_RATE = 0.70`, `cliff_drop`, `plateau_tol`, and the `ROBUST_II` floor are asserted bars**, not derived. Per-rung numbers ship so no reader is limited to the verdict. (Precedent: `MAX_UNRESOLVED_FRACTION = 0.10` is recorded as asserted.)
- **Taxonomic lineage provenance [unverified at plan time].** Step 1 requires each rung's lineage with a `source_url`; the exact NCBI Taxonomy endpoint shape was not confirmed during planning. Step 1 must verify the fetch before committing lineage strings and must NOT fall back to asserting ranks from memory. If no reachable source exists, the ordering becomes an explicitly-flagged asserted constant and the plan says so.
- **The phrasing floor is a FLOOR, not an estimate.** 57 genes survived an adversarial bar; the true mechanism-(ii) share sits somewhere in 0.10–0.25 and this instrument cannot narrow it. `coverage_lift_adjusted` therefore over-corrects less than the truth, in the conservative direction.
- **v2 has not itself been adversarially reviewed.** v1's review produced 3 critical issues; this plan changes the primary metric, deletes a step and adds a step in response. `/brainstorm` on v2 before `/execute-plan` is the repo's own recorded pattern.
- **No sentrux** → no architecture gate. The plan is additive (4 new modules, 1 new script, 3 edited files), so the gate's absence is low-impact.

## Open Questions

1. **Is the eukaryote rung yeast alone, or yeast AND a second eukaryote?** One eukaryote cannot separate "the eukaryote boundary" from "yeast specifically". A second (*S. pombe*, PomBase viability) would materially strengthen (B). Deferred: it doubles the fetch surface against an unverified source.
2. **Should a rung whose `coverage_lift` falls inside its own permutation p95 be reported as a number at all?** Currently yes, with the null beside it. Withholding it (the prospective layer's pattern) is arguably more honest but loses the decay shape, which is the measurement.
3. **Does the new regime row belong in `REGIMES` at all,** or is the essentiality cell deliberately not a regime (it is a curated-catalogue decoder, and `cell_regime.py` already records route→regime separately)? Step 9 assumes it belongs; this is a scope judgement.
4. **Whether to then extend `_CORE` with the ~57-gene phrasing vocabulary.** It is the named follow-on and it **ends the transfer experiment** for any rung it touches. Separate plan with its own before/after, not a tail of this one.

## Verification

1. `uv run pytest tests/ -q` → ≥ 5777 + new tests passed, **0 failed**. Capture `PYTEST_EXIT=${PIPESTATUS[0]}` explicitly — this harness reported exit 0 for a pytest run that exited 1 **four times** on 2026-10-05.
2. `uv run python scripts/essentiality_ladder_w0_probe.py --self-check` → offline, exit 0; then the live run writes the schema record every parser depends on.
3. `uv run python scripts/essentiality_transfer_ladder.py --self-check` → offline, exit 0, reconcile gate passes on the committed E. coli / human values (coverage_lift 0.3125 / 0.1633).
4. `uv run python scripts/essentiality_transfer_ladder.py` (live, D: + network) → artifact written; every declared rung present with a scored row or a named wall.
5. `uv run python scripts/essentiality_missed_vocabulary.py` → still exits 0 and its 5 guards pass after the Step-4 promotion (no behaviour change to the shipped probe).
6. `uv run python scripts/build_essentiality_report_card.py` → exit 0; before/after diff shows the 2 pre-existing rows unchanged in every field.
7. `uv run python scripts/regime_map.py` → exit 0; `organism_transfer_is_unmeasured()` length dropped by exactly 1.
8. `uv run python scripts/contract_number_audit.py` → still `NOTHING_TO_ADJUDICATE`.
9. `uv run python -m scripts.prospective_lock_validate` → `verify_lock ok=True`; `git diff --stat` shows zero changes to the five pinned frozen files. Use the REAL paths — `dna_decode/data/amr_rules.py` does not exist (it is `dna_decode/eval/amr_rules.py`), and a guard naming a nonexistent path passes vacuously.
10. `sha256sum dna_decode/essentiality/core_decoder.py` → unchanged from the plan-start value. The design is void if the decoder was retuned.

## Save-time amendments

Captured at: 2026-10-05
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work — amendments are NOT executable instructions.
If an amendment changed a Step contract (file lists, dependencies, structure), re-validation below
will fire and `/technical-plan` should be re-run before `/execute-plan`.

- v2: coverage_lift replaces AUROC as primary
- v1 Step 5 shipped so removed
- v1 Step 4 deleted on measurement and replaced by a provenance-gated 6-word phrasing floor
- W0 schema probe added as a serialising gate
- frozen bar relabelled a consistency lock
<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
