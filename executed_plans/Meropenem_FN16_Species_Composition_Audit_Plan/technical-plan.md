# Diagnose the 16 unexplained meropenem false negatives, and audit species composition across the provenance-disjoint arm

## Lens status

Repo-index: absent (`update-index` CLI not installed) — source-only reading, no degradation.
Graphify: absent. DESIGN.md: absent (no UI in scope). `project-rules.md`: absent → Project Rules Check omitted.
Tools: `codex` + `gh` present; **`sentrux` absent** → Architecture Gate records `n/a — sentrux not installed`.
Baseline tests: `5694 passed / 12 skipped / 0 failed` (full suite, verified 2026-10-05).

## Problem Statement

`Klebsiella × meropenem` is a **SCORED cell on the standing trust surface** (`wiki/decoder_validation_report_card.json`) reporting **acc 0.683 / sens 0.467 / spec 0.900, tp=14 fp=3 tn=27 fn=16, N=60**. Its 16 false negatives were attributed to the frozen rule being *"blind to porin-loss-mediated R (expected FN mode)"* — the rule's own `validated` string, verbatim in `dna_decode/eval/amr_rules.py::DRUG_RULE['meropenem']`.

**That attribution is measurably false.** The rule is `threshold=1, subclass_any={'CARBAPENEM'}`, and AMRFinder files `ompK35`/`ompK36` truncations under `Subclass=CARBAPENEM` — so the rule *counts* porin loss (10/10 porin-only isolates are called R, recorded 2026-10-04). The 16 FN therefore have **no explanation at all**, and a SCORED cell whose dominant error mode is unexplained is a trust-surface gap.

Reconnaissance done while scoping this plan narrows it sharply, and the finding is **not** the expected one:

- All **16/16** FN carry only `ampC-Kaer` (13) / `ampC_Kaer-1` (3), `Subclass=CEPHALOSPORIN`. **Zero `ompK*` rows in any FN** — so the porin story fails twice over: the rule can see porins, and these isolates have no porin call to see.
- `Kaer` is *Klebsiella **aerogenes***. Three sampled FN resolve to **`Klebsiella aerogenes`** from their own GenBank `ORGANISM` field; two sampled TP resolve to **`Klebsiella pneumoniae`**.

**Leading hypothesis:** the 16 FN are a **species-composition artifact** — *K. aerogenes* isolates scored with `-O Klebsiella_pneumoniae` under a rule validated on *K. pneumoniae*, whose meropenem resistance runs through AmpC hyperproduction plus porin loss rather than the acquired carbapenemases (`blaKPC`/`NDM`/`OXA-48`) the rule targets.

This is the **same defect class the repo recorded on 2026-10-04** for the CDC AR Bank panel (*"a genus-named reference panel is not a species cohort"*, 14 of 143 non-*pneumoniae*) — but found this time **inside a frozen-era provenance-disjoint cell**, which is more serious: a published SCORED number may be measured on a mixed-species cohort.

**Scope:** diagnose, generalize the check, and **disclose**. This plan does **not** re-score the cell. The 10 `provenance_disjoint_validation_*.json` artifacts are named frozen units of the reproducibility freeze (`wiki/reproducibility_freeze_2026-06-13.md`, Scored-cell artifacts row), so changing a published number is a **user authority call** — surfaced in Open Questions, deliberately not taken.

**Non-goals:** re-scoring any cell; editing any of the five sha256-pinned frozen files; changing the meropenem rule (that is the *other* authority fork, declined for this plan); the lone-porin over-call fix.

## Codebase Context

**The cell and its substrate (all verified locally, all offline):**

| thing | where | state |
|---|---|---|
| cohort labels | `data/raw/klebsiella_provdisjoint_meropenem/selected.tsv` | 60 rows, `GCA_x<TAB>R|S` |
| cached AMRFinder | `.../amrfinder_runs/<acc>/{main.tsv,mutations.tsv}` | **60/60 present** |
| cached assemblies | `.../refseq/<acc>/{genome.fna,annotations.gff3,annotations.gbk}` | **60/60 present** — `annotations.gbk` carries `  ORGANISM` |
| committed metrics | `wiki/provenance_disjoint_validation_klebsiella_merop_2026-06-10.json` | `_schema: provenance-disjoint-validation-v1`; `metrics{n_scored,tp,fp,tn,fn,abstain,acc,sens}` |

**Reconciliation is already proven feasible.** Re-running the deployed `amr_rules.call_resistance(main_tsv, 'meropenem', organism='Klebsiella')` over the 60 cached runs reproduces the committed matrix **exactly — tp=14 fp=3 tn=27 fn=16, 0 abstain** — so the attribution gate will pass and the 16 FN accessions are recoverable without any fetch.

**Generalization scope, measured per cohort** (`refseq` = species-auditable, `amrfinder_runs` = reconcilable):

| cohort | refseq | amrfinder | reconcilable? |
|---|---|---|---|
| campylobacter cipro | 40 | 40 | yes |
| e.coli ceftriaxone / cipro / gentamicin / tetracycline | 60 each | 60 each | yes |
| klebsiella cipro | 60 | 60 | yes |
| **klebsiella meropenem** | 60 | 60 | **yes (primary)** |
| klebsiella ceftriaxone | 60 | **54** | no — partial |
| klebsiella tetracycline | 60 | **33** | no — partial |
| klebsiella gentamicin | 60 | **3** | no — partial |

So **species composition is auditable for all 10** cells; **confusion-matrix reconciliation is possible for 7**. A cohort that cannot reconcile gets its species composition reported and **no outcome cross-tab** — an un-reconciled attribution is exactly what this plan's gate exists to refuse.

**The rule under test:** `amrfinder_classes_for('meropenem')` = `{BETA-LACTAM, CARBAPENEM, MULTIDRUG}`; `rule_for('meropenem')` = `{threshold: 1, subclass_any: {'CARBAPENEM'}}`; `calibrated_rule_for('Klebsiella','meropenem')` returns **None**, so the default `DRUG_RULE` path is the validated one (no opt-in registry entry fires).

**Independent cross-check available:** `dna-identify` shipped 2026-10-04. *K. aerogenes* is **not** one of its 14 supported organisms, and `identify/thresholds.py` records K. aerogenes landing at **0.10849+**, above `MAX_DISTANCE = 0.1036` — so the router **rejects** it. An `ABSTAIN` from `dna-identify` on an FN is therefore a second, method-independent signal that the isolate is outside the scored organism's set.

### Reusable-Code Survey

Searched: `dna_decode/eval/*.py`, `dna_decode/data/*.py`, `scripts/compute_lineage_metrics.py`, `scripts/build_identify_reference.py`, `dna_decode/data/annotations.py`, `dna_decode/data/organism_vocab.py` (graphify report absent; no `src/`-`lib/`-`utils/` dirs — the package is `dna_decode/`).

- **`scripts/compute_lineage_metrics.py::reconcile_raw_metrics` + `ReconcileMismatch` — REUSE, do not rebuild.** Already does precisely the needed gate: recompute the confusion from cached runs with the cohort's *original* rule-path args, raise rather than return on disagreement with the committed artifact, and return `(raw_conf, preds_by_acc)` — and `preds_by_acc` is exactly the per-accession prediction map the species cross-tab needs. Verified importable (`scripts/` has `__init__.py`; `from scripts.compute_lineage_metrics import reconcile_raw_metrics` succeeds).
- **`scripts/build_identify_reference.py::organism_from_genbank` — REUSE by relocating.** Reads only the GenBank header, returns `None` rather than raising, and its docstring already carries the right discipline (*"the caller MUST count the Nones, because a silently-shrinking corpus is the failure mode here"*). It currently lives in `scripts/`, so a package module importing it would invert the layering — Step 1 moves it to the module that already owns GenBank parsing.
- **`dna_decode/data/annotations.py` — the correct home.** Already the GenBank parsing module (`parse_genbank`, `load_annotation_table`) and has **no** organism extractor, so this is an additive gap-fill rather than a competing parser.
- **`dna_decode/data/cell_key.py::canonical_cell_key`** — the single join key shared by all report-card sidecars (M2 discipline). The new sidecar must key on it, not on a cohort directory name.
- **`scripts/provdisjoint_source_concentration.py` + `build_validation_report_card.py::load_source_concentration`** — the template for a namespace-separate, augment-only disclosure: own artifact, own loader, own rendered section, never merged into `load_scored()`.
- **`dna_decode/data/organism_vocab.py::canonical` / `supported_canonicals`** — normalizes an organism string to a canonical token and raises `UnknownOrganism` rather than defaulting. Reuse for classifying a resolved species, so *K. aerogenes* surfaces as unknown-to-the-vocabulary instead of being silently bucketed.

**Not reused, deliberately:** `dna_decode/eval/clonality.py` (lineage, a different axis); `scripts/salmserovar_*` species logic (serovar-specific).

## Pre-Change Baseline

- **The cell:** `Klebsiella × meropenem` SCORED — acc **0.683**, sens **0.467**, spec **0.900**, tp **14**, fp **3**, tn **27**, fn **16**, N **60**, abstain **0**.
- **The 16 FN have no valid explanation on record.** The only one offered (porin blindness) is false, and is stated inside a sha256-pinned file that this plan must not edit.
- **Species composition of every provenance-disjoint cohort: UNMEASURED.** No artifact, script, or test records which species any of the 10 cohorts actually contains. The cohort directory name (`klebsiella_*`) and the AMRFinder flag (`-O Klebsiella_pneumoniae`) disagree in granularity (genus vs species) and nothing checks that.
- **Report card:** 27 cells, 6 disclosure layers (`lineage`, `source_concentration`, `prospective`, `doubt_layer`, `organism_scope`, `error_rates`). No species-composition layer.
- **Suite:** 5694 passed / 12 skipped / 0 failed.

## Verification Signal

The plan succeeds iff all of these hold, each independently checkable:

1. **The reconcile gate fires before any attribution.** The audit reproduces `tp=14 fp=3 tn=27 fn=16` for the meropenem cell and **raises `ReconcileMismatch`** rather than reporting a cross-tab if any cohort disagrees with its committed artifact.
2. **A measured species × outcome cross-tab** for the meropenem cell, naming how many of the 16 FN are non-*pneumoniae* — a count, not an impression.
3. **Unresolved species are counted, never dropped.** The artifact carries `n_unresolved`, and the audit refuses a composition verdict for a cohort whose unresolved fraction exceeds a stated bar.
4. **The 9 other cohorts are audited for composition**, with the 3 non-reconcilable ones reporting composition and **explicitly no** outcome cross-tab.
5. **Augment-only, verified by diff:** report card has **27 cells before and after**, an identical cell-key set, and **zero** non-species field changes.
6. **Nothing frozen moved:** the five sha256-pinned files byte-unchanged, all 10 `provenance_disjoint_validation_*.json` byte-unchanged, `verify_lock` still `ok=True`.
7. **Suite green**, and the new tests fail when the defect is reintroduced (non-vacuity demonstrated, not asserted).

## Implementation Steps

### Step 1: Move the GenBank ORGANISM parser into the package as the single definition
Files: dna_decode/data/annotations.py, scripts/build_identify_reference.py, tests/test_genbank_organism.py
Depends on: none

**What changes:**
- Add `organism_from_genbank(gbk: Path, max_header_lines: int = 60) -> str | None` to `dna_decode/data/annotations.py` — the module that already owns GenBank parsing and currently has no organism extractor. Body moved verbatim from `scripts/build_identify_reference.py`, including the header-only read and the `None`-not-raise contract.
- `scripts/build_identify_reference.py` imports it and keeps a thin alias so its existing callers and tests are untouched. A package module must not import from `scripts/`, which is why the direction is scripts → package.
- Docstring records that the Nones are the caller's responsibility to count.

**Test strategy:**
- Unit: real `ORGANISM` line → value; absent line → `None`; unreadable/missing file → `None` not raise; a line beyond `max_header_lines` → `None` (the budget actually binds).
- **Behaviour-preservation pin:** the relocated function returns byte-identical output to the pre-move implementation on ≥3 real cached `annotations.gbk` files (one *K. pneumoniae*, one *K. aerogenes*, one from another cohort) — a refactor that changes a value is a defect, not a move.
- **The two real existing consumers must stay green WITHOUT being edited**, which is the whole reason Step 1 keeps a thin alias rather than rewriting call sites: the existing suite imports `organism_from_genbank` *from the script* at five assertion sites, and `scripts/identify_validate.py` imports it from the same place. Neither appears in this step's `Files:` because neither is modified — if either needs an edit, the alias is wrong and the step should stop rather than widen. Verified before planning: no other module imports it.

### Step 2: Pure species-composition logic, with an explicit refusal on thin resolution
Files: dna_decode/eval/cohort_species.py, tests/test_cohort_species.py
Depends on: Step 1

**What changes:**
- New non-frozen module. `resolve_cohort_species(refseq_root, accessions) -> dict[str, str | None]` using `organism_from_genbank`.
- `classify_species(resolved, expected_amrfinder_organism) -> dict` partitioning into `matches_expected` / `same_genus_other_species` / `other_genus` / `unresolved`. The genus/species split is the whole point: `Klebsiella aerogenes` is the same genus and a different species, which is the case a genus-named cohort hides.
- `compose(resolved, preds_by_acc, labels) -> dict` → species × (tp/fp/tn/fn) cross-tab. `preds_by_acc` is supplied by the caller (from `reconcile_raw_metrics`), so this module never calls `call_resistance` and stays pure/offline.
- `composition_verdict(...)` **REFUSES** (returns an explicit `INSUFFICIENT_RESOLUTION` state with no composition claim) when the unresolved fraction exceeds `MAX_UNRESOLVED_FRACTION`. A composition computed over a silently-shrunken denominator is the documented failure mode.
- Reuse `organism_vocab.canonical` to tag a resolved species as known/unknown to the supported set, catching `UnknownOrganism` so an unsupported species is *reported* rather than crashing the audit.

**Test strategy:**
- Unit on synthetic inputs: a pure-*pneumoniae* cohort → `other_genus`/`same_genus_other_species` both 0; a mixed cohort → exact counts; `Klebsiella aerogenes` classified `same_genus_other_species`, **not** `matches_expected` (the discriminating case).
- The refusal: unresolved above the bar → `INSUFFICIENT_RESOLUTION` **and no composition numbers in the payload** — asserted by key absence, so a future change cannot leak a number past the refusal.
- Denominator discipline: `matches + same_genus + other_genus + unresolved == len(accessions)`, asserted exactly.
- A cross-tab whose cells do not sum to the supplied confusion matrix raises rather than returning.

### Step 3: The audit runner — reconcile first, attribute second
Files: scripts/provdisjoint_species_audit.py, tests/test_provdisjoint_species_audit.py
Depends on: Step 2

**What changes:**
- New script. Per provenance-disjoint cohort: load `selected.tsv`, load the committed artifact, call the **reused** `reconcile_raw_metrics` (raising `ReconcileMismatch` on disagreement), then `cohort_species.compose`.
- A cohort whose cached AMRFinder runs are incomplete (kleb gentamicin 3/60, tetracycline 33/60, ceftriaxone 54/60) is **not** reconcilable: it reports species composition with `outcome_crosstab: null` and a stated reason. It must never fall back to attributing outcomes on a partial run set.
- Cohort → `(registry_organism, amrfinder_organism, drug)` is read from each committed artifact's own fields, never inferred from the directory name — the directory name being untrustworthy is the finding under test.
- Emits `wiki/provdisjoint_species_audit_<date>.{md,json}`, keyed by `canonical_cell_key`, schema `provdisjoint-species-audit-v1`. Date-stamped (the repo's literal-dated-filename lesson). Read-only; exit 0 always (a report, not a gate).

**Test strategy:**
- Offline unit tests with a tmp-dir fixture cohort (hand-written `selected.tsv` + minimal `main.tsv` + minimal `.gbk`): reconcile pass → cross-tab present; reconcile fail → raises and **writes nothing**.
- A partial-runs cohort → composition present, `outcome_crosstab` null, reason recorded.
- **Non-vacuity:** a fixture where the committed artifact disagrees by one cell must fail — demonstrated, not asserted.
- Real-data test (skipped when the cohort is absent): the meropenem cell reconciles to `14/3/27/16`.

### Step 4: Run the audit and write the diagnosis memo
Files: wiki/provdisjoint_species_audit_2026-10-05.json, wiki/provdisjoint_species_audit_2026-10-05.md, wiki/meropenem_fn16_diagnosis_2026-10-05.md
Depends on: Step 3

**What changes:**
- Run against all 10 cohorts; commit both artifacts.
- Write the memo for the primary question: what the 16 FN are, with the measured species counts, the `ampC-Kaer` determinant fingerprint (16/16), and the zero-`ompK*` finding that kills the porin account twice over.
- **Cross-check with `dna-identify`** on the FN set where Docker is available: an `ABSTAIN` is a second, method-independent signal that the isolate is outside the scored organism. If Docker is unavailable, record that the cross-check did not run — never imply it passed.
- The memo states the honest limits explicitly: a GenBank `ORGANISM` label is the submitter's, not a wet-lab identification; it is strong evidence of composition and is **not** itself a re-scoring; and it says plainly whether the hypothesis survived, including if the species split does **not** explain all 16.

**Test strategy:**
- `tests/test_provdisjoint_species_audit.py` gains a committed-artifact schema test (required keys, cell keys resolve, counts sum).
- The existing `test_every_fba_wiki_artifact_is_parseable_json`-style discipline: the new JSON is re-parsed after writing.
- Manual: read the memo against the artifact and confirm every number in prose appears in the JSON.

### Step 5: Render it as a namespace-separate, augment-only report-card disclosure
Files: scripts/build_validation_report_card.py, tests/test_build_validation_report_card.py
Depends on: Step 4

**What changes:**
- `load_species_composition()` + `build_species_block()` + a rendered **"Species-composition disclosure"** section, mirroring the source-concentration layer exactly: own artifact, own loader, own section.
- **Namespace-separate under its own `species_composition` key.** It must never feed `load_scored()` — a cell shares its `(organism, drug)` key with the scored row, and merging would overwrite a published metric with a composition field (the documented shared-key trap).
- The section **never changes a cell's state or metrics**; it qualifies them. A cell with clean composition renders nothing rather than a reassuring line, and a cohort in `INSUFFICIENT_RESOLUTION` renders the refusal with **no** numbers.
- Section text states the layer's scope: this describes the COHORT the number was measured on, not the rule.

**Test strategy:**
- **Augment-only by diff:** build the card before and after, assert 27 cells both times, identical key sets, and that every differing field lives under `species_composition`.
- Anti-overwrite test using deliberately **different** values for the scored and species blocks, so a merge bug cannot pass by coincidence (the pattern the prospective layer's test established).
- Non-vacuity: simulate the merge bug and confirm the test goes red.
- A clean-composition cell renders no line; an `INSUFFICIENT_RESOLUTION` cell renders no numbers.

### Step 6: Documentation, and prune the stale NEXT.md
Files: CLAUDE.md, LESSONS_LEARNED.md, NEXT.md
Depends on: Step 5

**What changes:**
- `CLAUDE.md`: one compact entry — the false porin attribution, the measured species finding, that the diagnosis is a disclosure and **not** a re-score, and the named authority fork. Compress within the existing word ceiling rather than raising it.
- `LESSONS_LEARNED.md`: the reusable bullet — *an error mode documented in a rule's own metadata is a claim, and a frozen string cannot be corrected by ordinary use*; plus the generalized form of the genus-vs-species trap now seen twice in two days, on two different cohorts.
- `NEXT.md`: stamped `2026-09-28` and describes the F1–F4 campaign as in-flight. Prune to current state and record the authority forks (re-score the cell; the lone-porin rule fix) as awaiting a user decision.

**Test strategy:**
- `tests/test_claude_md_growth_ceiling.py` + `tests/test_claude_md_citations.py` green (the ceiling is a hard bar and every cited artifact must exist).
- Full suite green; frozen five sha256 re-verified; `verify_lock` re-run.

## Execution Preview

| Wave | Steps | Parallel |
|---|---|---|
| 0 | Step 1 | 1 |
| 1 | Step 2 | 1 |
| 2 | Step 3 | 1 |
| 3 | Step 4 | 1 |
| 4 | Step 5 | 1 |
| 5 | Step 6 | 1 |

**Total waves: 6. Max parallelism: 1.** Fully serial, and that is the honest shape rather than a packaging failure: the parser must exist before the pure logic that calls it, the pure logic before the runner, the runner before its output, the output before the renderer that reads it, and the docs last so they describe what shipped. Splitting any step to manufacture parallelism would create intra-wave file overlap on `cohort_species.py` or the audit artifact. Confirmed by the toolkit: `overlaps: []`.

**Critical path:** Step 1 → 2 → 3 → 4 → 5 → 6 (the whole plan).

**Architecture Gate:** `n/a — sentrux not installed`.

## Risk Flags

- **The hypothesis may be wrong, or only partly right.** Three of 16 FN were sampled, not all 16. If some FN are genuine *K. pneumoniae*, the species split explains only part of the gap and the residual is a real catalog/detection question. Step 4 must report the residual rather than presenting a partial explanation as complete — the plan is a diagnosis, not a confirmation exercise.
- **A GenBank `ORGANISM` label is the submitter's claim, not a wet-lab identification.** It is the same evidence class the project's own G1 circular-label gate covers. It is strong for composition and must not be dressed up as independent species truth; `dna-identify` abstention is a *second* signal of the same broad kind (sequence comparison), not an orthogonal one.
- **Re-scoring is out of scope and must stay out.** The 10 provdisjoint artifacts are frozen units. If the audit shows the cell is mixed-species, the temptation is to immediately produce a "corrected" sens — that changes a published number and is a user authority call (Open Questions).
- **`reconcile_raw_metrics` lives in `scripts/`,** so this plan makes a second consumer of a script-level function. Acceptable (it is non-frozen, `scripts/` is a package, and the import is verified), but if a third consumer appears the function should move into `dna_decode/eval/`. Noted, not done here.
- **Three Klebsiella cohorts cannot be reconciled** (3/60, 33/60, 54/60 cached runs). Completing them means re-running AMRFinder via Docker — out of scope. The risk is that a reader treats "composition reported, no cross-tab" as a weaker finding when it is simply a different, honestly-scoped one.
- **`-O Klebsiella_pneumoniae` on a *K. aerogenes* genome drives species-specific point-mutation screening**, so the determinant calls for those isolates may be wrong in ways beyond the meropenem rule. This plan reports composition; it does not re-run AMRFinder under the right flag.
- **Adding a 7th entry to `DISCLOSURE_LAYERS` is deliberately NOT part of Step 5.** That tuple is enumerated against CLI renderers by test, so adding it would pull per-call CLI wiring into scope. Whether a cohort-composition caveat belongs on every call is a real design question (Open Questions), and shipping a renderer before answering it risks the `doubt_layer` mistake of a line that can never fire.

## Open Questions

1. **Should the `Klebsiella × meropenem` cell be re-scored on *K. pneumoniae* only?** A user authority call: it edits a frozen unit of the reproducibility freeze and changes a published SCORED number. The audit will say what the pneumoniae-only matrix *would* be; acting on it is the user's decision, not the executor's.
2. **Should species composition also surface at CALL time** as a 7th `DISCLOSURE_LAYERS` entry? Argument for: the CLI's `validation:` line quotes this cell's metric, so a reader of a live call inherits the mixed-cohort caveat. Argument against: it describes the cohort, not the call, and every existing call-time layer describes the call. Deferred rather than guessed.
3. **Does the same composition defect affect the *E. coli/Shigella* cohorts?** Those are genus-*and*-group named (`escherichia_coli_shigella`), where mixed composition may be intentional. Step 4 will report it; whether a *Shigella* genome inside an E. coli cell is a defect or the declared design needs the user's reading of the original cohort intent.
4. **Is `MAX_UNRESOLVED_FRACTION` the right bar, and at what value?** Proposed as an asserted threshold with its value recorded, not derived. All 60 meropenem assemblies have a `.gbk`, so it should not bind on the primary cell; it exists for the other nine.

## Verification

1. `uv run python -m pytest tests/ -q` — full suite green; compare against the 5694/12/0 baseline. Run with `python -u` and assert on the printed summary, **not** the harness exit code (a piped pytest reports the pipe's status, which has misreported twice in this repo).
2. `uv run python scripts/provdisjoint_species_audit.py` — exit 0; artifacts written; re-parse the JSON.
3. Reconcile proof: the meropenem cell reproduces `tp=14 fp=3 tn=27 fn=16`; a deliberately perturbed fixture raises `ReconcileMismatch`.
4. Augment-only diff: 27 cells before and after, identical key set, all differing fields under `species_composition`.
5. `sha256sum -c` the five frozen files; `git diff --stat` the 10 `provenance_disjoint_validation_*.json` (must be empty); re-run `verify_lock` → `ok=True`.
6. `uv run python scripts/contract_number_audit.py` — still clean, so no contract prose drifted from its artifact.
7. `uv run python scripts/project_status.py` — cell counts unchanged (this plan adds no decoder).
8. Manual: read `wiki/meropenem_fn16_diagnosis_2026-10-05.md` and confirm every number in prose appears in the JSON artifact.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,test-strategy-leak -->
