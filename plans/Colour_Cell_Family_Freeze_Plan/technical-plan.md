# Colour-Cell Family Freeze — Technical Plan

## Lens status

- **Correctness/feasibility:** applied — every file path, contract field and count below was verified live against the repo at `4e601d1`.
- **Risk/regression:** applied — no decoder call path changes; contracts, declared-gate strings and docs only. Frozen AMR surface untouched.
- **External-tool surface:** n/a — no new CLI/library/service is integrated. The only external surface touched is the in-repo `cell_registry` contract schema, verified by reading the dataclass.
- **Complexity/maintenance:** applied — the freeze roster is a literal guarded by a derived cross-check, mirroring `shipped_decoder_surface.py`.
- **Alternatives:** the curate branch was considered and deliberately excluded (see Open Questions).
- **Test gaps:** applied — the freeze guard must be proved non-vacuous by simulating a 20th cell.
- **sentrux gate:** n/a — sentrux not installed.

## Problem Statement

The animal colour/plumage family is **19 CLI cells / 65 loci**, every one shipping as `KNOWLEDGE_BASELINE`.
`scripts/colour_cell_substrate_screen.py` (shipped `4e601d1`) established, by derivation from the committed
catalogs, that the family is blocked twice over:

- **40 of 65 loci (62%) record no causal variant at all**, and **7 of 19 cells record none for any locus**
  (alpaca, cattle, mouse, pig, pigeon, rabbit, sheep). Those cells are **unvalidatable as written** — a locus
  whose causal variant is unspecified cannot be scored against any genotype file, so no cohort helps.
- Of the 25 loci that do record one, **14 (56%) are indel/structural**, off-panel for any SNP array. That is
  empirically what sank the one measured cell (dog: black 160/161 = 0.994, every other base colour unscorable).

Two concrete trust-surface defects follow, both verified live:

1. All 19 contracts carry `incoming_data_gate="n/a"` while gates **G9/G10** (added to the negative-results map
   in `4e601d1`) apply directly to them.
2. **7 cells promise a promotion path that does not exist.** `dna-rabbitcolor`'s `demotion_rule` reads
   *"KNOWLEDGE_BASELINE; MEASURED needs a free rabbit genotype+observed-colour cohort"* — but a cohort is
   necessary and **not sufficient**: with zero causal variants recorded, no cohort could score it.

**Goal:** freeze the family at 19, make each cell's contract state its real screen verdict, and make adding a
20th cell fail loudly rather than silently. **Non-goal:** curating the 40 unrecorded loci (separate plan —
see Open Questions). **Non-goal:** any change to decoder behaviour or output.

## Codebase Context

- `dna_decode/pigment/mammal_color.py` — `MAMMAL_CATALOGS: dict[str, MammalCatalog]`, 14 species; each
  `MammalCatalog.loci: dict[str, Locus]`. `Locus` is a dataclass whose free-text provenance lives in `source`
  and `note` (**singular** — a reader looking for `notes` silently under-reads).
- `dna_decode/pigment/{dog_coat,cat_coat,horse_coat,chicken_plumage,pigeon_plumage}.py` — 5 more catalogs,
  exposing module-level `LOCI` rather than `.loci`.
- `dna_decode/data/cell_registry.py` — `CellContract` dataclass; colour contracts defined around L446-L600.
  `cells()` is the accessor (there is **no** `CELL_CONTRACTS` symbol). Relevant fields: `evidence_tier`
  (`EvidenceTier` = INDEPENDENT_MEASURED / NEAR_INDEPENDENT / FAITHFUL_TO_TOOL / KNOWLEDGE_BASELINE /
  NO_FREE_SOURCE / NOT_CENSUSED), `incoming_data_gate` (documented **DECLARED, not executed**), `demotion_rule`.
  `cli_routable_manifest()` is the repo's canonical **derived-not-hardcoded** pattern to mirror.
- `dna_decode/cli.py` — `TRAITS` dict, 44 entries; each `{"summary","validation"}`. The 19 colour entries.
- `scripts/colour_cell_substrate_screen.py` — the screen. Pure logic currently lives in `scripts/`, which an
  in-package module **cannot import**; this forces the Step 1 move.
- `wiki/colour_cell_substrate_screen_2026-08-26.{md,json}` — the derived artifact and memo.
- `wiki/negative_results_map_2026-06-13.md` — G1–G10; G9/G10 added `4e601d1`.
- `tests/test_colour_cell_substrate_screen.py` (27 tests), `tests/test_cli_dispatch.py`,
  `tests/test_advertised_commands.py` — the guards that must stay green.

### Reusable-Code Survey

- **`dna_decode/data/cell_registry.py::cli_routable_manifest()`** — reuse its *derived-from-live-catalogs*
  shape for the freeze roster cross-check, per [[feedback_hardcoded_exclusion_list_undercovers]] (third
  instance in this repo).
- **`dna_decode/data/shipped_decoder_surface.py` + `tests/test_shipped_decoder_surface.py`** — reuse the
  "literal roster + test that nothing ships invisibly" pattern verbatim for the freeze guard.
- **`scripts/colour_cell_substrate_screen.py::{classify_variant, summarise, collect}`** — reuse as the single
  derivation; Step 1 relocates rather than reimplements.
- `graphify-out/GRAPH_REPORT.md` absent; no `src/utils`, `lib`, `common` dirs — searched:
  `graphify-out/GRAPH_REPORT.md`, `src/utils`, `src/lib`, `src/common`, `utils`, `lib`, `common`.

## Pre-Change Baseline

Captured live at `4e601d1`:

- 19 colour cells / 65 loci; SNV 11 · INDEL 11 · STRUCTURAL 3 · **UNRECORDED 40**.
- 7 cells `UNSCREENABLE_NO_CAUSAL_VARIANTS_RECORDED`; 9 `PARTIALLY_SNV_TRACTABLE`; 1 `NO_LOCUS_SNV_TRACTABLE`; 2 `FULLY_SNV_TRACTABLE`.
- All 19 contracts: `evidence_tier=KNOWLEDGE_BASELINE`, `incoming_data_gate="n/a"`.
- 7 `demotion_rule` strings state a cohort-only promotion path that cannot succeed.
- Suite: **3762 passed / 0 failed / 10 skipped** (excluding `tests/test_models_foundation.py`, whose single
  failure is the pre-existing `transformers` trust-remote-code drift).
- `uv run python scripts/colour_cell_substrate_screen.py --self-check` exits 0.

## Verification Signal

- Full suite still **3762+ passed / 0 failed**, with the new freeze tests added and none removed.
- `--self-check` still exits 0 and the three headline counts (65 / 40 / 14) still reproduce from the catalogs.
- The freeze guard is proved **non-vacuous**: a simulated 20th colour cell makes it fail.
- Every one of the 19 contracts reports a non-`n/a` `incoming_data_gate` naming G9 and/or G10 consistent with
  its own screen verdict.
- No cell whose screen verdict is `UNSCREENABLE_*` still claims a cohort alone would promote it.
- Behaviour unchanged: `dna-decode rabbitcolor --loci A=A/a,B=B/b,C=C/C,D=D/d,E=E/e` produces byte-identical
  output before and after.
- `amr_rules.py` + `calibrated_amr_rules.json` sha256 unchanged.

## Implementation Steps

### Step 1: Relocate the screen's pure logic into the package
Files: dna_decode/pigment/substrate_screen.py, scripts/colour_cell_substrate_screen.py, tests/test_colour_cell_substrate_screen.py
Depends on: none

**What changes:**
- New `dna_decode/pigment/substrate_screen.py` holding the pure logic moved verbatim from the script:
  `classify_variant`, `snv_panel_scorable`, `_loci_of`, `_source_of`, `collect`, `summarise`, plus a new
  `verdicts() -> dict[str, str]` returning `{cell: verdict}` as the single derivation both consumers use.
- **Forced by design, not gratuitous:** `dna_decode/data/colour_cell_freeze.py` (Step 2) cannot import from
  `scripts/`, so the logic must live in-package. The script becomes a thin CLI (`main`, `self_check`,
  `_DOG_TRUTH`, `_CATALOG_GAPS`, artifact writing) importing from the package.
- Regex ordering (STRUCTURAL → INDEL → SNV) and the `dataclasses.fields()` reader move **unchanged**; both are
  load-bearing and already pinned by tests.
- Update the test module's imports; keep every existing assertion identical.

**Test strategy:**
- All 27 existing tests in `tests/test_colour_cell_substrate_screen.py` pass unmodified apart from the import line.
- Re-run `--self-check` (must stay exit 0) and the full screen; assert the emitted JSON is byte-identical to the
  committed artifact — proving a pure move, not a behaviour change.

### Step 2: Add the freeze declaration module
Files: dna_decode/data/colour_cell_freeze.py, tests/test_colour_cell_freeze.py
Depends on: Step 1

**What changes:**
- `FROZEN_COLOUR_ROUTES: frozenset[str]` — the 19 routes, a literal roster (the deliberate friction point).
- `FREEZE_DATE = "2026-08-26"`, `FREEZE_RATIONALE` citing `wiki/colour_cell_substrate_screen_2026-08-26.md`
  and gates G9/G10.
- `freeze_status(route) -> {"frozen": bool, "screen_verdict": str, "gates": tuple[str, ...]}`, where
  `screen_verdict` comes from `substrate_screen.verdicts()` — **derived, never re-typed**.
- `gates_for_verdict(verdict) -> tuple[str, ...]`: `UNSCREENABLE_*` → `("G9",)`; `NO_LOCUS_SNV_TRACTABLE` →
  `("G10",)`; `PARTIALLY_SNV_TRACTABLE` → `("G9","G10")` when it has both unrecorded and off-panel loci, else
  the applicable one; `FULLY_SNV_TRACTABLE` → `()`.
- Docstring states plainly that this is an **attention/scope freeze, not enforcement** — it makes adding a
  cell fail a test, which a determined edit can still ratify.

**Test strategy:**
- Unit: `gates_for_verdict` over all five verdict values, including the empty tuple for fully-tractable.
- Unit: `freeze_status` on a frozen route and an unknown route.

### Step 3: Make the 19 cell contracts declare their real gate and a truthful promotion path
Files: dna_decode/data/cell_registry.py, tests/test_colour_cell_contracts.py
Depends on: Step 1

**What changes:**
- Replace `incoming_data_gate="n/a"` on each of the 19 colour contracts with the gate(s) that actually apply,
  matching that cell's screen verdict (e.g. rabbit → `"G9 — no causal variant recorded for any of its 5 loci"`;
  horse → `"G9+G10 — 3 of 5 loci off-panel (indel/structural)"`).
- **Correct the 7 false promotion paths.** For every `UNSCREENABLE_*` cell, the `demotion_rule` must stop
  claiming a cohort alone promotes it, and state the real precondition: *record the causal variants first; a
  cohort is necessary and not sufficient.*
- `dna-coatcolor`'s `demotion_rule` already records the measured substrate limitation — leave it, and only add
  the gate field.
- `evidence_tier` stays `KNOWLEDGE_BASELINE` for all 19 (see Open Questions — no lower tier fits).

**Test strategy:**
- Assert no colour contract has `incoming_data_gate == "n/a"`.
- Assert each declared gate string is consistent with that cell's live `verdicts()` entry — so a future catalog
  change that flips a verdict fails here rather than drifting.
- Regression: assert no `UNSCREENABLE_*` cell's `demotion_rule` matches a cohort-only promise (pin the exact
  rabbit wording that was wrong).

### Step 4: Surface the screen verdict in the 19 CLI trait contracts
Files: dna_decode/cli.py, tests/test_colour_cell_substrate_screen.py
Depends on: Step 1

**What changes:**
- Append the cell's screen verdict to each of the 19 colour `TRAITS[...]["validation"]` strings, in the shape
  already used by the corrected `coatcolor` entry: what the rule is, what the screen found, and the artifact
  reference.
- For the 7 unscreenable cells the wording must say **unvalidatable as written** — a user reading
  "deterministic curated OMIA epistatic rule" today learns nothing about that.
- `--help` text is generated from `summary`, not `validation`, so terminal help is unaffected; this is the
  trust-surface string only.

**Test strategy:**
- Assert every colour trait's `validation` names its verdict and cites
  `colour_cell_substrate_screen_2026-08-26`.
- Assert the 7 unscreenable cells' strings contain the phrase "unvalidatable as written".
- The CLI routing/advertised-command guards are untouched by construction (no `summary` or route
  changes); their greenness is confirmed by the full-suite run in Verification, not re-asserted here.

### Step 5: Add the freeze guard and prove it non-vacuous
Files: tests/test_colour_cell_freeze.py
Depends on: Step 2, Step 3, Step 4
(extends the module Step 2 created; Step 2 precedes it, so there is no intra-wave file conflict)

**What changes:**
- Guard: the colour-cell set **derived** from `substrate_screen.collect()` and from `cli.TRAITS` must equal
  `FROZEN_COLOUR_ROUTES`. A 20th cell fails with a message naming the freeze rationale and how to ratify.
- **Non-vacuity test (mandatory):** monkeypatch a synthetic 20th cell into the derived set and assert the guard
  *fails* — without this the guard could pass while checking nothing, the exact defect caught twice already in
  this session (the vacuous control fixture; the substring retraction guard).
- Pin the three headline counts (65 / 40 / 14) in one place so the memo and code cannot diverge silently.

**Test strategy:**
- The guard itself, plus the negative-control test above.
- Full-suite run to confirm no interaction with `test_cli_dispatch` / `test_shipped_decoder_surface`.

### Step 6: Record the freeze where a future session will actually read it
Files: CLAUDE.md, LESSONS_LEARNED.md
Depends on: Step 5

**What changes:**
- Add a CLAUDE.md gotcha: the colour family is **frozen at 19**; adding a 20th trips
  `tests/test_colour_cell_freeze.py`; screen a candidate against **G9/G10** first. CLAUDE.md is loaded every
  session, which is the only place a freeze can actually prevent cell #20.
- One LESSONS_LEARNED bullet: a contract that names a promotion path must state a **sufficient** one — 7 cells
  promised "get a cohort" when a cohort could never have worked.
- Verify the cited paths resolve — `tests/test_claude_md_citations.py` enforces this.

**Test strategy:**
- `uv run pytest tests/test_claude_md_citations.py -q` green (cited paths exist).
- Full suite green.

## Execution Preview

- **Wave 0:** Step 1
- **Wave 1:** Step 2, Step 3, Step 4 *(parallel — three distinct files: `colour_cell_freeze.py`, `cell_registry.py`, `cli.py`)*
- **Wave 2:** Step 5
- **Wave 3:** Step 6
- Total waves: 4 · Max parallelism: 3 · Critical path: Step 1 → Step 2 → Step 5 → Step 6 *(toolkit-computed)*

## Risk Flags

- **Step 1 moves an import path that a test committed ~1 hour ago depends on.** Mitigated by updating
  `tests/test_colour_cell_substrate_screen.py` in the same step and asserting the emitted JSON is byte-identical
  — a pure move must change no output.
- **`incoming_data_gate` is DECLARED, not executed** (per its own field comment). Steps 3's strings are
  audit evidence, not enforcement; the plan must not be described as "gating" anything at runtime.
- **The freeze is attention/scope, not enforcement.** A guard test is ratifiable by editing one literal. That
  is intentional (it forces a deliberate decision) but must be stated so nobody reads it as a hard invariant —
  the same honesty rail as the self-init cap.
- **Hand-listing risk.** `FROZEN_COLOUR_ROUTES` is a literal; the *derived* cross-check in Step 5 is what stops
  it drifting. Without that test this reintroduces
  [[feedback_hardcoded_exclusion_list_undercovers]] for a fourth time.
- **Wording-only steps are easy to fake green.** Steps 3/4 change prose; their tests must assert against the
  live `verdicts()` derivation, not against copies of the same prose.

## Open Questions

1. **Should the 7 unscreenable cells be demoted below `KNOWLEDGE_BASELINE`?** No existing tier fits —
   `NO_FREE_SOURCE` is about *labels* (there is no free rabbit cohort *and* no recorded variant), and
   `NOT_CENSUSED` means "never scored". The plan keeps `KNOWLEDGE_BASELINE` and carries the distinction in
   `incoming_data_gate`. Adding a tier (e.g. `UNSCOREABLE_AS_WRITTEN`) touches the report-card state machine
   and is deliberately out of scope — **user call**.
2. **The curate branch is deliberately excluded.** Recording causal variants for the 40 unrecorded loci would
   convert 7 cells from unvalidatable to screenable, but it means writing curated biological facts into shipped
   catalogs — a fabrication hazard unless every locus is sourced from OMIA/literature and verified. It needs its
   own plan with that sourcing discipline pinned; this plan does not foreclose it.
3. **Should `dna-donkeycolor` / `dna-roedeercolor` (fully SNV-tractable) be exempted from the freeze** as
   genuinely promotable if a cohort appears? They pass G9 and G10 cleanly. Current plan freezes all 19 uniformly.

## Verification

1. `uv run python scripts/colour_cell_substrate_screen.py --self-check` → exit 0.
2. `git diff --stat wiki/colour_cell_substrate_screen_2026-08-26.json` → empty after Step 1 (pure move).
3. `uv run pytest tests/test_colour_cell_substrate_screen.py tests/test_colour_cell_freeze.py -q` → all pass.
4. Non-vacuity: the synthetic-20th-cell test fails the guard when the freeze roster is not updated.
5. `uv run python -c "from dna_decode.data import cell_registry as cr; assert not [c for c in cr.cells() if c.route.removeprefix('dna-') in COLOUR and c.incoming_data_gate=='n/a']"`.
6. `uv run pytest tests/ -q --ignore=tests/test_models_foundation.py` → 3762+ passed, 0 failed.
7. `sha256sum dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json` unchanged vs `4e601d1`.
8. Behaviour spot-check: `dna-decode rabbitcolor --loci A=A/a,B=B/b,C=C/C,D=D/d,E=E/e` byte-identical pre/post.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
