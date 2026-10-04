## Lens status
Execution mode target: parallel (5 waves, max parallelism 2)
Repo-index: absent (`update-index` CLI not installed) — Step 1 research done by direct file reading
sentrux: not installed — Architecture Gate records `n/a`
gh: available | codex: available
project-rules.md: absent — `## Project Rules Check` omitted
External surfaces verified LIVE this session (not inherited): `mash` 2.3 subcommand list + `mash dist` / `mash screen` output fields from the pinned biocontainer; `amrfinder -l` organism allow-list from the pinned ncbi/amr image

## Problem Statement

Step 1 of the project's north star does not exist: **given a DNA sequence, identify what organism it is.** All 44 CLI traits REQUIRE the caller to name the organism (`--organism`, `-O`), so the pipeline's first step is currently a human typing the answer in. A taxonomy-keyword scan of `dna_decode.cli.TRAITS` returns `NONE`.

Two measured obstacles shape the design, both found during Step 1 research:

1. **There is no canonical organism vocabulary.** Across the 115 registered cells, E. coli appears as three distinct tokens — `Escherichia_coli_Shigella`, `Escherichia_coli`, `escherichia_coli`.
2. **There are TWO consumer vocabularies and they disagree by construction.** `amrfinder -l` on the pinned image accepts its own 32-value list (`Escherichia`, `Klebsiella_pneumoniae`, `Campylobacter`, `Salmonella`, ...), while the registry/rule layer takes different tokens (`Escherichia_coli_Shigella`, `Klebsiella`). `scripts/external_cohort_revalidate.py` already carries them as two separate parameters (`amrfinder_organism`, `registry_organism`). A router emitting one string satisfies neither consumer.

**Scope decision, the load-bearing one:** this is a **CLOSED-SET router over the organisms this tool actually supports, with abstention** — NOT open-world taxonomy. The supported set is finite and enumerable from `cell_registry` (41 organism tokens across 115 cells). Open-world identification would need a roughly 1 GB RefSeq reference sketch and would answer a question no consumer asks; closed-set routing with a hard ABSTAIN on anything out-of-set is both smaller and the honest shape, because abstention is the safe failure for a router whose output selects a downstream rule.

**Non-goals:** strain or subspecies resolution; metagenomic mixtures; open-world taxonomy; read (FASTQ) input; changing any decoder's behaviour; any edit to the five sha256-pinned frozen files.

## Codebase Context

**The four-place trait contract** (CLAUDE.md, verified in Step 1): a new CLI trait must touch (a) a `pyproject.toml` console entry (pattern `dna-plasmid = "dna_decode.plasmid.cli:main"`, line 87), (b) `dna_decode/cli.py::TRAITS` (dict at line 29; each value is `{summary, validation}`), (c) dispatch, (d) a `cell_registry` `CellContract`. Four existing guards fire on a new trait and MUST pass without being edited:

- `tests/test_cli_dispatch.py::test_traits_registry_matches_console_entries` (hand-maintained pin)
- `tests/test_cli_dispatch.py::test_every_trait_has_an_evidence_contract`
- `tests/test_advertised_commands.py` — every advertised command resolved through its REAL argparse parser
- `dna_decode/data/cell_regime.py::regime_for_route` **RAISES `UndeclaredRegime`** on an unknown route (line 84), so the new route must be declared or registration is a hard failure

**Module layout to match** (`dna_decode/mlst/`): `__init__.py`, `__main__.py`, `cli.py`, `core.py`, `runner.py`.

**Mash external surface, verified live rather than inherited:**

- image pinned at `quay.io/biocontainers/mash:2.3--hb105d93_10` (`dna_decode/eval/phylogeny.py:30`)
- subcommands present: `bounds`, `dist`, `info`, `paste`, `screen`
- `mash dist [options] <reference> <query>` emits fields `[reference-ID, query-ID, distance, p-value, shared-hashes]`
- `mash screen [options] <queries>.msh <mixture>` emits `[identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]`; `-w` is winner-takes-all; documented for read or contig **mixtures**
- **`dist` is the correct subcommand for an assembled genome; `screen` targets mixtures.** v0 takes assemblies, so `dist` against a reference sketch.

**AMRFinder external surface, verified live:** `amrfinder -l` prints `Available --organism options:` as a comma-separated list (32 values on DB `2026-03-24.1`). This is the AUTHORITATIVE allow-list and must be parsed and pinned, never hand-copied.

**Reference material is LOCAL, so v0 needs no download.** `D:/dna_decode_cache/refseq` holds 1,129 cached genomes and `data/raw/` holds 80 cohort directories covering acinetobacter, campylobacter, C. auris, E. coli/Shigella, Enterobacter cloacae, E. faecium, N. gonorrhoeae, Klebsiella, P. aeruginosa, Salmonella, S. aureus, S. pneumoniae and M. tuberculosis. No taxonomy reference sketch exists yet; the two `.msh` files on D: are TB-cohort artifacts, not references.

**Standing host rules:** every heavy or extra file goes on **D:**, never C: (C: runs 85-99 percent full). Direct `docker run` from Git Bash REQUIRES `MSYS_NO_PATHCONV=1`; `subprocess.run` through `tools/docker_runner.run` is unaffected because there is no shell.

### Reusable-Code Survey

- **`dna_decode/eval/phylogeny.py::compute_mash_distances`** — REUSE the Docker-routing pattern and the staging logic but NOT the function: it computes an all-pairs self-distance matrix (`mash sketch` then `mash dist sketch.msh sketch.msh`), whereas identification needs `mash dist <reference.msh> <query.fna>`. New thin runner, same image pin, same `use_docker` posture.
- **`tools/docker_runner.py::run`** — REUSE directly as the container boundary; it wraps FileNotFoundError and TimeoutExpired as `DockerRunnerError` and runs without a shell, so it is immune to the MSYS path trap.
- **`dna_decode/data/cell_registry.py`** and **`dna_decode/data/cell_regime.py`** — REUSE as the registration surface; both are NON-frozen.
- **`dna_decode/eval/clonality.py::greedy_representative_clusters_from_matrix`** — REUSE for reference de-duplication if the sketch is dominated by near-identical cohort genomes; it is chaining-resistant, unlike `phylogeny.cluster_by_ani`.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `dna_decode/eval/`, `tools/`, `dna_decode/{mlst,serotype,salmserovar}/`, `dna_decode/data/`, `scripts/`.

## Pre-Change Baseline

- **Organism identification does not exist.** A taxonomy-keyword scan of `dna_decode.cli.TRAITS` returns `NONE`: 44 traits, 0 that identify the organism.
- Every AMR and typing call path requires an explicit organism; `call_resistance(organism=...)` and AMRFinder `-O` are both supplied by the caller.
- 115 registered cells across 41 distinct organism tokens, including three spellings of E. coli.
- Full suite GREEN per CLAUDE.md (3875 passed / 11 skipped / 0 failed); re-measure before the first edit.
- No reference sketch, no canonical-vocabulary module, no `dna-identify` console entry.

## Verification Signal

The plan succeeds iff ALL of:

1. `dna-identify --genome-fasta X.fna` returns a canonical organism **plus both consumer tokens** (`amrfinder_organism`, `registry_organism`) or an explicit `ABSTAIN` — never a bare species string.
2. Every emitted `amrfinder_organism` is a member of the live `amrfinder -l` list, asserted by a test that PARSES that list rather than hard-coding it.
3. Every emitted `registry_organism` resolves to at least one `cell_registry` cell.
4. **Leave-one-out accuracy is reported PER ORGANISM, never pooled alone** — the reference set is dominated by E. coli and Klebsiella, so a pooled number measures cohort composition.
5. **Out-of-set control:** genomes from a species absent from the reference sketch ABSTAIN at a measured rate, and a router that never abstains is rejected regardless of accuracy.
6. The decision thresholds are frozen in a committed file **after** the distance-distribution analysis and **before** the validation run.
7. The five frozen files are byte-unchanged, `verify_lock` still returns ok, and the full suite is green.

## Implementation Steps

### Step 1: Canonical organism vocabulary with mappings to both consumer vocabularies
Files: dna_decode/data/organism_vocab.py, tests/test_organism_vocab.py
Depends on: none

**What changes:**
- New NON-frozen module defining one canonical organism token per supported organism, each carrying `amrfinder_organism` (or `None` where AMRFinder has no entry, for example C. auris and M. tuberculosis) and `registry_organism`.
- `canonical(token)` normalizes the three observed E. coli spellings onto one value; `consumers(canonical)` returns both downstream tokens.
- The `amrfinder -l` allow-list is PARSED at test time and every non-None `amrfinder_organism` asserted a member; the list is never hand-copied into source.
- Deliberately does NOT rename anything in `cell_registry` — this is an additive translation layer, so no existing cell changes behaviour.

**Test strategy:**
- Unit: all three E. coli spellings round-trip to one canonical value; an unknown token RAISES rather than defaulting (a silent default is how a wrong organism would reach a rule).
- Live-gated: parse `amrfinder -l` from the pinned image and assert membership; SKIP with an explicit reason when Docker is unavailable, never silently pass.
- Non-vacuity: assert the map is not constant (at least two distinct `amrfinder_organism` values) so the module cannot degenerate into ceremony.

### Step 2: Reference-sketch builder over local cached genomes
Files: scripts/build_identify_reference.py, tests/test_build_identify_reference.py
Depends on: none

**What changes:**
- Enumerate local genomes carrying a known organism label (cohort directory to organism), de-duplicate with `clonality.greedy_representative_clusters_from_matrix` so near-identical cohort genomes do not dominate, and emit one `mash sketch` to `D:/dna_decode_cache/identify/reference.msh` per the standing D: rule.
- Write a committed manifest sidecar `wiki/identify_reference_manifest_<date>.json` recording per-organism genome counts, the sketch k and s parameters, and the source cohort of every reference genome. The sketch itself is a regenerable D: artifact and is NOT committed.
- REFUSE to write a sketch covering fewer than 2 organisms, since a one-organism reference cannot discriminate.
- Record the composition imbalance explicitly in the manifest; it is a known caveat, not a defect to hide.

**Test strategy:**
- Unit on synthetic FASTAs: label derivation, the fewer-than-2-organism refusal, manifest shape.
- Integration, Docker-gated and skipped with a reason: build a 3-organism toy sketch and assert `mash info` reports the expected genome count.
- Assert the sketch path resolves under D: and never under the repository.

### Step 3: Identifier core, pure parsing plus the abstention decision
Files: dna_decode/identify/__init__.py, dna_decode/identify/core.py, tests/test_identify_core.py
Depends on: Step 1

**What changes:**
- `parse_mash_dist(stdout)` as a PURE function over the verified five-field format `[reference-ID, query-ID, distance, p-value, shared-hashes]`.
- `decide(hits, thresholds)` returns the canonical organism plus both consumer tokens, or `ABSTAIN` with a machine-readable reason drawn from `no_hit`, `above_max_distance`, `ambiguous_top2`, `reference_unavailable`.
- Thresholds are **PARAMETERS**, not literals, so the module is complete and fully testable before any value is chosen — the pattern the ECOFF seam established for a blocked constant.
- The returned record carries `top_hits` so a call is auditable, and never carries a strain or subspecies claim.

**Test strategy:**
- Unit: synthetic `mash dist` stdout reaching each abstention reason; a tie inside the ambiguity margin ABSTAINS rather than being resolved by sort order, which is the documented tie-breaking trap.
- Assert `decide` can never return a canonical organism absent from Step 1's vocabulary.
- Assert empty or garbage stdout ABSTAINS rather than raising, so a broken container cannot read as a confident call.

### Step 4: Distance-distribution analysis, and only then freeze the thresholds
Files: scripts/identify_threshold_analysis.py, dna_decode/identify/thresholds.py, tests/test_identify_thresholds.py
Depends on: Step 2, Step 3

**What changes:**
- Measure within-organism versus between-organism `mash dist` distributions over the reference set and emit `wiki/identify_distance_distribution_<date>.md` plus its JSON sidecar.
- **Only then** write the frozen thresholds `MAX_DISTANCE` and `AMBIGUITY_MARGIN` into `thresholds.py`, citing the analysis artifact inline and stating the mechanism each threshold polices.
- The ordering is deliberate. This repo has two recorded cases where the pre-registered bar was itself the error: a NET ceiling used to police a DIRECTIONAL failure, and a `+4 hits` bar that was unachievable by construction. The lesson on record is **do the mechanism analysis, THEN freeze** — not "do not freeze".
- Record the separation, or its absence, rather than tuning until it looks clean. Overlapping distributions are a finding, and the cell then ships abstention-heavy.

**Test strategy:**
- Unit: the constants are in range and the module exposes no alternative tunable path.
- Pin that `core.decide` is called WITH these constants at the single CLI call site, so a future constant change cannot be silently overridden — the documented shipped-default trap where a CLI `default=80.0` overrode a constant lowered to 40.
- Assert the analysis artifact cited inside `thresholds.py` actually exists.

### Step 5: CLI entry point and the four registration points
Files: pyproject.toml, dna_decode/cli.py, dna_decode/identify/cli.py, dna_decode/identify/__main__.py, dna_decode/identify/runner.py, dna_decode/data/cell_regime.py, tests/test_identify_cli.py
Depends on: Step 1, Step 3

**What changes:**
- `runner.py` invokes `mash dist <reference.msh> <query.fna>` through `tools/docker_runner.run`, with a native-binary path mirroring the posture in `phylogeny.py`.
- `cli.py` accepts `--genome-fasta` (required), `--reference` (defaulting to the D: sketch) and `--json`; argparse defaults DERIVE from the `thresholds.py` constants rather than restating them.
- Register all four places: the `pyproject.toml` console entry `dna-identify`, a `cli.py::TRAITS["identify"]` row, dispatch, and a `cell_regime` declaration for route `dna-identify`, without which `regime_for_route` raises.
- A missing reference sketch exits non-zero with an actionable message naming Step 2's builder, never a silent "unidentified".

**Test strategy:**
- No edit to the advertised-command guard is required and that is verified, not assumed: its `_console_scripts()` helper DISCOVERS routes by reading `pyproject.toml [project.scripts]`, so adding the console entry automatically brings `dna-identify` under it and a non-parsing `--help` example fails without anyone wiring it up. New tests therefore cover only this cell's own surface.
- In `tests/test_identify_cli.py`, resolve each advertised example through the cell's REAL argparse parser (patch `parse_args` to raise on a SUCCESSFUL parse, and patch `parse_args` ONLY — also patching `parse_known_args` skips the unrecognized-argument error and lets a bad example pass).
- Assert `regime_for_route("dna-identify")` no longer raises.
- Assert the four pre-existing guards pass UNEDITED.
- Assert each argparse default equals its `thresholds.py` constant and is not a hardcoded literal.

### Step 6: Leave-one-out validation with a mandatory out-of-set abstention control
Files: scripts/identify_validate.py, tests/test_identify_validate.py
Depends on: Step 2, Step 3, Step 4

**What changes:**
- Leave-one-genome-out over the reference set: rebuild the sketch without the held-out genome, identify it, score. Holding out is mandatory, because scoring a genome against a sketch that contains it is circular and returns a perfect match by construction.
- Report accuracy **PER ORGANISM** with each per-organism n, and emit the pooled figure only beside the composition table that explains it.
- **Out-of-set control:** identify genomes from a species deliberately excluded from the sketch and report the abstention rate. The script REFUSES a verdict when the out-of-set group is empty, and FAILS the run when out-of-set abstention is zero percent, because a router that never abstains is not usable as a router.
- Emit `wiki/identify_validation_<date>.md` plus its JSON sidecar, counting every abstention reason separately from errors.

**Test strategy:**
- Unit: the circularity guard (a held-out genome is absent from its own reference), the empty-out-of-set refusal, the zero-percent-abstention failure, and per-organism aggregation arithmetic reconciling with the pooled total.
- Assert an identification ERROR is counted separately from an ABSTAIN, so a crash cannot read as a cautious call.

### Step 7: Evidence contract written from the validation artifact
Files: dna_decode/data/cell_registry.py, tests/test_identify_contract.py
Depends on: Step 6

**What changes:**
- Add the `CellContract` for the identify cell with every number read FROM `wiki/identify_validation_<date>.json` and never from memory. The `forward` cell shipped for a day registered `NOT_CENSUSED` while being DMS-validated, and under-claiming is as much a trust-surface falsehood as over-claiming.
- Set the tier to `FAITHFUL_TO_TOOL` because the labels are cohort-directory-derived provenance rather than wet-lab identification, and record in `honest_limits` the closed-set scope, the reference-composition imbalance, species-level-only resolution and assembly-only input.
- Add a `metric_bindings` entry so `scripts/contract_number_audit.py` can verify each cited number is present in the artifact it names.

**Test strategy:**
- Assert every number in the contract prose appears in the cited artifact, importing the existing audit helpers rather than reimplementing the matcher (whole-number matching on BOTH the cited and decoy haystacks; a raw substring test is the retired defect that made `2.1` match inside `72.1`).
- Assert the cited artifact exists, via the existing `test_every_wiki_artifact_a_contract_cites_actually_exists` guard.
- The standing contract-number audit is read-only and exits 0 always, so it is INVOKED as a verification action (see `## Verification`) and is NOT modified by this step — no new tests are added to it.

## Execution Preview

| Wave | Steps | Parallelism |
|---|---|---|
| 0 | Step 1, Step 2 | 2 |
| 1 | Step 3 | 1 |
| 2 | Step 4, Step 5 | 2 |
| 3 | Step 6 | 1 |
| 4 | Step 7 | 1 |

Total waves: 5. Max parallelism: 2. Critical path: Step 1, Step 3, Step 4, Step 6, Step 7 (five steps).

Intra-wave file overlap: none. `cell_registry.py` is touched only by Step 7 and `cell_regime.py` only by Step 5, which sit in different waves from each other and from every sibling.

## Risk Flags

- **The reference set is composition-biased and that bounds every number it produces.** It is built from AMR and typing cohorts, so it is dominated by E. coli and Klebsiella and over-represents clinical isolates. Per-organism reporting in Step 6 is the mitigation; a pooled accuracy from this set would measure cohort composition, the documented "pooled ranking cannot answer a within-group question" failure.
- **Cohort-directory labels are PROVENANCE, not wet-lab identification.** A genome under `data/raw/klebsiella_*` is *recorded* as Klebsiella. The cell therefore cannot reach `INDEPENDENT_MEASURED` on this substrate, which is why Step 7 fixes the tier at `FAITHFUL_TO_TOOL`. Promotion would need an independently-identified cohort.
- **A mislabelled reference genome propagates silently.** Mitigation: Step 6's per-organism table surfaces a systematically misidentified organism as a block rather than as noise. No automatic relabelling is proposed.
- **Docker dependency.** Mash runs in a container, and a wedged Docker Desktop D: mount is a known failure on this host with `wsl --shutdown` as the recorded fix. The CLI must exit non-zero on an unavailable reference or runner rather than returning "unidentified".
- **Threshold overlap is a possible outcome, not a bug.** If within- and between-organism distances overlap for closely related taxa — Klebsiella species, or E. coli versus Shigella, which AMRFinder itself merges into a single `-O Escherichia` value — the honest result is a cell that abstains on those pairs. Do not tune until the overlap disappears.
- **`mash screen` is deliberately NOT used in v0** because it targets mixtures. The FASTQ and read path is a named non-goal, so a later read-input request is new work rather than a gap in this plan.
- Architecture Gate: `n/a — sentrux not installed`.

## Open Questions

1. **Reference scope:** local cohorts only (v0 as planned, roughly 13 taxa, zero download) versus adding Mash's prebuilt RefSeq sketch (about 1 GB onto D:, open-world capable). The closed-set choice is deliberate and reversible, and adding RefSeq later widens coverage without changing the interface.
2. **Abstention policy at the E. coli and Shigella boundary:** AMRFinder collapses them into one `-O Escherichia`, so the router arguably SHOULD NOT try to separate them. Proposal: adopt the AMRFinder allow-list granularity as the router's target granularity, which dissolves the question. Needs ratification.
3. **Viral and human targets are outside v0 scope.** A Mash sketch over bacterial and fungal genomes will not route HIV, SARS-CoV-2 or HCMV (26 plus 3 plus 5 cells) nor the human clinical cells; those need a different method, and per-cell reference-CDS BLAST already exists. Confirm that a bacterial-and-fungal-only v0 is the right first slice.
4. Whether the canonical vocabulary should eventually REPLACE the three E. coli spellings inside `cell_registry` (a rename touching many cells) or remain a translation layer permanently. v0 assumes the translation layer.

## Verification

1. `uv run pytest tests/ -q` — full suite green with zero failures, new tests included.
2. `uv run python scripts/build_identify_reference.py` produces the sketch on D: plus the committed manifest; assert the sketch is NOT inside the repository.
3. `uv run python scripts/identify_threshold_analysis.py` produces the distribution artifact, and the thresholds are committed citing it.
4. `uv run python scripts/identify_validate.py` produces per-organism accuracy plus the out-of-set abstention rate, and the run FAILS when out-of-set abstention is zero.
5. `uv run python scripts/project_status.py` shows the trait count moving 44 to 45 and the new cell present in the live registry.
6. `uv run python scripts/contract_number_audit.py` shows the new cell present, not LOW-discrimination, and the residual unchanged.
7. Frozen-surface check: `git diff --quiet` over the five pinned files AND an assertion that each exists and is non-empty, because `git diff --quiet` on a missing path exits 0 and tests nothing (the documented vacuous-guard trap); then re-run `verify_lock`.
8. `dna-identify --help` plus every advertised example resolved through the real parser.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
