# Execution Log — Closed_Set_Organism_Identification_Router_V0_Plan

Date: 2026-10-04
Waves: 5 planned (max parallelism 2) — executed **SEQUENTIALLY**; parallel refused on a measured basis
  (`data/` is a directory junction into D:, so a git worktree gets a plain `data/` with 59 tracked
  files and none of the 88 cohort dirs that Steps 2/4/6 need)
Files changed: `dna_decode/data/organism_vocab.py`, `dna_decode/identify/{__init__,core,runner,cli,
  thresholds,__main__}.py`, `scripts/{build_identify_reference,identify_threshold_analysis,
  identify_validate}.py`, `scripts/contract_number_semantics.py`, `dna_decode/data/{cell_registry,
  cell_regime}.py`, `dna_decode/cli.py`, `pyproject.toml`, `CLAUDE.md`,
  `tests/{test_organism_vocab,test_build_identify_reference,test_identify_core,test_identify_cli,
  test_cli_dispatch}.py`,
  `wiki/identify_{reference_manifest,distance_distribution,outofset_probe,validation}_2026-10-04.json`,
  `wiki/certification_capstone.{md,json}`
Sentrux verdict: n/a — sentrux not installed
Commit: `b57a092` (final) — steps at `131901a`, `51b90b8`, `28f18e5`, `ffffa6a`, `f3200dc`

## Result

`dna-identify` ships: 10 identifiable organisms, leave-one-out **211/211**, out-of-set abstention
**44/48 = 0.9167**. Registry 115 → 116 cells, 44 → 45 traits. Thresholds `max_distance=0.1036` /
`ambiguity_margin=0.015`, frozen AFTER the distance analysis. Frozen five byte-unchanged; prospective
lock re-verified `ok=True`.

## Deviations (6)

1. **Sequential, not parallel** — the `data/` junction makes worktrees unusable for the steps that read
   the genome corpus.
2. **ONE routing token, not two.** The plan's central claim — *"a router emitting one string satisfies
   neither consumer"* — was FALSE. `calibrated_rule_for`'s genus-prefix fallback makes the AMRFinder
   `-O` value resolve the registry key as well; verified identical across 6 organism/drug pairs in both
   the rule-resolved and None directions. The second field was dropped and the equivalence pinned.
3. **Label source = each assembly's own GenBank `ORGANISM` field**, not the cohort directory name.
   Per-genome, so a mislabelled cohort cannot propagate to every genome under it.
4. **De-duplication moved from reference CONSTRUCTION to VALIDATION** — duplicates change no call at
   deployment; they only inflate leave-one-out, which is a validation concern.
5. **Step 5 SPLIT**: runner early (Steps 4 and 6 need it), registration late (with Step 7's contract),
   so the suite never goes red between them and no `NOT_CENSUSED` placeholder ever shipped.
6. **Reference capped at 25/organism** after system memory pressure killed the 2,026-genome build.
   Also the better design: largest-organism share 0.6807 → 0.1179.

## Plan flaws found by executing it

- **Step 1 as specified silently lost 160 in-set genomes.** AMRFinder's `Campylobacter`/`Salmonella`
  are GENUS-level `-O` values, so species-strict aliases could not see `Campylobacter_jejuni` (66),
  `Salmonella_enterica` (60), `Campylobacter_coli` (34) — while correctly excluding 48 genuinely
  different species. Fixing it discharged the plan's **Open Question 2** by measurement: adopt
  AMRFinder's own granularity as the router's.
- **Step 6 as specified rebuilt the sketch per held-out genome** (~212 container spin-ups) for
  information one cached all-pairs `mash dist` already contains.
- **Tier is KNOWLEDGE_BASELINE, not the plan's FAITHFUL_TO_TOOL.** There is no tool comparator, and the
  NCBI `ORGANISM` label may itself have been assigned BY sequence comparison — the same evidence class
  this router uses — so the measurement is IN-DISTRIBUTION (gate **G1**).

## Defects found in my own work by running it (4)

- `--list-organisms` advertised **14** organisms while the reference covers **10**, promising coverage
  that does not exist (the `--help` over-claim pattern). Now splits IDENTIFIABLE from
  DECLARED-but-abstain-only.
- I invented four `cell_registry` field semantics (prose into `falsifier_ref` and
  `incoming_data_gate`, a bool into `native_abstention`, a tuple into `abstention_vocab`). The existing
  guards caught every one; the real shapes are a PATH, a `G1..G10` list, a term string, and a single
  enum member.
- Two of three `metric_bindings` named quantities `LABEL_VOCAB` cannot express (a raw count, a
  distance). Dropped rather than inventing vocabulary to make them typecheck; both numbers remain
  verified present (6/6, 0 of 208 decoys).
- A CRLF bug spent a whole build cycle **disguised as a mount failure**: Windows `write_text` turned
  `\n` into `\r\n`, mash tried to open `genome.fna\r`, and the stray CR mangled mash's own error text.
  Confirmed by the fix removing exactly 2,026 bytes — one per genome.

## Latent defect fixed in a shared checker

`scripts/contract_number_semantics.py::path_is_consistent` did `LABEL_VOCAB[canon]`, so ONE cell
binding a novel quantity `KeyError`'d the **entire** semantics audit. It now reports instead of raising.

## Two measured limits, both in the contract, neither fixable by tuning

1. **211/211 is OPTIMISTIC** — held out by GENOME, not by LINEAGE, and **146 of 211 (69.2%)** have a
   ~99.5%-ANI-or-closer twin still in the reference (median nearest same-organism distance 0.00165).
2. **4 out-of-set false calls are Klebsiella CONGENERS** (K. variicola ×3 at 0.05095, K. michiganensis
   ×1 at 0.06746) sitting BELOW the in-set ceiling 0.09880, so no threshold rejects them without
   rejecting legitimate in-set genomes. Lowering `max_distance` is REFUSED: it would trade a named
   blind spot for silent coverage loss.

## Open falsifier (named, not closed)

A **LINEAGE-disjoint hold-out** — the stronger test of limit (1). The cached all-pairs matrix at
`D:/dna_decode_cache/identify/allpairs.tsv` makes it cheap: cluster by distance and hold out whole
clusters rather than single genomes.

## Also still open

Four vocabulary organisms have **no reference coverage** and can therefore only abstain:
`candida_auris`, `enterococcus_faecium`, `mycobacterium_tuberculosis`, `streptococcus_pneumoniae`.
`--list-organisms` names them explicitly rather than implying coverage.
