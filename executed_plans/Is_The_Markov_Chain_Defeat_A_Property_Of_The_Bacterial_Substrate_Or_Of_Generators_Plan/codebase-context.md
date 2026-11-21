- `dna_decode/glm/baselines.py` — `MarkovGenerator(order, seed).fit(seqs: list[str])` then
  `.generate(n, length)`, and `UniformGenerator(seed).generate(n, length)`. **Both take a plain
  sequence list**, so neither needs a genome, a GFF3 or a reference download.
- `dna_decode/glm/falsifier.py` — `heldout_auroc(a: list[str], b: list[str], *, k=4, seed=0, mode,
  bins)` and `self_control_ceiling(natural: list[str], *, k=4, seeds=(0,1,2,3,4), mode, bins)`
  returning `(max_observed, all_values)`, plus `alphabet_and_length_ok(seqs, length)`. All pure,
  list-in, CPU-only. `self_control_ceiling` takes the MAXIMUM over seeds *"so a lucky split cannot
  manufacture a finding"*, and is re-derived per feature mode because *"a richer feature set has more
  capacity to fit split noise."*
- `scripts/glm_generator_falsifier.py` — the CLI. `build_natural(length, limit, seed, kind)` at line
  38 accepts `kind` in `{promoter, genomic}` but BOTH branches go through `load_fasta(genome)`, so the
  natural set is organism-locked by the genome file rather than by `kind`. A human natural set needs a
  third path that takes sequences directly.
- `dna_decode/glm/corpus.py` — builds the bacterial natural corpus from *"the committed MG1655
  reference plus its GFF3"*. Not reusable for human; not modified by this plan.
- `dna_decode/glm/lentimpra.py` — `load_cell(cell, data_dir)`, `load_design`, `load_activity`,
  `usable_cells()`. Sequences are **230 bp** with the shape *"15nt 5' adaptor - 200nt element - 15nt
  3' adaptor"*, and the authors' **10-fold assignment** is carried through. Usable: HepG2 61,463 and
  K562 98,332 joined pairs (join rate 1.000); **WTC11 REFUSES**, and the live reason is a join-rate
  bar rather than the gene-keying an earlier draft paraphrased: *"WTC11 joins at 0.0000 (0 of 23064
  measured rows), below the 0.95 bar."*
- **THE SUBSTRATE IS ALREADY ON DISK, and an earlier draft of this plan said otherwise.**
  `data/mpra_human/` holds 6 files, ~72 MB, dated 2026-11-16: `Table S3 - large-scale lib
  design.xlsx` (56.2 MB), `Table S6 - folds and performance.xlsx` (6.1 MB), `Table S7`, `Table S9`,
  plus `_fetch_manifest.json` and `_zenodo_record.json`. `usable_cells('data/mpra_human')` returns
  **HepG2 usable=True, K562 usable=True**.
- **How the earlier draft got it wrong, recorded because the mistake recurred three times in one
  session.** `usable_cells()` called with NO argument defaults to `DEFAULT_DIR =
  D:/dna_decode_cache/mpra_human`, which is on the detached drive, so it reports every cell
  unusable. Passing the real path works. The lesson is to check where data actually LIVES rather
  than where a default points; `CLAUDE.md` prescribes `scripts/data_inventory.py` for exactly this,
  and the committed `wiki/data_inventory.md` is stale (dated Oct 9, lists only `mpra` at 104.7 MB,
  has no `mpra_human` row).
- `scripts/fetch_lentimpra.py` — the MD5-verifying fetcher, already run on 2026-11-16; its manifest
  is on disk and is what Step 1 now verifies rather than re-running.
- **THERE IS NO HUMAN CEILING, measured 2026-11-16** (`wiki/human_frame_shape_2026-11-16.json`,
  status MEASURED): the deposit has no per-replicate columns, so no Spearman-Brown step-up is
  derivable. That artifact reports it as `null` *"rather than filled with the BACTERIAL 0.941, which
  both wiki/eukaryotic_arm_result_*.json currently carry as context"*. This plan's margin needs only
  the SELF-CONTROL ceiling, which is derivable, so the design survives — but the human side has no
  0.941 equivalent and none may be borrowed.
- `openpyxl>=3.1` is in the DEFAULT dependencies (declared 2026-11-19), which the xlsx loader needs.
- Bacterial anchors for the contrast, from `wiki/glm_generator_paired_sweep_2026-11-20.json`:
  Markov null mode-both **0.6679** at seed 7 (it is refit per seed and spans 0.5674-0.7024 across the
  grid), self-control ceiling mode-both **0.5865**, kmer **0.5447**, uniform control **0.9641**.

### Reusable-Code Survey
- `MarkovGenerator` / `UniformGenerator` — reused verbatim; they are already substrate-agnostic.
- `falsifier.heldout_auroc` + `self_control_ceiling` + `alphabet_and_length_ok` — reused verbatim as
  the measurement; no new statistics are written.
- `lentimpra.load_cell` + the authors' 10-fold assignment — reused as the natural set and as the
  fit/eval split, so the split is PUBLISHED rather than invented here.
- `scripts/fetch_lentimpra.py` — reused for the fetch; no new download code.
- The bacterial anchors are reused as committed numbers, not recomputed.
- Searched: `dna_decode/glm/`, `scripts/glm_generator_falsifier.py`, `scripts/fetch_lentimpra.py`,
  `tests/`; no `graphify-out/GRAPH_REPORT.md` on this host.
