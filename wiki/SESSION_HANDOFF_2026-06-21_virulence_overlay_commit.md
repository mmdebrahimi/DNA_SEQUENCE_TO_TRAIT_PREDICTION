# Session handoff — 2026-06-21: COMMIT + PUSH the finished virulence-overlay (fresh session)

Everything a cold session needs to land the genome-map virulence-determinant overlay. Repo:
`C:\Users\Farshad\PythonProjects\dna_decode` (origin `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`).

## ⇒ START HERE (the one thing that's owed)

**The virulence-overlay feature is FULLY IMPLEMENTED, TESTED, and VERIFIED — but 100% UNCOMMITTED.** The
prior session ran `/execute-plan` to completion (all 5 steps + post-exec docs) but never committed. The
work sits in the working tree on branch **`feat/genome-map-virulence-overlay`**, which is at the SAME
commit as `main` (`9b60b91`; `git rev-list --left-right --count main...HEAD` = `0 0` — no commits on the
branch, just uncommitted changes).

1. `cd C:\Users\Farshad\PythonProjects\dna_decode`
2. Verify it's sound (should be green — see **Verification** below; ~90s):
   ```
   .venv/Scripts/python.exe -m pytest tests/test_genome_map_virulence_overlay.py tests/test_vf_runner_coords.py tests/test_pathotype_vf_diff.py tests/test_genome_map_cli.py -q
   ```
   (use `.venv/Scripts/python.exe -m pytest`, NOT `uv run` — uv's cache has been wedging with `os error
   183` this week; the venv works directly.)
3. **Commit to `main`** (the project's standing convention = commit straight to main, the cross-machine
   sync channel). Cleanest path, since the branch == main: carry the working tree to main and commit there.
   ```bash
   git checkout main                 # working-tree changes carry over cleanly (branch==main commit)
   git add <the FEATURE files only — see the list below; do NOT git add -A>
   git commit -m "feat(genome-map): virulence-determinant overlay tier (v2, brainstorm-hardened) ..."
   git push origin main
   git branch -d feat/genome-map-virulence-overlay   # delete the now-redundant branch
   ```
   (If `git checkout main` complains about an overwrite, the branch and main have diverged after all —
   fall back to: commit on the feature branch, then `git checkout main && git merge --ff-only
   feat/genome-map-virulence-overlay && git push`.)

## ⚠️ Stage ONLY the feature files — NOT `git add -A`

The working tree has pre-existing dirty files that must be LEFT (same set the project has carried for
weeks). Stage exactly these (the virulence-overlay work + its docs):

**Modified:** `CLAUDE.md` · `dna_decode/genome_map/__init__.py` · `dna_decode/genome_map/build_map.py` ·
`dna_decode/pathotype/vf_runner.py` · `scripts/genome_map.py` · `scripts/genome_map_spike.py` ·
`tests/test_genome_map_cli.py` · `tests/test_pathotype_vf_diff.py` · `project_state/dna-decode-2026-05-11.md` ·
`wiki/decisions-log.md` · `wiki/genome_map_usage.md` · `wiki/plans-index.md` · `.claude/testable-modules.md`
**New (untracked):** `dna_decode/genome_map/virulence_overlay.py` · `tests/test_genome_map_virulence_overlay.py` ·
`tests/test_vf_runner_coords.py` · `executed_plans/Genome_Map_Virulence_Overlay_Plan/execution-log.md`
**Rename (already staged as R):** `plans/Genome_Map_Virulence_Overlay_Plan/technical-plan.md → executed_plans/Genome_Map_Virulence_Overlay_Plan/technical-plan.md`
  - ⚠️ The rename is the plan being archived to `executed_plans/`. `git add executed_plans/Genome_Map_Virulence_Overlay_Plan/` + `git rm --cached plans/Genome_Map_Virulence_Overlay_Plan/technical-plan.md` if the rename isn't already staged. (`git status` showed it as `R ` = staged-rename, so it should travel with the checkout.)

**LEAVE dirty (do NOT stage):** `uv.lock` · `wiki/ciprofloxacin_mechanism_audit_2026-06-05.{json,md}` ·
`bash.exe.stackdump` · `research_outputs/eukaryotic-...unsupported.md` · the untracked stale
`plans/TB_AMR_Decoder_CRyPTIC_Technical_Plan.md` (superseded; safe to `rm` if you like — NOT this feature's).

## What shipped (1 paragraph)

A **5th genome-map overlay tier — `virulence-determinant`** (`dna_decode/genome_map/virulence_overlay.py`).
Where a curated VirulenceFinder (VF) allele is PRESENT in ONE E. coli/Shigella genome, it's surfaced
behind the SAME coordinate-join integrity gate + a presence-only wall as the AMR `determinant-phenotype`
tier (presence of a curated determinant, NEVER a learned pathogenicity claim). The deterministic
pathotype-resolver call is shown SEPARATELY as a genome-level overlay (`genome_pathotype_call`, the
virulence analog of the AMR R/S call), QC-gated so a low-QC genome yields `AMBIGUOUS_LOW_QC` not a
confident commensal. E. coli/Shigella-scoped v1; offline-degrades to `UNAVAILABLE_NO_BLASTN` without
`blastn`. Plan: `executed_plans/Genome_Map_Virulence_Overlay_Plan/`. **Frozen AMR surface byte-unchanged.**

## What the prior session VERIFIED (per execution-log.md — do not re-litigate)
- **41 targeted tests pass** (virulence_overlay + vf_runner_coords + vf_diff + cli); execution log records
  baseline **1472 → 1500** (28 new), **0 regressions** (excl. `tests/test_models_foundation.py` host-torch limit).
- **Frozen AMR surface byte-unchanged** (verified empty `git diff` on `amr_rules.py` + `calibrated_amr_rules.json`).
- **LIVE run** (native blastn + committed VF DB, D:-free for the VF part) on cached E. coli ST131
  `GCA_002180195.1`: `virulence_status=FULL`, **27 virulence-determinant features** (2690/2691 high-confidence
  coord joins, 0 symbol-fallback), all-called-allele coverage incl. unclustered + tandem copies, genome
  pathotype call `ExPEC_COMPATIBLE/LOW_CONFIDENCE` (separate), DB sha `e94e6c6d4dae1ca6`; AMR side UNCHANGED
  (23 determinant-phenotype features — identical to the AMR-only spike).
- **Open Questions ratified before execution:** Q5/A = **YES** include the genome-level pathotype call (with
  the C1 QC-gated honesty contract); Q B = AMR `determinant-phenotype` wins tier precedence.

## The 5 design facts the implementation enshrined (so a reviewer doesn't "fix" them back)
1. **VF caller** `vf_runner.run_canonical_vf(all_hits=True)` (NON-frozen) — raises `-max_target_seqs` so
   tandem/multi-copy alleles survive; adds a coord-retaining `per_hit` list + `db_sha`. `per_gene`/`per_cluster`
   stay best-hit byte-identical → `build_vf_diff` UNCHANGED.
2. **no-`-parse_seqids` pin (empirically verified):** `makeblastdb` WITHOUT `-parse_seqids` keeps `sseqid` =
   the exact FASTA header first-token = the same token AMRFinder reports → the shared
   `phenotype_overlay.build_contig_name_map` reconciles BOTH overlays. WITH it, `sseqid` mangles to `gb|…|`
   and the contig-name map silently breaks. **Do not add `-parse_seqids`.**
3. **All VF metrics are ISOLATED** under separate keys (`virulence_join_quality` /
   `all_virulence_joins_symbol_fallback` / `virulence_determinant_feature_count` / `genome_pathotype_call`) —
   they NEVER feed the AMR `all_joins_symbol_fallback` nor the AMR GO/NO-GO spike gate (`gate.py` reads only
   the AMR `join_quality` keys — UNCHANGED).
4. **`genome_pathotype_call` mirrors the FULL deployed `resolve_call` contract** — `qc_pass` from
   `detect.assembly_qc.qc_verdict` (key is `qc_verdict`, NOT `pass`) + `support_gene_count` AND
   `meets_cross_axis_support` from a per-gene-coverage map. `status="insufficient_context"` ONLY when FASTA/VF
   unavailable (never when support is merely absent).
5. **Tier precedence:** high-confidence AMR join → `determinant-phenotype` > high-confidence VF join →
   `virulence-determinant` (its own `virulence` presence-only wall field) > `classify_feature_tier`; a
   symbol-fallback VF hit → `secondary_evidence` only (excluded from the tier).

## Verification (run before committing; all should be green)
1. `.venv/Scripts/python.exe -m pytest tests/test_genome_map_virulence_overlay.py tests/test_vf_runner_coords.py tests/test_pathotype_vf_diff.py tests/test_genome_map_cli.py -q` → 41 passed.
2. `git status --short -- dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json` → EMPTY (frozen surface untouched).
3. Optional full suite: `.venv/Scripts/python.exe -m pytest tests/ -q --ignore=tests/test_models_foundation.py -p no:cacheprovider` → expect ~1500 passed, 0 regressions (~90s).
4. Optional live re-run (D: is back online): `MSYS_NO_PATHCONV=1 uv run python -m scripts.genome_map --genome-fasta D:/dna_decode_cache/refseq/GCA_002180195.1/genome.fna --gff D:/dna_decode_cache/refseq/GCA_002180195.1/bakta/GCA_002180195.1.gff3 --organism Escherichia --sample-id GCA_002180195.1 --out-dir D:/dna_decode_cache/vir_smoke` → expect `virulence_status=FULL`, 27 virulence-determinant features. (Native blastn at `C:/Users/Farshad/ncbi-blast/bin`; no Docker needed for the VF part.)

## After this lands — the forward options (NONE auto-started; user's call)
The project is at a **terminal honest state on free public data** (`wiki/project_frontier_map_2026-06-19.md`).
After the virulence overlay commits, the genuinely-open moves are:
1. **Genome-map next tier / polish** — non-E.coli VF DBs, a virulence GO/NO-GO spike gate, or richer
   homology tiers (hmmer/Pfam/eggNOG). All deferred in v1; user investment call. The virulence spike stays
   AMR-only (`run_genome_map_for(virulence=True)` is opt-in, default False).
2. **TB ~1.6 TB regeno SCORED run** — now feasible (D: back, 4.4 TB free). Produces the v1b SCORED TB number.
   Multi-hour-to-days fetch with a documented USB-hiccup crash history → stage as a restartable detached batch
   (`scripts/stage_tb_vcf_subset.py` + a regeno fetch). Plan: `plans/TB_AMR_Decoder_RIF_INH_On_CRyPTIC_Plan/`.
3. **Evo 2 zero-shot variant-effect probe** (parked; `wiki/evo2_zeroshot_vep_lead_2026-06-19.md`) — label-free
   likelihood-delta on the cipro QRDR mutations. DIFFERENT method from the CLOSED embedding arm; LOW expected
   value (conservation≠resistance caveat); needs a REAL GPU (NOT the GTX 860M — hosted BioNeMo API or the
   Precision 7780). Not a pivot; does NOT reopen the embedding bet.
4. **Label ACQUISITION** (the only path that reopens the learned arm honestly) — a USER sourcing decision;
   anchor `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`.

## Working conventions (this project)
- `.venv/Scripts/python.exe -m pytest` (uv cache wedging this week, `os error 183`; `uv cache clean` when convenient).
- Commit straight to `main` (= the cross-machine sync channel; user syncs ~weekly). Footer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- Per-genome 18 MB maps under `wiki/genome_map_spike_*/` are gitignored (regenerable).
- D: is back online (4.4 TB free). Native blastn + the committed VF DB make the virulence overlay D:-free;
  only the refseq genome cache + Bakta/AMRFinder outputs live on `D:/dna_decode_cache/`.
- After committing, OPTIONAL: a `/project-state --append-action` row (ledger row ~132) recording the
  commit+push of the virulence overlay, and a `/retrospective` pass (the prior session's docs already updated
  CLAUDE.md / decisions-log / genome_map_usage / plans-index — `/documentation` + `/retrospective` ran as part
  of execute-plan, so this is light).

## Parked / not-this-thread
- The learned/embedding decoder is a CLOSED NEGATIVE on free data (3 de-confounded failures incl. Arabidopsis
  2026-06-12). Do NOT reopen without a label-clean substrate.
- Public-label AMR expansion is banked (`wiki/negative_results_map_2026-06-13.md`, the 8 rejection gates).
