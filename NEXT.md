## LIVE STATE 2026-11-25 — the finetune arm is RESOLVED; the slugs that must NEVER be re-pushed are still live



## DECISIONS NEEDED 2026-11-28 — read `wiki/decisions_pending_2026-11-28.md`

Three authority calls, with evidence and limits in that memo:

1. `frame_verdict.CHECK_ORDER` has 9 outcomes, its frozen bar names 8 — the grader can return a
   verdict the bar gives no reading for. Both committed; predates today.
2. Ratify or reject the 2-of-5 baseline-panel definability split. Nothing was written to the
   budget-4 gap artifact: that path is also the north-star step-4 expiry and a test fixture.
3. Do NOT commit the working-tree `wiki/data_inventory.{json,md}` — degraded census, restore
   from HEAD.

Leak DiD re-measured at 500 windows/half (**arms ledger unchanged at 0.6535**):
**LEAK_NOT_IDENTIFIABLE**, pre-registered -0.0054 now excluded. `run2_reseed_seed1` is the
CONTROL, which a first pass misread — retracted in `800de831`.
## HUMAN LOCUS-CAPABILITY RESULT 2026-11-26 — redirect the regulatory arm

The real two-cell lentiMPRA substrate probe ran on 2026-11-26 and returned
`WITHIN_LOCUS_SUBSTRATE_UNAVAILABLE`.

- Real shared elements: **7,321**, not the earlier 7,819 estimate.
- HepG2 and K562 observations per shared element: **exactly 2**.
- Loci meeting the `min_n=30` locus floor: **0**.
- Cross-tissue HepG2/K562 Spearman: **0.8595**.
- Artifact: `wiki/human_locus_capability_probe_2026-11-26.json`
- Runner: `scripts/human_locus_capability_probe.py`
- Pure and real-data regressions: 53 lentiMPRA/locus tests pass.

The HepG2/K562 panel cannot answer within-locus across-tissue conditioning.
It can answer per-element cross-cell association and per-cell sequence→activity
prediction, but it cannot support the instrument's locus floor. Do not launch a
generator or claim a locus-deconfounded result from this substrate.

Immediate redirect:

- stop expanding the lentiMPRA within-locus/tissue program,
- preserve the existing human oracle validity result,
- use lentiMPRA only for per-cell oracle and design-to-spec work,
- acquire or find a true multi-tissue/multi-biosample assay before revisiting
  the within-locus across-tissue metric.

**NEVER re-push these Kaggle slugs — every one holds retrieved output** (verified file-by-file against
`wiki/finetune_arm_kaggle_logs/`): `pa-finetune-arm-v1` / `v2` / `v3`, `pa-finetune-ckpt-probe` /
`probe2` / `probe3`, `pa-finetune-chain-s1` / `s2` / `s3`. The next link needs a FRESH slug.

**This supersedes the 2026-11-14 section that said to push to `pa-finetune-arm-v1`, which was
DANGEROUS.** That slug holds `attempt1_pa-finetune-arm-v1_2026-11-14.log.json`, so following the
instruction would have overwritten retrieved output — and `plan_run` would have refused the launch
anyway, for the reason below. Recorded rather than quietly deleted, because the instruction was written
in good faith from a route decision that a later measurement killed.

**THE ARM IS DONE — this section's old body said otherwise and was stale by two days.** Corrected
2026-11-25 against `soraya_runs/arms-ledger.json` + `wiki/finetune_arm_verdict_2026-11-17.json`.

- **`finetuned-promoteratlas-encoder-150bp` is `measured`, metric 0.6405, resolved 2026-11-18T03:28:09Z,
  verdict `FINETUNE_ADDS`**, all five of the frozen bar's `check_order` gates confirmed.
- **All ten seeds completed** (median 0.6405; min 0.6102, max 0.7070) against an incumbent median of
  0.5076 — delta **0.1329** against a required margin of 0.02.
- **All three shuffled-label nulls RAN** (0.1799 / 0.1997 / 0.1725). `verdict_secured` is **true** in
  `wiki/finetune_arm_session4_2026-11-16.json`, whose `status` is `COMPLETE`.
- Against the pooled 0.941 ceiling, fine-tuning moves the fraction of achievable signal from **0.3437 to
  0.6807** — roughly doubling it. Wall time 7.93 h against the 9.5 h ceiling, so it terminated on
  completion, not on the clock.

**What the old body got wrong, recorded rather than deleted:** it said 9 of 10 seeds, zero nulls run,
`verdict_secured` false, and that the arm "correctly reads `planned` / `blocked:untried-arm`". Every one of
those is superseded. It also presented a four-option authority fork (chain / re-size the bar / re-route /
accept a gradeable prefix) as *"unchanged and still the user's"*. **That fork is CLOSED — the chain was run
and it completed.** Nothing is owed on it.

**The honest limit that DOES survive, and is not in the old body:** the bar specifies a 38-epoch diagnostic
horizon and the run executed **7**. The verdict is unaffected, because it reads only the frozen 7-epoch
prefix while the horizon governs a pre-declared INELIGIBLE tail — but the tail existed to answer
*"would more budget have helped, empirically"*, and **that question is NOT closed.** This is a
frozen-budget result, not a converged one. Second surviving caveat: the nulls clear the ratio gate but are
not the benign case — they sit ~0.29 ABOVE an untrained baseline that is negative (-0.1123 / -0.0481 /
-0.1530), i.e. training on destroyed labels still reaches ~0.18 where an untrained encoder reaches ~-0.11.

**Two platform questions the older resume checkpoint lists as UNVERIFIED are RESOLVED — do not
re-litigate them:** the checkpoint-mount directory name
(`ckpt_mount_probe3_CONFIRMED_2026-11-15.log.json`) and cross-session resume, proven in all three
chained links by the reconciler AND independently by the kernel log — two routes sharing no code,
insisted on because the reconciler has a recorded vacuous pass on exactly that measurement.

## CORRECTION 2026-11-01 — CLAUDE.md's "14 failures, none a code defect" is wrong on the COUNT, and the one code defect is now FIXED

**Live count, grounded on a full run 2026-11-01 (9m48s): 13 failed / 7,238 passed / 180 skipped.**
CLAUDE.md still says *"14 failed"* and *"all 14 remaining are absent DATA, NONE a code defect"*. The
count is stale and **the second clause is wrong twice over** — once for most of today (a real
regression, fixed below), and once standing: **12 of the 13 are absent data and the 13th is not.**

> `tests/test_tile_chain_constraint_registered.py::test_the_eleven_pre_existing_constraints_are_still_present`
> ends in a hand-maintained `assert len(rows) == 12` while the constraint registry now holds **14**.
> Its eleven NAMED constraints are all still present — only the total drifted. **Pre-existing**, not
> from today's work (which registers no constraint; `git status dna_decode/constraints/` is empty).
> **Deliberately NOT re-pinned** — a count rail is guard structure, and bumping `12` → `14` is both
> the one-character fix and precisely the tidying the rail exists to interrupt. Open question
> recorded in `TODOS.md`: whether a total-count assertion earns its keep beside the named-membership
> loop that already carries the invariant.

CLAUDE.md **cannot be corrected in place** — it sits at 35,459 of its 35,460-word ceiling (room for
~one word) and additions there are USER AUTHORITY. This entry is the correction of record;
`README.md`'s suite paragraph carries the same figures and HAS been updated.

**RESOLVED — the code regression is gone.**
`tests/test_substrate_loader_ownership.py::test_each_substrate_accession_is_read_by_exactly_ONE_module`
fired on *GSE144621 is named in `['glm/tile_controls.py']` as well as its owner
`glm/genomewide.py`*, because `tile_controls._rows` opened the peak-tile file itself with a
hardcoded tab split — correct for that file and wrong by exactly the mechanism the guard names (the
fragment siblings are SPACE-delimited; `sniff_delimiter` in the owner is what decides per file).
`load_control_tiles` now calls `genomewide.load_peak_tiles` and filters the result, the parser is
gone, and that file is **8 passed** (was 1 failed / 7 passed). Routing changed **no data**: 936
`random` + 470 `neg_control` at 150 bp, identical. The signature narrowed to a cache **DIRECTORY**,
because the owner resolves the filename itself.

**One instruction in the superseded version of this entry was WRONG and was deliberately not
followed as written.** It said *"Do **not** add an `ALLOWED_MENTIONS` entry; that list is for
prose-only mentions"* — and after the fix the surviving reference IS prose-only: a single docstring
sentence *explaining that the owner parses the file*. The guard matches `\bGSE\d+` over source text
and cannot tell a docstring from an import, so a row was added, **but only after measuring that with
the parser deleted the guard still failed on that one mention** — which is the evidence the
exemption covers prose rather than a live reader. The structural fact is asserted at the structure's
level instead: `tests/test_tile_controls.py` uses `ast.parse` to check no import node names `gzip`
and no attribute access targets it, because a substring search for `gzip` would match the very
docstring that explains the delimiter trap. That substring form was specifically rejected and must
not be reinstated.

Full review:
`plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/brainstorm.md`.
Residuals (ordering divergence bounded in shape not magnitude; the duplicated corpus-window replica;
`pos_control` unavailable at 150 bp) are in `TODOS.md` under *Frame-conditioning / corpus-window
residuals (2026-11-01)*.

---

## DECISION NEEDED 2026-11-01 — no ledger in this project has a machine-readable MVP bar, and on FOUR of them that reads as permission to write a new one

**Do NOT run `/soraya --until-mvp` on `glm-condition-conditioning`, `glm-distribution-alignment`,
`glm-generator-finetune`, or `glm-self-distillation`.** CLAUDE.md names only the first. The
prohibition is four times wider than recorded, measured by sweeping every ledger:

| ledger | shape | `parse_diagnosis` | hazard |
|---|---|---|---|
| glm-condition-conditioning | TABLE | `empty-section` | **invites a replacement bar** |
| glm-distribution-alignment | TABLE | `empty-section` | **invites a replacement bar** |
| glm-generator-finetune | TABLE | `empty-section` | **invites a replacement bar** |
| glm-self-distillation | TABLE | `empty-section` | **invites a replacement bar** |
| glm-learned-representation | TABLE | `unparseable-items:4` | refuses — safe BY ACCIDENT |
| 6 non-GLM ledgers | list | `unparseable-items:N` | refuse — safe |

**11 ledgers, 0 parseable criteria.** `mvp_criteria._LIST_ITEM` matches markdown LISTS; these bars
are markdown TABLES. The four `empty-section` ledgers are the dangerous ones, because
`empty-section` is one of the two diagnoses that **licenses proposing a bar** — so an `--until-mvp`
run concludes there is no bar and offers to write one, silently replacing a frozen bar a human can
read perfectly well. G-A is safe only because its table happens to trip the *refusing* diagnosis
instead; that is luck, not design.

**A second, nested problem on G-D specifically.** Its frozen bar's criteria 1 and 2 name
`dna_decode/glm/finetune_corpus.py` and `tests/test_glm_finetune_corpus.py` — **both ABSENT**. The
work shipped as `scripts/glm_finetune_window_split.py` + `tests/test_glm_finetune_window_split.py` +
`tests/test_glm_finetune_kernel_builder.py`. The bar was written at project-init as a guess at the
implementation shape and the implementation chose another. So the *substance* of 1 and 2 is
satisfied (disjoint split, overlap 0, re-asserted in-kernel) while the *predicates* are false
forever, and the only mechanical repair is to create files at those paths — bar-gaming. The two
problems compound: fix the table-blindness and the newly-visible bar immediately fails on the path
mismatch.

**Why this is yours, not mine.** `mvp_criteria.py` lives at `C:/my_skills/skills/soraya/scripts/`,
is shared across every project, and is an acceptance-bar surface — that is the standing §2.4
authority call. Repointing G-D's criteria 1 and 2 at the real paths is an acceptance-bar edit after
seeing the result, which is also yours. **Recorded here and not in CLAUDE.md** because CLAUDE.md is
at 35,459 of its 35,460-word ceiling.

---

## D2 AUTHORISED 2026-10-19 — sourmash is an admissible G8 distance source (user)

The AUTHORITY gate is cleared; the EVIDENCE gate is not, and Step 6 must not fire on the
first alone. `plans/.../Mission_Reconciliation_A_Native_Clonality_Route_Plan` pre-committed
Step 6 (amend the CABBAGE G8 bar's distance source) to **Step 5 exit 0**. Step 5 returned
**exit 2 `INSUFFICIENT_DATA`** — all 10 cohorts blocked, `data/raw/` empty, D: offline — so
**zero agreement numbers exist**. Amending the bar now would record sourmash as validated
against a bar it has never been measured against: bar-shopping in the honest direction,
which this plan already declined once (see the `clause_logic` annotation).

**To fire it:** reattach D: -> `uv run python -m scripts.sourmash_vs_mash_crosscheck`
(needs `uv sync --extra lineage`). Exit 0 = all three clauses met -> Step 6 is unblocked
and pre-authorised, no further decision needed. Exit 1 = sourmash does NOT substitute for
Mash, and that is a real finding, not a failure to retry. Exit 2 = still drive-blocked.

The scorer is now EOL-safe (`verify_derived_from`), so the frozen bar no longer refuses to
load on an LF checkout — previously it would have exited 4 claiming the derived-from table
had changed. Nothing about the bar's three requirements moved.

---

# Open threads — what a fresh session should pick up
## TILE-RUNG SELECTION VALUE + THE HYBRID ARM — both ran 2026-10-14

Two measurements, and the first one **refutes the recommendation I had just given**.

**1. `SELECTION_IS_THE_BOTTLENECK`** (`wiki/glm_tile_selection_2026-10-14.md`). I argued the scoring
half was probably adequate and the prediction criterion should likely be retired for generator work,
citing the designed grid's 0.7937 of achievable selection value at k=10. On the tile rung the best
scorer realises **0.4807 median (2 of 10 seeds clear the frozen 0.50 bar)** — so the 0.7937 does not
transfer and prediction work on this rung IS justified. **It is a near-miss, not a collapse:** the k=10
null band is stdev 0.0485, so 0.48 sits ~10 stdevs above chance, and bar B passes 10/10 (0.485 vs gc
0.295). **The k-trend is the actionable part and it is monotone — 0.388 (k=2) → 0.485 (k=10) → 0.537
(k=20) — so the current oracle DOES clear the bar at k=20.** The answer is conditional on how many
candidates the loop can afford to assay, not a flat verdict on the scorer.

**2. `HYBRID_ADDS`** (`wiki/glm_tile_hybrid_arm_2026-10-14.md`) — **the first learned component to
measurably contribute on this substrate.** conv ⊕ max_tx rank-average: paired **+0.0339 over the
scalar, 10/10 seeds, CI95 [0.0238, 0.0409]** excluding zero; 0.4254 vs 0.3968, moving the rung from
42.2% to 45.2% of its 0.941 ceiling. The published GEMME pattern reproduced on DNA: **conv alone LOSES
to max_tx (−0.0767, 0/10) and still lifts it.** Control clean and decisive — max_tx with a *shuffled*
conv scores −0.1232, 0/10, so a content-free vector hurts and the gain is the conv's content.

**This does NOT overturn the G-A encoder-gate negative.** That gate asked "does the conv beat one GC
number?" — still no, CI includes zero. This asked whether the conv carries anything the *biophysical
scalar* lacks — yes. Both hold.

### The authority fork these two create

**Does a learned component contributing only IN COMBINATION satisfy the self-distillation ban's
precondition** (*"forbidden until a learned representation lands"*)? Not mine to declare. Note the two
results pull in opposite directions on the criterion question: #1 says keep improving prediction, #2
says a learned component now helps. Both support continued work on the *prediction* side, which is the
opposite of the pivot I proposed.

### Named follow-up, deliberately NOT run

**Re-run the selection measurement with the HYBRID as primary scorer.** It is now the best measured arm
(0.4254) and Step 3's bar was frozen with `primary_scorer: max_tx` before any number was read. Asking
whether the new best scorer clears the *same unchanged* bar is legitimate and cheap — but it needs its
own pre-registered run, not a retrofit into a committed artifact, and adding it mid-execution is the
scope drift the plan's non-goals exist to prevent.

Ledger: 10 arms, 6 measured / 4 planned; `ceiling_verdict` still **`blocked:untried-arm`** — a cheap arm
winning does not settle the ceiling question (R5 holding as designed).

### THE DRIVE LETTER IS PINNED — `D:` is now the 1 TB volume, permanently (2026-10-21, user-authorised)

The card reader that had been squatting on `D:` was evicted (`mountvol D: /D`) and the
**External1TB** volume (ST1000LM024, 931.5 GB, 907.5 free) reassigned to it. `MountedDevices` now
holds a **24-byte disk-signature binding** for `D:`, so it survives reboots and replugs — the
letter will not drift again for this disk. A junction `D:\dna_decode_cache ->
D:\PythonProjects\dna_decode_cache` reconciles the nested layout.

> **SCOPED 2026-10-23 — the binding holds, the ATTACHMENT does not, and the "measured consequences"
> below are all conditional on it.** Measured today: `wmic logicaldisk` lists **only `C:` and `S:`** —
> neither `D:` nor `E:` resolves, and `D:/` is not even enumerable. A disk-signature binding fixes
> **which letter a disk gets when present**; it cannot keep a removable disk plugged in, so "will not
> drift again" was the right claim about the wrong quantity. **Every `D:/…` bullet below should be read
> as "when that disk is attached".** What made the overnight promcalc run possible instead was a
> **repo-local mirror on `C:`**: `data/mpra` (96 MB) and `data/promcalc` (159 MB) are **real
> directories, not junctions** (`dir /AL data` finds no reparse point), so the whole stage-A → B → C
> chain ran with no removable volume at all. That is the FUTURE_FEATURES "mirror heavy artifacts off
> removable media" proposal, realised in practice for this one substrate — and it is also why the 12th
> expiry constraint, keyed `file-exists` to a `D:` path, can no longer fire.

**Measured consequences, not inferred (all conditional on that disk being attached):**

- `scripts/glm_tile_predictions.py` **with no arguments** prints `tiles 44106 | ceiling 0.9410`.
  The hardcoded default resolves; so does `scripts/promcalc_decomposition_arm.py`'s **bare**
  `load_peak_tiles()` (it has no `--cache-dir` flag and is now unblocked anyway).
- `D:/hf_cache` (4.4 GB, the GENERator-1.2B weights) and `D:/uv_cache` resolve at the documented
  paths. `run_glm_anchor_detached.bat` no longer hardcodes `E:` — **so the sweep's "`E:` attachment
  is per-session and never inherited" hazard is retired for this disk.**
- The 12th constraint's expiry test is now **live**: it names `D:/dna_decode_cache/promcalc/...`,
  whose parent resolves, so `file-exists` can actually fire. The drive-letter half of its stated
  limitation is discharged (the rationale keeps the original wording by design — that field is the
  reason *as written at the time*).

**STILL ABSENT and the wall classification is UNCHANGED:** `refseq`, `tools/` and
`dna_decode_data` live on a **different physical disk** that is not attached, so
`D:/dna_decode_cache` is `UNFIT` under `CACHE_REQUIRED_SUBDIRS` and the probe says so.
**The failure mode shifted and that cost was accepted up front:** those paths now resolve and fail
as *dataset* missing rather than *drive* missing, which `external_root.py` calls the worse
direction — the mitigation (`CACHE_REQUIRED_SUBDIRS` + `data_inventory`'s per-report probe) shipped
before the risk did.

**Two CLAUDE.md claims are now false** and were not edited (the file has **one word** of headroom
against its enforced ceiling): *"`D:` is a card-reader slot with no media"* and *"`E:` is a 931.5 GB
ST1000LM024"*. Also false and pre-existing: *"`HF_HOME`/`UV_CACHE_DIR` set persistently"* — both are
unset in User, Machine **and** process scope; they exist only inside the `.bat`.

### The tile chain's real blocker — DIAGNOSED 2026-10-21, and it is two walls not one

`wiki/tile_chain_substrate_wall_2026-10-21.{md,json}`. The follow-up above ("re-run the selection
measurement with the HYBRID as primary scorer") got its pre-registration on 2026-10-17 and has not run
because nobody had located its substrate. **I said in conversation that the substrate was present and
that was wrong** — I had found the raw MPRA tables on `E:` and taken them for the input the chain
consumes. Measured, there are **two** dependencies with **different** wall classes:

- **MPRA peak-tile table — `code-closable`, and now CLOSED.** `glm_tile_predictions.py:282` called
  `load_peak_tiles()` with no arguments, taking the hardcoded `D:` default and dying on an unhandled
  `WinError 21`. It now takes `--mpra-cache`; with the reachable copy the chain prints
  **`tiles 44106 | ceiling 0.9410`** — both the published anchors, exactly — and gets past the tile load.
- **promcalc per-tile decomposition — was `data-absent`, now CLOSED (2026-10-23).** The install decision
  was taken: an isolated env at `C:/Users/Laptop Portable/tools/promcalc_env` (sklearn 1.9.1, project env
  still 1.8.0 — the probe's `isolation_required` held, though the env is on `C:`, not the documented
  `D:/tools/promcalc_env`). Stage B ran 770 min / 7 workers to full 44,106-tile coverage, 0 errors.

**RESOLVED 2026-10-23 — the decomposition question is ANSWERED and the remaining items are narrower.**

- **Verdict `NO_GAIN`** (`wiki/promcalc_decomposition_arm_2026-10-23.json`): 9 dG components over 55
  features score **0.3925** vs **0.3968** for the single `max_tx` scalar reduced from the same table
  (paired −0.0043, 2/10 seeds, bar was margin 0.02 and ≥8). VALID not indeterminate — both reconcile
  anchors passed and the permuted control is clean at 0.0342. **Read it narrowly: both beat GC 0.3112 by
  ~+0.08 at 10/10 seeds, so promoter thermodynamics IS informative; only the nine-part breakdown over one
  number adds nothing (~730 discardable numbers per tile).** Not "thermodynamics doesn't help".
- **The hybrid-primary bar (2026-10-17) is STILL UNEXECUTED** — a different frozen artifact, now
  substrate-unblocked. Running it needs its own pre-registered run, not a retrofit.
- **MEASURED 2026-10-26 — the `orthogonal-modality-stack-promcalc-plus-N-channels` arm RAN** and is
  now `measured` in `soraya_runs/arms-ledger.json` at **0.4404**. Holding it `planned` through the
  decomposition's `NO_GAIN` was correct per R6; it now has its OWN evidence. **`WEAK_STACK`:**
  `decomposition`-as-third gains **+0.0145** over the 2-channel hybrid on **10/10 seeds, CI95
  [0.0108, 0.0172] excluding zero** — short only of the INHERITED 0.02 margin, which is not revisited.
  My pre-registered prediction (`NO_STACK` for both new combos) was **FALSIFIED**, and the
  orthogonality argument is **backwards**: `decomposition` is the MOST redundant with the comparator
  (+0.7812) and ADDS, while the less-redundant `gc` (+0.5508) SUBTRACTS (−0.0039, 4/10).
  `wiki/glm_tile_stack_arm_2026-10-26.md`; ceiling still blocked by 3 `planned` siblings.
- **The 12th expiry constraint cannot fire and needs retiring by hand.** It is `file-exists` on
  `D:/dna_decode_cache/promcalc/decomposition.csv.gz`; the data landed repo-local in `data/promcalc/`.
  Its own rationale named this case (the R9 pattern) — the stated limitation was the one that bit.

**The frozen bar was not edited either.** Its `not_runnable_on` note is half stale (`no uv` is false —
uv ran the full suite here today) and half right (`no D: data tree`), but amending a frozen
pre-registration is an authority call. Audited in the memo, not changed.

**Other verified-runnable high-VOI work** is resuming Steps 4–6 of
`plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/` — its substrate is byte-verified
present (GLM weights on `E:/hf_cache`, MG1655 reference in `data/cache/refseq`).

## THE MISSION IS NOW DEFINED — `project_state/trait-to-edit/big-idea.md` (authored 2026-10-09)

A validated mission definition for the north star, authored at the user's request from everything
established to date. **It parses and validates against the real `big_idea.py`**, and
`run_minister.preflight(...)` returns `(True, [], [])` — so the `minister` lane is genuinely available
for it, which it was not before (there was no mission big-idea.md; only two stale EP4 ones).

`mission-shape: greenfield` · 4 gaps · 10 checkable endpoint predicates · **Φ = 12** (finite by
construction, so the mission cannot generate unboundedly).

| gap | budget | what closes it |
|---|---|---|
| `amr-label-substrate-validated` | 3 | CABBAGE's G8 + BioProject-level diversity + one leakage-gated scored cell |
| `deconfounder-calibrated` | 2 | the F3 calibration actually run on 1001G + the IBS cross-check |
| `generator-beats-its-baseline-panel` | 4 | a generator beating the mandatory baseline panel on biologically-motivated splits |
| `trait-to-edit-end-to-end` | 3 | one bacterial trait → a ranked edit set, oracle-scored, rails audited |

It also names what must **not** be re-proposed (zero-shot embeddings on natural populations; ESM2 scale;
absolute cross-assay magnitude; the ETD route) so a generative lane cannot re-raise a closed negative.
**Step 4 (climb to animals) is deliberately absent** — premature until the generator gap has a result.

**Launching a minister mission on it is a separate decision and was NOT taken.** Authoring the definition
is not the same act as running it.



Short by design. **Transient** state only: what is in flight, what is waiting on the user, what was just
learned that isn't durable yet. Durable findings belong in `CLAUDE.md` / `wiki/`; scope facts belong in
`scripts/project_status.py` (derived, never written).

Prune aggressively. A stale entry here is worse than an empty file.

_Last updated: 2026-10-21._

> **Pruned 2026-10-05.** This file had sat at 2026-09-28 / 428 lines, and its own "known-stale" section
> had itself gone stale in the opposite direction (it warned that 46 traits was really 44; the live count
> is **46** again after `dna-identify` and `dna-tb` shipped). Everything removed was already durable in
> `CLAUDE.md` (gentamicin v2 lock, HCMV contracting, the declined NNRTI curation, PEAR's reclassification,
> the doubt-layer firing rate) or in `project_state/` (the five project families, F-C/F-D). Checked before
> deleting, not assumed.

---

## IN FLIGHT — deep-research campaign F1 / F3 / F4

Plan `plans/Deep_Research_Missing_Part_2026-09-28_Plan.md`. User ratified full breadth (F1–F4).
Screening spine is built and green (`scripts/deep_research_screen.py`, 13 tests; runs both
`rejection_gates` G1–G10 and `regime.screen_proposal`, and refuses to call anything promising until they
have run).

| family | state |
|---|---|
| **F2** condition coverage | **ANSWERED NEGATIVELY.** "Public K-12" adds samples, not condition breadth — glucose share 60% → **82.7%** of annotated rows, **42%** carry no carbon annotation, only **3** carbon sources reach 10 samples. Computed from the repo's own metadata; the "contact us" wall covers only `log_tpm`. |
| **F1** label substrates | **SEARCHED 2026-10-09 and it produced the campaign's strongest candidate.** `CABBAGE` / EMBL-EBI AMR Portal (release 2026-07) — 1.7M genome-phenotype pairs, lab-measured, CC-BY 4.0, **8.3 MB parquet, MD5-verified on the real FTP surface**. **9 of 10 gates pass**; only **G8 (`min_effective_lineages`)** is open. Through the repo's OWN `tier_for_isolate` it clears N=150 per-class at the STRICT bar on **four cells** — against the 2026-05-18 census verdict that **no** drug clears it. `wiki/f1_cabbage_label_substrate_2026-10-09.md`. Still unsearched: Zenodo/Figshare/Dryad + ETD aggregators; ATLAS is Vivli-request-gated. |
| **F3** de-confounded population designs | **The INSTRUMENT is CALIBRATED ON ONE PANEL PAIR — `VALID_DIFFERENTIAL_RESPONSE`** (2026-10-10, `wiki/f3_deconfound_calibration_1001g_2026-10-10.md`), superseding the 10-09 `blocked:positive_control_not_fold_compatible` whose substrate blocker was RESOLVED by using 1001 Genomes (1,135 natural Arabidopsis accessions + AraPheno phenotypes). It responds differentially to real population structure (**+0.0441** natural vs **0.0000** bi-parental cross). **What is NOT established is TRANSFER** — one positive and one negative panel only; it bounds the INSTRUMENT and validates no downstream organism-level prediction. `eval/regime.py::DECONFOUNDER_CALIBRATION` is the live status — read it rather than this row, which was stale for two days and contradicted this file's own `deconfounder-calibrated` row above. **The literature SEARCH half is still unrun.** Population DESIGN, not organism complexity, remains the measured barrier to the north star's step 4. |
| **F4** cross-assay dosage transfer | **NEVER SEARCHED.** |

**Relaunch SEQUENTIALLY (or 2 at a time) and forbid sub-agent spawning.** The first attempt spawned 12
concurrent agents; all died on a session rate limit with **zero** final reports. First-party search worked
fine immediately after — the fan-out was the constraint, not the session.

## IN FLIGHT — the cross-organism transfer ladder (drafted 2026-10-05, NOT saved)

A 9-step plan to turn the project's **single** step-4 data point (E. coli-tuned essentiality decoder applied
UNCHANGED to human, AUROC 0.5805) into a 5-rung ladder, to separate **(A)** decay-with-distance from
**(B)** a eukaryote-boundary cliff. Drafted + adversarially reviewed in conversation; **never written to
`plans/`** — deliberately, because the review mutated it (the save-once rule).

**Three findings that change it, all measured, none in the plan as drafted:**

1. **AUROC is the WRONG primary metric — it is tie-dominated.** 83.1% of human essentials score EXACTLY 0
   (E. coli 62.1%), and 11 of 13 patterns have P(hit|non-essential)=0.0000 — the decoder is **silent**, not
   wrong. Coverage and precision-where-it-fires move in OPPOSITE directions across the two rungs
   (coverage 0.3789→0.1689; precision 0.3684→**0.9583**) and AUROC blends them. **Use coverage.**
2. **The miss is ~2.5× host-specific vs phrasing, and the ladder survives** —
   `wiki/essentiality_missed_vocabulary_2026-10-05.md`. Mechanism (iii) floor 141/566 (24.9%) vs
   mechanism (ii) floor 57/566 (10.1%). So a cheap pre-test does **not** pre-empt the ladder, but ~10% is a
   catalogue PHRASING gap unrelated to phylogeny and must be separated out before any rung is scored.
3. **The FBA walls do not apply.** `fba/essentiality_labels.py` marks S. aureus `LABEL_WALLED` and
   P. aeruginosa `MODEL_WALLED`, but both are **model-gene-join** walls; this decoder joins to a genome
   annotation, so two "blocked" rungs are reachable. Reachability is still a HYPOTHESIS until fetched.

**Open AUTHORITY forks:** re-run `/technical-plan` around coverage, or `/save-plan` as-is? (The review also
wants a W0 schema probe before any parser, and the functional-class step demoted to a provenance-gated
audit — it is otherwise a second adjudication system.)

## DONE — GLM G-A + G-C v2 build plan (all 9 steps, 2026-10-08)

Pointer only — this section always said its derivations lived elsewhere, so it is reduced to the
pointers themselves. Plan `plans/GLM_G_A_And_G_C_V2_Build_Plan/technical-plan.md`, commits `e2aaf54`
to `f015cff`. Both families answered, both **VALID negatives** (clean nulls + registered protocols,
stamped machine-readably so a broken run could not satisfy the bar):
`wiki/glm_encoder_gate_result_2026-10-08.md` and `wiki/glm_condition_conditioning_result_2026-10-08.md`.
Pruned 2026-11-16 to bring NEXT.md back under its declared transient-store ceiling.

## The system design — drafted, awaiting ratification

`plans/Hybrid_Decoder_Architecture_Plan.md`. The hybrid is **not** "catalog + ML predictor" (that framing
scored 0 survivors). Measured shape: **CALL / DOUBT / EVIDENCE**.

| layer | status |
|---|---|
| **L1 CALL** — deterministic curated rules | shipped, 128 cells |
| **L2 DOUBT** — "this call may be incomplete, and why"; never a competing call | **shipped and firing** (`dna_decode/eval/doubt.py`); nobody in the field ships this |
| **L3 EVIDENCE** — de-confounding, nulls, denominators, leakage, provenance | built; **four** cohort-evidence layers now render (lineage, source-concentration, prospective, species-composition) |
| **L4 LEARNED** — forward/inverse, orthogonal modalities; molecular + constructed ONLY | shipped, bounded |

## Waiting on the user — authority calls, not executor tasks

1. **Re-score `Klebsiella × meropenem` on *K. pneumoniae* only?** (new 2026-10-05)
   Measured: the cohort is **38% *K. aerogenes*** and **all 16 false negatives are those off-species
   isolates** — on *K. pneumoniae* there are zero misses (sens 12/12). The pneumoniae-only matrix is
   **sens 1.000 / spec 0.880 / acc 0.919, N=37** against the published 0.467 / 0.900 / 0.683.
   Publishing it edits a **frozen unit** of the reproducibility freeze and changes a SCORED number on the
   trust surface. Current answer is **disclose** (`species_composition` layer, augment-only).
   `wiki/meropenem_fn16_diagnosis_2026-10-05.md`.

2. **Drop lone-porin counting from the meropenem rule?** *(evidence added 2026-10-05, pointing the
   other way — see the closed lever below)*
   All 3 false positives on the AR Bank cohort are porin-only, and the rule's `threshold=1` on a lone
   `ompK35`/`ompK36` truncation over-calls. But it edits sha256-pinned `amr_rules.py` → retires the
   active 2026-08-31 v2 lock and restarts the prospective clock for every cell, motivated by **3
   isolates** on a cohort whose own powering gate HARD-FAILED (n=5 susceptibles).

3. **Whether a single-source cell warrants more than disclosure.** 3 of 10 SCORED AMR cells rest on one
   BioProject. Current answer is *disclose*, namespace-separate. Demoting them is a scope decision.

4. **Do the two valid negatives CLOSE G-A and G-C, or defer them?** (new 2026-10-08; replaces the
   now-RESOLVED ratification fork, which the user answered *build the conv encoder anyway*)
   Both families now have a valid recorded negative, so the remaining question is purely an
   acceptance-bar one and it is the same shape for both. **G-C:** a `NO_GAIN` at n=2 media rules out
   conditioning *as tested on two growth media with 3-mer interactions* — not conditioning in general,
   and specifically not a model whose **representation** is modulated by condition. **G-A:** `WEAK`-but-
   not-distinguishable-from-zero rules out *this* encoder at *this* size on *this* substrate; a larger
   encoder is a different experiment and the result says nothing about G-B. Closing either is a scope
   decision; the executor default is to leave both as recorded negatives and start nothing new.

5. **~~Acquire a natural-population genotype panel so the de-confounder's positive control can run?~~
   ANSWERED BY MEASUREMENT 2026-10-09 — the substrate is resolved and it cost 330 MB.**
   `wiki/f3_substrate_resolved_2026-10-09.md`. The 1001-Genomes imputed SNP matrix is **0.33 GB**
   (HEAD-measured; a literature estimate said "a few GB"), HTTP 200, no auth, and **both MD5s match the
   publisher's own `md5sum.txt`**. It ships a **pre-computed IBS kinship matrix (10.3 MB)**, so the
   decisive question was answerable without touching the 919 MB genotypes. Through the SAME pipeline and
   guards the census ran on `dgrp`: **a fold-compatible grouping EXISTS at p5 — 108 clusters, largest
   9.5%, singletons 30.6%** — where `dgrp` passed at no rung. Phenotypes were **already on disk**: FT10
   overlaps **1,003 accessions (86.3%)**. And FT10 is the right control, not merely an available one —
   it is exactly the trait/panel of the closed Arabidopsis negative, where structure is MEASURED to
   dominate. **What remains is executor work, not an authority call:** run the calibration (needs the
   919 MB matrix for the ridge), and cross-check the grouping against this repo's own
   `ibs_distance_matrix` rather than the project's converted kinship — a disagreement there is a finding
   to diagnose, not a number to choose between. Caveats: genotypes are IMPUTED; p5 is a narrow window
   (p2 and p10 fail in opposite directions); 159 phenotyped accessions have no genotype.

6. **Whether to compress the other long CLAUDE.md bullets — NOW BINDING, not optional (escalated
   2026-10-08).** The growth-ceiling guard fired on an ordinary 60-word finding, because the file had
   arrived at **12 words of headroom**. The ceiling was raised 35,400 → 35,460 on the record (which is
   what the guard's docstring prescribes), but it is now constraining real findings rather than drift, and
   every future addition hits the same wall. `scripts/claude_md_weight.py` reports **43 bullets / 24,892
   words — 70% of the file — PROVABLY stored elsewhere in `wiki/`**, so the diet is available and would
   buy years of headroom. Compressing other sessions' bullets is a scope decision, which is why it is
   still here rather than done. Original note follows. The file loads every session; several long
   bullets cite a resolvable memo, so their derivations could become pointers. Measure first:
   `uv run python scripts/claude_md_weight.py`. Two long bullets have no external store and must stay
   whole — the tool already protects them.

7. **~~Does the CABBAGE label substrate reopen the AMR track?~~ RESOLVED BY THE USER 2026-10-09 — YES,
   it reopens.** *"CABBAGE should reopen the banked AMR track --> i think it should. use your best
   judgement."* The scope of that decision is narrow and was stated when it was made: what CABBAGE
   contradicts is the premise that **no usable label substrate exists** — it does **NOT** reopen the
   embedding negatives (closed on de-confounded signal-vs-structure grounds, which more labels cannot
   address, and independently corroborated at 24k-genome scale by Yu/Barquist 2025) and it does **NOT**
   change the generative north star, whose measured bottleneck is the generator.
   **Now executor work, no further authority needed:** close G8 (`min_effective_lineages`, needs the
   assemblies clustered) and measure BioProject-level source diversity (the `database` column is an
   overlapping aggregator label, and the screen flags its own G2 as "weak").
   `wiki/f1_cabbage_label_substrate_2026-10-09.md`.

8. **How should the GLM umbrella's `ceiling:` criterion be resolved?** (new 2026-10-12)
   Two arms were declared 2026-10-09. The 2026-10-11 run measured a **THIRD** (`encoder-pyramid`, 1.70
   decades) and closed **NEITHER** declared arm, so the criterion is still `blocked:untried-arm`. The
   accounting is now correct and the remaining question is an authority call, not a measurement: one
   declared arm is **SUBSTRATE-GATED** -- >=4 decades of span needs a >=10,000 bp coarsest branch against a
   150 bp tile, so **no budget reaches it** (pinned by arithmetic in `tests/test_glm_encoder.py`) -- which
   means resolving the criterion requires either changing the substrate or retiring the arm, and the
   arms-ledger resolution step was deliberately NOT run. `wiki/glm_pyramid_arm_2026-10-11.md`;
   negative already filed as `neg-5ad65fa5` with its revival condition.
   **Substrate fact for this fork, measured 2026-10-13 and now printed by `scripts/project_status.py`:**
   the "changing the substrate" branch has **cached** candidates — two fragment rungs sit at ~9% used
   (30,000 of 321,123 LB / 318,457 M9 rows, 48–499 bp), while `peak_tiles` and `sigma70_designed` are
   genuinely spent. **This does NOT un-gate the span arm:** 499 bp still cannot host a 10,000 bp coarsest
   branch, so that arm stays substrate-gated and the fork remains an authority call. What changed is only
   that "no data for longer sequences" was false — it came from `size_lever` being absent on three of the
   four rungs. `wiki/substrate_symmetry_2026-10-13.json`; `## The GLM` in `CLAUDE.md`.

## Cheap untried levers (executor work, no authority needed)

- ~~**`dna-identify` cross-check on the 16 meropenem FN.**~~ **DONE 2026-10-05** — ran it, and it
  CORROBORATES: 16/16 FN abstain (`above_max_distance`, 0.1129–0.1339, inside the 0.10849–0.13386 band
  `thresholds.py` records for K. aerogenes), 0 called `klebsiella_pneumoniae`, controls 12/14 called at
  0.0011–0.0090. **Partition agreement with the GenBank labels: 30/30, zero disagreements.** Scope: a
  CLOSED-SET router establishes NOT-pneumoniae, not *aerogenes* specifically.
  `wiki/meropenem_fn_identify_crosscheck_2026-10-05.json`.
- ~~**Re-run AMRFinder under the correct `-O` for the 23 off-species isolates.**~~ **ILL-POSED — CLOSED
  2026-10-05 without running it.** Measured from the pinned image (`amrfinder -l`): the organism list
  offers `Enterobacter_asburiae` / `Enterobacter_cloacae` / `Klebsiella_oxytoca` /
  `Klebsiella_pneumoniae` and **no `Klebsiella_aerogenes`**. So there is no "correct flag" to re-run
  under, and the species contamination in that cohort is **not remediable by re-running** — a structural
  limit of the tool's supported set, not a to-do. (The nearest available options sit under *Enterobacter*,
  which is where this organism's former name sits; choosing one would be a judgement about an unsupported
  species, not a correction.) This is also why the 23 isolates cannot simply be re-scored "properly".
- ~~**The 2 *K. aerogenes* true positives are unexplained.**~~ **DONE 2026-10-05, and they are NOT one
  class.** `GCA_003951185.1` is called R by `ompK35_G41TfsTer32` — a **lone porin truncation**, its only
  counted carbapenem determinant. `GCA_003951605.1` is called R by `blaNMC-A`, a real class-A
  carbapenemase. **This bears on authority fork 2 and points the other way:** a lone-porin call that is
  RIGHT, so dropping lone-porin counting converts a true positive into a false negative. n=1 against the
  3 AR Bank FPs — not settling, but the fix now has a measured cost on the same arm.
- **The 3 Klebsiella cohorts with incomplete cached AMRFinder runs** (gentamicin 3/60, tetracycline 33/60,
  ceftriaxone 54/60) have their species cross-tab withheld for that reason. Completing them needs Docker.

## The FBA switch cell, as it now stands

Direction is fine, **silence** is the problem: 61–76% of genes emit one identical ratio for every
condition. Do **not** build the relative ranking rule (oracle ceiling +4 genes, and a deployable version
must infer k = the original problem restated). **Do not quote a raw exact-set rate** — a conditionally
essential gene is two-sided by construction, so a constant call can never match one. The measured
bottleneck is fixed in both datasets: PRECISE-1K ∩ Keio carbon = **11 of 28**, 621 of 1,035 samples
glucose.

## THE ANCHOR REPRODUCED — the GLM sampling sweep, 3 of 6 steps shipped (updated 2026-10-21)

`plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/`. Steps 1-2 landed: the sampling
knobs are live (`wiki/glm_temperature_knob_2026-10-20.md`), the comparison mode is now a **contract
refusal** rather than two comments nobody enforced (`knob_contract.validate_comparison_mode`), and
`wiki/glm_sweep_invalidation_2026-10-20.md` records that the 2026-10-19 sweep was **doubly invalid**
(inert knob AND `conditional=False` against the aligned promoter set).

**Step 3 then RAN — `E:` was reattached and the anchor reproduced the recorded 2026-10-07 baseline
exactly**, all six rows bit-identical at `positional_bins=6` (`wiki/glm_sweep_anchor_2026-10-20.json`
— **committed; the earlier "untracked, commit it" note here was stale, verified 2026-10-21**):
`kmer` 0.7237 vs null 0.5990, `both` 0.7706 vs 0.6679, both
`WORSE_THAN_MARKOV_NULL`. So the harness, the device (CPU) and the corpus build are all vindicated;
what is NOT done is the sweep.

**G-D's unblock trigger is STILL UNDISCHARGED, and it is now ONE run away** — the anchor is the T=1.0
comparator, not a swept point. Three facts worth carrying: the sweep is **CPU inference** (the 8 GB GPU
blocker gates fine-tuning, not this); `E:` attachment is **per-session and never inherited**, so the
4.33 GiB weight file under `E:/hf_cache` can vanish between sessions; and launching must go through
`run_glm_anchor_detached.bat` / `glm_run_lock.py run`, since two concurrent torch jobs is this host's
access-violation hazard and `acquire`+`release` does NOT hold the lock from a launcher (measured — the
second job starts). Pass `--bins 6` explicitly: the runner default is 10 and one anchor attempt was
already voided by it (`..._bins10_UNCOMPARABLE_...json`, which also shows the `both` gap shrinking
0.1027 → 0.0612 from bins 6→10, so prefer `kmer` as the primary comparator). Full resume conditions +
the carried honest limits are in `TODOS.md`.

**Do not quote the 2026-10-19 headline** (GLM 0.8141 vs null 0.3876). Valid baseline is the committed
`wiki/glm_generator_falsifier_REAL_2026-10-07.json`: 0.7237 vs null 0.5990 (`kmer`),
`WORSE_THAN_MARKOV_NULL` — which is the project's **recorded premise**, not a new finding.

## User-authority calls open 2026-10-20

- **Install BLAST+ on this host?** Measured absent (`find_blastn()` -> `None`;
  `wiki/host_toolchain_and_crash_2026-10-20.md`). CLAUDE.md asserted it was installed for months —
  true of the host it was written on, never synced here. Until it is installed, every blastn path takes
  its offline branch: `vf_diff` / virulence overlay, fungal ERG11, HIV RT/PR/IN/CA and SARS-CoV-2 Mpro
  genome-FASTA modes, and the whole typing suite. **Not claimed: that any committed typing validation is
  invalid** — only that none can be RE-RUN here. `winget install --id NCBI.BLAST`, or the NCBI win64
  installer into `C:/Users/Laptop Portable/ncbi-blast/bin` (satisfies the resolver's 3rd branch, no env var).
- **Mash on server2?** Two blockers, neither executor-closable: Docker Desktop will not start from
  session 0 over SSH, and ~5.37 GB free RAM against a ~8 GB WSL2 default. Unblocks via a `.wslconfig`
  cap (`[wsl2] memory=2GB processors=2` — profile scope, so the user's file) or closing BlueStacks
  (~3.3 GB of the pagefile). Reversible route once capped: `wsl --install -d Ubuntu --no-launch` then
  `apt-get install -y mash`; undo with `wsl --unregister Ubuntu`.
- **The lineage metrics stay blocked regardless of both.** `amrfinder_runs/` is absent in 10 of 10
  cohorts and the committed artifacts carry zero per-isolate keys, so the chain ends at AMRFinderPlus ->
  Docker. Sourmash solved distances and does not touch this.

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**, including in this file.
  `uv run python scripts/project_status.py` is the authority — it is derived, never written.
- `wiki/project_distillation_2026-08-29.md` says 46 CLI traits. That was wrong when written (it was 44)
  and is now right again by coincidence, which is the cleanest possible illustration of why a written
  count is not evidence.

- **`torch` IS NOT IN THE DEFAULT DEPS, and no doc said so until now.** It was deliberately moved to
  the extras; `[forward]` carries `torch>=2.2` + `transformers>=5` and the pyproject records it as
  verified on torch 2.7.1+cu118. A plain `uv sync` therefore breaks EVERY GPU path (GLM encoder,
  generator, `dna-forward`) with a bare `ModuleNotFoundError`. Measured on a fresh clone 2026-10-29,
  where it cost a cross-machine round-trip: `uv sync --extra forward`, and run with
  `uv run --extra forward`. **Do NOT reach for `[ml]` just to get torch** — it pins `transformers<5`
  for the closed-negative DNABERT-2/NT track and the two are mutually exclusive in
  `[tool.uv].conflicts`, so uv refuses the pair rather than silently mixing them.
- **CLAUDE.md says `openpyxl` is "deliberately NOT installed". It IS installed, at 3.1.5**
  (measured 2026-10-28; it arrived transitively, `uv.lock` md5 `892c08339f82bd1a23f898821adb9735`
  unchanged). The sentence also has the gate backwards: the 18 `test_lentimpra.py` skips are
  DATA-gated, not import-gated — an import gate closes by installing, a data gate does not. Four of
  that file's tests now build a synthetic workbook with it and pass with no `D:` attached.
- **The fix for both of the above is blocked on a compression decision, not on knowing them.**
  CLAUDE.md measured **35,459 words against its own 35,460-word ceiling** on 2026-10-29 — ONE word of
  headroom — so any correction to it now has to pay for itself. `uv run python scripts/claude_md_weight.py`
  certifies **43 bullets / 24,106 words** as long AND provably stored elsewhere in `wiki/`, so the
  headroom exists; choosing which to reduce to a pointer changes what every future session reads first
  and is a user call. `tests/test_claude_md_growth_ceiling.py` is the rail, and raising the number
  rather than compressing is the thing it exists to prevent.

## FOUR AUTHORITY CALLS (2026-10-30) — THREE ANSWERED, ONE STILL OPEN; call 3 was later REVERSED

**ANSWERED later the same day. Read this before the four items below, which are kept verbatim as the
evidence each decision rested on.**

| # | call | decision | state |
|---|---|---|---|
| 1 | fragment-arm null bar | **spend the seeds (~270)** | **RESOLVED at 270 seeds, pre-registered Branch A: median +0.0092, exact `median_ci95` `[-0.0133, +0.0211]` (achieved coverage 0.9556). The validity gate PASSES, `INVALID_NULL_SCORED` is retired, and the arm closes `NO_GAIN`.** **IT IS A WEAK PASS and Branch A requires labelling it one** -- the interval still CONTAINS the 0.02 bar. **CORRECTION: the prior version of this row said the interval was "clear of the bar" at 215 seeds and "no longer a WEAK pass". b3's final 55 seeds moved the upper bound +0.0189 -> +0.0211, back across the bar -- exactly what that row's own caveat warned about, so the caveat was load-bearing rather than hedging.** Forecasts: **P1 falsified a SECOND time** (~80% that the point estimate would still be >= 0.02; measured +0.0092). **P2's claim holds** (270 buys a point estimate, not a confident resolution) but its forecast interval [+0.0144, +0.0340] sat entirely above zero while the measured one contains zero, so the sizing model mispredicted the LOCATION. Current number: `wiki/fragment_null_270seed_result_2026-10-31.md`; the 160-seed and 215-seed figures are SUPERSEDED snapshots. Pyramid-vs-kmer4 is unchanged -- it loses to a 4-mer ridge either way |
| 2 | G-B dependency predicate | **fix it** | DONE `0750ec63` — G-B is `blocked:predecessor` again; do NOT tidy the free-text clause back to a bare slug |
| 3 | G-D 8 GB hardware edge | ~~hold it~~ → **RE-OPEN IT (user, 2026-10-31)** | **EDGE LIFTED and the LoRA fine-tune is RUNNING on a Kaggle T4.** The sweep half was discharged first (3 of 3 points plan-eligible, `f9c47d7e`) and its answer — no temperature setting beats the Markov null — **falsifies GD-H4**, so sampling had not already done the work. The "8 GB GPU" clause was never missing hardware (free T4s are 16 GB). Bar FROZEN before any number in `wiki/glm_generator_finetune_acceptance_bar.json`: mode-`both` below **0.6679**, four outcomes named, **do NOT re-specify it**. `blocked_by` is a bare `none` on purpose — prose in that cell makes `advance_ranker` fail closed, which is what keeps G-B blocked. Umbrella row 21 / G-D row 4 |
| 4 | `mvp_criteria` table blindness | **STILL OPEN** | unchanged — do NOT run `/soraya --until-mvp glm-condition-conditioning` |

**Call 1 came with a correction to its own framing.** It was recorded as "re-specify the bar, or
accept no verdict — both legitimate". A third option existed: the 50-seed bootstrap CI
`[-0.0126, +0.0466]` **contains zero**, so the null is not measurably dirty and the bar is failing on
a point estimate narrower than its own SE; re-specifying it as a *leakage* bar (set where leakage
lands, not where noise lives) would have passed validity and let the arm close as **`NO_GAIN`** —
i.e. anti-flattering, which is what made it safe to propose. The user chose seeds on the same
evidence. Pre-registration, written before any of those seeds existed:
`wiki/fragment_null_270seed_prereg_2026-10-30.md` — it forecasts ~80% that the gate STILL fails, and
records that 270 seeds buys a point estimate rather than a confident resolution (the CI is forecast
to straddle 0.02 even at 270; excluding it would need ~1,463 seeds ≈ 121 h of T4).

**None of this moves the pyramid-vs-kmer4 result**, whose paired CI excludes zero on the negative
side at both seed ranges. The gate decides what the ledger may record, not what is true.

---

Original record, kept verbatim as the evidence the decisions rested on. None was taken AT THE TIME.

1. **Is the fragment arm's `|median| < 0.02` null bar satisfiable at all?** THE REAL FRONTIER. Pooled
   over 50 seeds the null median is **+0.0242** (seeds 0-9 +0.0517, the never-drawn seeds 10-49 +0.0191),
   so the gate FAILS, bootstrap P(clear)=0.317. The spread is INITIALISATION-driven and stable
   (sd 0.0661/0.0649/0.0654), so SE(median) shrinks only as 1/sqrt(n) — SE~0.005 needs ~270 seeds —
   while the point estimate sits ABOVE the bar. **More seeds would more likely CONFIRM a small real bias
   than clear it.** Either re-specify the bar or accept the arm gets no valid verdict. Do NOT spend more
   GPU expecting a clearance. Artifacts: `wiki/fragment_null_ext_summary_2026-10-29.json`, commit 87727f49.
2. **G-B's `blocked_by` edge encodes the wrong predicate.** Self-distillation's only edge is G-A, so marking
   G-A terminal-met would make it eligible — but the standing rule is *forbidden until a learned
   representation LANDS* and G-A closed NEGATIVE. The edge means "wait for G-A to FINISH" where the rule
   means "wait for G-A to SUCCEED". An umbrella flow-down change. First surfaced 2026-10-21, re-confirmed.
3. **Is G-D's 8 GB hardware half discharged?** server1 is a GTX 1070 with exactly 8 GiB and the Kaggle T4
   (16 GiB) ran the sweep, so the 2026-10-07 "server PC is mid-move" rationale is stale. Discharging it
   opens a ~3-day LoRA fine-tune. NOTE the sweep itself is NOT plan-eligible: `--verify-knob` was never
   passed, so 2 of 3 points are `excluded:requested-but-unverified`.
4. **`mvp_criteria.py` cannot read table-format MVP bars**, and the failure is dangerous in one direction.
   `_LIST_ITEM` matches markdown LIST items; G-A and G-C write their six criteria as TABLES, so both parse
   to 0 of 6. G-A gets `unparseable-items` (the rail STOPS — safe), but G-C's section is TIDY so it gets
   `empty-section`, which LICENSES proposing a NEW bar over a FROZEN one. **Do not run `--until-mvp` on
   G-C** until constraint `mvp-bar-parser-blind-to-table-format` expires. The file is shared machine-wide
   and decides what counts as a frozen bar, so fixing it is an acceptance-bar decision.

**And a ROADMAP TENSION worth a decision, not a fix.** `plans/Trait_Decoding_Roadmap.md` orders its seven
phases by ORGANISM COMPLEXITY (E. coli -> multi-organism -> ... -> eukaryotes), while
`dna_decode/eval/regime.py` records the measured discriminator as **POPULATION DESIGN** (constructed vs
natural), and the roadmap predates that correction. The evidence: the yeast segregant cross — a EUKARYOTE,
with splicing and chromatin fully present — decoded **12/12 traits at r 0.46-0.80**, while zero-shot
embeddings on natural populations are 0-for-5. So eukaryotic machinery is NOT what blocks the climb;
natural-population confounding is. The rung after bacteria may therefore be a CONSTRUCTED EUKARYOTE rather
than "animals". Not changed here — re-ordering a roadmap is a user call.


## CORRECTION 2026-11-05 - CLAUDE.md says `.claude/execute-plan-state/` is gone from this host. It is not.

CLAUDE.md's Project-workflow section reads: *"Execution state was tracked at
`.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json` (that directory no longer
exists on this host - the ledger + git history are the live record)"*. **The directory exists and
holds nine state files, the newest dated 2026-11-02.** It is gitignored (`.gitignore:89`), which is
why it is invisible to `git status` and easy to believe absent.

**This cost a real step.** `/execute-plan`'s Step 0 is "check for in-flight execution", and I
skipped it on the strength of that sentence. There WAS a state file for the frame-conditioning plan
(opened 2026-11-01, all steps pending, stopped before Step 1 because the review found three
blockers). My run happened to resume it correctly - same sequential mode, same recorded order
[1,2,4,3,5,6,7] - by independent agreement rather than by reading the file. An orientation claim
that is false in the "nothing is there" direction silently converts a resume into a fresh start.
CLAUDE.md cannot be edited (35,459 of 35,460 words; USER AUTHORITY), so this is the correction of
record.

## OPEN - four user-authority calls from the frame-conditioning run (2026-11-05, measured 11-14)

All four are in the first screen of `wiki/promoteratlas_frame_result_2026-11-01.md`, whose
MEASURED UPDATE 2026-11-14 block carries the detail and the revised recommendations. Still four
authority calls; two of them now rest on measurement rather than judgement. Nothing frozen was
edited: the pinned five, `uv.lock`, the bar, the conditioning artifact and the arm artifact are
byte-unchanged and `FRAME_PARTIAL` still reads `resolved: false`.

1. **The frozen bar grades a POINT MEDIAN where the statistic is unstable.** Unchanged as a fact:
   C3's weaker stratum margin splits 5-5 across the 0.02 bar with four NEGATIVE seeds, range
   [-0.0492, +0.1081], median 0.0004 above it. **My earlier "require an interval" recommendation is
   WITHDRAWN as stated - it is unachievable.** The 95% t-interval is [-0.0149, 0.0517], half-width
   0.0333, WIDER than the 0.02 threshold it must discriminate; ~21 seeds to get the half-width under
   the bar, ~3,253 to separate the estimate from it. The mean is 0.0184, BELOW the bar, while both
   reported estimators (0.0230, 0.0204) are above - the defensible point estimates already straddle
   it. Sign test 5/10, p = 1.000. Also under-specified: the grader reads a DIFFERENCE OF MEDIANS
   (0.0204) while the stability block reports PER-SEED PAIRED margins (0.0230). **Revised: edit no
   frozen artifact. Freeze a NEW bar for the successor run** (23 per-arm bars already live in
   `wiki/`) naming one statistic, grading it with a bootstrap interval, making `UNRESOLVED`
   first-class, and adding the gate in item 2. Drafted, NOT frozen - the criteria are yours.
2. **Ratify or reject the grader's eighth outcome**, `INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE`.
   **Still ratify, but it was INERT here** - re-grading from the artifact gives
   `absorbs_unattributably = False`; with seven outcomes the grader would still have said
   `FRAME_PARTIAL`, not `FRAME_EXPLAINS`. Ratify on prospective grounds (it precedes every
   substantive outcome in `CHECK_ORDER`, so it can only demote), not as having caught this confound.
   **The gate that would have caught it does not exist:** the artifact records
   `control_is_at_least_as_strong_as_the_frame = False` and grep proves it feeds NOTHING in
   `frame_verdict.py` or `freeze_frame_bar.py`. The missing one is
   `INDETERMINATE_CONTROL_PREMISE_UNMET`, reading a boolean already computed.
3. **ANSWERED BY MEASUREMENT - the shuffled-window control is DEAD and the replacement is built,
   run and passing.** `frame_features` takes windows only as the union mask, so a window shuffle
   leaves the two `cds_rows`-derived columns byte-identical; those two ALONE score 0.3183 vs the
   full table's 0.3170 (at least the full table on 6/10 seeds; the window columns below it on
   10/10) and absorb +0.0705 vs the frame's +0.0687 - the "control" would have absorbed MORE than
   the treatment and fired a spurious indeterminate
   (`wiki/promoteratlas_frame_column_split_2026-11-05.json`). Replacement is a CIRCULAR ROTATION of
   the whole annotation (`dna_decode/glm/annotation_shift.py`), matched exactly on gene count,
   gene-length multiset, strand composition, covered bases and density. Run over 5 offsets
   (`wiki/promoteratlas_frame_shift_null_2026-11-14.json`): rotated absorption -0.0009 vs true
   +0.0687, so **attributable_absorption = +0.0696 against a control whose premise holds by
   construction**; rotated frame table predictively dead (0.0114 vs 0.3170); and the in/out collapse
   does NOT reproduce (+0.2459 true vs -0.0164 rotated, both rotated strata ~0.19, agreeing with
   the committed size-matched null's 0.1860). REPORTED, NOT GRADED - the bar does not name it, so it
   settles ATTRIBUTION and does not resolve PARTIAL-vs-EXPLAINS. Your call is now between a
   measured-dead option and a measured-working one.
4. **Candidate 8 is stale** - the G-A table says pending, the arms ledger says its arm is measured.
   Untouched, because it is another session's work; that reasoning is unchanged and correct. **The
   real gap is that nothing DETECTS the divergence:** 7 test files reference `arms_ledger` and none
   reconciles a candidate row against it, while `advance_ranker`'s `next_action` is `table[0]`
   unconditionally - so a stale `pending` row actively misdirects the next `--advance` run, which
   already happened once in this arc. A report-only reconciliation check is the addition.

Unchanged and still blocking Steps 6 and 9 of the vacation plan: **server2 is not on this L2
network** (no ARP entry, four independent routes). Recovery stays armed via
`scripts/retrieve_server2_finetune.py --apply --then-score`.
