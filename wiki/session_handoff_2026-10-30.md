# Session handoff — 2026-10-30

<!-- doc-status: PARTLY SUPERSEDED as of 2026-10-30. Three of the four authority calls below
     were answered later the same day; the banner immediately under this comment is the
     current state and outranks any section text that still reads as open. The one open call
     is 2.4 (the MVP-bar parser). Sections 1 and 4-7 remain LIVE. Re-derive every number
     from the cited artifact before repeating it. No repo-wide doc-status convention exists
     in code (`doc-status` returns 0 hits in dna_decode/ and scripts/), so this follows the
     global first-30-lines rule rather than a local schema. -->

> **SUPERSEDED IN PART — three of the four authority calls were ANSWERED on 2026-10-30 (later the
> same day). Read this block before anything below it.**
>
> | call | decision | state |
> |---|---|---|
> | §2.1 null bar | **spend the seeds (~270)** | **GATE PASSES at 160 seeds — median +0.0112, branch A. See `wiki/fragment_null_160seed_interim_2026-10-31.md`.** b3/b4 both RUNNING to 270 (b3 FAILED on Kaggle once and was retried 2026-10-31) |
> | §2.2 G-B edge | **fix the predicate** | DONE `0750ec63` — G-B is blocked again |
> | §2.3 G-D hardware | **hold the edge** | HELD. **Sweep half now DISCHARGED 2026-10-31** — 3 of 3 points plan-eligible, `wiki/glm_generator_sampling_sweep.json` emitted, constraint EXPIRED -> RE-OPEN, result still negative (`f9c47d7e`). **This supersedes §4's "do not create that file" trap** — see the note there. G-D is still blocked BY THIS CALL, not by hardware |
> | §2.4 MVP parser | *still open* | unchanged — do not run `--until-mvp glm-condition-conditioning` |
>
> §2.1 carried a **correction**: this handoff framed it as "re-specify or accept no verdict, both
> legitimate". The 50-seed bootstrap CI `[-0.0126, +0.0466]` **contains zero**, so the null is not
> measurably dirty and the bar is failing on a point estimate narrower than its own SE. A third
> option existed — re-specify it as a *leakage* bar, set where leakage actually lands — and
> unblocking the verdict makes the arm **fail** (`NO_GAIN`), so that option was anti-flattering
> rather than motivated. The user chose seeds anyway, a legitimate call on the same evidence.
>
> **The pyramid-vs-kmer4 science does not depend on any of this** and is unchanged: the paired CI
> excludes zero on the negative side at both seed ranges. The gate decides what the LEDGER may
> record, not what is true.

**Read in this order:** (1) the DO-NOT-REDO list, because it protects ~7 hours of completed compute;
(2) the four authority calls, because work is waiting on them and none is an executor task; (3) the
measured frontier. Everything else is reference.

Prior durable state, all current as of this file:
`~/.claude/soraya-resume/dna_decode-77502f867d.json` (rewritten 2026-10-30 17:18) ·
`NEXT.md` §FOUR OPEN AUTHORITY CALLS (commit `effac868`) · `wiki/decisions-log.md` tail.

---

## 1. DO NOT REDO — all finished and committed

| work | cost | where the result lives |
|---|---|---|
| 40-seed null extension (seeds 10–49) | 197.8 min Kaggle T4 | `wiki/glm_fragment_pyramid_arm_kaggle_ext_2026-10-29.json`, `wiki/fragment_null_ext_summary_2026-10-29.json` — commit `87727f49` |
| G-D sampling sweep (T = 1.0 / 0.8 / 0.6) | ~8 min Kaggle T4 | `wiki/glm_sampling_sweep_2026-10-29_*.json` — commit `1733d8ae` |
| promcalc decomposition rebuild | ~5 h equivalent | `D:/dna_decode_cache/promcalc/decomposition.rebuilt.csv.gz` — validated twice |
| 10-seed fragment arm + server1 replication | 54 + 81 min | `wiki/glm_fragment_pyramid_arm_kaggle_full_2026-10-28.json` |

**Also do not create `wiki/glm_generator_sampling_sweep.json` by hand.** That exact filename is the
expiry test for constraint `gd-finetune-gated-on-sampling-sweep-not-compute`. Writing it would fire a
LoRA fine-tune's unblock trigger on a sweep the project's own gate excludes (see §4).

> **RESOLVED 2026-10-31 — and the rule above is UNCHANGED, not lifted.** That file now EXISTS, and
> legitimately: the sweep was re-run with `--verify-knob` (`glm-sampling-sweep-vk`), all 3 grid points
> came back plan-eligible, and the artifact was emitted **BY CODE** from the real outputs via
> `scripts/emit_glm_sweep_unblock_artifact.py` — which asserts the discharge condition and refuses
> otherwise. "By hand" is still forbidden and is still exactly what that emitter exists to prevent;
> two of its refusals are pinned by `tests/test_glm_sweep_unblock_emitter.py`. **Do not delete the
> file thinking this trap was violated, and do not re-emit it by hand.** The gate no longer excludes
> the sweep, so the premise of the trap ("a sweep the project's own gate excludes") no longer holds.
> The RESULT is unchanged and still negative: no grid point beats the Markov null.

---

## 2. The four authority calls — work is waiting on these

Recorded verbatim in `NEXT.md`. None was taken; none is an executor task.

### 2.1 Is the fragment arm's null bar satisfiable at all? — THE REAL FRONTIER

Not a compute question. Measured:

| sample | n | null median | sd |
|---|---|---|---|
| seeds 0–9 (twice-published) | 10 | +0.0517 | 0.0661 |
| seeds 10–49 (never drawn before) | 40 | +0.0191 | 0.0649 |
| **pooled 0–49** | **50** | **+0.0242** | 0.0654 |

Bar is `|median| < 0.02`. Pooled **fails** (bootstrap P(clear) = 0.317, CI `[-0.0126, +0.0466]`).
The fresh 40 clear it as a point estimate only — P = 0.490, a coin flip, with a 0.0009 margin about
one fourteenth of its own SE (0.0129). **Quote the pooled 0.0242, never the 0.0191**; the latter is
the flattering subsample.

**Why more seeds will probably not help:** the spread is driven by random *initialisation*, not
sample size, and is stable across all three samples (0.0661 / 0.0649 / 0.0654). SE(median) shrinks
only as 1/√n — SE ≈ 0.005 needs ~270 seeds — while the point estimate sits **above** the bar. More
seeds would more likely **confirm a small real positive bias** than clear it.

**DECIDED 2026-10-30: spend the seeds.** Target n = 270, where SE(median) reaches 0.0050.

- Pre-registration, written BEFORE any of these seeds existed and fixing all three outcome
  readings: `C:/Users/Laptop Portable/pythonprojects/dna_decode/wiki/fragment_null_270seed_prereg_2026-10-30.md`
- Kernel: `scripts/kaggle/fragment_null_batch_kernel.py` — four batches of 55 tiling seeds 50..269,
  `PINNED_SHA` UNCHANGED so all 270 seeds come from byte-identical code.
- **Resume with the driver. It is idempotent and safe to re-run from any session:**
  `uv run python scripts/fragment_null_batch_driver.py step` — pulls whatever finished, launches the
  next batch when a GPU slot frees. Also `status` and `pool`.
- Slugs `glm-fragment-null-b1..b4`. b1 (seeds 50..104) and b2 (105..159) launched 2026-10-30 and
  were RUNNING at handoff; b3 and b4 wait on the GPU cap. ~4.5 h per batch, two at a time, ~9 h all.

**THE FIRST FORECAST IS NOW FALSIFIED BY MEASUREMENT — do not repeat it.** (1) It gave ~80%
that the gate would still fail at 270. It PASSED at 160 (median +0.0112): the first fifty
seeds were the high samples, and 110 fresh ones land at ~0. Quote the measurement, not the
forecast. (2) 270 seeds buys a **point estimate, not a confident
resolution** — the CI is forecast to straddle 0.02 even at 270, and excluding the bar from above
would need ~1,463 seeds, about 121 h of T4. **Normal theory is NOT the decision statistic**: it and
the measured bootstrap disagree about zero at n=50, so the verdict comes from a bootstrap recomputed
on the pooled per-seed values.

### 2.2 G-B's dependency edge encodes the wrong predicate

`glm-self-distillation`'s only `blocked_by` edge is G-A. Marking G-A `terminal-met` would make
self-distillation eligible — but the standing rule is *forbidden until a learned representation
**LANDS***, and G-A closed **negative** (+0.0122 against a 0.05 bar). The edge means "wait for G-A to
FINISH" where the rule means "wait for G-A to SUCCEED" — identical except in the case that occurred.
Fixing it is an umbrella flow-down change. First surfaced 2026-10-21; re-confirmed 2026-10-29.

**FIXED 2026-10-30 (`0750ec63`).** The `blocked_by` cell now carries a free-text predicate beside
the slug: `glm-learned-representation-2026-10-07, a learned representation LANDS (G-A terminal-met
is NOT sufficient)`. An unmatched predecessor is fail-closed in `advance_ranker` — the same
mechanism that holds G-D — so the block can only be cleared by a human deleting the clause, which is
the authority the rule requires. The slug is RETAINED, so the graph edge, cycle detection and the
critical path are unchanged. **Do not "tidy" the free text back to a bare slug;** the umbrella's
Requirements Flow-down section now says why. Verified against a baseline taken before the edit:
G-B moved from `+ eligible` to `x blocked:predecessor`, with G-A, G-C, G-D and G-E unchanged.

### 2.3 Is G-D's 8 GB hardware half discharged?

server1 is a GTX 1070 with exactly 8 GiB; the Kaggle T4 (16 GiB) ran the sweep. The 2026-10-07
*"server PC is mid-move"* rationale is three weeks stale. Discharging it opens a ~3-day LoRA
fine-tune, which is why it is not an executor call. The frozen plan
`plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/` (still in `plans/`, not archived) Step 6 explicitly says
**do NOT change G-D's `blocked_by` edge for the 8 GB GPU** — honour that until the user
overrides it.

**HELD 2026-10-30 — the edge is untouched, as the frozen plan requires.** Its *other* half, the
`--verify-knob` attestation blocker named in §3, is FIXED (`c6a49923`): `run_cell`'s argv now passes
the flag unconditionally, in the builder and in the tracked emitted kernel, with three guard tests
(the builder had none at all, which is why it shipped). The sweep is re-pushed as
**`glm-sampling-sweep-vk`** — a NEW slug, because `glm-sampling-sweep` still holds the retrieved
10-29 results and Kaggle keeps one output per kernel. **It is NOT yet launched:** Kaggle refuses a
third concurrent batch GPU session, so it waits for b1 or b2 to finish. Re-push it then; ~10-15 min.

### 2.4 `mvp_criteria.py` cannot read table-format MVP bars

`mvp_criteria._LIST_ITEM` matches markdown **list** items; G-A and G-C write their six criteria as
markdown **tables**, so both parse to **0 of 6**. The failure is asymmetric and the tidier ledger is
the exposed one:

- G-A → `unparseable-items:4` → the R18 rail **stops**. Safe, but only because its MVP section happens
  to contain four unrelated prose bullets (the GA-H1..H4 lines).
- G-C → `empty-section` → that verdict **licenses proposing a new bar**, over a frozen six-row bar
  carrying a documented D1 authority decision.

**Do not run `/soraya --until-mvp glm-condition-conditioning` until this is fixed.** Registered as
constraint `mvp-bar-parser-blind-to-table-format` with a `test-exit-0` expiry that passes when G-C
parses to 6. Not fixed here: the file lives at `C:/my_skills/skills/soraya/scripts/mvp_criteria.py`,
is shared by every project on this machine, and decides what counts as a frozen bar — an
acceptance-bar surface.

---

## 3. The measured frontier — where the GLM actually stands

**The pre-registered comparison loses, and more sharply than before.** Pyramid vs kmer4:

| seeds | delta median | wins | ci95 | clears margin |
|---|---|---|---|---|
| 0–9 | −0.0127 | 2/10 | `[-0.0205, -0.0040]` | false |
| 10–49 | **−0.0061** | **13/40** | **`[-0.0133, -0.0037]`** | false |

Half the magnitude, but the CI now **excludes zero** — a smaller, more definitive loss to a 4-mer
baseline. `status: PARTIAL_SEED_RANGE`, `verdict: null`, correctly: the arm refuses a pre-registered
verdict on a shifted seed range, so none of this overturns the formal verdict.

**The encoder is NOT broken, and this matters because it was published otherwise.** kmer4 holds
ρ 0.1009–0.1241 across all ten seeds while the pyramid holds 0.0839–0.1249. The comparator therefore
validates the task, the splits, the labels and the scoring, and the pyramid is learning substantially.
What is **not** validated is the torch path specifically — kmer4 goes through `kmer4_features` +
closed-form `_ridge_fit_predict`, the pyramid through `build_masked_pyramid` + torch.

**`learned_contribution` is NOT independent good news.** It reads +0.0965 with 37/40 seeds positive
(was +0.0419, 9/10) — but it is defined as `pyramid_spearman − null_spearman` on the same init and
split, so a lower null *mechanically* raises it. Its rise and the null's fall are the **same fact**.
It stays DIAGNOSTIC, never a bar.

**The G-D sweep result.** Nulls reproduce the recorded 0.5990 (`kmer`) / 0.6679 (`both`) **exactly**
at all three points, both modes (`sweep_gates.null_reproduced` ok, deltas 0.0). Temperature is a real
lever — the kmer gap over the null falls 0.1365 → 0.0401 → 0.0733 across T = 1.0 → 0.8 → 0.6, best at
T = 0.8, a ~3.4× reduction. But **no point beats the null in either mode**, and no `both` value clears
0.6679. So G-D's own `:5` note does **not** fire: the family's premise is not weakened, no re-rank due.

**Trigger not discharged, and the reason is a kernel defect not a result.** `--verify-knob` appears
**zero** times in `scripts/kaggle/build_glm_sampling_sweep_kernel.py` or the emitted kernel, so
`sweep_gates.verdict_eligibility` returns `excluded:requested-but-unverified` for T = 0.8 and 0.6 and
only `eligible:default-sampling` for T = 1.0. Two of three points are not plan-eligible. The knob
demonstrably moved the output (0.6391–0.7355 with the null pinned), so this is a missing
**attestation**, not an inert knob. Re-run with `--verify-knob` to discharge.

---

## 4. Compute roster — MEASURED, not assumed

| box | GPU | verdict |
|---|---|---|
| Kaggle | T4, 16 GiB, CC 7.5 | **use this for all GPU work** |
| server1 | GTX 1070, 8 GiB, Pascal CC 6.1 | only local box worth a real GPU job |
| server2 (`DESKTOP-SJG795O`) | **GTX 860M, 4 GiB** — same as the laptop | CPU/parallel work only |
| this laptop | GTX 860M, 4 GiB, CC 5.0 | too weak |

**Standing user rule:** *"use kaggle, forget server2 for heavy usecase, only for parallelising simpler
work"*. server2's GPU was long recorded as "fairly good" and is not — measured 2026-10-29 as the same
chip as the laptop.

**server2 is reachable by `ssh server2` even when all its Claude sessions read `offline`.** A dead
Claude session is not a dead machine: check `ping DESKTOP-SJG795O` → `ssh server2 hostname` before ever
reporting it unavailable. Its `D:` is the project's 4.6 TB Seagate drive and holds `dna_decode_cache`
(refseq, promcalc, embeddings), `dna_decode_data`, `hf_cache` (including the GENERator weights),
`uv_cache` and `tools`.

Mechanics that cost real time: the remote shell is **cmd.exe**; pipe Python to stdin
(`ssh server2 "python - args" < script.py`) rather than fighting four layers of quoting; Git-Bash
`/tmp` is not Windows `C:\tmp`; filter the SSH banner and PowerShell's `#< CLIXML` blob.

---

## 5. Traps that cost time this session — read before writing checks

1. **A verification criterion written against the PRE-FIX state.** Five instances across three
   sessions, three in this one. The worst: a plan step asserted that `NO_GAIN` / `+0.0001` /
   `null_clean=true` keep "unchanged counts" — impossible for a row whose job is to quote those
   strings to say they stand; the assertion fired *after* the write landed. **Use `git diff --numstat`
   insertions-only and "pre-existing occurrences intact", never "count unchanged".** Memory:
   `feedback_read_the_verification_criterion_against_the_fixed_state`.
2. **Never equality-compare worktree text to a git blob.** These files are **stored as LF** and the
   worktree is **CRLF**, so the comparison fails and reads like data loss. Measure endings in Python
   per file and assert the *predicted direction and magnitude*, not equality.
3. **Grep the owner's vocabulary, not the artifact's name.** `fragment_pyramid` returns 0 hits across
   `project_state/` and `features/` and nearly became "no family owns this arm". The family exists —
   `glm-learned-representation-2026-10-07` — and calls it "pyramid arm" / "multi-resolution". Memory:
   `feedback_search_the_owners_vocabulary_not_the_artifacts_name`.
4. **An instrument answers a narrower question than the claim.** `ListAgents` measures *session*
   state, not machine state; `external_root.py` probes *this host's* drives, not the world. Name the
   instrument's scope, then run the cheap refutation. Memory:
   `feedback_instrument_scope_is_not_claim_scope`.
5. **Counts cannot validate a promcalc table.** `decompose.json` reports `coverage: 1.0`, 44,106
   tiles, `sites_per_tile_mean 73.0`, zero errors — on a gzip that **cannot be read to the end**.
   Stream the artifact; never trust the manifest.
6. **Kaggle OAuth expires ~12 h** and the error misdirects ("most likely a wrong kernel slug"). The
   tell is HTTP 401 on `kernels_list` *plus* the same failure on a known-good slug. Fix is
   `uv run kaggle auth login --force` **in the user's own terminal** — credentials never route through
   a session.
7. **Kaggle retrieval for repo-cloning kernels.** The output listing is ~3,577 files, so a plain pull
   429s — that was always self-inflicted. Walk `kernels_list_files` at `page_size=200` (1000 → HTTP
   400), artifacts are on the **last** page, then one `kernels_output` call with that `page_token`.
   Matched files download from signed storage URLs, not the throttled API.
8. **Kaggle allows exactly TWO concurrent batch GPU sessions, and a refused push used to report
   SUCCESS.** Measured 2026-10-30: `kernels_push` does not raise on refusal -- it returns
   `ref=''`, `url=''`, `versionNumber=0` and the reason in **`.error`**
   (`"Maximum batch GPU session count of 2 reached."`). The helper printed the two empty fields and
   returned 0, and the follow-up `status` then 404'd, which is the recorded signature of a *wrong
   slug* -- so the diagnosis went to the wrong place until the response was dumped field by field.
   **Fixed `1f63590e`:** a refusal now exits 1 and names the cause. Nothing is created on a
   refusal, so a retry is always safe and cannot duplicate a run.
9. **A kernel that does not exist raises a MISDIRECTING `ValueError`, identical in shape to expired
   auth.** `kernels_status` on an absent kernel says *"Permission 'kernels.get' was denied ... the
   most likely cause is a wrong kernel slug"* -- not an HTTP 404. An expired token produces the
   same thing, so mapping it straight to "absent" makes an auth outage look like "nothing launched
   yet". Discriminate by repeating the call on a **known-good slug** first; that is what
   `fragment_null_batch_driver._api_healthy` does, and it refuses rather than guessing.
10. **`kernels_output(file_pattern=...)` is a REGEX, not a glob, and it does NOT avoid the 429.**
   It filters what is DOWNLOADED, not what is LISTED, so the call still pages the full ~3,600-file
   listing server-side and throttles anyway; a loose pattern also matched a false positive *inside
   the cloned repo*. Trap 7 above therefore stands, and is now encoded in
   `scripts/fragment_null_batch_driver.py::_pull` rather than living in prose.
11. **This harness collapses doubled backslashes inside Bash heredocs.** It bit THREE times in one
   session, including once while writing this very entry: a doubled-backslash-n inside a patch
   string became a REAL newline and wrote a file with an unterminated string literal, and
   `C:\\Users` collapsed to `C:\Users` and died as a truncated-escape
   error. Write Windows paths with forward slashes, as the rest of this repo does, or use the
   dedicated file-write tool for anything containing a backslash. Recovery is
   `git restore --source=HEAD -- <one file>`, which is targeted and not on the forbidden list.
12. **`python3` is the Windows Store stub here. Use `python`.** It silently prints an install nag
   instead of running.

---

## 6. Hard constraints

- **CLAUDE.md is at 35,459 of its 35,460-word ceiling — one word.** Nothing can be added until
  something is compressed out. `uv run python scripts/claude_md_weight.py` certifies 43 bullets /
  24,106 words as long AND provably stored in `wiki/`; choosing which to reduce to a pointer changes
  what every future session reads first, so it is a user call. `tests/test_claude_md_growth_ceiling.py`
  is the rail; raising the number instead of compressing is what it exists to prevent.
- **Every heavy file lives on `D:`, never `C:`.** `D:` attachment is per-session and never inherited —
  check it, and note the laptop's `D:` (External1TB, 827 GB free) is a *different disk* from server2's
  `D:` (the 4.6 TB Seagate).
- **The frozen five must never be modified:** `dna_decode/eval/amr_rules.py`,
  `dna_decode/data/calibrated_amr_rules.json`, `dna_decode/data/mic_tiers.py`,
  `dna_decode/data/shipped_decoder_surface.py`, `dna_decode/eval/cohort_manifest.py`. `uv.lock` is
  also a frozen reproducibility unit.
- **Forbidden:** `git reset --hard`, `git clean -fd`, `git checkout .`
- **`torch` is NOT in the default deps.** It lives in the `[forward]` extra; a plain `uv sync` breaks
  every GPU path with a bare `ModuleNotFoundError`. Use `uv sync --extra forward`. `[forward]` and
  `[ml]` are mutually exclusive in `[tool.uv].conflicts`.
- **Never run the full pytest suite while a torch job holds a model** — it dies with an access
  violation that reads like a test failure.

---

## 7. A roadmap tension worth a decision

`plans/Trait_Decoding_Roadmap.md` orders its seven phases by **organism complexity** (E. coli →
multi-organism → … → eukaryotes). `dna_decode/eval/regime.py` records the *measured* discriminator as
**population design** (constructed vs natural). The roadmap predates that correction.

The evidence: the yeast segregant cross — a **eukaryote**, with splicing and chromatin fully present —
decoded **12/12 traits at r 0.46–0.80**, while zero-shot embeddings on natural populations are
0-for-5 de-confounded. So eukaryotic machinery is **not** what blocks the climb; natural-population
confounding is. The rung after bacteria may therefore be a **constructed eukaryote** rather than
"animals".

What *does* transfer upward: the method (pre-registered bars, null controls, baseline panels), the
generative loop (propose → score with a validated oracle → select), and the molecular layer — the
`forward` cell is validated across TEM-1, CTX-M-14, PTEN, HIV RT, SARS-CoV-2 Mpro and a fungal enzyme.
What does **not**: the regulatory layer. A bacterial promoter is a −10 box, a −35 box and a spacer
(~100–300 bp); a eukaryotic one is TATA + ~6 general factors + Mediator + enhancers up to ~1 Mb away +
chromatin. Different object, not a scaled one — which is why `promcalc` does not transfer up, the
mirror of the recorded finding that AlphaGenome's architecture does not transfer down. Recorded in
`scripts/promcalc_decomposition_arm.py` beside the existing sigma70 statement.

---

## 8. Recent commits

```
effac868  record four open authority calls and a roadmap tension
03ff0cd6  fill the execution log's commit field with the SHA it records
1313cae0  record the fragment rung in its family ledger, and write down two negatives' boundaries
87727f49  resolve the null at 40 fresh seeds: the +0.05 bias was the ten shared inits, gate still fails
c9009335  record the sweep in G-D, and register an MVP-bar parser defect
1733d8ae  land the G-D sampling sweep: lower temperature helps 3.4x and still loses to a 3-mer Markov chain
726fc0f5  correct the sweep grid against its frozen plan, and dispatch both GPU arms to Kaggle
67b1891f  discharge two G-D constraints and record why the third stays binding
```

Working tree is clean except one stray `plans/Mission_Reconciliation_A_Native_Clonality_Route_Plan/technical-plan.md.v1bak`
that predates this session.

---

## 9. If you want a next action rather than a decision

Everything high-value is gated on §2. The honest reading is that **the frontier is a decision, not a
task** — so the right move is to answer §2.1 rather than to find work around it.

If you want a task anyway, the two cheapest that are genuinely useful and need no authority:

- **Re-run the G-D sweep with `--verify-knob`** at every non-default grid point and confirm
  `knob_verified_end_to_end == true`. ~10 min of Kaggle T4. That alone discharges §2.3's blocker half
  and makes the existing numbers plan-eligible.
- **Add a positive control for the torch path.** The comparator already proves the task is learnable;
  nothing proves this particular conv encoder can hit a known-achievable number on *any* task. That
  separates "our torch path has a subtle bug" from "conv-on-fragments is the wrong inductive bias".
  Prefer an **in-domain** control — the splice-site task an external source suggested is eukaryotic and
  its labels are annotation-derived, which trips this repo's own gate G1.
