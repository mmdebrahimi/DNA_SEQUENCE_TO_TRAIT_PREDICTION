# Plans Index
<!-- Auto-maintained by /save-plan. Do not edit manually. -->

## [plan_file: plans/Two_Machine_Operating_Contract.md] 2026-05-26
**Status:** active (in force; first production test = cef audit-aware closeout integration)
**Summary:** Durable operating contract between Codex (Precision 7780) + Claude (GTX 860M). Ratifies Codex's 2026-05-26 division-of-labor proposal + 5 amendments from same-day /brainstorm round.
**Key decisions:**
- 2 lanes: Codex = execution; Claude = discovery + planning + contract-lock conversion
- Handoff gate = workflow STATE (not a permanent 3rd lane) with 5 mandatory checks (push status / pull / sync check / locked-parameter coverage / contract tests run)
- `Contract Locks` section mandatory in every locked-design memo; each parameter MUST point to a regression test OR `KNOWN_DIVERGENCE_TARGETS` marker; empty target = unacceptable
- Single-commit rule: lock + enforcement target land in the SAME commit (historical splits grandfathered)
- Source-of-truth artifact ownership (NOT author): Claude owns planning contracts; Codex owns execution artifacts; cross-lane proposals land as `provisional`
- 5 falsification triggers in §6 (the contract's own falsification test)
- Pinned by 8 self-tests at `tests/test_two_machine_operating_contract.py`

---

## [plan_file: plans/EP_4_Pathotype_Idea_Anchor_Draft.md] 2026-05-26
**Status:** candidate (pre-staged; trigger = cef audit-aware closeout pushes to origin)
**Summary:** Ready-to-paste `/idea-anchor` candidate sentence for E. coli pathotype prediction (Rank 1 EP-4 entry per `plans/EP_4_Non_AMR_Phenotype_Candidates.md`). Bridges the scoping memo to the actual planning-chain invocation.
**Key decisions:**
- Sentence shape mirrors v1 "Honest AMR resistance predictor" framing: genome FASTA → multiclass pathotype call (EPEC/EHEC/ETEC/UPEC/EAEC/commensal) + audit-grade provenance + CGE VirulenceFinder side-by-side comparison
- CGE VirulenceFinder benchmark choice is the direct analog of AMRFinderPlus comparison for AMR
- Concentrated mechanism (pathotype = acquired-gene clusters) → reuses NT mean-pool + XGBoost architecture per 2026-05-17 cross-drug finding
- 3 pre-emptive risks named: label inconsistency on EnteroBase, multiclass audit-framework extension, external-benchmark-discipline-lift decision

---

## [plan_file: plans/Cef_Audit_Aware_Packet_Design.md] 2026-05-26
**Status:** active (IN FLIGHT on Precision 7780 per Codex 2026-05-26 ask)
**Summary:** Smallest-credible design for the cef audit-aware closeout that matches cipro release discipline without reopening broader scope. 5 artifacts; ~3.5-4 hr Precision 7780 compute; zero MODEL architecture change + one new drug-agnostic merge component (~150-300 LOC).
**Key decisions (post-/brainstorm revision 2026-05-26):**
- D1: `signal_quality = clean_count / len(merged)` (matches cipro; not heuristic clean/(clean+opacity+suspect))
- D2: SUSPEND threshold LOCKED at 0.40 (matches cipro 2026-05-17 calibrated value; supersedes heuristic 0.50)
- D3: 4 canonical example strains: R + S + 562.28389 (FP miss) + 562.7695 (FN miss) — both misses surfaced for release honesty
- D4: Pre-committed release-wording branches per verdict (RUN_FULL_AND_CLEAN / MIXED / SUSPEND_CONDITION_4)
- D5: JSON schema contract pinned to `pipeline._load_audit_verdict`'s expectations: top-level `gate_verdict` + per-strain `noise_class` + `mic_tier` + `primary_mechanisms` (plural list) + `co_resistance_modifiers`
- D6: After this ships → EP-4 pathotype prediction is the next jump per the trait-decoding roadmap

---

## [plan_file: plans/Trait_Decoding_Roadmap.md] 2026-05-26
**Status:** candidate (long-horizon spine; doesn't fire as a single execution)
**Summary:** Connective tissue between current v0/v0.1 ship state and original 2026-05-11 goal ("DNA input → phenotype/trait identification at gene level"). 7 phases with per-phase terminal claims + dataset prerequisites + architecture-transferability bets + phase-transition gates.
**Key decisions:**
- Phase 0 (cipro cached) ✓ → Phase 1 (cipro genome-input) ✓ → Phase 2 (multi-drug E. coli; cef in flight) → Phase 3 (multi-organism Klebsiella) → Phase 4 (non-AMR bacterial phenotypes; pathotype prediction recommended) → Phase 5 (multimodal DNA + image/RNA) → Phase 6 (eukaryotic; fungal → plant → human GWAS)
- Each phase requires `/idea-anchor + /project-init` when it fires
- Phase ∞ (decode arbitrary trait from DNA) is aspirational; not a literal target

---

## [plan_file: plans/EP_4_Non_AMR_Phenotype_Candidates.md] 2026-05-26
**Status:** candidate (gated on cef audit-aware closeout)
**Summary:** First conceptual jump out of AMR per the trait-decoding roadmap Phase 4. Surveys 6 candidate bacterial non-AMR phenotypes (growth rate / pathotype / biofilm / plasmid stability / lactose / metabolic auxotrophy) + ranks.
**Key decisions:**
- Rank 1: **pathotype prediction** (EnteroBase substrate; multiclass EPEC/EHEC/ETEC/UPEC/EAEC/commensal) — concentrated mechanism + huge public dataset + clinical relevance + minimal architectural change
- Rank 2: biofilm formation (mixed mechanism; defer until EP-1.5 architecture decision)
- Rank 3-6: skipped — too AMR-adjacent OR too easy OR no buyer signal

---

## [plan_file: plans/Cef_V0_1_Promotion_Slice_Plan.md] 2026-05-26
**Status:** partially superseded (by plans/Cef_Audit_Aware_Packet_Design.md Steps 2-5)
**Summary:** Pre-execution cef promotion plan. PARTIALLY superseded by Codex's overnight ship (Codex went bigger than mini N=12 reuse — built fresh 67-strain cohort) + the higher-fidelity Cef_Audit_Aware_Packet_Design.md.

---

## [plan_file: plans/DNA_Decoder_v0_1_Genome_Input_Cipro_Plan.md] 2026-05-25
**Status:** executed (Codex on Precision 7780; cross-path concordance 4/4)
**Summary:** Codex's authoring of the v0.1 cipro genome-input slice 1 plan; shipped 2026-05-25 with full validation.

---

## [plan_file: plans/Post_V0_EP_Ladder_Plan.md] 2026-05-25
**Summary:** LOCKED EP ladder from cipro v0 (step 0) toward the long-term goal "DNA input → phenotype/trait identification at gene level." Convergent recommendation across 3 /brainstorm + 1 /probe + 1 /review rounds. Supersedes the misframed "Commercial Discovery Prerequisites" draft from earlier this session (NOT saved).
**Key decisions:**
- D1: Use Evidence Packet labels (EP-0, EP-1A/B/C, EP-1.5, EP-2, EP-2.5, EP-3) per the 2026-05-15 framing reset; "Phase" labels stay retrospective-only.
- D2: Each EP has an explicit terminal claim to bound scope (resolves the 2026-05-11 unbounded "decode any DNA" problem per-EP).
- D3: EP-2 (multi-drug) forks by mechanism class — cef (concentrated) reuses v0 arch; tet (distributed) requires EP-1.5 architecture decision first.
- D4: EP-1C honest-output framework generalization is its own EP, not background.
- D5: EP-1B external benchmark is its own EP (10 public AST-labeled E. coli genomes vs AMRFinderPlus + RGI).
- D6: EP-4+ (non-AMR / eukaryotic / multimodal) = separate `/idea-anchor + /project-init` cycles per paired dataset; NOT roadmap commitments.
- D7: "Commercialize" REMOVED from ladder — it's a downstream outcome, not a development phase.

---

## [plan_file: plans/v1_Horizon_Framing_Plan.md] 2026-05-24
**Summary:** Drafts the "AI DNA decoder tool at maturity" framing as input for a future `/idea-anchor` + `/project-init` cycle. Resolves the long-overdue Pending Decisions row from 2026-05-17 ("/idea-anchor + /project-init proper entry"). NOT auto-invoking those skills; producing the candidate sentences + scoping doc the user can run them against. 3 candidate framings (1=Honest AMR predictor, 2=Bacterial G2P platform, 3=Multimodal decoder) + 7 open questions + ready-to-paste /idea-anchor invocation.
**Key decisions:**
- D1: 3 framings cover narrow / medium / broad scope; user picks one before /idea-anchor fires
- D2: Recommended Framing 1 (Honest AMR resistance predictor) — tightest alignment with v0 + v0.1 paths, honest-output discipline as moat
- D3: v1 success criteria (6 hard gates) + v1 explicit non-goals + v2+ horizon candidates
- D4: Pending answers from user on framing / timeline / compute / honest-output gate / ship vehicle / license / co-authorship

---

## [plan_file: plans/v0.1_Ingestion_Contract_Plan.md] 2026-05-24
**Summary:** v0.1 ingestion contract covering BOTH candidate paths (Path G = real-genome-input cipro decode; Path C = cef-cached expansion). Planning-only. 7 design decisions + 5 risk flags + 5-wave sequencing.
**Key decisions:**
- D1: Both paths share the v0 JSON output schema (additive, not breaking)
- D2: leave_one_accession_out CV mandatory for any new training
- D3: Top-K=10 mechanism recovery is target, NOT hard gate (matches v0 relock)
- D4: Path G compute defaults to Precision 7780 (GTX 860M can't fit NT v2 100M)
- D5: Path C gated on BV-BRC cef categorical-MIC feasibility check
- D6: --audit-merge-json required for canonical reporting (carries from v0)
- D7: Per-genome scoring; batch mode is v0.2+

---

## [plan_file: plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md] 2026-05-22 (EXECUTED 2026-05-23)
**Summary:** Verdict-conditional ship-path plan locked the response to all 4 falsifier outcomes BEFORE results landed. Codex on Precision 7780 ran the falsifier 2026-05-23 — verdict FAIL. Step F (FAIL branch) fired: v0 shipped 2026-05-24 with documented scope-limit. v0 spec RELOCKED to match the cached-strain implementation. Pre-commitment discipline (verdict-vs-budget LESSON 2026-05-14) prevented motivated-reasoning around the FAIL result.
**Key decisions (executed):**
- D1: Pre-commit verdict-conditional response before results land
- D2: PASS + FAIL both ship v0 — only attribution_scope_confidence default differs (FAIL fired)
- D3: Saturation gate + lineage-confound are co-causes, not either/or
- D4: Mash-cluster gated on PASS — DID NOT FIRE (verdict was FAIL)
- D5: Pre-written commit-message templates per branch
- D6: Schema lock between Codex's runner + Claude's consumer

---

## [plan_file: wiki/cipro_bounded_falsifier_coordination_plan_2026-05-22.md] 2026-05-22
**Summary:** Comprehensive coordination plan for the cipro bounded-falsifier experiment between Claude (GTX 860M laptop) + Codex CLI (Precision 7780, RTX 3500 Ada). Codex executes the falsifier; Claude owns subset selection + leakage check + scope-limit doc template. 12 sections covering roles, naming convention, 12-strain subset (4 ERS control / 4 ELX-family failure / 4 all-negative-Δ), parallel leakage check (<5s), diagnostic exports spec (saturation_flag via baseline_proba_R + max_abs_delta), Mash-cluster gating, v0 ship-or-document-scope-limit fork, verdict matrix.
**Key decisions:**
- D1: 12-strain subset = 3 buckets × 4 strains, each testing a distinct failure mode (control / failure / negative-delta)
- D2: Leakage check on `GCA_025200635.1` runs BEFORE/PARALLEL to falsifier, blocks interpretation if LOSO same-genome leakage present
- D3: Diagnostic exports add logit deltas + max-abs-Δ across all genes — disambiguates saturation from lineage confound (the gap /brainstorm caught)
- D4: Mash-cluster work GATED on PASS verdict only — saturation diagnosis supersedes
- D5: FAIL path ships v0 WITH documented `attribution_scope_confidence` field — symmetric success outcome per north star
- D6: Codex owns runner mechanics; Claude owns subset + leakage check + v0 scope-limit doc; disagreements default to Codex's spec on runner, Claude's plan on subset

---


## [plan_file: Phase2_Framing_Brainstorm_Technical_Plan.md] 2026-05-17
**Summary:** Apply 5 /review edits to `plans/Phase2_Framing_Brainstorm_Plan.md` to convert it from "framing brainstorm with implicit preference" → "framing brainstorm with explicit 3-candidate slate + decision gates." User-confirmed scope: add classifier-tier Candidate 3; Candidate 1 merge sentence permits Months-including-Databricks budget. 2 steps, documentation-only.
**Key decisions:**
- D1: Documentation-only plan — no code changes; downstream code work after /idea-anchor + /project-init in separate session
- D2: Third candidate is classifier-tier (post-feasibility-census A1), not exit-tier — user-confirmed via AskUserQuestion
- D3: Candidate 1 merge sentence permits full-budget scope (Months including Databricks) — user-confirmed
- D4: Strip `(preferred)` label; framing brainstorm exposes choice, doesn't pre-rank
- D5: Open Questions Q1-Q3 promoted to gate (heading + intro paragraph)
- D6: Honest Read on User Intent reframed as hypothesis-to-confirm, not premise

---

## [plan_file: Phase2_Framing_Brainstorm_Plan.md] 2026-05-17
**Summary:** Pre-/idea-anchor ideation pass capturing the Phase 2 anchor search space, axis critique (A1-D3), missing anchors, and 2 candidate /idea-anchor sentences after Phase 1 closeout 2026-05-17.
**Key decisions:**
- D1: Axis menu structured by tier (Product / Research / Horizon / Meta), not by topic — makes timeframe + artifact + success-metric trade-offs visible
- D2: Generative-ideation v2.1 payload pattern — output is critique + generation in one round (axis critique + missing anchors + cross-domain analogs + sub-problem decomposition + 2 candidate anchors + user-intent read)
- D3: Stopping discipline preserved — brainstorm produces framing options NOT commitments; no Phase 2 experiment fires; fresh /idea-anchor + /project-init required before action
- Preferred candidate: audit-first AMR evidence-packet system absorbing A2 (audit-as-product) + B2 (4th-mechanism falsifier) + label-quality cartography
- Honest user-intent read: epistemic infrastructure (knowing when evidence is valid), not AUROC chasing

---

## [plan_file: EP1_EP2_Cross_Drug_Synthesis_Plan.md] 2026-05-17
**Summary:** Single-deliverable plan to write the cross-drug architectural-finding synthesis from EP1 (cipro) + EP2 (cef + tet). Closes Phase 1 evidence collection internally. External publication deferred per EP1 closeout discipline.
**Key decisions:**
- D1: One deliverable (`wiki/ep1_ep2_cross_drug_architectural_finding_*.md`); no new code
- D2: Internal scope only (PC1=`internal_closeout` unchanged)
- D3: Defined stopping rule — synthesis closes Phase 1 evidence collection
- D4: Residual uncertainty explicitly listed (cef mechanism question, tet failure-mode disambiguation, Stage 1 architectural verdict, BV-BRC cef feasibility, per-gene NT windows, Bakta annotation, mean+max preflight v3)
- D5: Cef cross-tab from `Cef_Mechanism_Audit_Plan.md` cited as corroborating evidence ONLY; synthesis stands without it

---

## [plan_file: Cef_Mechanism_Audit_Plan.md] 2026-05-17 (revised post-brainstorm round 3 — reduced to 3 steps)
**Summary:** Smoke-runner bug fix + per-strain JSON sidecar + cheap cef NT-vs-mechanism cross-tab diagnostic. Reduced 2026-05-17 from 8 steps to 3 after /brainstorm round-3 framing critique caught scope re-inflation, threshold substitution, Mash-clade-at-N=12 degeneracy, and "rigor theater" evidence-object pattern. The cef audit is now appendix-tier corroborating diagnostic, NOT a decision-bearing framework — the cross-drug architectural finding is captured in EP1+EP2 verdicts; primary deliverable moved to `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md`.
**Key decisions:**
- D1: Cef audit is appendix-tier, NOT decision-bearing (round-3 reframing)
- D2: Smoke-runner bug fix + JSON sidecar is the only load-bearing step (serves future EPs beyond cef)
- D3: Dropped framework discipline disproportionate to N=12 — PC pre-conditions, Mash-clade baseline, evidence-object schema, attribution preflight, 4-cell decision matrix all DROPPED
- D4: External publication still deferred per EP1 closeout
- 3-step structure: Step 1 (bug fix + JSON sidecar) → Step 2 (re-fire cef smoke) → Step 3 (NT-vs-mechanism cross-tab; appendix-tier interpretation only)
- Plan revision history: v1 (rejected 11 steps) → v2 (3 steps) → v3 (re-inflated 8 steps after round-2 detail accumulation) → v4 (this revision, 3 steps after round-3 framing critique)

---

## [plan_file: Cipro_Decision_Bundle_Technical_Plan.md] 2026-05-17
**Summary:** Implementation blueprint for the post-SUSPEND_CONDITION_4 decision bundle, scoped down per /review reductions: collapse Tier 0 Steps 1-2 into one script, drop Bakta + circular variants, defer mean+max preflight v3.
**Key decisions:**
- 12 steps across 5 waves: Wave 0 (3 parallel foundation refactors) → Wave 1 (5 parallel test+consumer) → Wave 2 (2 parallel tests) → Wave 3 (manual census+manifest run) → Wave 4 (conditional label-sensitivity LOSO)
- Critical path: Step 1 → Step 7 → Step 9 → Step 11 → Step 12 (5 waves; 3 code + 2 runtime); max parallelism = 5 agents (Wave 1)
- /review reductions absorbed: collapse census + manifest into one script (Step 7), drop Bakta smoke + mechanism completeness + 3 of 5 manifest variants + --ignore-gate flag, defer mean+max preflight v3 to a conditional follow-up
- Pre-conditions PC1 (D9 framing lock) + PC2 (D4 numeric threshold) are user-locked BEFORE Wave 3 runs (not enforced by code)
- Source plan drift flagged: `plans/Cipro_Decision_Bundle_Plan.md` still contains un-reduced spec; should be edited in-place before /execute-plan per 2026-05-14 HIGH-salience lesson

---

## [plan_file: Cipro_Decision_Bundle_Plan.md] 2026-05-17
**Summary:** Tier-0 cheap decision bundle on N=38 + BV-BRC MIC census, fired BEFORE any Databricks burst or per-gene NT diagnostic. Replaces the rejected binary Path A vs Path B framing after SUSPEND_CONDITION_4 verdict landed 2026-05-17.
**Key decisions:**
- D1: Replace binary A/B with 4-tier decision bundle (Tier 0 cheap → Tier 3 large spend gated)
- D2: BV-BRC census must cover multiple phenotype policies (HIGH_R/HIGH_S + CLSI-strict + EUCAST-strict)
- D3: Frozen label-policy manifest before any relabeling experiment (single source of truth, no hand-coded overrides)
- D4: Primary estimand for relabeling = per-strain error concentration + rank-order stability, NOT max AUROC
- D5: Mean+max attribution preflight v3 is a closeout falsifier (Tier 1), not a fork-decider
- D6: Curated baseline informational run needs 3rd verdict field `given_suspended_gate: INFORMATIONAL_ONLY`
- D7: Bakta 4-strain smoke selection requires a negative-control strain (borderline_S no-mech)
- D8: Mechanism completeness — AMRFinder differential test first, manual blastn only on discordant rows
- D9 (open tradeoff): Phase 1 EP1 deliverable framing ("publish" vs "ship working classifier") affects Tier 2/3 weights — user must lock before Tier 3 fires

---

## [plan_file: Return_Decision_Tree_Patch_Plan.md] 2026-05-16
**Summary:** Apply the /review synthesis's 3 surgical correctness fixes + structural restructure to `wiki/return_decision_tree_2026-05-16.md`, plus the one-line `-u` flag fix to `run_stage1b_detached.bat`. Implementation-ready under 4 steps; max parallelism 2 (doc + bat are independent).
**Key decisions:**
- D1: Eng correctness patches come FIRST, restructure SECOND (sequenced Step 1 → Step 2 because both modify the same doc file)
- D2: `-u` flag (not `PYTHONUNBUFFERED=1` env var) for the .bat fix — explicit + tighter blast radius
- D3: Merge enumerated 8 sub-steps into 4 actual steps; adjacent same-file ops serialize anyway
- D4: Don't touch the running Stage 1b process; Step 3 affects FUTURE relaunches only

---

## [plan_file: Sidework_Sequence_Ship_Path_Plan.md] 2026-05-13 (decisions locked 2026-05-14)
**Summary:** Delta from `Sidework_Sequence_Plan.md` after `/review` (2026-05-13) + post-save `/brainstorm` (2026-05-14). Resolves the load-bearing B-scope problem (B-B locked: drop clade-only from smoke; 12 unique MLST → singleton clades make clade-only degenerate). Fixes deterministic-hashing reproducibility bug. Narrows ARCHITECTURE.md to per-line judgment. Locks --per-class=20. Excludes wiki/GATE_A_REPORT.md from C scope.
**Key decisions (all locked 2026-05-14):**
- B-resolution = B-B (drop clade-only from smoke; 4 variants) (D1)
- Step C scope = 13 files (11 originals + LESSONS_LEARNED + docs/ARCHITECTURE.md with per-line judgment); wiki/GATE_A_REPORT.md NOT in scope (D2)
- Step C is edit-then-stage, not stage-only (D3)
- `mlst_to_clade_id` helper lives in `dna_decode/data/cohort.py` (D4)
- Test scope = 1 parametrized in tests/test_data_cohort.py + 1 integration in tests/test_pipeline_cli.py (both files exist; extend) (D5)
- Drop `fallback_counter: dict`; use deterministic hashing (D6)
- Multi-scheme MLST + scheme-collision: scheme-aware tuple hash via `zlib.crc32` (NOT Python `hash()` — process-salted) (D7)
- Step M demoted out of wave graph (D8)
- Time-box C to 30 min (D9)
- B commit message includes numerical before/after of clade-only AUROC (D10)
- Post-populate slow tests gate the smoke gate (D11)
- Step D `--per-class 20` locked (N=40 with 5R margin) (D12)
- Register `slow` pytest marker (D13)

---

## [plan_file: Sidework_Sequence_Plan.md] 2026-05-13
**Summary:** Ordered work to do while NT v2 100M populate runs in background (~45 min remaining). Brainstorm-revised after Codex critique surfaced under-scoped C, mixed-scope TODOS hunks, and B's deeper-than-15-min reality.
**Key decisions:**
- Sequence A → C → E → B; D optional last (D1)
- C scope = current-state files only; skip archived plans + project_state snapshots (D2)
- TODOS.md hunk staging via `git add -p` across A, C, E (D3)
- B = narrow scope (helper extraction + unit test); defer per_clade_baseline strain-keying (D4)
- pytest discipline during populate = CPU-only target tests via `-m "not slow"` (D5)
- Auto-memory user_environment.md update is separate from C (D6)

---

## [plan_file: Phase2_Decision_Gate_Plan.md] 2026-05-13 (Step 1/2 updated 2026-05-14 for B-B lock)
**Summary:** Split the originally-conflated "12-strain decision gate" into a 12-strain smoke/falsification gate (4 variants, clade-only dropped per B-B lock 2026-05-14) + a tiered N=50 → N=150 staged decision gate (local screen → Databricks burst). Stage 2 acceptance: NT ≥ best classical + 5 pp AUROC AND top-10 NT attribution includes gyrA / parC / parE (biological-plausibility check).
**Key decisions:**
- 12-strain = smoke gate, tiered N=50 → N=150 = staged decision gate (D1)
- Fix clade-only `hash(mlst) % 10` placeholder BEFORE running smoke gate (D2)
- Keep 12-strain smoke at existing 5 variants — no RF/TabPFN/SNP-table additions yet (D3)
- Document the smoke result as smoke, not decision (D4)
- Tiered Option-C threshold: N=50 local screen "any positive lift" → N=150 Databricks "5 pp + biology check" (D5)
- TODOS additions for 4 deferred /research items; SNP-table scope corrected to "parse AMRFinderPlus POINT* rows" per Codex 2026-05-13 (D6)

---

## [plan_file: NT_Deferral_Docs_Cleanup_Ship_Path_Plan.md] 2026-05-13
**Summary:** Scope-tightened delta from `NT_Deferral_Docs_Cleanup_Plan.md` after `/review` synthesis — drops D1's `/save-plan` hedge (already disproven on disk), trims the over-prose `[BLOCKED]` bullet, sharpens the deferral annotation, attempts cheap NT revision retrieval before recording the gap, and adds an explicit untracked-file staging note.
**Key decisions:**
- Drop the `/save-plan` hedge — direct manual edit only (D1)
- Trim `[BLOCKED]` bullet to ≤3 sentences + separate Environment line (D2)
- Attempt NT revision retrieval before recording the gap (D3)
- Sharpen "gate failed" → "equivalence test failed at model load" (D4)
- Step 3 explicit `git add` reminder for untracked plan file (D5)

---

## [plan_file: NT_Deferral_Docs_Cleanup_Plan.md] 2026-05-13
**Summary:** Docs-only follow-up to commit `d4a4652` — apply the three issues surfaced by `/brainstorm` against the just-shipped NT AutoModel refactor deferral: stale plans-index, conflated TODOS scope, missing reproducibility metadata.
**Key decisions:**
- Re-run `/save-plan` before manual index edit (D1)
- Split the TODOS entry into specific (BLOCKED) + general (OPEN) (D2)
- Lean reproducibility metadata, one line (D3)
- Diagnostic spike deferred, NOT killed (D4)

---

## [plan_file: Audit_Calibration_NT_AutoModel_Ship_Path_Plan.md] 2026-05-13
**Summary:** Scope-reduced delta from `Audit_Calibration_NT_AutoModel_Plan.md` after `/review` synthesis — 5 steps / 4 waves → 2 commits, 2 waves. Drops dual-verdict columns (institutionalization risk), drops wiki update (no text to replace), splits NT refactor into a separate gated commit (equivalence test required).
**Key decisions:**
- Asymmetric warning banner replaces dual-verdict columns (D1)
- Drop wiki/GATE_B_REPORT.md update entirely (D2)
- Split NT refactor into a separate, gated commit (D3)
- Keep default-semantics test as the regression lock (D4)
- thresholds_block(rules) helper, not inline string list (D5)
**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.

---

## [plan_file: Audit_Calibration_NT_AutoModel_Plan.md] 2026-05-13
**Summary:** Fix a credibility bug in the just-shipped audit cohort generator (`scripts/audit_cohort.py`) — the "GO" verdict was emitted under silently-relaxed thresholds; defaults produce "WARN" — AND simultaneously replace `NucleotideTransformerModel`'s `AutoModelForMaskedLM` with `AutoModel` to eliminate the `output_hidden_states=True` workaround.
**Key decisions:**
- Audit report header MUST surface threshold values (D1)
- Two verdict columns — Phase 1 production + Gate B infra-only (D2)
- Pin default semantics in tests (D3)
- NT switches to `AutoModel`, not `AutoModelForMaskedLM` (D4)
- Pooling-strategy tag stays "single_seq_mean" (D5)
**Status:** Commit 1 shipped (473b8eb); Commit 2 deferred 2026-05-13 — equivalence test failed at model load (AutoModel.from_pretrained state_dict mismatch on NT v2 100M trust_remote_code checkpoint). See plans/NT_Deferral_Docs_Cleanup_Plan.md and TODOS.md [BLOCKED] NT AutoModel swap.

---

## [plan_file: BVBRC_Genome_Metadata_Adapter_Plan.md] 2026-05-12
**Summary:** Wire `BVBRC_genome.csv` (BV-BRC Genomes-tab export) into the cohort path as a new adapter module, bypassing the wrong-contract `pilot.fetch_ncbi_assembly_quality` scaffold and feeding the existing `--assembly-metadata` wire that `cohort.candidates_from_bvbrc_ast` already accepts.
**Key decisions:**
- Bypass the scaffold instead of implementing it (D1)
- New CLI flag rather than overloading existing `--assembly-metadata` (D2)
- Coverage-log line surfaces ID-namespace mismatches early (D3)
- `fetch_ncbi_assembly_quality` stays scaffolded (D4)

---

## [plan_file: Ecoli_G2P_Phase1_Closeout_Plan.md] 2026-05-12
**Summary:** Wrap up the stalled `/execute-plan` epilogue for `Ecoli_G2P_Phase1_Ship_Path_Plan.md` — toolchain restore, doc reconciliation, first authoritative test pass, archive, state cleanup, push, final report.
**Key decisions:**
- Selective expansion over hold scope (D1)
- Real-data validation = Phase 2 entry criterion, not Phase 1 closeout (D2 — pending)
- Toolchain restore approach: uv vs pip (D3 — pending)
- Archival convention: status-header + git tag recommended (D4 — pending)
- `/documentation` before commit (D5)
- Retrospective re-derivation, not skip (D6)
- Test outcome recorded, not gated (D7)
- Delete both state files at end (D8)
- `.claude/execute-plan-state/` added to `.gitignore` (D9)

---

## [plan_file: Ecoli_G2P_Phase1_Ship_Path_Plan.md] 2026-05-12
**Summary:** Contracted path to ship Phase 1 of `Ecoli_G2P_Platform_Technical_Plan.md`. Captures the `/review` synthesis verdict (HOLD scope + selective contraction within remaining steps) plus the deferred Wave 3.5 hardening fixes from the post-Wave-3 `/brainstorm`. Estimated remaining work: ~700 LOC across 5 implementation steps + 4 hardening edits.
**Key decisions:**
- HOLD scope, do not expand (D1)
- Reorder — Step 15 (smoke + fixtures) BEFORE Step 14 (CLI) (D2)
- Step 14 collapses to one `scripts/pipeline.py` with subcommands (D3)
- Step 13 visualization uses matplotlib + TSV export, NOT pygenometracks (D4)
- Step 17 leaderboard collapses to a shell loop over `pipeline.py train` (D5)
- Step 16 docs trimmed to README + ARCHITECTURE.md only (D6)
- Apply Wave 3.5 hardening BEFORE Step 14 wiring fires (D7)
- Add quantization-fidelity micro-step (selective addition) (D8)

---

## [plan_file: Gene_Presence_AUROC_Bug_Fix_Plan.md] 2026-05-14
**Summary:** Strengthen the diagnostic, confirm the strain-unique-identifier-domination hypothesis on real data, then add a `gene_symbol` column to `AnnotationTable` so the gene-presence smoke variant returns a non-degenerate AUROC at N=12.
**Key decisions:**
- Add `gene_symbol` column to AnnotationTable; do NOT rewrite `gene_id` (preserves embedding cache key + fixes parse_gff3/parse_genbank asymmetry) (D1)
- Strengthen diagnostic before mounting F: drive (absolute counts + per-prefix namespace breakdown + side-by-side dual-extractor AUROC) (D2)
- Strengthen synthetic falsifier with strain-unique-blocks + shared-core LOSO + all-zero held-out row case (D3)
- Add `INDETERMINATE_IDENTIFIER_OOV` smoke verdict as defense-in-depth guardrail (D4)

---

## [plan_file: Stage1_N40_Cipro_Engineering_Screen_Plan.md] 2026-05-14
**Summary:** Run a 4-experiment matrix (NT-XGBoost gate + NT-logreg sanity + k-mer-XGB classical + NT+k-mer-fusion-logreg diagnostic) under LOSO on the N=40 cipro cohort (effective N=38) with paired bootstrap CI, MLST diagnostic appendix, and a 3-bucket verdict to decide whether to spend Stage 2 N=150 Databricks burst budget.
**Key decisions:**
- Restore NT-XGBoost as primary gate-bearing head; add NT-logreg as sanity-check baseline; fusion is diagnostic-only NOT gate-bearing (D1)
- All variants run with `calibrate=False` for primary AUROC (uniform calibration discipline matching smoke-gate; calibration is small-N footgun) (D2)
- Add diagnostic appendix: MLST distribution + per-strain LOSO predictions + paired bootstrap CI (B=1000) + 3-bucket verdict (≥5 pp CLEAN / 3-5 pp NOISY / <3 pp FAIL) (D3)
- Gene-presence + AMRFinderPlus POINT* baselines explicitly out of scope; result packet notes 'best classical' is bounded (D4)

---

## [plan_file: Stage1_Refactor_And_Test_Hardening_Plan.md] 2026-05-14
**Summary:** Convert the /review synthesis into a 3-step refactor: pre-commit decision rules in the Stage 1 plan, reduce `scripts/stage1_n40_cipro.py` to thin orchestration over existing infrastructure, and pin the two critical untested behaviors (fusion-exclusion from gate + `calibrate=False` discipline).
**Key decisions:**
- Treat Stage 2 burst as atomic; pre-commit deterministic per-bucket actions (CI-lower-bound rule converts borderline NOISY PASS → FAIL) (D1)
- Refactor runner to reuse existing infrastructure: `leave_one_strain_out_cv` for NT variants, factored `dna_decode/eval/loso_kmer.py` for k-mer + fusion, `_train_baseline_logreg(calibrate=False)` for logreg path (D2)
- Replace silent mean-fallback on `ClassifierTrainingError` with re-raise; eliminate the failure-masking pattern (D3)
- Add fusion-exclusion + `calibrate=False` discipline regression tests; pin /brainstorm-flagged failure modes (D4)
- Loud MLST handling (raise on None instead of "unknown" fallback) + bootstrap-skip-count reporting (D5, D6)

---

## [plan_file: Stage2_N150_Prep_Plan.md] 2026-05-14
**Summary:** Resolve the three deferred Stage 2 decisions (annotation source, AMRFinderPlus integration, Databricks vs local) and ship the infrastructure needed for a Stage 2 N=150 cipro decision-gate run, so the gate runs cleanly once Stage 1 PASSes.
**Key decisions:**
- Annotation source = Bakta re-annotation for cross-strain stable gene symbols (defer Roary; accept-degenerate rejected) (D1)
- AMRFinderPlus POINT* SNP-table baseline IS in scope for Stage 2 (gyrA/parC/parE textbook signal; load-bearing comparator) (D2)
- Compute = Databricks burst for N=150 NT populate (~3-5 hr A100); local CPU for everything else (Bakta + AMRFinder + analysis) (D3)
- Stage 2 cohort = N=150 expanded from gate_b_cohort.parquet (67 strains) via audit-cohort pipeline with relaxed assembly-quality thresholds (D4)

---

## [plan_file: Stage2_Docker_Tools_Install_Plan.md] 2026-05-14
**Summary:** After user starts Docker Desktop, install Mash + Bakta + AMRFinderPlus via pinned Docker images, write a single tools/docker_runner.py Python orchestration module (NOT .sh wrappers), and smoke-validate on K-12 + one cipro-R strain. Resolves the Phase A.1 / A.2 / A.5 install steps from Stage2_N150_Prep_Plan.md.
**Key decisions:**
- Docker (containers) — NOT WSL2 Ubuntu (D1)
- One Python tools/docker_runner.py module — NOT three .sh wrappers (D2)
- Pin Docker image tags — NOT :latest (D3)
- Correct AMRFinderPlus invocation: amrfinder_update for DB, --database (not --database_path), --mutation_all <file> takes a path (D4)
- Stage Bakta DB on C: drive if room; verify before flagging install complete; record versions+SHA digests (D5)
- Smoke-test on TWO strains: K-12 (binary works) + one cipro-R (POINT-row parsing actually exercised) (D6)

---

## [plan_file: EP_4_Pathotype_Discovery_Closeout_Memos_Plan.md] 2026-05-27
**Status:** executed (2026-05-27; commits 272e90d + bdc4c1e on origin/main)
**Summary:** Two-step sequential execution of the remaining discovery-machine actions for the EP-4 pathotype project: NCBI Pathogen Detection `host_disease` facet audit + EP-1 SUSPEND-gate reuse feasibility memo. User overrode the prior "deferred until Gate A signal" recommendation.
**Outcome:** Step 1 → HONEST-GAP (programmatic facet access blocked; 3 manual paths documented at `research_outputs/ncbi-pathogen-detection-host-disease-facet-2026-05-27.md`). Step 2 → ADAPTED_REUSE (falsifier not triggered; 7/7 noise-class symbols have direct pathotype analogs; workhorse-side reuse recommendation at `research_outputs/ep1-suspend-gate-pathotype-reuse-feasibility-2026-05-27.md`). Step 4 conditional ledger update fired the Step 2 side (Decisions Made +1, Action 5 retired); Step 1 Pending Decision NOT retired (HONEST-GAP path).
**Key decisions:**
- D1: Override the "deferred until Gate A" recommendation — execute now; accept the risk that if Gate A fails, the memos become reference docs rather than direct inputs
- D2: Discovery memos are feasibility-grade (1-page verdict + analogy table), NOT implementation specs — the v0 contract belongs on the workhorse per the handoff doc
- D3: Single commit covering both memos — narrow blast radius
- D4: Conditional project-ledger update only if either memo produces a clean go/no-go verdict that retires a Pending Decision or updates a Hypothesis status

---
