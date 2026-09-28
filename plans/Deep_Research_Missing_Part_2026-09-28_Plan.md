# Deep-research campaign — can recent theses/papers build the missing part? (PLAN, 2026-09-28)

**Status:** PLAN ONLY. Nothing searched yet. Awaiting scope ratification.
**Mode requested:** `--until-mvp minister … decompose & --plan`
**Preconditions verified (not assumed):** money hook hardened + wired (real-surface: armed → money exit 2,
destructive exit 2, benign exit 0; dormant unarmed), `Agent` present, cwd = project root, self-init 25/25
free, no live peer session. **Minister preflight returns `(True, [], [])` — the autonomous lane is open.**

---

## 0. The finding that reorders the whole search

Before choosing sources I derived the target from the live registries rather than from memory.
`wiki/negative_results_map_2026-06-13.md` line 7 states the project's own measured conclusion:

> **"decoding is bounded by LABELS, not models."**

So a campaign that hunts **architectures** in recent theses searches the axis this project has already
measured as *not binding*. That is the R2 framing error, and it has bitten here before: one run named
"score gLM2-650M vs the curated baseline" the highest-VOI move, and `screen_proposal` put it in
`natural × molecular × zero_shot = LOSES_TO_CATALOG` — a **recorded negative**, retired this session.

**Consequence for this plan:** datasets/labels rank above methods, and every methods candidate must clear
the regime screen *before* it is called promising.

## 1. What "the missing part" actually is (derived from `eval/regime.py` + the negative-results map)

| # | Missing capability | Current measured state | Can new literature move it? |
|---|---|---|---|
| **M1** | Organism-level trait decoding from **natural** populations | zero-shot **CLOSED_NEGATIVE** (0-for-5, de-confounded); **supervised = REQUIRES_DECONFOUNDING** (open, conditions attached) | Yes — but look for *population DESIGNS*, not models |
| **M2** | **Condition-switch** prediction (same gene, different environment) | the **one genuinely OPEN** regime cell; within-gene AUROC 0.71–0.81 but the model emits one identical ratio for 61–76% of genes. Bottleneck **measured**: 11 of 28 Keio conditions ∩ PRECISE-1K, 6 of 11 rest on ≤5 samples | **Yes — highest-leverage.** A richer condition-resolved expression/fitness corpus directly unblocks it |
| **M3** | **Dosage** on the inverse (effect→edit) | it RANKS, cannot DOSE; calibrators provably cannot transfer (assays share no scale: CcdB's whole range lies below TEM-1's minimum) | Maybe — calibration transfer / conformal under shift is an active methods area |
| **M4** | **The label supply itself** | the repeatedly-measured binding constraint; G1–G10 kill most candidates | **Yes — this is what literature most reliably adds** |
| **M5** | Catalog **completeness** (not accuracy) | the doubt layer's measured failure mode (`rmt` gap; HIV V179F) | Partly — already partly covered by the 2026-09-03 AMRrules scan |

## 2. Ground already covered — do NOT re-run

| Artifact | Covers | Implication |
|---|---|---|
| `wiki/prior_art_decoder_landscape_2026-09-03.md` | AMRrules/ESGEM-AMR (our rule layer, being built by an ESCMID working group); Hu et al. 2024 *Brief Bioinform* bbae206 (rule-based beats ML under phylogeny-aware splits) | Start from these, don't rediscover them |
| `wiki/prior_art_genomic_language_models_2026-08-31.md` | the gLM landscape | gLM-as-zero-shot-scorer is a closed direction |
| `wiki/prior_art_conservative_resistance_blosum_2026-08-31.md` | conservative-resistance / BLOSUM prior art (Friedman 2013 + EGFR/ALK counterexample class) | The mechanism is published; scope claims accordingly |

## 3. Where I will look (concrete)

**A. Theses — the user's explicit ask.** Aggregators: **OATD.org**, **NDLTD Global ETD**, **DART-Europe**,
**CORE.ac.uk**, ProQuest (mostly paywalled — noted, not relied on), EThOS (verify post-2023-incident status).
National/institutional: DiVA (SE), Theseus (FI), theses.fr/HAL (FR), DNB (DE), TU Delft / Wageningen / Utrecht
(NL), Oxford ORA / Cambridge Apollo / Edinburgh (UK), **DTU Findit (DK)**, ETH Research Collection (CH).
**Group-targeted, because these groups own the substrates we already consume:** DTU **CGE** (ResFinder /
PointFinder / MLST — our finder cells wrap their DBs), Oxford **Modernising Medical Microbiology** (CRyPTIC +
the Oxford cohort), **Sanger** (AllTheBacteria / 661k), **EMBL-EBI** (AMR Portal), **Stanford** (HIVDB /
CoV-RDB), **Institut Pasteur** (BIGSdb / `wzi` — our ktype), Fowler & Rubin (MaveDB / **Atlas of Variant
Effects**), Marks lab (EVE), Rost lab (VespaG).

**B. Papers / preprints.** bioRxiv + medRxiv, arXiv `q-bio.GN` + `cs.LG`, **Europe PMC** (full text +
preprints, good API), **OpenAlex / Semantic Scholar** for **citation-graph expansion from our own anchors**
(Rhee 2003 PhenoSense · Bloom 2013 segregants · Orth 2011 Table S1 · Napier barcode · WHO catalogue v2 ·
ProteinGym · Hu 2024 bbae206 · Friedman 2013). Venues: *Bioinformatics*, *NAR*, *Genome Biology*,
*Microbial Genomics*, *Nature Methods*, *eLife*, *AAC/JAC*, *Lancet Microbe*, ISMB/RECOMB/MLCB.

**C. Dataset registries — ranked FIRST per §0.** MaveDB · ProteinGym · **Atlas of Variant Effects** ·
CRyPTIC · BV-BRC · NCBI Pathogen Detection · EBI AMR Portal · CARD Resistomes · **AllTheBacteria (2024)** ·
661k · **Fitness Browser / RB-TnSeq** · **PRECISE-1K / iModulon** (the direct M2 lever) · BacDive ·
EnteroBase · Pathogenwatch · FDA-ARGOS · NARMS · EuSCAPE · and **Zenodo / Figshare / Dryad**, because a
thesis usually deposits its data there even when the PDF is embargoed.

## 4. Skills, and what each one buys

| Skill | Why it improves *this* search |
|---|---|
| **`/research-department`** (FIRST) | Retrieves what's already tried/foreclosed and accrues negatives. With 3 prior-art scans + a negative-results map + a regime map on disk, this is the correct entry point and is what stops re-covering ground |
| **`/research`** | The orchestrator: web research → 13-column audit table → auto-runs intake + follow-up queue. Used **per question**, not once for everything |
| **`Agent` fan-out** | One sub-agent per source family (theses / preprints / datasets / methods) — independent, parallel, no shared state. The accelerator |
| **`/research-leads`** | **Thesis defences and conference talks are on YouTube**, and this host has a *verified working* transcript fetcher (`fetch_transcript.py` + `search_yt.py`). A defence talk states the contribution plainly and is reachable when the PDF is embargoed — a real edge for the thesis half |
| **`/research-verify`** | Re-fetches every cited URL and checks quote-verbatim + DOI cross-check. Non-negotiable here: writing a biological constant from memory is this repo's named fabrication hazard |
| **`/research-intake`** | Splits supported vs unsupported claims; enforces per-row locators |
| **`/graphify`** | Once the corpus passes ~50 items, build the knowledge graph to surface cross-paper patterns |
| **`/sme-panel`** | Multi-domain expert lenses + devil's-advocate pass on the shortlist |
| **`/probe`** (by hand) | Repo-grounded feasibility for each survivor |

**The spine of the plan — two screens, applied before anything is called promising:**
- `scripts/screen_candidate_gates.py` → **G1–G10** (does a usable LABEL exist / is the rule scoreable?)
- `dna_decode/eval/regime.py::screen_proposal` → **regime verdict** (is this a recorded negative?)

**A literature find is a LEAD, not a finding, until it clears both.** Without this the campaign reproduces
the `/innovate` sweep pattern, where most survivors fell as soon as they were acted on.

## 5. Decomposition — 4 families (+1 optional)

| Family | Question | Gap | Why separable |
|---|---|---|---|
| **F1 · label-substrate-scan** | What *new* public phenotype/label substrates (2024-09 → 2026-09) clear G1–G10? | M4 | Highest VOI per §0; pure dataset work, no method claims |
| **F2 · condition-coverage-scan** | Is there a condition-resolved expression/fitness corpus deeper than PRECISE-1K ∩ Keio? | M2 | Targets the ONE open regime cell; a single dataset could unblock it |
| **F3 · deconfounded-design-scan** | Do recent theses/papers demonstrate a *population design* that de-confounds organism-level g→p? | M1 | Designs, not models — the discriminating variable is population design |
| **F4 · dosage-transfer-scan** | Any calibration-transfer / conformal-under-shift method that could give the inverse a magnitude? | M3 | Methods-only; must clear the regime screen |
| *F5 (optional)* · **thesis-harvest-infra** | A reusable harvester + `/graphify` corpus over ETD aggregators | — | Only if F1–F4 show the manual search is the bottleneck |

**Flow-down:** F1 and F2 are independent and run first (both are dataset questions). F3 depends on F1 (a
de-confounded design is worthless without a substrate to run it on). F4 is independent but lowest-ranked.
**Critical path: F1 → F3.**

## 6. Draft MVP bar (checkable predicates only — for ratification)

1. `file-exists wiki/deep_research_missing_part_2026-09-28.md` — the synthesis memo
2. `file-exists wiki/deep_research_candidates_2026-09-28.json` — one row per candidate with **both** screen
   verdicts (G1–G10 + regime), recorded even when the verdict is REJECT
3. `test-exit-0 uv run python scripts/screen_candidate_gates.py --verify` — every committed screen still reproduces
4. `test-exit-0 uv run pytest tests/test_deep_research_candidates.py -q` — the harvest is machine-checked
   (schema + every cited URL/DOI syntactically resolvable + no candidate admitted without a screen verdict)
5. `project-state-row` — an Action-Log row per family with its verdict

**Deliberately NOT in the bar:** "find something usable." That is not checkable and not in my gift — the
honest deliverable is *a screened, verified, reproducible answer*, including "nothing clears the gates,"
which on this project's record is the likelier outcome and is still worth having.

## 7. Budget + honest expectations

Estimated **35–60 steps** of the ~100-step `--plan` self-budget. Gates unchanged: paywalled theses
(ProQuest) and any paid API are a **money pause** — I will use only free/open access and say plainly where
a paywall blocks a lead rather than paying or guessing.

**Base rate, stated up front:** of the last ~6 candidate substrates screened here, **one** cleared
(PEAR/G6) and it still wasn't buildable as specified. Expect mostly REJECTs with named gates. The value is
in the screened record, not in a hoped-for win.
