# Can recent theses/papers build the missing part? — INTERIM (2026-09-28)

**STATUS: INTERIM, AND THE SEARCH IS ~15% DONE.** The campaign was planned for four families (F1–F4,
user-ratified full breadth) and cut off by a **session rate limit** before any family agent returned a
final report. F3 and F4 were **never searched**. Do not read the small candidate count as "the literature
is empty" — read it as "the literature was barely looked at".

Plan: `plans/Deep_Research_Missing_Part_2026-09-28_Plan.md` · Harvest:
`wiki/deep_research_candidates_2026-09-28.json` · Screen: `wiki/deep_research_screen_2026-09-28.json`

## The planning finding, which outranks every search result here

Derived from the repo's own `wiki/negative_results_map_2026-06-13.md`:

> **"decoding is bounded by LABELS, not models."**

A campaign that hunts **architectures** in recent theses therefore searches the axis this project has
already measured as *not binding*. That error was made here days ago — one run named "score gLM2-650M vs
the curated baseline" the highest-VOI move when the regime screen puts it in
`natural × molecular × zero_shot = LOSES_TO_CATALOG`, a **recorded negative**. So this campaign ranks
datasets above methods, and every method candidate is regime-screened *before* being called promising.

## What the screen says (3 candidates)

| class | n | meaning |
|---|---|---|
| **INCONCLUSIVE** | 2 | the screen RAN, condemned nothing, and awaits **one named measurement** |
| **LEAD_ONLY** | 1 | recorded as a negative for our axis; deliberately not label-screened |
| REJECTED | 0 | nothing was condemned |

### F2 (the one genuinely OPEN regime cell) — two candidates, both INCONCLUSIVE, both `regime=OPEN`

**1. "Public K-12" (PRECISE-1K + public SRA) — `https://github.com/SBRG/precise1k`.** The only candidate
**verified by direct fetch and quotation**, not from a search summary. Repo's own wording: PRECISE-1K
**1,035** samples + **1,675** high-quality public samples = a **2,710-sample** resource.

Why it bears on M2: the measured bottleneck is the *intersection* of Keio conditions with PRECISE-1K
expression — **11 of 28** carbon conditions, **621 of 1,035** samples on glucose, **6 of 11** usable
conditions resting on ≤5 samples — and PRECISE-1K's condition space is only **4 temperatures / 9 media /
39 supplements / 76 knockouts**. The 1,675 added samples come from *many* labs, so they plausibly span a
wider condition space, which is precisely the binding side.

**Three things verification found that the search summary did not say:**
- **The decisive number is absent.** The repo states **no** distinct-condition/media/carbon count for the
  combined set, so whether it widens the Keio intersection is **UNVERIFIED**. This single number decides
  the candidate, and it is what the G6 gate is waiting on.
- **Access wall**, verbatim: *"`log_tpm` files for the public dataset are not provided as they are too
  large for GitHub—please contact us for these."* An author contact is an **external** wall…
- **…but there is a code-closable path:** the samples are public SRA data and the repo ships the raw
  metadata (`Escherichia_coli_20220127.tsv`), so the compendium could in principle be **re-derived from
  SRA**, converting the external wall into a code wall. Cost unestimated, not attempted.
- **NOT NEW:** SRA metadata pulled **2022-01-27**, well outside the stated window. It does not answer
  "what is new"; it may still answer the user's actual question, "is there anything we can use."

**2. Fitness Browser Feb-2024 release** (`figshare.com/articles/dataset/…/25236931`) — **search-reported,
UNVERIFIED**: 7,552 fitness experiments across 46 bacteria + 2 archaea (up from 6,570 / 42 in Nov 2020).
No 2025/2026 release located. **Widening the fitness side alone does not fix M2** — the bottleneck is the
intersection with condition-matched *expression*, and for a non-E.-coli organism that pairing is
unestablished. **Known data defect landing on our exact panel:** the authors report bad sucrose and
D-mannitol stocks, such that *E. coli* BW25113 appeared to grow on sucrose in the original media but not
with fresh stock — any carbon-panel work must account for it.

### F1 — one candidate, and it is a NEGATIVE with useful shape

**MaveDB / Atlas of Variant Effects growth** (search-reported, UNVERIFIED): 7M+ variant-effect
measurements over 700+ genes. But **no new β-lactamase / AMR-target DMS surfaced for 2025–26** — the only
β-lactamase entry found is the long-known Stiffler 2015 TEM-1 study. The Alliance's growth is oriented to
**human clinical** variants (a Clinical Atlas), not the pathogen targets our forward/inverse cell needs.
And its headline infrastructure win is coordinate **mapping** (98.61% of ~2.5M variants), which **does not
touch M3**: knowing *where* a variant sits says nothing about putting two assays' magnitudes on one scale.

## A defect found in my own tooling by inspecting its output

The first screen run reported both F2 candidates as **REJECTED**. They were not. The gate screen had
returned `INCOMPLETE` ("no measurement supplied for G6") and my classifier folded that into `REJECTED` —
condemning two candidates that had passed **9 of 10 gates** and carried **`regime=OPEN`**, the one
genuinely open cell in the project. "Needs one more number" is not "rejected", and the conflation errs in
the **discarding** direction. Fixed with a distinct `INCONCLUSIVE` class that names the awaited gate;
pinned by two tests, one of them a non-vacuity control asserting a genuinely tripped gate is *still*
rejected.

## Process findings (route around these next time)

- **Do not fan out 4 agents × sub-agents.** That 12-way parallel burn is what hit the session limit. Run
  families sequentially or two at a time, and forbid sub-agent spawning.
- Agent-reported, **UNVERIFIED**: CORE.ac.uk's search page is a JS shell and its API 429s without a key;
  Edinburgh ERA reportedly has a 2024–25 indexing gap.
- `research-leads` YouTube search works (real PhD defences returned) but carries **no `upload_date`**, so
  the 2-year window needs a per-video lookup.
- **Screen coverage limit, found by using it:** `rejection_gates` supports exactly two intended layers
  (`L1_AMR_RS`, `L4_forward_continuous`). A condition-resolved expression or fitness corpus fits **neither**
  cleanly; both F2 candidates are filed L4 as the closer fit. That is a limitation of the screen for this
  campaign, not a property of the candidates.

## THE DECISIVE MEASUREMENT — taken, and it is a NEGATIVE

The "single highest-value next action" identified above was *get the distinct-condition count*. **It was
then taken, in this same run, and it kills the candidate for M2.**

The access wall turned out to be **only on `log_tpm`** (the expression matrix) — **the metadata was fully
fetchable**. `gh` located `data/k12_modulome/metadata_qc_with_p1k.csv` (1.18 MB, 2,738 curated samples with
`Carbon Source`, `Nitrogen Source`, `Base Media`, `Temperature`, `pH`, `Supplement`), and it was downloaded
and computed over directly. **No author contact was needed for the count.**

| | PRECISE-1K (published) | **Public K-12 (measured here)** |
|---|---|---|
| glucose share | 621/1,035 = **60%** | **1,313/1,587 annotated = 82.7%** |
| distinct carbon sources | 9 media | **18**, but only **10** with ≥5 samples and **3** with ≥10 |
| carbon unannotated | — | **1,151 / 2,738 = 42% blank** |

**Verdict: `FAILS_THE_M2_REQUIREMENT`.** The 1,675 added public samples are mostly *more glucose* and
*un-annotated rows*. Glucose dominance gets **worse** (60% → 82.7% among annotated), and the long tail is
thin — only three carbon sources reach ten samples. The expansion adds SAMPLES, not the CONDITION BREADTH
that M2's bottleneck requires.

Honest caveat: PRECISE-1K's published bottleneck was an *intersection with Keio's 28 carbon conditions*,
while the figure computed here is *distinct carbon sources present in the compendium*. Those are related
but not identical, so the comparison is directional. The 42%-unannotated and 82.7%-glucose figures stand
on their own regardless.

Incidental but real: sucrose appears on 11 samples here, and the sibling Fitness Browser candidate carries
a **known bad-sucrose-stock defect** on *E. coli* BW25113 — any sucrose-based cross-use inherits it.

## A SCHEMA GAP this exposed, and the temptation not to take

**The 10 rejection gates have no gate for CONDITION COVERAGE** — which is precisely M2's blocker. G6
measures *phenotype-value* degeneracy (mode share / distinct values), not how many growth conditions a
corpus spans. So the standing screen **structurally cannot adjudicate** these candidates' real problem,
and both correctly remain `INCONCLUSIVE` on the gate screen while failing the campaign criterion recorded
beside it.

**The tempting fix — feed the condition count into G6 — is refused.** That would silently redefine a
committed gate to mean something it does not mean, and every prior screened candidate's G6 verdict would
become incomparable. Whether the gate family should GAIN a condition-coverage gate is a schema decision,
i.e. an authority call, not a patch to make today's candidate resolve.

## Where that leaves the four families

- **F2 — answered, negatively, on the strongest lead.** Public K-12 fails the condition-breadth
  requirement. Fitness Browser Feb-2024 remains unverified and, structurally, widens only the *fitness*
  side while M2 binds on the *intersection with condition-matched expression*.
- **F1 — one shallow search.** The DMS axis returns a clean negative (no new AMR-target DMS; the Alliance's
  growth is human-clinical and its mapping work does not touch M3). The AMR measured-AST axis,
  AllTheBacteria, the EBI AMR Portal and Zenodo/Figshare thesis deposits were **not searched**.
- **F3, F4 — not searched at all.**

## Next actions, in order

1. **F3 and F4 from scratch** — they are wholly unsearched and F3 sits on the project's critical path.
2. **F1's non-DMS axes** — measured-AST collections, AllTheBacteria, EBI AMR Portal, thesis data deposits.
3. **Re-verify the Fitness Browser counts** against the Figshare record (currently search-reported only),
   and, if pursued, check whether any of its 48 organisms has a condition-matched *expression* partner —
   that pairing, not the fitness breadth, is the thing M2 needs.
4. **Relaunch sequentially, not 12-way parallel** — that fan-out is what hit the session limit.
