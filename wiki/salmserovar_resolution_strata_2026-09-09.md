# The serovar cell's 0.8400, stratified by how each call was resolved (2026-09-09)

**The headline survives, and the measurement found two things the headline could not say.**

`typing:Salmonella:salmserovar` publishes **0.8400** (168 hit / 11 miss / 21 no-call) on 200 wet-lab
slide-agglutination-labelled isolates. That number pooled three interpretability classes, and the
proportions were unmeasured — a rate whose interpretable fraction is unknown is this project's own
*"a rate is a claim about its denominator"* error in a different costume. This measures the split.

Reproduce: `uv run python scripts/salmserovar_resolution_strata.py`
(offline, seconds, no blastn — it reads two committed checkpoints and the shipped table).

## Result

| route | hit | miss | no-call | what it means |
|---|---|---|---|---|
| **exact** | **165** | 6 | — | the `(O,H1,H2)` key names one serovar under *both* ambiguity policies |
| **arbitrary_winner** | **3** | 5 | — | the key is contested; the table records an arbitrary winner |
| **fallback** | **0** | 0 | — | the `O+H1` phase-incomplete route produced nothing (see below) |
| unresolved | — | — | 21 | the formula did not resolve at all |

**165 of 168 hits (98.2%) are uniquely resolved.** Accuracy counting *only* those: **0.8250** vs the
published 0.8400 — a difference of **3 isolates**.

**So the headline holds.** The brainstorm's concern that an unknown fraction of 0.8400 was
"right by luck" is answered: the fraction is 3/168 = **1.8%**, and it is now recorded rather than
assumed small. Quote **0.8400** with the stratum disclosed, or **0.8250** for the
tie-break-independent figure — both are defensible; what was not defensible was not knowing.

## Method — and why it needs no rebuild

The contested set was **measured, not reconstructed**. The two ambiguity policies were each scored on
the same 200 isolates, and the policy affects only the name lookup. The script **REFUSES (exit 2) unless
all 200 antigenic formulas are byte-identical across the two checkpoints** — that gate is what makes a
serovar diff attributable to the policy alone, and it is pinned by a test that feeds it a deliberate
mismatch. All 200 matched. An isolate named under `first` and absent under `omit` is exactly one whose
call depends on the arbitrary winner: there are **8**, with **0** disagreements in the other direction.

Docker being down did not block this — rebuilding the candidate table from `Initial_Conditions` was the
obvious route and the unnecessary one.

## Finding 1 — the arbitrary-winner policy trades in both directions, and only one was reported

The name-table memo records `first` as **strictly dominant**: +3 hits over `omit`, zero regressions.
That is true *on the hit count* and it is not the whole trade. All 8 of these isolates would **abstain**
under `omit`, so the policy also converts **5 abstentions into confidently-wrong serovar calls**. Its
hit rate within the stratum is **3/8 = 0.375**, far below the cell's own accuracy.

The two framings disagree about the same 8 isolates:

- **hit-count reading** (the shipped one): net **+3**, adopt.
- **asymmetric reading**: net **+5 wrong calls**. This cell's *own* threshold work committed to exactly
  this asymmetry, in these words — *"an abstention is recoverable by a human, a confident wrong serovar
  is not"* — and pre-registered a bar that was deliberately harsher on new errors than on abstentions.

**By the standard this cell already adopted, the shipped policy is the losing side of the trade.** Both
numbers now ship in the artifact. **Which reading governs is an acceptance-bar question and is
deliberately NOT resolved here** — it is the user's call, and it is the third time in one session that
an acceptance bar has turned out to be the load-bearing question rather than the code.

The 8 isolates are concrete: `13:z29:-` → `Agoueve` **three times** (truth Johannesburg, Cubana, Cubana)
— the known genuine within-subspecies-I collision, and one more than the 2 previously recorded.

## Finding 2 — the fallback contributes zero, and that is a measurement, not an artifact

**All 179 calls came through an exact key.** The `O+H1` phase-incomplete fallback produced **no calls at
all**. This was checked rather than inferred, because a stratum of exactly zero is also the shape a
broken classifier produces: the fallback **was reached** on the 6 isolates carrying a complete formula
whose exact key missed, and **failed to collapse to a single name on each**. Not dead code — exercised
and unsuccessful.

**Consequence for the open work:** the brainstorm's C3 (ambiguity metadata must cover the fallback
route, not just exact keys) has **measured headroom of zero on this cohort**. That is a reason not to
build it *yet* — not evidence it can never matter.

## Finding 3 — the largest no-call cluster is one formula, and it is monophasic Typhimurium

**4 of the 21 no-calls (19%) are the single formula `4:i:-`.** The caller resolves the antigens
correctly and abstains only because **8 serovars share O=4 / H1=i**, so `O+H1` cannot collapse. Their
wet-lab labels are `4,[5],12:i:-`, `I 4,5,12:i:-`, and `Salmonella Typhimurium - monophasic` — i.e. the
labels are themselves **antigenic formulas**, and the caller already emits the matching formula.

This is the **largest single addressable no-call cluster** and the clearest next lever. It is **named,
not taken**: it needs its own pre-registered bar, and two bar overrides in one session is already not a
pattern to normalise. The remaining 2 complete-formula no-calls (`22-gene3:z:1,6` → Poona,
`23-gene2:z:l,w` → Worthington) are the previously-recorded malformed-DB-name class.

## Honest limits

- This measures the **composition** of an already-published number. It cannot and does not change it.
- The `arbitrary_winner` stratum is a **LOWER BOUND** on contested-key involvement: the diff captures
  every isolate whose *call* depends on the winner, but not a contested key whose winner happens to
  equal what would have been returned anyway.
- The `fallback` stratum cannot separate arbitrary winners inherited through the `O+H1` candidate set —
  moot at zero calls here, live if that ever changes.
- Per-serovar cap of 12 flattens prevalence: every figure is **per-isolate accuracy on a diverse mix**,
  not population-weighted.
- One cohort, one DB build. The zero-fallback and `4:i:-` findings are properties of **this cohort**.

Artifact: `wiki/salmserovar_resolution_strata_2026-09-09.json` (carries the full per-isolate table).
Tests: `tests/test_salmserovar_resolution_strata.py` (10).
