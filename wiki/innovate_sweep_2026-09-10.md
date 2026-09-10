# /innovate sweep — 2026-09-10

Four forced operators (G1 reframe · G9 abduce · G8 recombine · G5 relax-a-constraint), 7 deduped
candidates, every kill-test **executed by the falsification engine**, not by hand.

**Engine verdict (verbatim):**

```
python skills/soraya/scripts/falsification.py --run-ledger <ledger> --cwd .
OK (executed+ranked): 7 candidate(s) ->
  6 survived ['AR-unwired','DOUBT-bacterial-gap','ARBANK-source-conc',
              'COUPLING-not-popdesign','ERRORS-isolate-clustered','POOLING-dilution']
  | 1 killed ['AR-onesided-novel'] | 0 unfalsified []
exit 0 (validator clean)
```

Kill-tests are pytest nodes in `tests/test_innovate_killtests_2026_09_10.py` with **inverted polarity**:
each PASSES exactly when disproving evidence is found (engine reads that as `killed`). Every failure below
was verified to be a substantive `AssertionError`, not a crash — a test that errors also "fails" and would
otherwise manufacture a false survivor.

`survived` means **passed one executed falsification attempt**. It does not mean true.

## Killed

**`AR-onesided-novel`** — *claim:* reading the single-class AR Bank cohorts one-sided is an unexploited
move. **Disproved.** `wiki/ar_isolate_bank_external_validation_2026-07-18.md` already states *"the curated
bank is resistance-enriched → per drug, ONE class is powered … Read one-sided, it's a coherent, powered,
independent test."* The work ran in July 2026 and produced 26 committed artifacts.

**This candidate existed because the data inventory shipped hours earlier produced a FALSE LEAD.** It
flags `ar_bank_kleb_extval_ceftriaxone` as referenced by nothing — correctly, since no file names that
*directory* — while the work was done and reported under a different naming convention
(`wiki/external_validation_ar_bank_kleb_ceftriaxone_*.json`). `orphan` means "no file names this
directory", never "unused"; the artifact's own `known_limits` says so, and this is the first live instance.

## Survivors, ranked

### 1. `POOLING-dilution` — mean-pooling, not population design, may be the mechanism
Every one of the 0-for-5 used a **mean-pooled whole-genome** representation; every success resolves
individual loci (segregant-cross markers, per-gene GPRs, per-variant scoring). Pooling weights each causal
gene ~1/N_genes. **Measured: PC1 explains 0.807 of pooled-embedding variance** on the committed cache — the
leading direction is ancestry by construction, for any model. Explains why scale never helped: more
parameters do not change 1/N weighting.
*Kill-test:* PC1 fraction < 0.5 would disprove. **0.807 — survived.**
*Limit:* correlational and on one cache; the decisive separation (a pooled representation on a
**constructed** cross) was built at `scripts/yeast_bloom_fm_prep.py` and never ran — blocked on coordinate
liftover.

> **DEPLOYABLE HALF EXECUTED AND REFUTED, 2026-09-10** (`wiki/pooling_vs_locus_resolution_2026-09-10.md`).
> The survivor's actionable prediction — *restricting to candidate causal loci before pooling recovers
> signal pooling destroys* — was run against a bar frozen beforehand and returned the pre-registered
> **`REFUTED_LOCUS_IDENTITY`** branch. QRDR-restricted (0.8386) does beat whole-genome pooling (0.8029),
> but **a random draw of 4 arbitrary single-copy genes reached 0.9020** — higher than the four genes that
> actually cause the resistance. Narrowing helps; **locus IDENTITY does not**, which is what
> ancestry-correlated resistance looks like. The PC1 = 0.807 observation stands; the conclusion drawn
> from it does not. **The random-gene control was the load-bearing arm and it is what killed the claim.**

### 2. `COUPLING-not-popdesign` — the discriminator may be determinant↔clade coupling, not population design
HGT is nature's own randomization: a horizontally mobile determinant is decoupled from clade the way
meiosis decouples a locus in a cross. **Measured on `wiki/provdisjoint_lineage_metrics.json`: acquired-gene
cells' mean sensitivity drop under lineage collapse is −0.057 (they *gain*) while chromosomal ciprofloxacin
loses +0.051.** Population design and inheritance-coupling are perfectly correlated across all five
failures and three successes, so the current attribution in `dna_decode/eval/regime.py` is under-determined.
*Kill-test:* acquired collapsing ≥ chromosomal would disprove. **Split holds — survived.**
*Limit:* n=10 cells, one drug in the chromosomal arm. Suggestive, not established.

### 3. `DOUBT-bacterial-gap` — **RETRACTED AS WRITTEN 2026-09-10, and the correction found two real defects**

**What was claimed:** "a false clean bill on bacterial target-site cells."
**Why that was an over-claim:** the kill-test called `target_site_doubt("ciprofloxacin", …)` **directly**,
and no shipped surface makes that call. `--observed` is fungal/viral-only (`ERROR: --observed is
fungal-only; bacterial drugs use --amrfinder-run`), and the bacterial record is built separately, carrying
`validation: trust_block(...)` — which **does** include a scored `doubt_layer`
(`arm: determinant_completeness`, verified live for ciprofloxacin and gentamicin). Bacterial cells are
covered; they are covered by a *different* function. **A function reachable only by direct call is not a
shipped disclosure gap** — the same "does it reach a call?" discipline recorded on 2026-09-02.

**What the correction found instead — both live on the CLI, both now fixed:**

1. **A genuinely MUTANT-LEVEL cell was told its catalog was position-based.** `lenacapavir` (CAI) ships
   the CAPELLA emergent-substitution set *precisely because* capsid polymorphisms (K70R/A105T) made a
   position-based rule call 140/140-R — yet it was absent from `_MUTANT_LEVEL_CELLS`, so the CLI printed
   `doubt: n/a -- this catalog is position-based` **two lines above its own `MUTANT-LEVEL v0` caveat**.
   A self-contradicting human-facing disclosure, the same class as the `signals[0]` bug of 2026-09-02.
2. **NOT-MEASURED rendered as SILENCE on four cells.** `doubt_one_line` knew `applicable is False` and
   `assessed is False` but not `measured is False`, so every UNMEASURED cell fell through to the
   honest-silence return — `sarscov2-mpro`, `fungal-fluconazole-erg11`, `fungal-voriconazole-erg11`
   **since 2026-09-02**, plus `lenacapavir`. The record said "has NOT been measured"; the human output
   said nothing. That is exactly the failure `doubt_one_line`'s own docstring exists to prevent.

Fixed in `dna_decode/eval/{doubt,position_novelty}.py` + `dna_decode/data/target_site_completeness.py`;
15 guards in `tests/test_doubt_not_measured_rendering.py`, proven non-vacuous (4 fail when the fix is
reverted). Registering `hiv-cai` in the doubt map alone **crashed** a real lenacapavir call with a bare
`KeyError` — a cell must be registered in the position-novelty catalog too, now pinned by a parity test.

**Net: the survivor's headline was wrong and its investigation was worth more than the headline.**

### 4. `ARBANK-source-conc` — our own standard, never pointed at our own highest numbers
**0 of 26 AR Bank artifacts carry any source-diversity field.** The 0.60 largest-source bar refused to
report 6 of 10 frozen SCORED cells; the AR Bank arm — single curated deposits by construction, reporting
several 1.00 accuracies — was never screened by it. `[[feedback_apply_your_own_standard_to_your_own_cohort]]`
with the roles reversed.
*Kill-test:* any such field present would disprove. **Zero — survived.**
**MEASURED 2026-09-10 (`scripts/ar_bank_source_concentration.py`) — and the survivor's predicted
magnitude is REFUTED.** The hypothesis was that the profiler "returns largest_share > 0.60 for every
cell". It does not: **7 of 13**. The six gonorrhoeae cohorts are **single-panel (share 1.000)** — the arm
reporting several 1.00 accuracies — while **Klebsiella (5 panels, 0.333–0.408) and E. coli (6 panels,
0.263–0.286) PASS the same bar**; S. aureus levofloxacin marginally fails at 0.615. The gap is real and
**CONFINED**, not blanket. Grouped by CDC `panel_id`, a PROXY: the cohorts carry no BioProject field and
the 0.60 bar was calibrated on BioProjects, so read a failure as "draws on one curated panel", never as a
like-for-like BioProject number.

### 5. `ERRORS-isolate-clustered` — error mass may be isolate quality, not catalog gaps
`SAMN04014846` is mis-called under **both** kleb-ciprofloxacin and kleb-gentamicin — mechanistically
unrelated rules. It is also the single false-R behind the published Klebsiella gentamicin spec 0.889.
Two independent catalog gaps is not the parsimonious explanation.
*Kill-test:* no multi-cohort error isolate would disprove. **One exists — survived.**
*Limit:* n=1 isolate. This is a lead, not a finding.

### 6. `AR-unwired` — computed evidence that reaches no consumer
**0 AR Bank cells on `wiki/decoder_validation_report_card.json`.** The one-sided numbers exist
(kleb gentamicin spec 0.889, kleb ceftriaxone sens 1.00) stamped `hard_fail` + `run_degraded`, because the
powering gate requires both classes. Real gap, lowest novelty — wiring, not insight.

## Honest caveats on this sweep

- **6-of-7 survival is a high rate and partly an artifact of test shape.** Four kill-tests
  (`AR-unwired`, `DOUBT-bacterial-gap`, `ARBANK-source-conc`, `ERRORS-isolate-clustered`) test whether a
  gap *currently exists* — easy to survive, and surviving establishes the gap, **not** that closing it is
  worth anything. Only #1 and #2 make substantive mechanism claims where the number carries weight.
- **One survivor's headline was WRONG and a `survived` verdict did not catch it (2026-09-10).**
  `DOUBT-bacterial-gap` passed an executed, valid, gated kill-test and was still an over-claim, because
  the kill-test exercised a function the shipped CLI never calls for that drug. **A kill-test inherits the
  reachability of whatever it invokes** — calling a library function directly proves a property of the
  function, not of the product. The engine cannot see that distinction; only reading the call path can.
  Retracted in place above, with the two real defects the retraction uncovered.
- **No discrimination controls (H1).** They are supported only for `file-exists` / `project-state-row`;
  every kill-test here is `test-exit-0`, for which controls are deferred in v0 and would be marked
  controlled-attempt-invalid. Stated rather than omitted silently.
- **Survivor #1's deployable half is now REFUTED BY EXECUTION** (see the inset above), which is the
  second of six survivors to fall once acted on. Both fell the same way: the kill-test that let them
  survive was weaker than the claim they were carrying — #3's exercised a function no shipped surface
  calls, #1's measured a correlation without the control that could contradict it. **A `survived` verdict
  bounds the test that ran, never the claim.**
- **Deeper adversaries not yet run.** Per protocol the expensive pass belongs on survivors only:

  ```
  /brainstorm wiki/innovate_sweep_2026-09-10.md
  /idea-validation-council wiki/innovate_sweep_2026-09-10.md
  ```
  Until then #1 and #2 are `unfalsified-until-council`, **not** ratified findings.
- **G1's reframe was not shared with G9/G8/G5** — they were launched in parallel for cost, so each ran on
  the original framing only. A deviation from the ideal sequencing, recorded rather than hidden.
- Nothing here is auto-executed. The top survivor goes to `/technical-plan`, not to a build.
