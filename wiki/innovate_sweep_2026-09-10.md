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

### 2. `COUPLING-not-popdesign` — the discriminator may be determinant↔clade coupling, not population design
HGT is nature's own randomization: a horizontally mobile determinant is decoupled from clade the way
meiosis decouples a locus in a cross. **Measured on `wiki/provdisjoint_lineage_metrics.json`: acquired-gene
cells' mean sensitivity drop under lineage collapse is −0.057 (they *gain*) while chromosomal ciprofloxacin
loses +0.051.** Population design and inheritance-coupling are perfectly correlated across all five
failures and three successes, so the current attribution in `dna_decode/eval/regime.py` is under-determined.
*Kill-test:* acquired collapsing ≥ chromosomal would disprove. **Split holds — survived.**
*Limit:* n=10 cells, one drug in the chromosomal arm. Suggestive, not established.

### 3. `DOUBT-bacterial-gap` — a false clean bill on bacterial target-site cells
`target_site_doubt("ciprofloxacin", {"gyrA": {"D86N"}})` returns **only `position_novelty`**, `any_doubt:
false`. The L2 completeness screen is never constructed when `doubt_cell_for(drug)` is None, so every
bacterial target-site cell is a fourth, undeclared state the module's own three-state discipline does not
admit. The early-return rationale is sound for position-*novelty* and a non-sequitur for *completeness*:
a position-based catalog can be incomplete at the position level — the V179F shape.
*Kill-test:* a completeness signal being produced would disprove. **None produced — survived.**

### 4. `ARBANK-source-conc` — our own standard, never pointed at our own highest numbers
**0 of 26 AR Bank artifacts carry any source-diversity field.** The 0.60 largest-source bar refused to
report 6 of 10 frozen SCORED cells; the AR Bank arm — single curated deposits by construction, reporting
several 1.00 accuracies — was never screened by it. `[[feedback_apply_your_own_standard_to_your_own_cohort]]`
with the roles reversed.
*Kill-test:* any such field present would disprove. **Zero — survived.**

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
- **No discrimination controls (H1).** They are supported only for `file-exists` / `project-state-row`;
  every kill-test here is `test-exit-0`, for which controls are deferred in v0 and would be marked
  controlled-attempt-invalid. Stated rather than omitted silently.
- **Deeper adversaries not yet run.** Per protocol the expensive pass belongs on survivors only:

  ```
  /brainstorm wiki/innovate_sweep_2026-09-10.md
  /idea-validation-council wiki/innovate_sweep_2026-09-10.md
  ```
  Until then #1 and #2 are `unfalsified-until-council`, **not** ratified findings.
- **G1's reframe was not shared with G9/G8/G5** — they were launched in parallel for cost, so each ran on
  the original framing only. A deviation from the ideal sequencing, recorded rather than hidden.
- Nothing here is auto-executed. The top survivor goes to `/technical-plan`, not to a build.
