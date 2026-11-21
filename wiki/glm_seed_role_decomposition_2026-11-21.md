# The instrument CAN resolve the generator question — in one feature mode, and the seed lever I planned is dead

**Status:** `DIAGNOSTIC`. Artifacts: `wiki/glm_seed_role_ceilings_2026-11-21.json` (per-seed noise floor) and
`wiki/glm_seed_role_variance_2026-11-21.json` (variance decomposition). Both read-only, CPU-only, no torch,
no GPU, no network, no D: volume. Reproduce with
`uv run python scripts/glm_seed_role_ceilings.py` and `uv run python scripts/glm_seed_role_variance.py`.
Plan: `plans/Can_This_Instrument_Resolve_The_Generator_Question_At_All_Plan/`, Steps 1–3 of 7 executed.

## DECISIONS NEEDED

**One. `kmer` mode has no winnable region. Retire it, or keep it as a reported-but-not-scored diagnostic?**
Across the committed grid's own ten seeds the `kmer` Markov null averages **0.5302** against a noise floor of
**0.5422** — margin **−0.0121**, negative at 6 of 10 seeds. The target a generator is asked to beat sits
inside the measurement noise, so no `kmer` result can register an improvement. This changes what every
future generator artifact reports, which is why it is yours and not mine. My recommendation is to keep
reporting it and stop *scoring* it, because a mode that reads near-zero is informative about the instrument
and dangerous as a headline.

**Two. Steps 4 and 6 of the saved plan have lost their rationale — strike them?** Both exist to decouple the
single `--seed` into named roles and re-run the grid on the decoupled harness, on the premise that the shared
seed inflates the error bar. That premise is now measured **false, with the opposite sign** (see CORRECTIONS).
Step 6 also needs a gated Kaggle session, so striking it saves a real cost. I have already done the one piece
that survives on its own merits: the artifact now records its seed. I have **not** edited the plan — that is a
scope call.

**Three. The generator ceiling criterion is still open, and it now has a denominator.** In `both` mode the
zero-shot generator sits **+0.2150** above the noise floor, of which **+0.1185** merely reaches parity with a
3-mer Markov chain. That is a concrete target an arm set could be scored against. Whether the criterion
should optimise *distinguishability from natural sequence* (cheap, oracle-free, now with a measured winnable
region) or *oracle-scored design value* (what the north star names, but it makes the predictor half
load-bearing out of distribution) remains the authority decision it was. Nothing has been written to
`soraya_runs/arms-ledger.json`.

## CORRECTIONS

**My "0.0543-wide kmer band" was wrong twice over, and the error flattered the arm.** I derived it as
`0.5990 − 0.5447`. Both inputs are unrepresentative: **0.5990 is seed 7's kmer null, the MAXIMUM of the ten**
(range 0.4623–0.5990), and 0.5447 is a single-seed ceiling. Measured per seed, the kmer margin is **−0.0121**,
not +0.0543 — the sign is opposite, so the mode I described as having a narrow winnable region has none.

**"The model-vs-null gap is 0.1186" was the ZEROSHOT figure quoted without its qualifier.** Recomputed:
zeroshot **+0.1185**, epoch1 +0.1201, epoch2 +0.1388, epoch3 +0.1643, all-model mean **+0.1355**. That series
looks like monotone degradation under fine-tuning and **I am not claiming it** — it is the same delta that
returned `DELTA_NOT_ESTABLISHED` against a paired seed spread of 0.2072.

**The plan's own central hypothesis is measured false, and this is the load-bearing correction.** I planned
the seed-role split as a free precision lever, reasoning that one `--seed` driving corpus shuffle, generation
RNG and discriminator split together mixes nuisance variance into every reported spread. Measured on the
Markov null across 3 arms × 10 replicates × 2 modes, **the coupling SUPPRESSES the spread rather than
inflating it.** In `both` mode: all-coupled SD **0.0140**, generation-only SD **0.0280**, corpus-only SD
**0.0249**. Pinning two of the three nuisance sources *doubles* the spread. So decoupling would not narrow
the error bar — it would widen it.

**Seed 7 is unrepresentative in both modes at once, in the two directions that each mattered.** Its `kmer`
null 0.5990 is the maximum of ten; its `both` ceiling 0.5865 sits well above the 0.5419 mean. Every prior
reading of this arm that leaned on seed 7 was leaning on an outlier twice.

## The result

**The two feature modes are not two views of one question.**

| mode | null (mean) | noise floor (mean) | margin | seeds positive | verdict |
|---|---|---|---|---|---|
| `kmer` | 0.5302 | 0.5422 | **−0.0121** | 4/10 | `NO_WINNABLE_REGION_ON_AVERAGE` |
| `both` | 0.6383 | 0.5419 | **+0.0964** | 10/10 | `WINNABLE_REGION_AT_EVERY_SEED` |

**This number had never existed, for a structural reason.** The committed grid's cell schema is
`{model, epoch, seed, both, kmer, null{kmer,both}}` with **no ceiling field**, so the per-seed margin could
not be computed from it. The one committed ceiling pair came from a different artifact at a single seed and
was read as though it applied at all ten. The probe refuses to write unless it reproduces that pair
(0.5447 / 0.5865) at seed 7, and the test proves the gate is non-vacuous by perturbing the *measurement*
rather than the constant.

**So the answer to the plan's title question is yes — in `both` mode.** There is a 0.0964-wide region,
unanimous across seeds, and the generator is 0.2150 away from the floor. That is a real, measurable target.
In `kmer` mode the answer is no, and no amount of seeds fixes it.

**The variance result has a measured mechanism, not a shrug.** `Var(B) + Var(C)` is **7.18×** `Var(A)`;
independent contributors would make those roughly equal, so the roles are negatively correlated and cancel
inside a single seed. The concrete route: the Markov null is **fitted on the fit half of the very corpus draw
it is then scored against**, so a draw that happens to be easier to discriminate also hands the null a
matched fit set. Coupling preserves that matching; pinning the corpus breaks it.

**The obvious rival is refuted, not argued away.** Arm B pins the corpus at one seed, so its inflated spread
could be a quirk of that corpus. Re-run at anchors 0, 3, 7 and 9 the generation-only SD is **0.0220–0.0288**,
above the coupled 0.0140 at *every* anchor.

**This corroborates the existing paired-delta discipline from a new direction.** The coupling is a feature of
a matched within-seed comparison, which is what `scripts/glm_paired_sweep_verdict.py`'s docstring already
said: *"within one seed both models see the same windows, the same split and the same null, so the confound
cancels inside each delta."* I had read that as a defence of the paired delta and planned to remove the
coupling anyway. The coupling is what makes the delta tight.

## Honest limits

- The variance decomposition is measured on the **Markov null**, so it bounds the *instrument's* nuisance
  variance. The HF generator's own sampling variance is a different quantity and may be larger or smaller;
  nothing here measures it.
- The cancellation mechanism is supported by the additivity ratio and the anchor sweep but is **not isolated**
  by an experiment that removes the fit-set matching directly.
- Arm C varies corpus **and** split together, so it bounds their combined contribution without separating
  them from each other.
- The nulls are **read** from the committed grid rather than re-run, so this inherits that grid's generation.
  The ceilings are freshly computed.
- The ceiling is the **MAX** over 5 internal self-control seeds (the falsifier's own rule), so it is
  conservative by construction — a mean-based floor would be lower and every margin more positive.
- 10 replicates per arm, so each SD is itself imprecise; a ratio of two SDs at n=10 is a coarse reading.
- The verdict thresholds (0.5 / 0.8 on the SD ratio) are **stated, not derived**. Raw ratios ship beside them.
- One geometry (150 bp promoter windows, limit 1200), one organism, one null. The margin is a property of a
  (substrate, geometry, feature-mode) triple, not of discriminators in general.
- **A winnable region is necessary, not sufficient.** A k-mer-plus-positional discriminator sees LOCAL
  composition, so clearing it would not establish that generated sequence is functionally realistic.
- Steps 1–3 of 7 ran. Step 5 (freeze a bar) is unexecuted; Steps 4 and 6 are the ones this memo recommends
  striking, and that recommendation is **not** self-executed.
