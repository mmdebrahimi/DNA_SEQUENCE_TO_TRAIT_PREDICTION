# The completeness screen rediscovers `rmt` from first principles — L2, step 1

First build under the CALL / DOUBT / EVIDENCE plan (`plans/Hybrid_Decoder_Architecture_Plan.md`), F-A step
1. Today's framing sweep verified this object **did not exist** (`F4-same-class-as-rmt` survived its
kill-test).

## What it is

The `rmt` and HIV blind spots are one shape: **a determinant family present in the data but
unrepresentable by the rule**. `scripts/determinant_completeness_screen.py` detects that shape generically
from cached AMRFinder output — no model, no network, no Docker, frozen surface untouched.

**It never emits a call.** An uncounted determinant is a *candidate* for human review. Many exclusions are
deliberate and correct (`blaTEM-1` is not ceftriaxone-R; `aph`/`aadA` are not gentamicin-R), which is why
the ranking matters more than the detection.

**It asks the deployed rule, it does not reimplement it.** For each distinct determinant the screen writes
a table carrying the original header and the row *verbatim*, then asks `call_resistance`. Re-deriving the
rule's logic here would drift from `DRUG_RULE` the moment either changed.

## The headline: the known gap, re-found blind

| signature | symbol | subclass | genomes | R | S |
|---|---|---|---:|---:|---:|
| **`rmt_like`** | **rmtE1** | AMINOGLYCOSIDE | 36 | **36** | **0** |
| `mixed` | aph(6)-Id | STREPTOMYCIN | 97 | 62 | 28 |
| `mixed` | aph(3'')-Ib | STREPTOMYCIN | 95 | 62 | 29 |
| `mixed` | aph(3')-Ia | KANAMYCIN | 60 | 43 | 10 |
| `mixed` | aadA5 | STREPTOMYCIN | 71 | 38 | 31 |

A screen that knows nothing about gentamicin specifically ranks the **known** blind family first — 36 R
carriers, **zero** S — while the deliberate exclusions sort below it as `mixed`. That is the validation:
the general shape reproduces the specific finding.

Across the other deployed drugs the screen flags **no `rmt_like` family with meaningful support**:
meropenem and tetracycline return only `mixed` (blaTEM-1 9R/31S, oqxA 19R/17S — correct exclusions),
ciprofloxacin returns `qnrA1` (4R/0S) and `oqxA10` (3R/0S), both small-n and both already documented
deliberate exclusions.

## Two defects, both caught by reading the output rather than trusting it

**1. The ranking buried the answer.** Sorting by raw R-count put `rmtE1` (36R/0S) **fifth**, beneath
`aph`/`aadA` at 62R/28S. Purity is what separates a gap from a deliberate exclusion, so the signature now
leads and volume breaks ties.

**2. A one-row probe can never satisfy a multi-hit rule.** Ciprofloxacin's rule requires **two** QRDR
hits, so the single-determinant probe reported *0 of 51 determinants counted* and flagged every QRDR point
mutation as a gap — `parC_S80I` at 60R/0S on top — when the rule represents them perfectly. That
conflates *"the rule cannot represent this"* with *"the rule needs more than one"*, and only the first is a
completeness gap. The probe now repeats the row to the rule's threshold, and the QRDR false positives
disappear.

Both produced plausible-looking output. Neither would have been caught by a green suite.

*(A third, smaller one: the first smoke run capped at 250 genomes hit only **unlabelled** accessions, so
the ranking — which is driven by R/S carriers — degenerated to raw prevalence. The scan is now
labelled-first.)*

## Honest limits

- **Ranking depends on labels**, and only 200 of 1,818 cached genomes carry NCBI-PD calls. Families with
  no labelled carrier are reported `unlabelled` and sort last — unassessable, not innocent.
- **`rmt_like` is a heuristic** (≥3 R carriers, 0 S), not a test. It is a triage signature for review.
- **Oxacillin is uninformative here** — no oxacillin labels in the census, and that cell is
  `LABEL_CONFOUNDED` anyway.
- This run scanned 220 genomes per drug. A full-index run is cheap and unrun.

## What it does not do

No resistance prediction. No change to any existing call. No touch to the frozen surface. The next step
(F-A step 3, wiring a `doubt` block into the record) is separate and not taken.

Reproduce: `uv run python scripts/determinant_completeness_screen.py [--drug X] [--limit N]`
