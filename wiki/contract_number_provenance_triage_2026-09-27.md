# Contract-number provenance triage — all 34 unverifiable numbers, by kind (2026-09-27)

Research step (Step 3) of `plans/Contract_Number_Provenance_Sound_Decoy_Control_Typed_Provenance_Kinds_Plan/`.
No production code. Every number the standing audit reported as UNVERIFIABLE at commit `bcfcbc2` appears
below exactly once, with its kind, the quoted prose that is the evidence for that kind, and — for the
`artifact` kind — the specific artifact plus the **full-pool discrimination rate** measured by the Step 1
control.

## The rule, fixed before looking

> **Cite only when the artifact resolves AND every one of the cell's numbers is present AND full-pool
> discrimination is better than LOW (rate < 0.50).** A citation that cannot discriminate is recorded as
> such rather than added, because adding it would restate the `4d18399` over-claim this plan exists to
> correct. A citation that would surface drift is NOT added — it is reported as a finding.

Applied as written. It admitted one cell the plan expected it to reject (`essentiality`) and rejected two
the plan expected it to admit (`pneumoserotype`, and the single-file `hla` starting point). Those
disagreements are recorded rather than resolved by moving the rule.

## Headline finding — the audit's artifact-side matcher accepts a substring inside a longer number

Found while checking whether the candidate citations were genuine. The prose side is carefully
boundary-guarded (`DECIMAL_RE` uses `(?<![\w.])` and `(?![\d.])`), but the artifact side is a raw
substring test — `any(v in haystack for v in _number_variants(tok))`. So a short token matches inside a
longer number: `2.1` "matches" inside `72.1`.

**Measured exposure: 4 of the 157 currently-clean numbers rest on nothing else.**

| cell | number | prose says | what it actually matched | real status |
|---|---|---|---|---|
| `typing:bacteriophage:phage` | `0.862` | `25/29 called = 0.862` | `0.86` inside `"similarity": 0.867` | **derived** (25/29 = 0.8621) |
| `typing:bacteriophage:phage` | `0.291` | `25/86=0.291` | `0.3` inside `0.36363636…` | **derived** (25/86 = 0.2907) |
| `typing:Salmonella:salmserovar` | `0.900` | `Coverage 0.705->0.900` | `0.9` inside `0.9217877…` | **derived** (180/200) |
| `finder:any:forward` | `0.5115` | `CcdB 0.5115` | `0.51` inside `0.5189` / `0.517` | **unresolved** |

Three are legitimately derived, which this script's own posture explicitly allows; their *evidence* was
nonetheless coincidental. The fourth is genuinely unresolved and is the one that needs adjudication.

**Not fixed in this run, deliberately.** A boundary-aware matcher would flip the standing verdict from
`NOTHING_TO_ADJUDICATE` to `ADJUDICATION_REQUIRED` with 4 candidates. That is a change to the project's
standing trust surface and a user-authority call, the same class as editing a frozen file — and this
plan's own discipline for a drift-surfacing discovery is to *report it*, not to fix it silently. It is the
recommended next step. Until then **the clean `0/157` is four numbers weaker than it reads**, and Step 5
records that caveat rather than the bare figure.

## Decisions

| cell | numbers | decision | artifact(s) | discrimination |
|---|---|---|---|---|
| `typing:klebsiella:kleb` | 3 | **ADMIT** | `wiki/klebsiella_topk_ksweep_2026-07-25.{json,md}` | **8/599 = 1.3% HIGH** |
| `essentiality:any:essentiality` | 7 | **ADMIT** | `essentiality_report_card.json` + `essentiality_e3_learned_2026-07-28.json` + `essentiality_e3_human_2026-07-28.json` | 39/199 = 19.6% MODERATE |
| `hla:human:b5701` | 5 | **ADMIT** | `wiki/hla_validation_2026-07-06.md` | 264/600 = 44.0% MODERATE *(weak — see note)* |
| `typing:Streptococcus_pneumoniae:pneumoserotype` | 5 | **DO NOT ADD** | — | would surface drift |
| `pgx:human:cyp2d6` | 2 | **DO NOT ADD** | — | 251/299 = 84.0% LOW |
| `typing:human:pigment` | 2 | **DO NOT ADD** | — | 254/299 = 85.0% LOW |

Unverifiable residual after Step 4: **24 → 9**.

### ADMIT — `typing:klebsiella:kleb` (the one strong citation)

All 6 numbers present with exact, boundary-checked matches. Discrimination **1.3%** — 591 of 599
unrelated haystacks fail to contain this number set, so the citation genuinely could have failed.

The plan warned this had to be **narrowed from the 33-file bundle** used in the planning probe, whose
6/6 result was confounded by haystack size. Confirmed and narrowed: the source is the single
`klebsiella_topk_ksweep_2026-07-25` pair, and its HIGH grade is measured against equally-sized decoys.

| number | prose | evidence |
|---|---|---|
| `0.45` | `top-1 ~0.45 / top-5 ~0.60 over 147-165 KL-types` | exact match in artifact |
| `0.60` | same clause | variant match in artifact |
| `0.49` | `lift +0.49 over a 0.10 prior null` | exact match in artifact |

(The cell's other three numbers were already typed in Step 2: `0.90` and `0.10` are `threshold`,
`10.5281` is a `structural-non-measurement` Zenodo DOI prefix.)

### ADMIT — `essentiality:any:essentiality` (against the plan's expectation)

The plan expected this to land as **not worth citing**, because its seven numbers match many unrelated
artifacts. The rule as fixed admits it: 19.6% is MODERATE, which is better than LOW. Recording the
disagreement rather than moving the rule — a pre-registered rule that is overridden once it produces an
inconvenient answer was never a rule.

It needed **three** artifacts, and that is provenance rather than citation-shopping: the report card holds
the genome-wide figures, `e3_learned` holds the E. coli E3 number and `e3_human` holds the human one,
because the prose spans three separate results.

| number | prose | evidence |
|---|---|---|
| `0.695` | `offline E. coli AUROC 0.695 genome-wide` | rounds to artifact `0.6952` |
| `9.3` | `base rate 9.3%` | rounds to artifact `0.0928` (percent/fraction pairing) |
| `0.373` | `sens 0.373 / spec 0.984` | rounds to artifact `0.3732` |
| `0.984` | same clause | exact match |
| `0.37` | the same sensitivity, quoted at 2dp | rounds to artifact `0.3732` |
| `0.795` | `The learned E3 complement lifts it (E. coli 0.795 / human 0.911)` | **only in `essentiality_e3_learned`** |
| `0.911` | same clause | **only in `essentiality_e3_human`** |

**The two-file version of this citation was rejected.** `essentiality_report_card` alone is missing
`0.795` and `0.911` under boundary-checked matching — they had appeared to match only through the
substring defect above. Had that not been checked, this step would have shipped a citation that a correct
matcher immediately reports as drift.

### ADMIT — `hla:human:b5701`, but the plan's starting point was wrong

The plan's located candidate was `wiki/hla_b5701_validation_2026-07-06.json`, flagged "unverified". It is
**not** the right artifact: the cell's prose spans three alleles (B\*57:01 sens 0.979, B\*58:01 sens 0.61,
A\*31:01 sens 0.0) and there are three separate per-allele artifacts. `0.61` is absent from the b5701
file, so citing it alone would surface drift.

Two alternatives were measured. The three per-allele JSONs together contain all five numbers but grade
**LOW (161/199 = 81%)** — a three-file haystack matches almost anything, which is the size confound Step 1
neutralised, here showing up as a reason not to cite a bundle. The single roll-up
`wiki/hla_validation_2026-07-06.md` contains all five with four exact matches, and grades MODERATE.

| number | prose | evidence |
|---|---|---|
| `0.979` | `sens 0.979 / spec 0.992 / PPV 0.855` | exact |
| `0.992` | same clause | exact |
| `0.855` | same clause | exact |
| `0.61` | `B*58:01 rs9263726 sens 0.61 weak` | exact — artifact reads `(sens 0.61, PPV 0.18), unsafe for an SJS/TEN screen` |
| `0.0` | `A*31:01 rs1061235 not-paneled sens 0.0` | variant |

**Note the weakness, and do not over-read it.** 44.0% means 264 of 600 unrelated artifacts also contain
all five numbers; it clears the fixed bar by 6 points. Two of the five are `0.0` and `0.61`, which are
common values. A second asymmetry: the citation is a `.md`, which is not in the json decoy pool, so its
haystack is one markdown file compared against single-JSON decoys. This citation records provenance; it
does not strongly verify it.

### DO NOT ADD — `typing:Streptococcus_pneumoniae:pneumoserotype` (would surface drift)

The plan required narrowing from "the 7-file glob used in the probe" to one artifact. Narrowed to
`wiki/pneumo_serotype_cohort_validation.json`, which holds 4 of the 5 numbers (serogroup 0.939, exact
0.661, and the Quellung-subset 0.952 / 0.690) — but **not `2.1`**.

`~2.1%` is the no-call rate, a **derived** figure the artifact stores as counts rather than a rate
(5 no-call against 230 scored = 2.17%, quoted as "~2.1%"). Citing this artifact would put a derived
number into `candidate_drift` and manufacture an adjudication item out of a correct claim.

**A two-file set appeared to fix it and was rejected as a false match.** Adding
`pneumo_selection_rule_probe_2026-09-04.json` made the drift disappear — but its only literal `2.1` is
inside **`72.1`**, a coverage metric. That is the substring defect in the headline finding, and adding the
file would have manufactured provenance. This is the cleanest demonstration in the run of why "no drift"
is not by itself evidence that a citation is sound.

The honest fix is a `derived` kind (a number correctly computed from counts the artifact does hold), not a
citation. Not built here — it is a new kind, and Step 2's taxonomy was closed before this was known.

### DO NOT ADD — `pgx:human:cyp2d6` and `typing:human:pigment` (vacuous by measurement)

Both fail the rule on discrimination, and both fail it badly.

| cell | numbers | artifacts holding all of them | discrimination |
|---|---|---|---|
| `pgx:human:cyp2d6` | `0.62`, `1.0` | **730** of ~1,200 wiki files | 251/299 = **84.0% LOW** |
| `typing:human:pigment` | `0.468`, `1.0` | **737** of ~1,200 wiki files | 254/299 = **85.0% LOW** |

Two numbers, one of them `1.0`, is not a fingerprint. `wiki/cyp2d6_structural_2026-07-06.json` and
`wiki/pigment_1000g_population_2026-07-29.json` are the semantically correct artifacts and both contain
their cell's numbers — citing them would be *true* and would still verify nothing, because five of every
six unrelated artifacts would satisfy the same check. Recorded as a measured property of the numbers, not
a defect of either cell.

## The 10 numbers Step 2 typed out of the residual

Already resolved in code; listed here so all 34 are accounted for.

| cell | number | kind | evidence from the prose |
|---|---|---|---|
| `typing:klebsiella:kleb` | `10.5281` | structural-non-measurement | `Zenodo 10.5281/zenodo.14065540` — a DOI prefix |
| `typing:cat:catcolor` | `5.1` | structural-non-measurement | `ARHGAP36 5.1-kb intron-1 deletion` |
| `typing:horse:horsecolor` | `4.6` | structural-non-measurement | `STX17 4.6-kb dup grey` |
| `typing:Escherichia_coli:pathotype` | `0.833` | enforced-by-test | `EXACT-equality asserts` in `tests/test_pathotype_expec_recall.py` |
| `typing:Escherichia_coli:pathotype` | `1.0` | enforced-by-test | same asserts (precision and EPEC recall) |
| `typing:Escherichia_coli:pathotype` | `0.917` | superseded-value | `the earlier flat-K=1 support rule scored 0.917 (11/12) but OVER-RESCUED` |
| `typing:Escherichia_coli:pathotype` | `0.80` | threshold | `each >=0.80` — the per-gene coverage floor |
| `typing:klebsiella:kleb` | `0.90` | threshold | `greedy-rep @0.90` — the clonality-correction cut |
| `typing:klebsiella:kleb` | `0.10` | threshold | `a 0.10 prior null` |
| `pgx:human:nudt15` | `9.5` | external-reference | `*3 EAS AF ~9.5%` — a published population frequency |

The DOI pattern also caught `doi:10.5061/dryad...` in **both dog cells**, which this triage had not
predicted; those two were previously counted as artifact matches by coincidence. A per-number list would
have shipped that miss.

## Reconciliation against the script

The check that matters: if the memo's totals disagree with the audit's, one of them is wrong and both get
re-read.

| quantity | memo | `wiki/contract_number_audit_2026-09-27.json` |
|---|---|---|
| numbers unverifiable at `bcfcbc2` | 34 | 34 |
| typed out of the residual by Step 2 | 10 | 10 (structural 3 + enforced 2 + superseded 1 + threshold 3 + external 1) |
| unverifiable after Step 2 | 24 | `n_numbers_unverifiable` = 24 |
| admitted for citation in Step 4 | 15 (3 + 7 + 5) | — |
| projected residual after Step 4 | 9 | to be confirmed by re-running |
| total extracted, all cells | — | `n_numbers_extracted` = 191 |

Step 2's structural pattern types 5 numbers in total; 3 came out of these 34 and 2 are in cited cells
(the dog DOIs), which is why the structural row above lists 3 and the script reports 5.

## Limits

- **Discrimination is measured against this repo's `wiki/*.json` pool as it stands today** (599–600
  files). A rate is a statement about that corpus, not a universal property of the numbers.
- **`MODERATE` is a weak pass.** Two of the three admitted citations are MODERATE, and `hla` clears the
  bar by 6 points. The bar `rate < 0.50` is asserted, not derived — as is the 5% HIGH boundary.
- **A citation records provenance; it does not make a number correct.** The audit checks that a cited
  number is present in the artifact, which is a much weaker claim than that the measurement was sound.
- **The substring defect is unfixed**, so every "present" verdict in the current artifact — including the
  ones this memo relies on — carries that looseness. The three admitted citations were verified with a
  boundary-aware matcher specifically so that this step does not add a citation resting on the defect it
  reports.
