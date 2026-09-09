# Porting SeqSero2's O procedure: +0.120 wet-lab accuracy, zero regressions — and the frozen bar still says REJECT

The serovar cell trailed the pinned reference tool by **−0.175**, and the entire gap was the O axis.
The cause was established as a **lost header qualifier**: our antigen DB is SeqSero2's, re-headered into
`O__<antigen>__<index>`, a schema with no field for the semantic qualifier SeqSero2 attaches to six of
its O entries. The caller therefore read a **differential marker** as a standalone positive allele and
flattened the **O9-vs-O9,46 discriminator**, so plain O9 was unemittable.

This run restores the qualifier and ports the decision procedure that consumes it. **Every number below
is measured on the same 200 wet-lab-labelled isolates, against a bar frozen before any of the code
existed** (`wiki/salmserovar_o_fix_acceptance_bar.json`, `REGISTERED_BEFORE_IMPLEMENTATION`).

## The result

| | before | after | SeqSero2 1.3.2 |
|---|---|---|---|
| wet-lab accuracy (n=200) | 0.7050 | **0.8250** | 0.8800 |
| hit / miss / no_call | 141 / 39 / 20 | **165 / 14 / 21** | 176 / 24 / 0 |
| accuracy on called | 0.783 | **0.922** | 0.880 |
| delta vs the reference tool | −0.1750 | **−0.0550** | — |
| O agreement (fixed 200 denom.) | 159 agree / 33 differ | **183 / 2** | — |
| H1 agreement | 200 / 0 | 200 / 0 | — |
| H2 agreement | 146 agree / 54 unresolved | 146 / 54 | — |

**The H axes are byte-identical, by construction rather than by luck.** Their headers were left exactly
as they were, so the H side of the BLAST database is the same file content it was before; the change
cannot reach them. Confirmed: H1 200/0/0 and H2 146/0/54, unchanged in every cell.

## What the bar said, and why it is a false positive

    [PASS] O agreement rises >= +15 over 159            (159 -> 183, +24)
    [PASS] H1 agreement / unresolved not worse          (200 -> 200, 0 -> 0)
    [PASS] H2 agreement / unresolved not worse          (146 -> 146, 54 -> 54)
    [PASS] wet-lab accuracy not below 0.7050            (0.7050 -> 0.8250)
    [FAIL] overall no_call not above 20                 (20 -> 21)
    VERDICT: REJECT

Six of seven clauses pass, several by wide margins. The seventh fails **by one isolate**, and the
per-isolate transition table shows the mechanism that clause exists to catch is **absent**:

    hit  -> hit       141
    miss -> hit        20
    no_call -> hit      4
    miss -> no_call     5
    miss -> miss       14
    no_call -> no_call 16
    hit  -> anything else:  0        <-- zero regressions

The bar's own `explicitly_not_a_success_signal` names the forbidden pattern: *"a rise in O agreement
achieved by converting wrong calls into abstentions."* That is not what happened. **Twenty misses became
hits.** The net `no_call` rose by one because **five wrong calls became honest abstentions while four
old abstentions were rescued into hits** (+5 − 4 = +1).

And the five new abstentions are the right answer, not a dodge: **SeqSero2 itself returns no O antigen
on four of the five** (`I -:z:1,5`, `I -:f,g:-`, `I -:d:1,7` ×2). All five fire
`C_guard_1319_marker_only` — the guard that stops the 130 bp differential marker from asserting
`1,3,19` on its own. Previously each of those five was a **confident wrong serovar**.

**A net count cannot distinguish laundering from rescue.** That is the defect in the clause: I wrote a
*net* ceiling to police a *directional* failure mode. The transitions are now emitted by the script for
exactly this reason, so the question is answered by measurement rather than by argument.

## The decision: ADOPTED, as a recorded OVERRIDE

The bar is **not edited** and the artifact's `verdict` field still reads `REJECT` — that is what
pre-registration means, and rewriting either would make the exercise theatre. The code ships anyway,
and this is an **override**, logged as one:

- it is strictly dominant on the measured evidence (**zero** hit→miss, **zero** hit→no_call);
- the failing clause's stated mechanism is measurably absent;
- the margin is one isolate against a +24-hit, −25-miss, +0.120-accuracy change.

**This repo has hit this exact class before and recorded it**: the 2026-09-04 coverage-threshold sweep
found that following its own pre-registered "most conservative clearing cut" rule *would have regressed
the cell by 19 correct calls*, and the lesson banked was that **an anti-overfitting rule can itself be a
source of error**. Same shape here — with the important difference that there the rule's *premise* was
wrong, and here the rule's *operationalisation* (a net count for a directional property) is.

If this override is judged wrong, the reversal is one commit and the evidence to judge it is all above.

## Two defects the first cut had, both caught on real data

The port did **not** work first time. Run 1 scored **0.7300** and failed the bar on two clauses; the
diagnosis produced both fixes below. Its checkpoint is kept at
`D:/dna_decode_cache/ss2_checkpoint/ofix_rows.run1_coverage_argmax_and_tyr.jsonl`.

**1 — the tyr O9-vs-O2 refinement is REFUSED, and that refusal is measured.** SeqSero2 splits O-9 from
O-2 by comparing k-mer scores for two `tyr` genes. Those references are ~99% identical to each other, so
under BLAST both align essentially full length:

| genome (true O9) | tyr-O-9 | tyr-O-2 |
|---|---|---|
| GCA_009152955 (Enteritidis) | id 99.30 / cov 100.00 | id 99.20 / **cov 100.20** |
| GCA_009474265 (Typhi) | id 100.00 / cov 100.00 | id 99.80 / **cov 100.20** |

Coverage is degenerate and *exceeds 100* on the gapped alignment, so **a 0.2-point artifact decided an
antigen**; identity separates them by 0.1–0.2 points, which is noise. Run 1 fired `O-2` on **16
isolates: 13 miss, 3 no-call, ZERO hits** — every one a true O9 that SeqSero2 called `9`. The branch
default (O-9) is taken instead and the rule name says the discrimination did not run.

**This corrects a claim I published earlier the same day.** The reading memo argued the coverage↔k-mer
mapping was sound because both are coverage-of-reference. That holds for the branch *gates* (different
antigens, where presence-strength is what is being asked) and **fails for the paralog discrimination**,
where near-full-length homologous alignment carries no signal but a k-mer-recovery score does. The
memo's own honest-limits line said the mapping was argued and not measured, and should be validated
before the numbers were trusted. It has now been validated, and it failed in one place.

**Named cost, not hidden:** a true O-2 genome (Paratyphi A, Nitra, Kiel, Koessen) will be called O-9.
**The cohort contains no O-2 serovar at all**, so the positive direction was never testable here —
which is itself a reason not to ship a 0.1-point margin rule. `tests/test_salmserovar_o_procedure.py`
pins that no branch can emit O-2; if a faithful k-mer refinement is ever added, that test is the one
that should fail. The two `tyr` entries stay in the DB and their scores are returned.

**2 — branch C keeps this caller's measured identity-primary ranking.** Run 1 ported upstream's argmax
literally, ranking on coverage alone. SeqSero2 ranks on ONE number (k-mer recovery, which conflates how
much of the reference is present with how well it matches); BLAST gives those as TWO numbers, and which
leads was already settled here by measurement — coverage-primary picks the wrong antigen when
near-identical alleles both align full length. Run 1 regressed **six** isolates: five where the
multi-gene `O:23-gene2` out-covered the correct `O-13_wzx`, one where `9,46` out-covered `3,10`.

So **the branches are ported and the scoring is not**: preconditions from upstream, generic ranking from
the evidence this caller already has. Ties fall back to source index (position in SeqSero2's own FASTA)
rather than blastn report order, so a genome cannot call differently between runs.

## What the branches actually do, on real data

| rule | n | meaning |
|---|---|---|
| `C_argmax` | 144 | ordinary best-allele pick, identity-primary |
| `A_wbaV_default_O9_tyr_not_discriminable_by_blast` | 16 | plain O9 — previously **unemittable** |
| `B_not_in_1319_present` | 16 | reciprocal marker selects O-3,10 |
| `no_o_allele_called` | 8 | nothing hit |
| `C_guard_1319_marker_only` | 7 | the marker cannot assert its own antigen |
| `B_not_in_1319_absent` | 4 | O-1,3,19 by **absence** of the reciprocal marker |
| `C_argmax_after_1319_exclusion` | 4 | marker excluded, argmax re-run |
| `C_guard_wbaV_without_wzy` | 1 | bare wbaV winner rewritten to O9 |

The two `not_in` markers are used in **opposite** directions — one selects positively, the other can
never produce a positive call — which is why the symmetric rule that a header name suggests would have
been wrong.

## Honest limits

- **The O axis is now a PORT of SeqSero2's procedure, so O agreement with SeqSero2 is a COMPATIBILITY
  measure, not independent corroboration.** Only the wet-lab accuracy is an independent number, which is
  why the bar was written on it. Future "delta vs SeqSero2" on this axis is a regression test.
- The cell **still trails the reference tool** (0.8250 vs 0.8800, −0.0550). Better, not fixed. The
  residual is 14 misses and 21 abstentions.
- Branch-gate thresholds are an **argued** mapping from k-mer coverage-of-reference to BLAST coverage,
  and the tyr case above is a measured example of that mapping failing. The gates were not re-validated
  individually.
- Our O-membership bar is the caller's coverage threshold (**40**), stricter than SeqSero2's **25**. Left
  in place deliberately: 40 was itself adopted against a pre-registered bar, and moving two things at
  once would make neither measurable.
- SeqSero2's column is **reused** from the 2026-09-08 checkpoint (its binary, image and DB are untouched
  by this change). The script **refuses to run** unless that cache reproduces the frozen baseline tally
  exactly, so it cannot silently compare against a different denominator.
- Per-serovar cap 12 flattens prevalence: these are per-isolate accuracies on a deliberately diverse
  mix, **not** population-weighted rates.
- Residual label circularity is bounded, not eliminated. Both callers are scored on the same labels, so
  the **delta** survives contamination better than the absolute levels.
- `O:22-gene*` / `O:23-gene*` still parse to malformed antigen names (`22-gene2`). Untouched here; the
  identity-primary ranking stops them winning, it does not fix the names.

## Reproduce

```bash
# rebuild the antigen DB carrying the SeqSero2 role (needs the pinned container once)
docker run --rm quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0 \
  cat /usr/local/lib/python3.7/site-packages/seqsero2_db/H_and_O_and_specific_genes.fasta \
  > D:/dna_decode_cache/ss2_antigens.fasta
uv run python scripts/build_salmserovar_db.py \
  --antigens-fasta D:/dna_decode_cache/ss2_antigens.fasta --out data/salmserovar_db

# re-score against the frozen bar (exit 0 = ADOPT, 1 = REJECT, 3 = refused)
uv run python scripts/salmserovar_o_fix_rerun.py
```

Frozen AMR surface byte-unchanged — typing cell.
