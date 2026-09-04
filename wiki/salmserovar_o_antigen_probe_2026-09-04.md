# Two causal claims about this cell, both asserted, both measured wrong

A short chain worth recording as a method note as much as a result.

| # | claim | how it was made | outcome |
|---|---|---|---|
| 1 | "phase-2 flagellin is the single largest defect" | counted formulas ending in `-` | **wrong** — that count conflates absent H2 with nothing resolving on the H axes |
| 2 | "the real priority is O-antigen **DB coverage** — data engineering" | asserted while correcting #1 | **wrong** — 14 of 21 O-unresolved isolates *do* hit an O allele, below threshold |
| 3 | "so lower the coverage threshold" | **not made** | measured cost unknown; deliberately left open |

Claim 2 was asserted one commit after correcting claim 1 for being asserted. That is the finding worth
carrying: the failure mode is not the specific wrong answer, it is reaching for a cause without a
measurement each time.

---

## What is measured

Of 21 O-unresolved abstentions, re-probed at identity ≥ 60 / coverage ≥ 30 (deployed: 90 / 80):

- **14 have a sub-threshold O hit** — the information is already in the DB
- **7 have no O hit at any threshold** — a genuine coverage gap

**All 14 sub-threshold hits name the *correct* O group.** Zero name a wrong one. So the rejected
alignments are correct matches being discarded, not noise being correctly filtered.

**The shape is sharper than "a threshold problem":**

| | value |
|---|---|
| identity (median) | **99.8** — near-perfect |
| coverage (median / max) | **58.4 / 78.9** — all below the deployed 80 cut |
| O group | **7 (×11)**, 40, 18, 3,10 |

Near-perfect identity at partial coverage, **concentrated on the O7 wzx/wzy allele** (11 of 14). That is
the signature of an allele-length or partial-alignment mismatch on one reference, not a uniform
threshold miscalibration across the DB.

## What is deliberately NOT concluded

**Nothing about what to change.** Having been wrong twice by asserting a cause, the obvious third
assertion — "relax the coverage cut" — gets no free pass:

- **Lowering the cut trades abstentions for wrong calls, and that cost is not measured here.**
  Abstention is the safer failure; this cell already abstains on 29.5% and is *not* deployed as a
  drop-in caller. Turning silence into confident error would be the wrong direction to be wrong in.
- **Whether the fix is the threshold, the O7 reference allele, or partial-alignment handling is a
  hypothesis.** The concentration on one O group argues for an allele-level cause, but that is exactly
  the kind of inference this memo exists to stop asserting.
- **`Infantis` appears in both buckets** — 3 sub-threshold, 4 no-hit-at-any-threshold. O7 detection is
  inconsistent *within a single serovar*, which neither hypothesis explains.

The decisive next measurement is not another diagnosis: it is running the relaxed threshold over the
whole 200-isolate cohort and counting **new wrong calls against new rescued calls**. Until that exists,
the honest state is "the information is present and correctly matched; the cost of admitting it is
unknown".

## Honest limits

- 21 isolates from one cohort and one antigen-DB build (49 O groups present, 62 alleles).
- Only the O-unresolved abstentions were probed — this says nothing about the H axes, nor about
  isolates that already resolve.
- A permissive hit demonstrates the information is present. It does **not** demonstrate that admitting
  it is safe.

## Reproduce

```bash
uv run python scripts/salmserovar_nocall_anatomy.py    # which axis fails (offline)
uv run python scripts/salmserovar_o_antigen_probe.py   # threshold vs coverage (needs blastn)
```

Frozen AMR surface byte-unchanged — typing cell.
