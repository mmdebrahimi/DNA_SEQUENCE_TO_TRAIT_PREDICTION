# PointFinder reproduces an independent caller on the QRDR — and the one real miss is an alignment failure, not a calling error

The last untested finder cell. `finder:Escherichia_coli:pointfinder` shipped `FAITHFUL_TO_TOOL` —
checked against the reference *method*, never against reality. It is the exact complement of the
ResFinder locus-collapse work, which compared **acquired genes** and could not touch point mutations.

**Why it matters beyond the cell:** the catalogued genes are gyrA / gyrB / parC / parE — the QRDR — and
the **frozen cipro decoder rule (`qrdr_point`) consumes exactly these determinants from AMRFinder**. So
this independently re-derives a determinant class the *deployed* surface depends on. Read-only; the
frozen surface is untouched.

## Where the comparator actually lives

AMRFinder's point screen is **not** in `main.tsv` — across all 1,818 committed runs that file's `Type`
is only `AMR`/`STRESS`, with **zero POINT rows**. The screen is in a separate `mutations.tsv`, in
`--mutation_all` form, so it lists every *screened* position including the clean ones marked
`[WILDTYPE]`. Keeping those would turn "screened and found nothing" into a called mutation; they are
filtered out.

## Result — 648 genomes

| | |
|---|---|
| scored (aligned ≥1 reference gene) | **646** |
| abstained (aligned nothing) | **2 (0.3%)** |
| AMRFinder calls, all | 4,951 |
| AMRFinder calls, catalogued positions | 1,527 — **restriction removed 3,424** |
| PointFinder calls | 1,527 |
| **agree / PF-only / AMR-only** | **1,524 / 3 / 3** |
| genomes with an exact set match | **641/646 = 0.9923** |
| mean per-genome Jaccard | **0.9967** |

Per-mutation, on the determinants the cipro rule cares about: **gyrA_S83L 418 both / 1 PF-only / 0
AMR-only**, **parC_S80I 353 / 0 / 0**, **gyrA_D87N 324 / 1 / 0**.

**Verdict `POINTFINDER_AGREES_WITH_AN_INDEPENDENT_CALLER`.**

## The six discordances, read rather than averaged away

**All three AMR-only calls are catalogue differences, not misses.** Each sits at a position PointFinder
catalogues but names a residue its catalogue does not list as resistance-conferring, so abstaining is
correct *by its own catalogue*:

| AMRFinder call | PointFinder's Res_codon set at that position | |
|---|---|---|
| `parE_S458W` | S458 → A, T | W not catalogued |
| `parE_L445F` | L445 → H | F not catalogued |
| `parE_E460A` | E460 → D | A not catalogued |

**So across 646 scored genomes the calling logic has zero genuine misses.**

The three PointFinder-only calls are `parE_I529L` in one genome and `gyrA_S83L` + `gyrA_D87N` in
another. The second is worth noting — those are the two commonest cipro determinants and AMRFinder did
not report them there — but it is **one genome**, and "PointFinder-extra" is not the same as
"PointFinder-correct". A tool comparison cannot adjudicate that.

## The single real false negative traces to abstention, not calling

Two genomes (0.3%) aligned **no** reference gene at all — one is a 137-contig draft where the genes are
presumably split across contigs. AMRFinder found 3 catalogued mutations in those two, and exactly one of
them, **`parC_S80I`, is a residue PointFinder does catalogue**. That is the only genuine miss in the
entire run, and its cause is the failed alignment rather than the calling logic.

**This is visible to a user, which is why it is not filed as a defect.** The shipped CLI prints
`(no catalogued point mutations; genes aligned: none)`, so a reader can tell an abstention from a clean
result. What is *not* distinguished is the machine-readable `status`, which stays `"ok"` either way —
the distinguishing information lives in `genes_aligned`. That is a weaker, real observation, and it is
recorded rather than dressed up as a bug.

## Two accounting traps this run had to fix

1. **Abstention is not agreement.** The first version scored a genome where PointFinder aligned nothing
   against an equally empty comparator set as a *perfect match* — counting a genome the caller did no
   work on as a success. Abstained genomes are now excluded from the agreement metrics and reported
   separately, together with how many real comparator calls fell inside them.
2. **A restriction that removes nothing is not a control.** The sibling ResFinder run shipped claiming a
   POINT-row exclusion that excluded nothing, and had to be corrected. So this run **reports** what the
   position restriction removes (3,424 of 4,951) and refuses to call it meaningful if it removes none.

## Honest limits

- **The comparator is a TOOL, not a wet-lab label.** This measures agreement with an independent
  implementation, not correctness — both callers could be wrong together, and both ultimately derive
  from the same published QRDR literature. **The cell stays `FAITHFUL_TO_TOOL`.**
- Scope is the four genes with a committed reference CDS. The PointFinder catalogue also lists rpoB /
  16S / 23S / ampC-promoter / pmrAB / folP positions with **no reference sequence on disk**; this says
  nothing about them.
- Agreement on a position *both* callers catalogue does not validate either **catalogue**. A resistance
  position missing from both is invisible here.
- Position-level restriction does not equalize the two catalogues — residue-level differences remain,
  which is exactly what the three AMR-only calls are.
- Epistasis (`Required_mut`) is recorded but **not enforced** by the v0 caller, so a mutation whose
  effect depends on a partner is still reported as conferring resistance.
- Genomes are AMR-cohort leftovers, **enriched for resistance**, not a random sample.

## Reproduce

```bash
uv run python scripts/pointfinder_amrfinder_concordance.py
```

Needs blastn + cached assemblies + committed AMRFinder runs. Frozen AMR surface byte-unchanged.
See [`pointfinder_amrfinder_concordance_2026-09-05.json`](pointfinder_amrfinder_concordance_2026-09-05.json).
