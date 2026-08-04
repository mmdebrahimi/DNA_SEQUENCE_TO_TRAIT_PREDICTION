# Motility catalog v0 — the first NON-metabolic trait catalog (2026-08-03)

`dna-decode motility` decodes **flagellar swimming motility** from gene presence — the AMR/metabolic
determinant→phenotype paradigm applied, for the first time, to a **non-metabolic physical behaviour**. It's
the ④ rung of the edit→cell-trait ladder: a curated determinant catalog for a trait FBA can't model.

```
dna-decode motility --genes flhD,flhC,fliA,fliC,motA,motB,fliF,fliG,flhA,fliI   -> MOTILE
dna-decode motility --genes fliC,motA,motB                                       -> NON-MOTILE (no master flhDC)
dna-decode motility --feature-table GCF_..._feature_table.txt.gz                 -> genome-wide
```

## The rule (5 modules; MOTILE iff all present)

Flagellar motility is a strict regulatory + structural cascade — missing any load-bearing module = no
working flagellum:

| module | genes | why it gates |
|---|---|---|
| master_regulator | flhD **AND** flhC | class-1 master; no flagellar gene transcribes without it |
| sigma28 | fliA | class-3 sigma; filament + motor aren't expressed without it |
| flagellin | fliC **OR** fljB | the filament protein |
| motor | motA **AND** motB | stator/torque generator (a flagellum that can't rotate is non-motile) |
| basal_export | fliF, fliG, flhA, fliI | MS-ring/switch + type-III export apparatus |

**Chemotaxis (cheA/W/Y/Z) is reported SEPARATELY and does NOT gate motility** — a che-mutant still *swims*
(it just tumbles randomly / can't chase a gradient). Gating swimming on chemotaxis would be a biology error;
the cell reports `chemotaxis_competent` as an independent field.

## Validation (KNOWLEDGE_BASELINE — curated catalog vs literature anchors)

| genome | call | anchor |
|---|---|---|
| E. coli K-12 MG1655 (full set) | **MOTILE** + chemotactic | motile ✓ |
| Salmonella Typhimurium | MOTILE | motile ✓ |
| Shigella flexneri (flagellar pseudogenes) | **NON-MOTILE** | non-motile ✓ |
| ∆flhDC / ∆fliC / ∆motAB | NON-MOTILE (that module) | knockout non-motile ✓ |

There is **no free, genome-keyed swim-plate motility cohort** fetchable (the recurring label wall), so this
is a KNOWLEDGE_BASELINE tier — a curated catalog validated against literature-known anchors, like the
`metabolic` cell — **not** a big measured-cohort claim.

## Honest scope (load-bearing)

- v0 = the flagellar **swim/no-swim DIRECTION** from gene PRESENCE. NOT swim speed / rate / gradient chase.
- **Presence-based:** it cannot see a gene that is annotated PRESENT but silently inactivated — the classic
  **K-12 flhD IS1-insertion** (many lab strains are non-motile with a present-but-dead flhD) would be
  mis-called MOTILE. A sequence-integrity genome-mode is the named v0.1 follow-on (not fabricated here).
- Type-IV-pilus **twitching**, gliding, and swarming-specific regulators are OUT of v0 (flagellar swim only).
- NON-frozen cell (like metabolic/flowering); the frozen AMR decoder surface is byte-unchanged. NOT clinical.

`dna_decode/motility/flagellar_catalog.py` + `cli.py`; 12 tests `tests/test_motility.py`. The pattern is
reusable for the next non-metabolic trait catalogs (biofilm: curli/cellulose/PGA; acid resistance; etc.).
