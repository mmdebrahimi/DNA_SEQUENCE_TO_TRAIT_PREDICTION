# PEAR's last open gate is closed — G6 PASSES, and the trap that nearly said otherwise

**G6 was the one rejection gate PEAR left open**, for a blunt reason: nobody had read a single PEAR
fitness value. The repo ships no plain text — only two `.RData` workspaces. Today they were read, and
the gate was **measured** rather than declared open.

**Verdict: G6 PASS. The ten-gate screen now returns `CLEARS` (was `INCOMPLETE`).**

---

## The wall was disk, and it was self-inflicted

The recorded blocker was *"needs R, which is not installed; C: is at 99% (3.0 GB free)"*. That was true
and it was also the wrong place to be looking. **D: has 4.1 TB free.** Installing R there costs nothing
and touches C: not at all:

```
micromamba (11 MB, no admin)  ->  D:/tools/micromamba.exe
R 4.6.1                       ->  D:/tools/r_env
PEAR repo + extracted TSVs    ->  D:/dna_decode_cache/pear/
```

The blocker was never "R is unavailable" — it was that the install was being aimed at the wrong drive.

## Why R, and not another Python attempt

Not a preference. `pyreadr` fails, and `rdata` fails **at the parser, not the converter**:
`rdata.parser.parse_file` raises `NotImplementedError: Type RObjectType.WEAKREF` on a weakref nested
inside serialized *bytecode* — a closure carried by the ggplot objects. Parsing is sequential, so the
failure cannot be stepped over without risking silent byte-desync, which would yield numbers that *look*
like data. R reads its own format exactly. That is the entire argument.

## What the workspaces actually contain

**No source data frames** — only plot objects. But a ggplot object carries its data in `$data`, and
these are genome-wide tile maps, so that slot is the full scan:

| object | shape | what it is |
|---|---|---|
| `Figure.2A` | 3,957 × 9 | ceftazidime, **nucleotide-level** (`pos` 1–792 × `gt` ∈ A/C/G/T + per-position mean) |
| `Figure.2B` | 3,957 × 9 | cefotaxime, same shape |
| `Figure3.A` | 2,114 × 6 | **per-variant** in `C648T` notation, with `CTX` and `CAZ` columns |
| `Figure.2F/2G`, `Figure3.B` | 12 / 5 / — | small summary panels; `3.B` carries no data slot |

**This answers U5, and not entirely favourably.** These are the **aggregated per-variant effect sizes
the authors plotted — NOT the ~23,000 raw barcoded strains.** The scan is complete for the positions
shown (792 nucleotide positions, all four bases), but the raw library is not here. Anyone quoting
"23,000 variants" off this extraction would be wrong.

## The measurement

`assay_degeneracy` from the shipped forward/inverse cell — same function, same bars (mode-share > 25%,
or fewer than 20 distinct levels).

| table | n | distinct | mode share | degenerate |
|---|---|---|---|---|
| Figure.2A ceftazidime, per-nt | 2,165 | 2,105 | 0.0051 | **no** |
| Figure.2B cefotaxime, per-nt | 2,163 | 2,151 | 0.0018 | **no** |
| Figure3.A cefotaxime, per-variant | 2,114 | 2,106 | 0.0019 | **no** |
| Figure3.A ceftazidime, per-variant | 2,114 | 2,066 | 0.0047 | **no** |
| *Figure.2A with WT rows left in (diagnostic)* | *2,957* | *2,106* | *0.2678* | ***yes*** |

That last row is the point of this memo.

## The trap: the assay looked censored, and it was the baseline

The first run flagged the nucleotide tables as **degenerate at mode-share 0.2678** — just over the 25%
bar. That would have rejected a usable substrate.

The number is exactly `792 / 2957`, which is one row per position. Checked directly: the modal value is
`effect_size == 1` (relative growth 1.0), it occurs 792 times, **all 792 are rows where `gt == isWt`,
and zero non-WT rows carry it.** It is the wild-type base at each position — the normalizer measured
against itself, 1.0 *by construction*. Not a measurement of anything, and certainly not evidence of
censoring.

Excluding it, mode-share falls from 0.2678 to 0.0051 and the verdict inverts.

**This is the same defect class as the NNRTI `L234L` / `K238K` entries** — a non-variant admitted as a
variant, where WT and MUT are the same symbol. It has now bitten twice in three days, in unrelated data,
and in **opposite directions**: there it inflated a catalog's apparent performance, here it would have
condemned a good assay. The diagnostic row is kept in the artifact permanently so the effect stays
visible rather than being quietly corrected away.

**It is also the exact mirror of CcdB.** There, a 79.3% tie at the assay ceiling was *real* censoring
that flattered every metric. Here a 26.8% tie was an artefact of including the reference. A degeneracy
screen has to ask *what the tied rows are*, not just how many there are.

## What this does and does not license

- **It closes G6 and completes PEAR's gate screen.** All ten gates now resolve; verdict `CLEARS`.
- **It does not say the forward cell will work on CTX-M-14.** A cleared gate bounds the *label*; transfer
  is an empirical question, and the honest prior is mixed — our forward path was validated on TEM-1, a
  different β-lactamase.
- **It does not deliver ~23,000 variants.** ~2,100 per drug, aggregated.
- `Figure3.A`'s `C648T` notation is nucleotide-level and maps directly onto the shipped genome-edit path
  (`predict_genome_edit`, `--cds-fasta` + HGVS `c.`), so the natural next step is a real comparison — not
  another screen.

## Reproduce

```bash
micromamba run -p D:/tools/r_env Rscript scripts/pear_extract_fitness.R   # -> D:/dna_decode_cache/pear/extracted
uv run python scripts/pear_g6_screen.py                                   # -> wiki/pear_g6_screen_2026-09-01.json
```
