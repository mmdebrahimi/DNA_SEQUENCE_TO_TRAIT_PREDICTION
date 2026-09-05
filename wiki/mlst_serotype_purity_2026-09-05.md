# The MLST caller recovers real lineages — measured against a wet-lab label, with no curated biology

`typing:bacteria:mlst` ships CLI-routable and was registered `FAITHFUL_TO_TOOL` — checked against the
reference *method*, never against reality. It became cheap to test the moment the serotype
lineage-disjoint run called sequence types on 400 real E. coli genomes that **also** carry a wet-lab O:H
serotype label. Zero new compute: this reads that run's committed checkpoint.

## The test, and why it needs no remembered biology

E. coli clonal lineages conserve their O:H antigens, so a working MLST caller must produce sequence types
that are **serotype-pure**. A broken one — wrong allele matching, a mis-joined profile table, an
off-by-one in the profile lookup — carves the cohort arbitrarily and lands at chance.

The obvious alternative would score against remembered associations ("ST131 is O25:H4"), which means
**asserting biology from memory as the scoring key**. This measures per-ST purity against a shuffle null
instead and lets the data speak.

**The null is the load-bearing part.** Purity alone is gameable in both directions: a caller giving every
genome its own ST scores 1.000, and one pooling everything into a single ST scores the cohort's modal
frequency. So the null **shuffles the wet-lab serotype labels while holding the observed ST partition
fixed** — identical group sizes, identical serotype frequencies, only the association destroyed.
Over-splitting inflates null and observation equally and cannot manufacture a result.

## Result

398 genomes, 145 sequence types, 27 reaching the size-3 floor and covering 271 genomes (68.1%), against
**170 distinct wet-lab serotypes**.

| statistic | observed | null mean | **null max** (1000 shuffles) |
|---|---|---|---|
| serotype (O:H) purity | **0.7860** | 0.2156 | 0.2583 |
| top-2 serotype purity | **0.9188** | 0.3699 | 0.4170 |
| O-antigen only | **0.8266** | 0.2235 | 0.2694 |

Every statistic is compared against **its own matched null** — scoring a top-2 figure against a top-1
null would manufacture significance. All three exceed the null's *maximum* over 1000 shuffles, roughly
three times over.

**Verdict `MLST_RECOVERS_REAL_LINEAGES`.**

## The corroboration I did not supply

The modal serotypes fell out of the data. They are the textbook lineage↔serotype pairings:

| ST | n | modal wet-lab serotype | |
|---|---|---|---|
| ST11 | 53 | **O157:H7** | 53/53 = 1.00 |
| ST21 | 50 | O26:H11 | 48/50 = 0.96 |
| ST17 | 23 | O103:H2 | 20/23 = 0.87 |
| ST655 | 7 | O121:H19 | 7/7 = 1.00 |
| ST678 | 6 | O104:H4 | 5/6 = 0.83 |

That ST11 comes back as O157:H7 on **53 of 53** genomes, without O157:H7 ever being named in the scoring
code, is the check the design was built to allow.

## The low-purity STs are real biology — and two of them strengthen the result

These were **measured**, not interpreted away. The composition of every ST below 0.70 purity ships in the
artifact.

- **ST131 — O25:H4 (11) + O16:H5 (9) + O25:H17 (2) + O62:H4 (1).** Top-1 purity 0.48 looks like a
  failure; the two leading serotypes are the two canonical ST131 clades, and together they are 20/23 =
  0.87. This is a genuinely **bimodal** lineage being recovered correctly, which is exactly what the
  **top-2 statistic (0.9188 against a null max of 0.4170)** measures rather than assumes.
- **ST16 — O111:H8 (8) + O111:H− (4).** Every genome is O111; the split is entirely on the H/motility
  axis. That is what the **O-antigen-only statistic** exists to capture, and it is higher (0.8266) than
  the combined figure (0.7860) as predicted.
- **ST10 — 11 distinct serotypes across 14 genomes**, purity 0.14. ST10 is the well-known generalist
  commensal lineage. **Low purity here is the correct answer**, and a caller that made ST10 look pure
  would be more suspicious, not less.

## Honest limits

- **This is a COHERENCE check, not a correctness check.** It shows the sequence types track real clonal
  structure. It does **not** show they carry the same *numbers* a reference MLST implementation would
  assign — a caller with a systematically shifted profile table would pass this and still report wrong ST
  numbers. That needs the reference tool installed and pinned locally, which was **not** done. **The cell
  stays `FAITHFUL_TO_TOOL`.**
- One organism (E. coli), one scheme (Achtman 7-locus), one cohort.
- Genomes within a sequence type are clonal, so the genome-level p-value overstates independence. **The
  claim does not rest on the p-value** — it rests on the observation exceeding the null's *maximum* over
  1000 shuffles.
- Serotype labels are NCBI-PD submitter strings; E. coli O:H typing is traditionally slide agglutination,
  but per-isolate method is not provable from the metadata.
- STs below the size-3 floor are excluded, so this says nothing about singleton or near-singleton
  lineages — and 31.9% of the cohort sits below that floor.

## Reproduce

```bash
uv run python scripts/mlst_serotype_purity.py
```

Offline — reads the committed lineage-disjoint checkpoint. No blastn, no network, no Docker.
See [`mlst_serotype_purity_2026-09-05.json`](mlst_serotype_purity_2026-09-05.json).
