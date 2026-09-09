# Scored against the real reference tool: the caller trails it, and the whole gap is one axis

The 2026-09-04 salmserovar validation measured our caller at 0.783 against a comparator at 0.925 and
called that a −0.142 delta "vs the in-silico incumbent". That comparator was NCBI Pathogen Detection's
`computed_types` — a **production call of undocumented tool, version and configuration**. The artifact
recorded the limit honestly: a delta against a production field cannot separate

- **(a)** our caller is worse than the reference **method**, from
- **(b)** our caller merely **diverges** from one undocumented NCBI pipeline.

Resolving it needed SeqSero2 installed and pinned locally. This is that run — **SeqSero2 1.3.2** from
the pinned biocontainer, the **same 200 isolates**, the **same wet-lab agglutination labels**, and the
same scoring function for all three callers.

## Result

| caller | accuracy | hit / miss / no-call | abstain |
|---|---|---|---|
| **ours** | **0.7050** | 141 / 39 / 20 | 10.0% |
| **SeqSero2 1.3.2** | **0.8800** | 176 / 24 / 0 | 0% |
| NCBI `computed_types` | 0.9200 | 184 / 15 / 1 | 0.5% |

**Delta vs the reference tool: −0.1750.** The ambiguity resolves toward reading **(a)**: our caller
genuinely trails the reference *method*, not merely an undocumented pipeline. That is the honest
outcome, and it is worse for us than the alternative reading would have been.

## The entire gap is the O antigen

SeqSero2 emits per-axis predictions, so the earlier O-axis suspicion could be tested directly rather
than inferred from abstention anatomy:

| axis | agreement where both resolved |
|---|---|
| **H1 (fliC)** | **200 / 200 = 1.000** |
| **H2 (fljB)** | **146 / 146 = 1.000** |
| **O** | 159 / 200 = 0.795 |

**Our H-antigen calls agree with the reference tool perfectly.** The formula lookup and the H axes are
not the problem; the O antigen is, entirely.

Splitting O further — abstention is not error:

| | |
|---|---|
| O agrees | 159 |
| O **resolved but different** | **33** |
| O unresolved (`O?`) | 8 |
| O resolution rate | 192/200 = 0.960 |
| O accuracy **where we resolved** | 159/192 = **0.828** |

So the O problem is mostly **confidently wrong calls (33)**, not abstentions (8) — a different emphasis
from the 2026-09-04 no-call anatomy, which was about where the caller stays silent.

## The O errors are two strings

| ours → SeqSero2 | n |
|---|---|
| **`9,46` → `9`** | **15** |
| **`1,3,19` → `-` / `4` / `9` / `8`** | **13** |
| `3,10` → `1,3,19` | 3 |
| `22-gene2` / `23-gene3` → `13` | 2 |

**28 of 33 O disagreements (85%) are two patterns**: `9,46` over-called where the reference says plain
`9`, and `1,3,19` emitted where the reference says something else entirely.

`1,3,19` **independently corroborates a prior observation**: the 2026-09-04 run, using a completely
different route (a comparator field, not a reference tool), recorded Typhi being called `1,3,19` rather
than `9,12`. Two unrelated methods landing on the same string is meaningful.

The malformed DB names `22-gene2` / `23-gene3` are the previously-recorded DB-hygiene defect, and are a
small contributor (2 of 33).

**This is a lead, not a diagnosis.** The concentration names *where* to look; it does not establish
*why* those strings are emitted. Fixing them without first measuring the cause is the error this project
keeps catching in itself.

## Two things not to over-read

**SeqSero2 itself scores 0.880 here, below its own published accuracy (~0.95).** That is expected and it
matters: the cohort caps each serovar at 12 isolates, deliberately flattening prevalence, so it
over-weights rarer and harder serovars. That depresses **every** caller. **Quote the delta, not the
levels** — the delta shares a cohort, labels and scoring; the levels do not transfer to a population.

**NCBI `computed_types` (0.920) beats the standalone reference tool (0.880) on these isolates.** Reported
as an observation, not explained: it could be a newer or differently-configured pipeline, a combination
of callers, or residual label correlation that the earlier circularity probe bounded but could not
eliminate. Nothing here distinguishes those.

## Honest limits

- **SeqSero2 is a TOOL, not the wet-lab assay.** It is scored against the same agglutination labels as
  our caller, so this compares two in-silico callers against a shared external truth — not a claim that
  SeqSero2 is correct.
- **Fairness is structural:** all three callers go through the same `equivalence.equivalent` (notation
  normalisation plus the committed Kauffmann-White table, no fuzzy near-misses), and abstentions are
  counted separately from misses for each. No caller gets leniency another is denied.
- Per-serovar cap 12 flattens prevalence, so these are per-isolate accuracies on a deliberately diverse
  mix, **not** population-weighted rates.
- Residual label circularity is **bounded, not eliminated** — per-isolate agglutination provenance is
  unprovable from metadata. Contamination would inflate every caller, which is another reason the delta
  is more trustworthy than the levels.
- SeqSero2 was run in **k-mer mode on assemblies** (`-m k -t 4`); its read-based or microassembly modes
  could score differently.
- 0 isolates were excluded — every one had an assembly and SeqSero2 errored on none — so all three
  callers share the full 200-isolate denominator.

## Reproduce

```bash
docker pull quay.io/biocontainers/seqsero2:1.3.2--pyhdfd78af_0
uv run python scripts/salmserovar_seqsero2_validate.py
```

Two traps encoded in the runner: SeqSero2 **basenames** its input and chdirs, so the file must be staged
into the mount with the container's workdir set to it; and the bind-mount source must be a **Windows-form
path** — a Git-Bash `/tmp/...` source silently produces an anonymous volume while the command still exits
0. Frozen AMR surface byte-unchanged.
See [`salmserovar_seqsero2_2026-09-08.json`](salmserovar_seqsero2_2026-09-08.json).
