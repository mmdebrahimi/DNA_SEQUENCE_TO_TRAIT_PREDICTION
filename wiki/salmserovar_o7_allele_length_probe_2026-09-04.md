# The O7 allele-length hypothesis is falsified — and my "concentrated on O7" claim was scoped wrong

I flagged this hypothesis in four consecutive runs as "untested". It is cheap, so it should have been
tested sooner. Testing it took one command and killed it.

## The hypothesis

The O-antigen sub-threshold hits were **11 of 14 on O group 7**, at near-perfect identity (median 99.8)
but partial coverage (median 58.4). The natural reading: **the O7 `wzx/wzy` reference allele has a
length or structure mismatch**, causing systematically partial alignments.

## Falsified

| | |
|---|---|
| O7 reference length | **1080 bp** |
| O-allele median | **1164 bp** (range 130 – 4477) |

O7 is **entirely unremarkable in length**. Length is not the explanation.

## And the "concentrated on O7" reading was scoped wrong

Looking at the *whole* alignment population rather than only the 14 hits that happened to drive
abstentions:

| O group | near-perfect-identity / partial-coverage hits |
|---|---|
| **3,10** | **13** |
| 7 | 12 |
| 4 | 7 |
| others | 7 combined |

**O group 3,10 has *more* such hits than O7.** The concentration I reported was a property of the small
abstention-driving subset, not of the underlying alignments. Both statements are now in the record with
their scopes visible, because the narrow one is still true of what it described — it just never
licensed the general claim I hung on it.

## What the pattern actually is

Near-perfect identity over roughly 60% of the reference (O7: ~678 bp aligned of 1080) is the signature
of **genuine partial homology** — the genomic `wzx/wzy` diverges from the reference over part of its
length, or the reference carries a segment these genomes lack.

That is a **DB-content question, not a length artifact** — and it is precisely why relaxing the
**coverage** cut rather than the identity cut was the change that paid.

## Honest limits

- This tests the **length** hypothesis only. It does not establish *what* the divergent region is —
  that needs aligning the reference against a genomic hit and inspecting where the alignment stops.
- One cohort (200 Salmonella), one DB build.
- A falsified hypothesis is not a closed question: the O-antigen abstentions still have a cause, and it
  is now known to be neither reference length nor a single O group.

## Reproduce

Offline from the committed antigen DB + the cached alignment pass
(`D:/dna_decode_cache/salm_asm/sweep_hits.jsonl`). See
[`salmserovar_o7_allele_length_probe_2026-09-04.json`](salmserovar_o7_allele_length_probe_2026-09-04.json).
