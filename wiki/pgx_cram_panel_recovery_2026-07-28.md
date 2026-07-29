# Read-level CRAM recovery of panel-limited PGx sites (Tier-1a, 2026-07-28)

**Result:** the "residual 3 silent" TPMT non-core samples are **proven, on real reads, to be a phased-VCF
panel limitation — NOT a decoder gap.** All 3 non-core alleles are present in the 1000G 30x CRAMs; the
sentinels would fire given full genotyping. New reusable tool: `scripts/pgx_cram_genotype.py`.

## The problem (recap)

After the sentinel layer + widening the TPMT fetch, 3 GeT-RM samples stayed silently mis-called (`*6`/`*12`/
`*40`). Root cause was **inferred**: those sites aren't genotyped in the NYGC 30x **phased VCF** — a documented
artifact (Star Allele Search, PMC10811916: *"several star-allele-defining variants present in the Phase 3 10x
dataset were absent from the NYGC phased 30x VCF files"*). But the **reads** (the CRAMs) are the full data.

## What we did (Tier-1a, $0, no permissions)

`scripts/pgx_cram_genotype.py` — read-level genotyping of arbitrary PGx sites directly from a 1000G 30x CRAM,
generalizing the CYP2D6 pileup path (Docker samtools `-B -q 0 -Q 0` + ENA reference auto-fetch, range-requests
over the remote CRAM, no full-reference download). Give it a gene's `SENTINELS` (or a raw site list); it
counts ALT vs REF reads per site → a read-level genotype call + an `alt_present` flag.

## Real-CRAM validation (the proof)

Ran it on the 3 residual-silent samples at their exact TPMT non-core sites (CRAM URLs via
`scripts/resolve_1000g_cram.py`):

| sample | GeT-RM truth | site (GRCh38) | ALT reads / depth | VAF | read-level call | ALT present? |
|---|---|---|---|---|---|---|
| NA18603 | TPMT `*1/*6` | chr6:18133845 T>A | 16 / 36 | 0.44 | 0/1 (het) | ✅ |
| NA12751 | TPMT `*1/*12` | chr6:18139710 G>A | 25 / 41 | 0.61 | 0/1 (het) | ✅ |
| HG01474 | TPMT `*1/*40` | chr6:18130729 C>T | 18 / 31 | 0.58 | 0/1 (het) | ✅ |

**3/3** — every non-core allele is present in the reads as a clean heterozygote matching the `*1/*N` truth.
So the phased-VCF simply doesn't emit these sites; the decoder's sentinels are correct and would withhold
these samples on a WGS/clinical VCF (or on read-level calls). The residual is **data-representation, not a
code gap** — now demonstrated, not inferred.

## Why this matters beyond the 3 samples (reusable infrastructure)

`pgx_cram_genotype.py` is the general form of the CYP2D6 pileup path: **any** panel-limited PGx site, **any**
gene, is now recoverable from the 1000G CRAMs we already fetch — the free, read-level substrate. It answers
"is this non-core allele's ALT actually in the reads, and at what genotype?" — the exact question that
separates a decoder blind-spot from a public-VCF filtering artifact.

## Scope / honesty

- This is a per-site ALT-presence + genotype caller (mpileup allele counting), NOT a full star-allele caller
  — sufficient to prove recoverability + to feed a read-level withhold, not to replace the VCF diplotype path.
- `UNCALLABLE` below depth 8; het band VAF 0.15–0.85 (documented, tunable).
- The GeT-RM *consensus* truth (`*1/*6` etc.) is the label; the reads confirm the ALT is real.

## Files
- `scripts/pgx_cram_genotype.py` — the general read-level site genotyper (+ `--gene`/`--sites`/`--json`)
- `tests/test_pgx_cram_genotype.py` — 7 offline tests (pure parse/call logic, mock docker_run)
- Companion to `wiki/pgx_validation_cohort_acquisition_memo_2026-07-28.md` (Tier-1a of the acquisition memo)

Frozen AMR/forward surfaces byte-unchanged (new tool + test only; no pgx-catalog edit needed — the sentinels
were already correct).
