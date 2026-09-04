# The E. coli fix does not transfer to pneumococcus — and the cached cohort was corrupt

Two results, one of which is a data-integrity finding that would have silently poisoned any reuse.

## 1. The selection rule is live here, but it never helps

The E. coli serotype fix (coverage-only → identity-primary allele selection) lifted H accuracy
0.770 → 0.926. A sweep found the **same coverage-first pattern** in `pneumoserotype` and `plasmid`, and
those were deliberately left alone pending measurement. This measures it.

On 25 freshly-fetched GPS assemblies:

- **1 of 25 calls flips** between the two orderings — so the rule is *not* inert (my initial guess that
  100/100 winners would make it irrelevant was **wrong**: only 29% of winners sit at 100/100).
- **The single flip is wrong under both rules.** `ERS629073`, measured serotype **1**, is called `19AF`
  under coverage-primary and `19A` under identity-primary. Both are serogroup-19 misses.
- **improved 0 · worsened 0 · neither-right 1.**

**Verdict `RULE_FLIPS_BUT_NEVER_IMPROVES_ON_SAMPLE`.** No evidence identity-primary helps this cell —
the opposite of the E. coli result, and consistent with pneumococcal misses having a **different cause**:
they are systematically *within-serogroup* (9A/9V, 6B/6E, 15B/15C), a single-best-reference v0 ceiling
rather than cross-hybridization between antigens.

**Counting flips alone would have been misleading.** The first version of this probe reported
`RULE_IS_LIVE_FULL_RUN_JUSTIFIED` on the strength of that one flip. A flip only matters if it turns a
wrong answer into a right one; the verdict now scores flips against the measured label.

## 2. The cached assembly cohort is 100% unusable

Of 23 cached `.fa.gz` files under genome filenames: **11 are 199-byte HTTP 403 pages, 2 are empty, and
10 are truncated gzip streams.** Zero are usable assemblies.

The previous fetcher wrote whatever the server returned. Fed to blastn these produce a generic
*"invocation failed"* — **indistinguishable from a genuine no-match**, i.e. corrupt input wearing the
costume of a biological result. It also explains the original run's 25 "assembly-unavailable".

A fresh fetch works fine: the same accession downloads at 643 KB compressed / 2.1 Mb uncompressed,
versus the 36 KB truncated file on disk. The probe now **validates every download by decompressing it in
full** and discards partial writes rather than caching them.

**A head-only check is the wrong guard, and I wrote it first.** A truncated gzip decompresses its opening
blocks perfectly and raises only at the end — so reading the first 2048 bytes passes exactly the files
that will later fail. Validation must decompress the whole stream.

## 3. The cell was under-claimed on the trust surface

`typing:Streptococcus_pneumoniae:pneumoserotype` was registered `FAITHFUL_TO_TOOL` while its own report
card recorded **independent validation against wet-lab Quellung** (n=230). It is now
`INDEPENDENT_MEASURED`, written from the artifact:

| | value |
|---|---|
| serogroup concordance | **0.939** (n=230) |
| exact serotype | 0.661 |
| explicit-Quellung subset (n=42) | 0.952 serogroup / 0.690 exact |

**Quote the serogroup number.** Exact-serotype is the within-serogroup-limited lower bound, not the
cell's accuracy. Under-claiming is as much a trust-surface falsehood as over-claiming.

## Honest limits

- **n=25, not 260.** A single flip on 25 assemblies is thin evidence. It argues against spending ~235
  ENA fetches on a full both-rules run, but it does not prove the rule harmless.
- This measures whether the **call** changes, and whether it changes for the better at **serogroup**
  level. Exact-serotype effects at n=25 are not powered.
- The label is wet-lab phenotypic serotype, but the cps DB and serotype universe are a shared reference
  system — **reference-coupled, not circular**.
- `plasmid` still carries the coverage-only pattern and remains **unmeasured**. Same discipline applies:
  do not propagate the fix without a cohort.

## Reproduce

```bash
uv run python scripts/pneumo_selection_rule_probe.py --fetch 25
```

Needs blastn + network. Frozen AMR surface byte-unchanged.
