# A third archive returns zero carriers — and zero is worth recording

**Oxford has no `rmt` at all.** Across **4,979** AMRFinder-scanned *E. coli* from the Oxford bacteraemia
cohort — 2,873 of them carrying a measured broth-microdilution gentamicin MIC — there is **not one
`rmt` or `npmA` carrier**. One `armA`, and nothing else in that gene family.

So Oxford **cannot** test the gentamicin v2 rescue's specificity. No carrier, no test. What it does
instead is put a number on how rare the target is, and that number changes how one standing limit reads.

---

## Why Oxford was the right third opinion

The [archive search](rmt_independent_archive_search_2026-09-03.md) established that the corpora richest in
`rmt` are the ones **structurally incapable** of answering this question: they ascertain on high-level
aminoglycoside resistance, so they contain zero susceptible carriers by construction. Oxford inverts every
one of those problems:

- **ascertained on bacteraemia**, not on resistance — the sampling is independent of the outcome;
- **independent on both axes** — Oxford's own MIC *and* Oxford's own AMRFinder run, not NCBI-PD re-served;
- **already on disk and already plumbed** (`scripts/oxford_score.py`), so asking cost nothing.

It is the same cohort on which the frozen decoder already scored cipro 0.960 / gent 0.990 accuracy
(2026-06-15) — a cohort we trust and have used.

## The zero is proven, not merely reported

A zero from a broken parse looks exactly like a zero from real absence, so the probe refuses to report
unless the scan demonstrably worked. It found **26 distinct aminoglycoside determinant symbols over 3,771
hits**, including six present **exactly once** — `aac(2')-IIa`, `aac(6')-Ib4`, `aadA15`, `aadA16`,
`aadA4`, `aph(3')-IIa`. **Genes present once are detected, so a gene present once would have been found.**
The guard is pinned in both directions by test: it refuses (exit 3) on a dead scan and leaves no artifact
behind, and it does not block a live one.

## What this changes, and what it does not

**It reframes the E. coli scope limit.** The standing caveat has been *"whether the E. coli scope is
genuinely safe or merely under-sampled — twelve carriers is not many."* Twelve is still not many, but the
reason is now measured: **`rmt` in E. coli is genuinely rare in unselected collections**, at under 1 in
4,979 here. The twelve are few because carriers are few, not because the search was shallow. That does not
make the E. coli scope safe — it makes "safe" harder to demonstrate for anyone, by any means.

**It closes Oxford as a candidate**, so the archive is not re-attempted. The remaining lever for the
over-call question is unchanged: a Klebsiella-rich, resistance-independent collection, which is where the
measured over-call (PPV 0.475) actually lives.

**It corroborates the v2 rescue's premise on a second AMRFinder version.** The single `armA` carrier is
filed under `Subclass=GENTAMICIN` by this older, independent run — confirming that the frozen rule always
counted `armA` and that the gap was the `rmt` family only. That distinction is load-bearing: the deployed
`symbol_rescue` deliberately excludes `armA`, and this is independent evidence the exclusion is right.
The carrier's MIC is 32 (resistant), consistent with the determinant.

## What it explicitly does not show

- **Nothing about the rule's specificity.** Zero carriers is neither support nor counter-evidence. A
  cohort containing none of a determinant is structurally incapable of detecting a rule keyed on it —
  the same blindness the [source-concentration layer](gentamicin_rmt_bvbrc_2026-09-03.md) exists to
  disclose, arriving here as absence rather than concentration.
- **No generalisation past this population.** One UK region, 2008–2018 bacteraemia. `rmt` is
  substantially more prevalent in Asia and in *Klebsiella*.
- **DB coverage is unverified.** This deposit's AMRFinder version is not recorded. Total blindness to 16S
  methyltransferases is unlikely — `armA`, the same family and curation era, *is* detected — but
  `rmt`-specific coverage in that version was not confirmed.

## Reproduce

```bash
uv run python scripts/oxford_rmt_prevalence_probe.py
```

Offline (reads the Oxford deposit on `D:`). No frozen file is touched.
