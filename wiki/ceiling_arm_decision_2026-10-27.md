# The two "data-blocked" ceiling arms are not one wall (2026-10-27)

**Why this exists.** `wiki/ceiling_arm_wall_classification_2026-10-26.md` classified three untried
ceiling arms and filed two of them under one label: **"blocked on data acquisition."** The user asked
for the problem in plain terms, for the path forward, and specifically *"why not run the whole thing on
server1?"* Answering that required separating three different things the one label had merged:
**assay physics**, **data availability**, and **compute**.

**Nothing was killed, re-scoped, downloaded, or run.** Both arms remain `planned` in
`soraya_runs/arms-ledger.json`, byte-untouched. Changing either changes what `blocking()` reports, and
R6 permits a kill only on the arm's own executed evidence — which a substrate fact is not. This memo
is the analysis the decision rests on, not the decision.

## Correction to the 2026-10-26 memo

That memo reported both arms as blocked on data acquisition. That was measured against **our attached
drives**, not against what exists. One of the two is not data-blocked at all.

| arm | 2026-10-26 | measured 2026-10-27 |
|---|---|---|
| `multiresolution-pyramid-4-decades-UNRUN` | needs data not on any attached drive | **no such dataset exists anywhere** — the wall is the ASSAY, not the disk |
| `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` | needs a panel not on any attached drive | **public and free** — a download, not an acquisition decision |

Filing an impossibility and a download under one "expensive/blocked" tier is why `blocked:untried-arm`
read as an undifferentiated wall for two consecutive `--advance` runs.

## Arm 1 — the 4-decade pyramid: blocked by ASSAY PHYSICS

### What it asks

A multiresolution encoder reads a sequence at several widths at once (`PYRAMID_WIDTHS = (3, 12, 48,
150)` in `dna_decode/glm/encoder.py`). "Decades" is **sequence length**, pinned by arithmetic rather
than inferred:

| quantity | value | check |
|---|---|---|
| `rungs.fragments_LB.max_span_decades_at_3bp_floor` | 2.22 | `log10(499/3) = 2.2210` ✓ |
| the measured sibling arm `pyramid-4scale-1.70decades-150bp` | 1.70 | `log10(150/3) = 1.6990` ✓ |

So the field is `log10(max_sequence_length / 3)`, and four decades needs **~10,000 bp** (1 bp floor) or
**30,000 bp** (the stated 3 bp floor). **The longest cached sequence is 499 bp — ~20× short** [grounded:
`wiki/mpra_substrate_census_2026-10-12.json`, `rungs.fragments_{LB,M9}.max = 499`].

### Why no dataset closes it

A pooled reporter assay measures activity by synthesizing or shearing many short elements and reading
them together. Published designs sit in the hundreds of bp: a CapSTARR-seq comparison reports fragments
with a median of **~464 bp** against 300 bp target regions [grounded: *Challenges and considerations for
reproducibility of STARR-seq assays*, Genome Res 33(4):479]. lentiMPRA elements are ~200–230 bp
[grounded: Agarwal 2025, below]. Oligo-pool synthesis caps at roughly 200–300 nt
[**inferred** — general knowledge; a targeted search did not return a vendor-sourced limit, so this is
labelled, not asserted].

A 10,000 bp element therefore means one construct at a time, which yields tens of measurements rather
than the 10^4–10^6 a model needs. **This is a property of the measurement technology, not of our disk,
our budget, or our hardware.** Acquiring more MPRA data cannot fix it; a different *assay class* or a
different *label type* would be required.

### The encoder's own comment asserts this, and its scope is wrong

`dna_decode/glm/encoder.py` carried, flagged **`(measured, not estimated)`**:

> every category in GSE144621 is <= 150 bp ... so the absolute ceiling is log10(150/1) = 2.18 decades
> even with a degenerate 1 bp branch

`grep -n "499\|fragment" dna_decode/glm/encoder.py` returned **nothing** before today, so that module
never knew the fragment rung exists. The claim is true of the **peak-tile file** and false of the
**accession**, which also ships 321,123 fragments at 48–499 bp from the same cache directory. The
reachable fragment span, `log10(499/3) = 2.22`, **exceeds the 2.18 the comment called absolute.**

Corrected in place as a scoped addition; the original reasoning is kept and **`PYRAMID_WIDTHS` is
untouched**, so no behaviour changed. Third instance of the forgotten-substrate error, and the loudest,
because `(measured, not estimated)` is exactly the marker that makes a reader stop checking.

### Recommended: RE-SCOPE to 2.22 decades on the fragment rung

- **+0.52 decades** over the measured sibling (1.70 → 2.22), i.e. **3.33× more span**, reaching a real
  biological scale at the top end rather than buying span at the bottom. The encoder's own comment
  already rejects the cheap alternative — *"a width-1 kernel cannot see a motif at all, so buying span
  by shrinking the finest branch to 1 bp inflates the decade metric without testing the hypothesis the
  metric stands for."* The fragment re-scope is the honest direction: every branch stays a scale
  someone would choose on its own.
- **Substrate is already local**: 321,123 rows (LB) / 318,457 (M9), of which the largest measurement
  has used **30,000 = 9.3%**. Ceiling 0.7928 (LB) / 0.8014 (M9).
- **Code-closable**: `PYRAMID_WIDTHS` + `dna_decode/glm/genomewide.py::load_fragments` both exist.
- **It merges with the size-lever question into ONE run** rather than two.

**Open design point for the pre-registration, not for improvisation:** fragments are VARIABLE length
(48–499 bp, 359 distinct lengths), and `load_fragments`' own docstring warns that *"k-mer features are
length-invariant and are therefore the honest default here."* A fixed-window pyramid needs a declared
pad/crop rule, and that rule must be frozen before the run, not chosen after seeing a number.

**Do NOT kill this arm.** Killing it discards the record that the multi-scale question is open at a
reachable scale. Re-scoping keeps the question and makes it answerable.

## Arm 2 — the eukaryotic panel: NOT data-blocked

- The entire local MPRA cache is **three files, one accession, GSE144621, *E. coli*** (U00096.2)
  [grounded, re-measured]. `wiki/data_inventory.md` contains **zero** matches for
  `human|eukaryot|hg38|hg19|lentiMPRA|sharpr|starr|GRCh`. So the absence here is confirmed.
- **The data exists publicly.** Agarwal et al., *Nature* 2025 lentiMPRA measured **>680,000 sequences**
  across HepG2 / K562 / WTC11, with ~60,000 elements tested in all three cell types; deposited on
  Zenodo [grounded, web]. **MPRAbase** aggregates ~130 MPRA experiments, ~17.7M tested elements, 35
  cell types, 4 organisms [grounded, web].
- So the wall is a **download plus a run**, not an acquisition decision. Re-label it
  `runnable: needs download`, which turns it from an authority question into scheduling.
- **It does not help Arm 1** — lentiMPRA elements are ~200–230 bp.
- It is the arm that tests the north star's step 4 ("climb to animals"), so it is the strategically
  interesting one and now also a reachable one.

## "Why not run the whole thing on server1?"

**Because server1 buys compute, and Arm 1 is not compute-bound.** Its input does not exist. A faster
GPU processing 499 bp sequences still receives 499 bp sequences; so does a cluster. The constraint is
which sequences were ever measured.

**But the question is right for everything else**, and that is the useful half — the answer differs per
arm, which is precisely what the single "blocked" label hid:

| work | bound by | server1 the right home? |
|---|---|---|
| 4-decade pyramid (literal) | assay physics | **No** — nothing closes it |
| 4-decade pyramid **re-scoped to 2.22** | compute (321k rows) | **Yes** |
| eukaryotic panel | download + GPU training on ~680k seqs | **Yes** |
| `promoteratlas` (the ready arm) | weights fetch + GPU | **Yes** |
| fragment size-lever at full N | compute | **Yes** |

`Server1 [82de27]` is listed **idle and reachable** over Remote Control. A capability request was sent
(GPU, free disk, internet reach, shell vs HTTP-only, existing genomics data). **Remote Control does not
confirm reads and nothing reports back**, so as of this memo there is **no reply and no verified
capability claim** — the row above says where server1 *would* help, not what it has been measured to
have. `memory/reference_server1_gpu_box.md` records a GTX 1070 reachable by Ollama HTTP with no shell,
which if still true would make it a *runner* only after a shell is confirmed.

## What this changes, and what it does not

- `blocked:untried-arm` should be read as **one actionable arm (`promoteratlas`), one re-scopable arm,
  and one download** — not three equivalent walls.
- Both arms stay `planned`. The ledger is byte-unchanged. The kill/re-scope/leave decision is the
  user's, and a re-scope is itself a declaration that should go through a pre-registered arm.

**Limits.** The oligo-synthesis cap is labelled `[inferred]`, not sourced. "No such dataset exists" is
a claim about every reporter-assay design this search surfaced, not a proof of non-existence — a
long-fragment or whole-locus method at scale would refute it, and that is the revival condition. The
lentiMPRA figures are from the paper's abstract and search results, **not** from a fetched file
manifest, so size-on-disk and licence terms are unverified. Nothing here measures whether the re-scoped
arm will beat 0.4404; the prior for a frozen-embedding arm on this rung is losing to promcalc.

## Sources

- [Challenges and considerations for reproducibility of STARR-seq assays](https://genome.cshlp.org/content/33/4/479.full)
- [Massively parallel characterization of transcriptional regulatory elements in three diverse human cell types](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10028905/)
- [Agarwal et al., Nature 2025](https://shendure-web.gs.washington.edu/documents/Agarwal_Nature_2025.pdf)
- [Zenodo deposit](https://zenodo.org/records/10558183)
- [MPRAbase](https://cmg.psu.edu/publication/2025-04-genome_res)
