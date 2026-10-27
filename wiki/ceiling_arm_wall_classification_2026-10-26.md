# The three untried ceiling arms — which wall each one actually sits behind (2026-10-26)

**Why this exists.** `soraya_runs/arms-ledger.json` holds three `planned` arms on the tile-ceiling
criterion. All three carry an **empty `evidence` string and no `rationale` key**, so nothing in the
repo distinguished *"expensive but runnable today"* from *"needs data that is on no attached drive."*
`blocking()` therefore reported `blocked:untried-arm` naming three arms with no way to act on any of
them, and two consecutive `--advance` runs read that as an undifferentiated wall.

**Nothing was written into the ledger.** `evidence` is reserved for EXECUTED evidence — the ledger's
own guard says so (`tests/test_arms_ledger_stack_arm.py:73`: `assert arm["evidence"] == ""`, *"a
planned arm has no evidence yet"*), and that is the R6/R7 discipline, not an accident. The three arms
also carry no `rationale` key, so recording a wall there would reshape the schema of a ledger with
**four currently-red rails that TODOS marks `[OPEN — USER AUTHORITY]`**. So this is a memo.

## What was searched

A substrate claim is a claim about the search, so here is the search:

- `D:/dna_decode_cache/` — present, contains exactly `amr_portal` + `mpra`
- `D:/dna_decode_cache/mpra/` — **exactly three files**, all one accession:
  `GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz`,
  `GSE144621_..._M9_expression.txt.gz`, `GSE144621_peak_tile_expression_formatted_std.txt.gz`
- `dna_decode/data/external_root.probe()` — two candidates, both `D:/dna_decode_cache` and
  `D:/PythonProjects/dna_decode_cache`, both with entries `('amr_portal', 'mpra')` and **no
  `refseq`**, so `CACHE_REQUIRED_SUBDIRS = ("refseq",)` refuses both
- `D:/hf_cache/hub` — **one model only**: `models--GenerTeam--GENERator-v2-prokaryote-1.2b-base`
- `data/promcalc/tile_predictions/` — **10 seed tables present locally** (`manifest.json` +
  `seed-00..09.json.gz`). Note the tile substrate is LOCAL, not on D:; `DEFAULT_OUT` points at
  `D:\dna_decode_cache\promcalc\tile_predictions`, which is **empty** — a runner without
  `--in-dir data/promcalc/tile_predictions` refuses rather than producing an empty arm.
- census span fields, per rung (`wiki/mpra_substrate_census_2026-10-12.json`)

## The three classifications

### 1. `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` — **CODE-CLOSABLE TODAY**

Cost tier `medium`, declared 2026-10-13. The only one of the three that is reachable now.

- **Substrate: PRESENT.** The arm scores tiles, and all 10 seed tables are on local disk.
- **Missing: the model weights only.** `D:/hf_cache/hub` holds exactly one model and it is
  GENERator-v2, not PromoterAtlas. So the gap is a public-weights fetch plus disk, not data
  acquisition. C: is 85–99% full, so the fetch belongs on `HF_HOME=D:/hf_cache` (already set).
- **It is also the arm that answers the user's standing position.** `[[learned-model-is-the-product]]`
  has been restated repeatedly; this is a LEARNED representation on a CONSTRUCTED panel, which
  `dna_decode/eval/regime.py` puts in the regime that WORKS — not the closed natural-population
  zero-shot negative. It is therefore the cheapest ambitious arm AND the on-position one.
- **Not started here.** Fetching weights and running it is the next move, not a record.

### 2. `multiresolution-pyramid-4-decades-UNRUN` — **EXTERNAL WALL, AND THE CENSUS ALREADY SAID SO**

Cost tier `expensive`, declared 2026-10-09. **`decades` means SEQUENCE LENGTH, not expression
dynamic range** — and the derivation is now pinned by arithmetic rather than inferred:

| quantity | value | check |
|---|---|---|
| `rungs.fragments_LB.max_span_decades_at_3bp_floor` | 2.22 | `log10(499/3) = 2.22` ✓ |
| the measured sibling `pyramid-4scale-1.70decades-**150bp**` | 1.70 | `log10(150/3) = 1.70` ✓ |

So the field is `log10(max_sequence_length / 3)`, both values reproduce exactly, and the sibling arm
reached 1.70 decades **because peak tiles are 150 bp**. A 4-decade arm needs ~10,000 bp (from a 1 bp
reference) or 30,000 bp (from the stated 3 bp floor). **The longest sequence in any cached rung is
499 bp.** So the arm is not expensive — it is **impossible on this substrate**, by a factor of ~20x.

**The census states this outright, in a field I did not read until after writing the first version
of this memo**, under `the_forgotten_substrate.what_remains_true`:

> the LITERAL 4-decade arm needs ~10,000 bp; fragments top out at 499 bp

**My own first classification of this arm was wrong** — I called it an unmeasured wall needing a
cheap span computation, having searched for the computing CODE (there is none; the field was derived
by hand) and concluded the definition was unavailable. The definition was recoverable two ways I
skipped: the arithmetic above, and a sibling key in the same file. The irony is exact and worth
keeping, because that same census field carries the generalisation
**"a substrate blocker is a claim about what you looked for"** — the rail I invoked while not having
looked at the field that answered the question. Recorded, not smoothed over.

**Consequence:** this arm should not sit in `planned` implying runnability. Whether to kill it
(R6: only on its own executed evidence — and "impossible on this substrate" is a substrate fact, not
executed evidence of the arm) or re-scope it to the reachable 2.22 decades is a **USER AUTHORITY**
call, because killing or re-scoping a declared arm changes what `blocking()` reports.

### 3. `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` — **EXTERNAL WALL, MEASURED**

Cost tier `expensive`, declared 2026-10-09.

- The entire MPRA cache is **one accession, GSE144621, which is *E. coli*** (U00096.2). There is no
  eukaryotic or human panel on any attached drive, and `external_root` refuses both cache candidates
  for want of `refseq`.
- The wall is **data acquisition** — a new panel has to be fetched or licensed. No amount of local
  work closes it. This is the one of the three that genuinely deserves `expensive` + deferred.
- It is also the arm that would test `[[substrate-is-not-the-subject]]` (the north star's step 4,
  "climb to animals"), so it is strategically the most interesting and the least reachable. That
  tension is the honest state, not a contradiction to resolve here.

## What this changes

`blocked:untried-arm` named three arms as if they were one wall. They are not one wall, and they are
not three equivalent ones:

| arm | cost label | real wall | next action |
|---|---|---|---|
| `promoteratlas-base-embeddings-frozen-150bp` | medium | **public weights fetch** — substrate present (10 local seed tables) | runnable today |
| `multiresolution-pyramid-4-decades` | expensive | **external data** — needs ~10,000 bp; cached max is 499 bp (~20x short) | kill or re-scope: USER AUTHORITY |
| `substrate-change-constructed-eukaryotic-or-human-panel` | expensive | **external data** — whole cache is one *E. coli* accession | acquisition: USER AUTHORITY |

So **exactly one of the three is reachable today**, and two are blocked on data acquisition rather
than on effort. That matters for how `blocked:untried-arm` should be read: it is not "three things
nobody got round to", it is "one actionable arm plus two acquisition decisions".

The ambitious lane was not starved for lack of candidates — it was starved because **no candidate
carried a reason**, so every `--advance` re-derived the same undifferentiated wall and moved on.

**Limits.** Cost tiers are the author's free-form labels, not measurements. "Runnable today" bounds
the SUBSTRATE and the WEIGHTS only — it says nothing about whether the arm will beat 0.4404, and a
frozen-embedding arm on this rung has a measured prior of losing to promcalc. Nothing here was run.
