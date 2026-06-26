# Tet functional-alphabet within-lineage — PARTIAL result (2026-06-26)

## Why partial
The external D: drive (Seagate Portable) disconnected mid-session, taking the genome FASTAs with it. The
probe's two arms have different inputs:
- **Functional-determinant arm** needs only the AMRFinder TSVs -> **all 118 are on C: (`data/amrfinder_runs/`)** -> RECOVERABLE offline.
- **k-mer comparator arm** needs the genome FASTAs -> on the disconnected D: -> BLOCKED.

So the head-to-head BEATS_KMER / TIES / FAILS verdict can't be computed until D: is reconnected. But the
functional half — the load-bearing one for the headroom question — was recovered.

## Result (functional arm only, tet, N=118, 20 shared lineages, 174 within-lineage pairs)
- **Functional-determinant within-lineage concordance = 0.963** (p < 0.0001 vs the in-MLST label-permutation null; null mean 0.502).
- Admitted 118 (58R / 60S), 0 dropped, 104 functional tokens (tetA/tetB/... efflux + tetM/tetO ribosomal-protection + acrB/acrR efflux/regulatory alleles).

## The honest read — this CONTRADICTS the expected outcome
The pre-committed expectation (user + Soraya) was that tet would **fail/tie** within-lineage because the
curated determinants are *incomplete* for mobile-element tet resistance — and that failure would be the
first concrete **headroom** signal for a learned model (#3). The data says the opposite:

- On tet, the **curated determinant alphabet DOES separate R/S within a lineage** (0.963), nearly as cleanly
  as cipro's QRDR (1.000). The "tet determinants are incomplete" hypothesis is **not supported** for the
  strains where AMRFinder makes calls — tetA/tetB/tetM presence tracks within-lineage tet resistance.
- This points **AWAY from DNA-LLM headroom**, consistent with the cipro finding: where the deterministic
  determinant alphabet already separates within-lineage, a learned representation has little room to add value.

## Caveat (do not over-read)
This is the functional arm ALONE. The full headroom verdict needs the k-mer comparison: if k-mer ALSO reaches
~0.96 within-lineage, neither alphabet has an edge; if functional >> k-mer, functional wins but it IS the
curated determinants = the shipped deterministic decoder (no learned-model headroom either way). Both
readings point to "no DNA-LLM headroom on tet." A BEATS_KMER result would NOT change that — it would just
re-confirm the determinant decoder.

## To finish (when D: reconnects — fully offline, no internet)
```bash
export UV_CACHE_DIR=C:/Users/Farshad/AppData/Local/uv_cache_c   # see env note below
uv run python scripts/functional_alphabet_probe.py \
  --cohort data/processed/shared_lineage_tetracycline_cohort.parquet \
  --drug tetracycline --refseq-cache D:/dna_decode_cache/refseq \
  --out wiki/functional_alphabet_probe_tet_n118_<date>
```

## ENV NOTE (separate, important): uv cache broke when D: disconnected
`uv`'s default cache (`C:/Users/Farshad/AppData/Local/uv/cache`) is a reparse-point into the now-dead D: ->
every `uv` command fails with `os error 183 (ALREADY_EXISTS)` until D: returns. Workaround used this session:
`export UV_CACHE_DIR=C:/Users/Farshad/AppData/Local/uv_cache_c`. Durable fix: set `UV_CACHE_DIR` to a C:
path in your environment, or reconnect D:.
