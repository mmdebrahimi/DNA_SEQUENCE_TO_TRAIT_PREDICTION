<!-- recommendation.md for 2026-06-02-1600-ep4-v0-cli-build -->
# Recommendation — next moves

## Stop: mvp-reached
The v0 compatibility-resolver CLI is shipped + tested. This advances ledger Goal 3 (v0 marker catalog + decision rules + parser) substantially: 23 clusters codified, 11-class table + abstention implemented, parser + provenance JSON working on real genomes.

## Ranked next actions (VOI)
1. **Full-cohort resolver eval (ledger Goal 5 / H4)** — run the CLI over all 24 cached genomes (+ the ETEC von Mentzer set when convenient), tabulate calls vs labels, and measure: confident-call precision, abstention rate, and abstention enrichment for known-ambiguous cases. This will quantify the ExPEC under-call (UPEC ≥2-strong rule) and tune thresholds. ~36 min compute (cache the per-genome k-mer set to speed reruns, like the NT arm did).
2. **Calibrate ExPEC sensitivity** — the conservative UPEC rule abstains on 1-strong-fimbrial ExPEC. Options: count iha/usp/vat/sat/hlyF as additional strong/support markers, or add a `UPEC_LOW_CONFIDENCE` tier at 1 strong + ≥2 support. Decide against the H4 abstention-rate target (≤15% on unambiguous).
3. **CGE VirulenceFinder side-by-side diff (ledger requirement)** — currently `caller_is_independent_baseline=false` and detection is k-mer-seed, not BLAST. For a true diff, run real VF (PyPI Virulencefinder 3.2.0) on a few genomes and diff gene-calls. Needs the VF software (a dep-install — gate it).
4. **Speed**: cache per-genome k-mer sets to disk so repeated CLI runs / the cohort eval don't rebuild (~90 s/genome today).

## EMIT (user-only `/project-state`)
```
/project-state ecoli-pathotype-prediction-cli-2026-05-26 --append-action \
  --class build "v0 compatibility-resolver CLI shipped: dna_decode/pathotype/ (markers 23 clusters + resolve 11-class table/abstention + detect k-mer-seed + cli FASTA->provenance JSON). 18/18 decision tests pass. Demoed: EPEC->tEPEC_COMPATIBLE (eae+bfpA provenance), ExPEC->AMBIGUOUS (on-spec conservative UPEC rule). Advances Goal 3. Next: full-cohort H4 abstention eval + ExPEC sensitivity calibration + real-VF side-by-side diff."
```

## Held / done this turn
- `git push` — user authorized push this turn; pushing the build commit too (consistent with "push … and then proceed"). Will report transparently.
