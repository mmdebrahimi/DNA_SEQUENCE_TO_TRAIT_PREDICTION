# PGx precision-leak audit — the non-core → `*1` silent-mis-call exposure (2026-07-28)

**Read-only verification (R4).** Turns the sentinel-layer handoff's *"no sentinel layer (v0)"* caveat into a
**number** DNA-11 can prioritize against. Computed OFFLINE from the committed GeT-RM truth set + the caller's
own core allele sets — no 1000G VCF, no network, no Docker, and **no pgx code touched**.
Script: `scripts/pgx_precision_leak_audit.py` · data: `wiki/pgx_precision_leak_audit_2026-07-28.json`.

> **UPDATE 2026-07-28 (post-population):** TPMT (10 sentinels, 6 withheld) **and CYP2B6** (3 distinctive-SNP
> sentinels, **18 withheld**, core 62/62 held — no false-withhold despite the `*6`-shared 785 SNP) are now
> populated (PharmCAT-sourced, Ensembl-verified). **CYP2C8 too** (4 PharmVar-sourced sentinels → 5 withheld,
> silent mis-calls 0/87, core 82/82 held). The live leak is now **0 samples** (CYP3A5 has 0 non-core carriers
> in its n=9 overlap → nothing to leak; its non-core alleles position-collide with core + are unvalidatable
> here → safely deferred).
> re-run `scripts/pgx_precision_leak_audit.py` for the current numbers (the `.json` auto-updates as genes are
> populated). The 36-across-4 figures below are the ORIGINAL pre-population audit that motivated the work.

## Headline

On the committed GeT-RM 1000G-overlap truth set, **36 real samples carry a non-core star-allele that the
v0 core-SNP proxy silently mis-calls `*1`** across the 4 quantifiable **leak** genes (SENTINELS = []) —
**10.1% of the 357 scored samples**. This is a **lower bound** (the ~87-sample-per-gene GeT-RM overlap is
small; non-core alleles are rarer than in a real clinical panel).

## Per-gene exposure

| gene | guard | non-core / scored | rate | non-core alleles observed | reading |
|---|---|---|---|---|---|
| **CYP2B6** | ❌ LEAK | **16 / 114** | **14.0%** | *18×7, *7×6, *2×4 | highest exposure — the single-SNP *6-proxy is blindest |
| **TPMT** | ❌ LEAK | **15 / 147** | **10.2%** | *8×4, *2, *16, *46 ×2 each; *40/*24/*32/*21/*12/*6/*33 | highest clinical stakes (thiopurine) |
| **CYP2C8** | ❌ LEAK | **5 / 87** | 5.7% | *17×2, *18, *15, *16 | moderate |
| **CYP3A5** | ❌ LEAK | 0 / 9 | 0.0% | — | **n=9 underpowered** — no non-core in this tiny overlap ≠ no leak |
| CYP2C19 | ✅ guarded | 8 / 87 | 9.2% | *8×2, *39×2, *4, *35, *13, *15 | **control** — these are correctly **WITHHELD** by the sentinel |
| CYP2C9 | ✅ guarded | 13 / 87 | 14.9% | *9×4, *8×3, *5/*6 ×2; *61, *11 | **control** — correctly WITHHELD |
| CYP2D6 | — SNP-surface | 7 / 47 | 13.0% | *21×2, *40×2, *14, *15, *46 | separate case; structural alleles (n excluded) = CRAM stack |

**Guarded control:** 21 non-core samples across CYP2C19 (8) + CYP2C9 (13) are **correctly withheld** — the
sentinel layer doing its job. That is exactly the behavior the 4 leak genes lack.

## Faithfulness proof (verify-in-batch)

My tier classification is faithful to the concordance script's own logic: my **CYP2D6 count = 7 non-core /
47 core** exactly matches DNA-11's committed `pgx_getrm_concordance_cyp2d6_2026-07-06.json`
(`noncore_snp_n=7, core_snp_n=47`). The audit imports the per-gene `core` sets + `ref_equiv` aliases from
`scripts.pgx_getrm_concordance.GENES` (single source of truth) rather than re-hardcoding them; the
`ref_equiv` fix (CYP2C19 `*38 == *1` → core, not non-core) was caught by this cross-check and applied — it
does not affect the leak headline (all 4 leak genes have `ref_equiv = {}`).

## What DNA-11 should do with this

- **Order the sentinel propagation by measured exposure × stakes:** **TPMT** (15 samples, thiopurine) and
  **CYP2B6** (16 samples, highest rate) first; then **CYP2C8** (5); **CYP3A5** last of the harnessed set
  (0/9 but underpowered — don't skip, just deprioritize).
- The non-core allele lists above ARE the guard target lists (source their defining rsIDs verbatim from
  PharmVar per the handoff's anti-fabrication rail): TPMT *8/*2/*16/*46/…, CYP2B6 *18/*7/*2, CYP2C8 *17/*18/*15/*16.

## Caveats (load-bearing)

- **Lower bound.** The GeT-RM 1000G overlap is ~87 samples/gene; a real clinical panel has more non-core
  carriers. 36 is a floor, not the field rate.
- **NUDT15 / UGT1A1 / DPYD are NOT in the GeT-RM concordance harness at all** (no committed truth column) —
  their leak is real but **un-sizeable offline**; size it after a v0.1 truth fetch (the handoff's v0.1 item).
- **CYP2D6 structural** alleles (*5/*13/*36/…) are excluded from this SNP count — a different surface,
  already handled by the CRAM structural stack.

## Provenance

Offline from `tests/data/pgx_getrm/star-allele-comparison_common.tsv` + `getrm_{cyp3a5,tpmt,cyp2b6}_consensus.tsv`
(committed) + `GENES` core sets (`scripts/pgx_getrm_concordance.py`), read 2026-07-28. No pgx code modified;
no VCF/network/Docker. Companion to `wiki/pgx_sentinel_layer_handoff_2026-07-28.md`.
