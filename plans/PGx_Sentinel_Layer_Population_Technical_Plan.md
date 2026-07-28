# Technical Plan — Populate the PGx SENTINELS lists (non-core → withhold, not silent `*1`)

**Date:** 2026-07-28 · **Author:** Soraya (new-phenotype-cells session) · **Class:** (d) cross-cutting — touches shared pgx contracts (though additive per-cell).
**Companions:** `wiki/pgx_sentinel_layer_handoff_2026-07-28.md` (the grounded finding) + `wiki/pgx_precision_leak_audit_2026-07-28.md` (the measured 36-sample exposure).
**Status:** DRAFT plan — NOT executed. Planning STOP in effect; awaits the execution-modality authority decision (§7) + a pre-exec `/brainstorm` before any `/save-plan`.

## ⚠️ Ownership decision precedes execution (R4) — read §7 FIRST

Populating `SENTINELS` writes **pgx code = DNA-11's owned lane** (I handed this off two turns ago *because* it's theirs). This plan is reversible/non-colliding; **executing it is not** without an ownership decision. §7 is a hard authority fork — do not start §5 until it's resolved.

## 1. Problem

The sentinel-withhold mechanism is generic + working in `dna_decode/pgx/caller.py`, but 7 star-allele catalogs ship `SENTINELS = []`. The measured leak: **36 real GeT-RM 1000G-overlap samples (10.1%) carry a non-core allele silently mis-called `*1`** (a lower bound). Populating each gene's `SENTINELS` with its non-core-allele defining sites makes those samples resolve to `phenotype_withheld` (correct abstention) instead of a confident wrong call.

## 2. Scope — validatability-bound v0 (R2 framing check; a reversal from the handoff's ordering)

The handoff said "TPMT/NUDT15/UGT1A1 first" (clinical stakes). The **audit refines this**: NUDT15/UGT1A1/DPYD are **NOT in the GeT-RM concordance harness** → a sentinel there is **un-validatable on real samples** (populating unverifiable guard data is exactly the fabrication risk). So the binding constraint for v0 is **validatability**, not clinical stakes.

- **v0 (this plan) = the 4 GeT-RM-harnessed leak genes:** **TPMT** (harnessed + top clinical stakes + 15 exposed) → **CYP2B6** (16 exposed, highest rate) → **CYP2C8** (5) → **CYP3A5** (0 in the tiny overlap — populate, but validation deferred; don't skip).
- **v0.1 (DEFERRED, blocked): NUDT15 / UGT1A1 / DPYD** — populate only after their GeT-RM truth is fetched (the handoff's v0.1 item), else the sentinels can't be proven to fire correctly.
- **OUT: CYP2D6** — SNP-surface non-core is a different caller path; structural alleles are the CRAM stack. Separate plan.

## 3. The load-bearing crux — accurate PharmVar GRCh38 coordinates (anti-fabrication)

Each `SentinelVariant(rsid, chrom, pos, ref, alt="*", implies, note)` must be a **real PharmVar allele-defining SNP at its correct GRCh38 coordinate** (the caller reads a GRCh38 VCF). There is **no fetch script** — the existing CYP2C9 sentinels were hand-curated (`rs28371686 10:94981301 C→* *5`). **A wrong coord = a silent no-op sentinel** (never fires → the leak persists invisibly) **or a false-withhold** (fires on the wrong site). This is the single highest-risk, un-fabricatable part.

**Mitigation (make the anti-fabrication rail EXECUTABLE):** a new `scripts/verify_sentinel_coords.py` that, for every proposed `SentinelVariant`, cross-checks `(rsid → GRCh38 chrom:pos, ref)` against **Ensembl REST** (`/variation/human/<rsid>`) and **fails closed** on any mismatch. This converts "hand-curated, fabrication-prone" into "sourced + machine-verified." (Network read — auto-classified; offline-degrades to SKIP with a loud `UNVERIFIED` stamp, never a silent pass.)

Per-gene chromosome (to be VERIFIED per-allele by the script, not asserted): TPMT chr6 · CYP2B6 chr19 · CYP2C8 chr10 · CYP3A5 chr7.

## 4. Per-gene guard targets (from the audit's observed non-core alleles)

Source each allele's defining SNP(s) verbatim from PharmVar; the audit's observed non-core alleles are the target list:
- **TPMT:** *8, *2, *16, *46, + *40/*24/*32/*21/*12/*6/*33 (observed); prioritize the no/low-function ones (CPIC).
- **CYP2B6:** *18, *7, *2.
- **CYP2C8:** *17, *18, *15, *16.
- **CYP3A5:** the known non-expressor non-core alleles (none observed in the n=9 overlap → source from PharmVar/CPIC, validation deferred).

## 5. Steps (waves)

**Wave 1 — sourcing discipline (shared, critical path):**
1. Build `scripts/verify_sentinel_coords.py` (Ensembl-REST rsID→GRCh38 verifier, fail-closed, offline-degrade-to-UNVERIFIED) + its unit test (mock the REST call).
2. For each v0 gene, source the candidate `SentinelVariant` rows from PharmVar allele definitions (rsID + GRCh38 coord + ref + implies + protein-change note). Run the verifier → every row must pass before use.

**Wave 2 — per-gene population (parallelizable across the 4 genes; each independent):**
3. Populate `dna_decode/pgx/<gene>_catalog.py::SENTINELS` with the verified rows (mirror the CYP2C9 `alt="*"` wildcard pattern).
4. Guard against a sentinel site coinciding with a `CORE_DEFINING` site (overlap check → error if so).
5. Add a per-gene withhold test (mirror `test_pgx_cyp2c19.py::test_sentinel_*_withholds`): synth a VCF with the sentinel ALT → assert `phenotype_status == "phenotype_withheld"`, `phenotype is None`, `sentinel_hits` implies the right `*N`.

**Wave 3 — validation + trust-surface (acceptance gate):**
6. Re-run `scripts/pgx_getrm_concordance.py` for each edited gene (real 1000G VCF — R3 real-surface; needs the cached VCFs / fetch). Assert **core concordance UNCHANGED** AND the previously-exposed non-core samples now resolve to **withheld** (not silent `*1`). Re-run `scripts/pgx_precision_leak_audit.py` → the leak count for v0 genes drops to 0 *guarded*.
7. Truth-up each edited cell's `cell_registry` `demotion_rule` ("mis-called *1 (no sentinel v0)" → the sentinel-guarded wording).
8. Full suite green + frozen AMR/forward surfaces byte-unchanged.

## 6. Critical path & verification (the MVP bar if this later runs `--until-mvp`)

Critical path: **Wave 1 (verifier + verified rows) → Wave 2 (per-gene, parallel) → Wave 3 (concordance gate)**.

Acceptance criteria (all checkable):
1. `verify_sentinel_coords.py` passes for every populated row (no `UNVERIFIED`, no mismatch).
2. Every v0 gene's withhold test passes (`test-exit-0`).
3. `pgx_getrm_concordance.py` core concordance **byte-identical** to the committed numbers for each edited gene (no regression) — the previously-exposed non-core samples now withheld.
4. `pgx_precision_leak_audit.py` v0-gene leak → guarded (exposure moves from LEAK to withheld).
5. registry `demotion_rule` truth-up per gene; frozen surfaces byte-unchanged; full suite green.

## 7. Execution-modality AUTHORITY FORK (decide before §5) — the one thing I won't self-resolve

This writes pgx code, DNA-11's lane. Three clean paths — **your call**:
- **(A) DNA-11 executes** (default per R4): this plan + the audit are the spec; hand it over. Zero collision; slowest.
- **(B) I execute on a non-colliding side branch** (`git worktree`, R4/C): I do Waves 1-3 on a side branch; DNA-11 / you cherry-pick to `main`. No `main` collision; you review before merge.
- **(C) You authorize me on `main`**: fastest; requires you to reassign the pgx lane for this task (DNA-11's board row is stale — last update 2026-07-07 — so live collision risk is low, but this is your call, not mine).

**My recommendation:** **(B)** — I have the full grounded context (audit + mechanism), a side branch keeps `main` clean and gives DNA-11/you a reviewable diff, and it doesn't require reassigning the lane. But (A) is the most conservative and (C) the fastest; the choice is an authority decision.

## 8. Risks (seed for the pre-exec `/brainstorm`)

- **Wrong GRCh38 coord → silent no-op sentinel** — the dominant risk; mitigated by the Wave-1 verifier (fail-closed).
- **Sentinel coincides with a core defining site** → double-count/wrong withhold — Wave-2 overlap check.
- **`alt="*"` wildcard fires on a benign ALT at the same coord** → over-conservative false-withhold (abstain > mis-call, but reduces yield) — accept + note; a specific `alt` narrows it where PharmVar gives one.
- **CYP3A5 0-exposed** → withhold un-validatable on this cohort — populate + mark validation-deferred.
- **Ownership collision** (§7) — the authority fork.

## 9. Not doing (scope discipline)

Full `decompose` → per-family `/project-init` is **overkill** for a ~4-gene additive change (would burn self-init slots for a task that is one plan + parallel per-gene units). This is planned as a single technical plan with parallelizable Wave-2 units, not a multi-family portfolio.

## Progress (2026-07-28 — execution begun after the pre-exec brainstorm)

The pre-exec `/brainstorm` reshaped the plan (4 grounded findings). Ownership resolved: this session IS
DNA-11 (user directive), so execution runs on `main`.

- ✅ **Wave 1.5 — shared withhold helper (commit 9831e15).** Extracted `apply_sentinel_withhold` +
  `_scan_sentinel_counts` (list-keyed) used by BOTH `call_diplotype` and `assemble_compound_diplotype`;
  generalized the CYP2C19 `*35` rule into `SentinelVariant.accounted_by_core`; the compound path (TPMT)
  now withholds. All 4 brainstorm findings fixed. 309 pgx/compound/trio/registry tests pass (CYP2C19/
  CYP2C9 byte-identical); 7 new tests. Compound path is a no-op while `SENTINELS=[]` (zero live change yet).
- ✅ **Coord verifier (commit ddc257f).** `scripts/verify_sentinel_coords.py` — Ensembl-GRCh38 fail-closed
  check; 11 offline tests + live-verified the 4 CYP2C9 sentinels OK. The anti-fabrication rail is now executable.
- 🔎 **Sourcing path CONFIRMED + de-risked (finding).** PharmVar's own API is now **401 key-gated** (requires
  a PharmVar account key). The SAME allele→rsID→GRCh38 definitions are **FREE** from **PharmCAT**
  (`raw.githubusercontent.com/PharmGKB/PharmCAT/.../alleles/<GENE>_translation.json`, source CLINPGX/CPIC,
  GRCh38) and the **CPIC API** (`api.cpicpgx.org/v1/allele_definition`). Both reachable this session.
- ✅ **TPMT POPULATED + validated (commit ab3d381).** 10 non-core sentinels sourced from PharmCAT
  `TPMT_translation.json` (GRCh38 chr6), all 10 Ensembl-verified (10/10 OK). Real-1000G GeT-RM concordance:
  **core 85/85 UNCHANGED + 6 non-core samples now WITHHELD** (were silent *1). Leak audit 36→21; TPMT
  guarded. 321 tests pass; registry truthed-up. THE FIRST REAL PRECISION GAIN.
- ✅ **CYP2B6 POPULATED + validated (commit 6670468).** 3 distinctive-SNP sentinels (`*2`/`*7`/`*18`),
  PharmCAT-sourced + Ensembl-verified, with the `*6`-haplotype-shared 785 SNP deliberately EXCLUDED (no
  false-withhold). Real-1000G GeT-RM: **core 62/62 UNCHANGED + 18 non-core WITHHELD** (silent 26→8). Leak
  21→5. Verify-in-batch caught 2 real issues: (a) sentinel sites outside the narrow *6-proxy VCF → re-fetched
  a wider chr19 region (40991000-41017000); (b) the GRCh37 lift map (`pgx_decode_pgp_uk.GRCH37_POS`) lacked
  the new sentinels AND `_all_variants()` never yielded `tp.SENTINELS` (TPMT sentinels silently un-lifted) →
  added 13 GRCh37 coords + wired `tp.SENTINELS`. `*4`/`*9` (absence-defined) = documented v0 gap.
- 🔎 **CYP2C8 sourcing wrinkle CONFIRMED:** not in PharmCAT (404) AND the CPIC API schema (`allele` /
  `sequence_location` joins) needs real exploration (my first-guess columns errored) → a focused sourcing
  sub-task (CPIC schema OR PharmVar download), not a quick continuation. 5 exposed samples.
- ⚠️ **Per-gene wrinkles found (the remaining 2 are NOT clean repeats — handle each carefully):**
  - **CYP2B6:** non-core `*7`/`*18` are COMPOUND with SNPs shared on the core `*6` haplotype (rs2279343/785
    rides on `*6`) → a naive sentinel there FALSE-WITHHOLDS a core `*6` call. Needs distinctive-SNP selection
    (`*7`→rs3211371, `*18`→rs28399499) + shared-SNP exclusion. Only `*2` (rs8192709) is a clean single-SNP.
  - **CYP2C8:** NOT at the PharmCAT allele-def path (HTTP 404) → source from the CPIC API
    (`api.cpicpgx.org/v1/allele_definition` + linked variant table) or a PharmVar download (different parser).
  - **CYP3A5:** 0 non-core in the n=9 GeT-RM 1000G overlap → populate-but-unvalidatable on this cohort
    (underpowered); lowest priority.
- ⏭ **REMAINING (CYP2B6 / CYP2C8 / CYP3A5).** Parse the source `<gene>` allele definitions
  format → each non-core allele's defining variant(s) → `SentinelVariant` rows → run `verify_sentinel_coords`
  on every row (must pass) → populate `<gene>_catalog.SENTINELS` → per-gene withhold test → re-run
  `pgx_getrm_concordance.py` (core concordance unchanged) + `pgx_precision_leak_audit.py` (leak → guarded)
  → registry `demotion_rule` truth-up. Order: TPMT → CYP2B6 → CYP2C8 → CYP3A5. This is a distinct,
  fabrication-sensitive ingestion pass — do it focused, not rushed; the infra above makes it SAFE (verifier
  gates every coord). NOT money/hardware gated (code-closable; PharmCAT is free).

## Pipeline next steps (Planning STOP — not executed)

Per the planning pipeline: **pre-exec `/brainstorm`** on this plan (class (d) cross-cutting → mandatory) → `/save-plan` → resolve §7 → `/execute-plan` (or hand to DNA-11). I stop here and await the §7 authority decision + your go-ahead.
