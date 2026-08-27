# Deep dive: where the label constraint actually binds (2026-08-27)

**Task:** a full deep dive on "label acquisition — the binding constraint".
**Result: the premise needs splitting, and one of the eleven walled cells is not walled at all.**

Everything below is verified live against this repo and the network today, not recalled. Where I did not
verify something, it says so.

---

## 1. "Labels are THE binding constraint" is REGIME-SPECIFIC, and the newer regime is already free

That framing is from the reproducibility freeze (2026-06-13) and it is true **for natural-population,
organism-level phenotype**. It does **not** hold for the regime this project's learned models actually win
in — and the repo's own later evidence says so:

| substrate | free? | confound | result |
|---|---|---|---|
| Bloom-2013 yeast segregant cross (1,008 segregants) | **yes** | removed BY CONSTRUCTION (one cross) | **12/12 traits decode**, r 0.46–0.80 vs null p95 0.03–0.12 |
| BXD mouse recombinant-inbred panel (198 strains) | **yes** | removed BY CONSTRUCTION | generalises fungi→mammal; brain weight r=0.574, 6/12 beat null |
| ProteinGym / MaveDB DMS | **yes** | designed libraries | ESM2+ProSST+GEMME beats best-single on 84–90% of proteins |
| GBM vs ridge on Bloom | — | — | nonlinear **beats** linear 7/12 (maltose 0.728→0.889) — real epistasis, not inflation |

So for **constructed variation**, labels are *not* binding — free substrates exist, are in hand, and work.
`wiki/design_epoch_plan_2026-08-07.md` already draws this line ("natural variation (confounded) vs
constructed variation (confound-free by construction)") and remains a **DRAFT for ratification**. Ratifying
or rejecting it is a live user decision that this deep dive did not resolve.

**Where labels genuinely still bind:** validating the *deployed deterministic decoder* against independent
measured phenotype — the report card's own walled cells, below.

## 2. The wall, quantified: 11 cells — and they are two different things

From `wiki/decoder_validation_report_card.json` (27 cells; 10 SCORED / 3 UNDERPOWERED / 2
ABSTAINS_BY_DESIGN / 1 LABEL_CONFOUNDED / **11 NO_FREE_PHENOTYPE_SOURCE** / 0 NOT_CENSUSED):

| group | cells | status of the `no_free_source` claim |
|---|---|---|
| **C. auris** × fluconazole, voriconazole, micafungin, caspofungin | 4 | **FACTUALLY WRONG — see §3** |
| **P. falciparum** × artemisinin, artesunate, chloroquine, dihydroartemisinin | 4 | right in substance, wrong (unrecorded) reason — see §4 |
| **Influenza A** × oseltamivir, peramivir, zanamivir | 3 | **UNVERIFIED** — never attempted; see §5 |

## 3. C. auris is NOT label-walled — a free source already produced a POWERED result

`dna_decode/data/shipped_decoder_surface.py:45-46` declares `no_free_source` for C. auris. The repo
contains the refutation:

- `wiki/ar_bank_caur_powered_result_2026-07-20.md` — **POWERED**, CDC AR Isolate Bank measured MICs:
  **12 isolates (5R/7S), sens 1.00 / spec 0.714 / acc 0.833**; the mechanism-attributable HIGH-confidence
  subset (n=9) is **1.00 / 1.00**.
- `data/raw/ar_bank_caur_extval_*/selected_strict.tsv` — real per-BioSample R/S labels on disk.
- The validation JSON records `label_source: ar_bank_MIC_cdc_tentative_breakpoint`,
  `independence_tier: CDC AR Isolate Bank measured MIC; provenance-disjoint (0 overlap vs the fungal G1
  tuning cohort)`.

**This is NOT lazy staleness, and the artifact says why:** it carries `not_in_shipped_surface: true` and
`rule_status: CURATED_NONFROZEN` / `rule_scope: scorer_local`. The validation was earned by a NON-FROZEN
scorer-local rule, deliberately kept outside the shipped surface — the same discipline as the TMP-SMX
experimental overlay.

**The defect is that the trust surface conflates two different claims:**

> "no free isolate-level phenotype source **exists**"  ≠  "no **shipped-surface** validation exists"

The first is false for C. auris. The second is true. The report card renders the first.
Under-claiming is as much a trust-surface falsehood as over-claiming — the project's own standing rule.

### Why I did not just fix it — a hard constraint, verified

`dna_decode/data/shipped_decoder_surface.py` is **sha256-pinned by the prospective lock**
(`wiki/prospective_lock_manifest_2026-06-22.json` → `surface_sha256`, alongside `amr_rules.py`,
`calibrated_amr_rules.json`, `mic_tiers.py`, `cohort_manifest.py`). `verify_lock` returns
**`ok=True, drifted=[]`** right now. Editing that file breaks the lock and invalidates the
leakage-free-by-construction guarantee the 2026-08-24 prospective accrual rests on. Same wall as the
gentamicin `rmt` catalog gap: **a fix is an unfrozen revision needing its own validation and a NEW lock —
a user authority decision, not a technical one.**

**The non-frozen path exists and needs no authority:** an *experimental-validation disclosure layer* in
`scripts/build_validation_report_card.py`, exactly mirroring the prospective-lock and lineage-disclosure
layers already there — it AUGMENTS a cell (never changes its state, never touches the pinned file) with
"a free measured source exists; validated under a non-frozen scorer-local rule". That is buildable,
reversible, and would stop the surface under-stating itself. **Not built in this run** — named as the
next move.

## 4. P. falciparum: the free source is real, rich, live — and trips G1

MalariaGEN Pf7 is open and reachable **right now** (Sanger FTP `226`, listing verified):
`Pf7_crt_haplotypes.txt`, `Pf7_csp_c_terminal_haplotypes.txt`, `Pf7_genetic_distance_matrix.npy`,
`Pf7_vcf/`, `Pf7.zarr.zip` (591 GB), `Pfalciparum.genome.fasta`, GFF.

That is a **GENOTYPE** resource. Its resistance classifications are **marker-derived** — produced by the
same class of genomic rule the decoder implements — so they trip **G1 circular-label**, the project's own
first gate. `no_free_source` is therefore correct *in substance* for malaria, but the recorded reason is
absent, and "no source exists" is the wrong description of "a large free source exists whose labels are
circular". WWARN (clinical efficacy) is reachable and is the non-circular candidate; **its per-isolate
paired genome+phenotype availability is NOT verified here.**

## 5. Influenza: the claim is unverified, not established

No `data/raw` dir, no wiki artifact, no script — the three influenza cells were **never attempted**.
Reachability today: BV-BRC `200`, NCBI Influenza Virus Resource `200`, GISAID `200` (registration).
**I verified reachability only — NOT whether any of them carries per-isolate NA-inhibition IC50 paired to a
genome.** Calling this "no free source" is an assertion; calling it "verified absent" would be false.

## 6. The transient-block pattern held again

`wiki/label_wall_data_sources_2026-08-03.md` recorded **DEG 15 (tubic)** as "timed out from this host".
Probed today: **`200`**, and its download endpoint also `200`. That file's own 2026-08-07 correction had
already caught the same thing for the P. aeruginosa SI ("the bot-block was transient"). OGEE v3 is
genuinely down (`000`). **Two of three "blocked" verdicts in that doc were transient** — consistent with
the standing rule: *try the public fetch before declaring a data wall.*

## What I checked, and what I did not

- **Verified:** the report-card split; the C. auris artifacts + their `not_in_shipped_surface` flag; the
  lock pinning + live `verify_lock`; Pf7 FTP contents; reachability of 11 endpoints.
- **NOT verified:** whether WWARN or any influenza source carries per-isolate measured phenotype paired to
  a genome (reachability only); whether the AR Bank holds enough *additional* C. auris isolates to power
  the other three antifungal cells.

## Recommendation

1. **Build the experimental-validation disclosure layer** (no authority, no frozen-file edit) so the
   surface stops under-stating C. auris. This is the concrete next move.
2. **Decide the design-epoch draft** (`design_epoch_plan_2026-08-07.md`) — it has been awaiting
   ratification for three weeks and it governs whether the constructed-variation regime becomes the
   project's centre of gravity. **User authority.**
3. **Only then** consider biobank acquisition (`label_acquisition_anchor_2026-07-04.md`, UK Biobank path).
   It remains the lever for *human-cell* external validation — but it is no longer accurate to call it the
   single unlock, because the constructed-variation regime is already open and free.
