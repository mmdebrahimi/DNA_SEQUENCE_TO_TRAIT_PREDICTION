# `j3-abo-embedding-v1` — artifact contract (workhorse ⇄ laptop)

**Purpose.** Fix the exact shape of the two embedding artifacts the **workhorse (GPU)** produces so the
**laptop** can run the J3 downstream falsifier (JEPA `learned_rep` vs frozen-FM vs the deterministic
LD/mechanism baseline) with **no index drift**. Staged 2026-07-07 by the laptop session; the embeddings
themselves are GPU-gated (a genomic foundation model does not run on the laptop's GTX 860M — this is why the
J-track lives on the workhorse).

Locus: **ABO** (chr9q34.2). Trait: **ABO blood group**, self-reported (near-independent, non-circular,
~15% noisy — the recurring label-noise + non-deletional-O caveat).

---

## 1. The substrate (INPUT — already produced, laptop, CPU-only)

`scripts/abo_opensnp_ingest.py` → `data/j3_abo/j3_abo_substrate.{json,tsv}` (schema `j3-abo-substrate-v1`).
Source: the verified local OpenSNP 2017-12-08 dump (`D:/dna_decode_cache/opensnp/…zip`, 21.4 GB, 3570
entries, valid). **N = 395 samples** (users with an ABO self-report AND `rs8176719`); 563 ABO-labelled users,
168 lacked `rs8176719` on their array.

Per-sample columns (the TSV is the canonical **sample index** — row order = embedding row order):

| column | meaning |
|---|---|
| `user_id` | OpenSNP public user id (the join key; **stable, ordered — DO NOT re-sort**) |
| `raw_label` / `abo_group` / `label_o_vs_nonO` | self-report → {O,A,B,AB} → {O, non_O} |
| `rs8176719` | O-deletion (DD=O / DI,II=non-O) — 23andMe D/I calls |
| `rs8176746`, `rs8176747` | A/B-distinguishing tag SNPs |
| `deterministic_o_call` | rs8176719 DD→O else non_O — **the mechanism baseline** |
| `det_vs_label_concordant` | deterministic call == self-report O-vs-non-O |
| `genotype_member` | the source genotype file inside the zip |

**Deterministic O-vs-non-O baseline: concordance 0.924 (363/393).** This is the number the learned rep must
BEAT on a **de-confounded** slice to claim an embedding niche — consistent with the project's standing
finding that embeddings capture population structure, not mechanism (do NOT re-litigate that here; J3 is the
fair test, not a foregone conclusion).

## 2. The two embedding artifacts (OUTPUT — workhorse, GPU)

Both are produced by embedding a fixed genomic **window around ABO** for **each substrate sample in TSV row
order**. Emit as `.npz` (portable, laptop-loadable with numpy only):

### `frozen_fm` → `data/j3_abo/j3_abo_frozen_fm.npz`
- `emb`: `float32[N, D_fm]` — a FROZEN genomic-FM embedding (no fine-tuning), mean-pooled over the ABO
  window, one row per substrate sample, **row i ↔ substrate TSV row i**.
- `user_id`: `str[N]` — MUST equal the substrate `user_id` column, same order (the alignment assertion).
- `meta`: json string — `{model, revision, window_chrom, window_start, window_end, assembly, pool, D_fm}`.

### `learned_rep` → `data/j3_abo/j3_abo_learned_rep.npz`
- `emb`: `float32[N, D_rep]` — the **JEPA learned representation** (J2 encoder) over the SAME window, SAME
  row order.
- `user_id`: `str[N]` — identical index to `frozen_fm` + substrate.
- `meta`: `{encoder_ckpt, jepa_run_id, window_*, assembly, D_rep}`.

**Hard invariants (the laptop falsifier asserts these before scoring):**
1. `frozen_fm.user_id == learned_rep.user_id == substrate.user_id` (elementwise, same order). Any mismatch →
   the falsifier HARD-FAILS (no silent re-align — index drift is the classic embedding-eval footgun).
2. `emb.shape[0] == N == 395`; no NaN/Inf rows (a half-flushed embedding must not look like a valid mean).
3. `assembly` recorded; the window is the SAME for both arms (same coordinates, same pooling).

## 3. The falsifier (laptop, CPU — once §2 lands)

Mirrors the eukaryotic embedding tests (the fair-test protocol, `wiki/embedding_niche_cross_domain_synthesis`):
small head (logistic/ridge) on each rep, **de-confounded** (within-ancestry / within-lineage slice; ABO
freq varies by ancestry, so an ancestry-blind AUROC is confounded), compared to (a) the deterministic
`deterministic_o_call` baseline (0.924) and (b) a plain LD/genotype baseline on the 3 ABO SNPs.
**Verdict:** the rep has a niche ONLY if `within-group` metric > the mechanism baseline AND > structure-only.
Reuse `dna_decode/eval` + the `clade`/`within-group` machinery; no new confound-handling invented here.

## 4. Transfer

Git is the cross-machine channel. The **substrate TSV** (`data/j3_abo/`, small, public open-consent OpenSNP
data) is force-committed so both machines share the EXACT sample index. The 21 GB zip is NOT git-transferred
(the workhorse either has its own copy at the same `D:` path or the user transfers it out-of-band). The
`.npz` embeddings are workhorse outputs → committed back (or synced) for the laptop falsifier.

## 5. Status

- ✅ OpenSNP zip: real, local, verified (21.4 GB, valid, 3570 entries).
- ✅ Substrate `j3-abo-substrate-v1`: 395 samples staged + verified; deterministic baseline 0.924.
- ✅ Contract (this doc) + parser tests (`tests/test_abo_opensnp_ingest.py`, 5 offline).
- ⏳ `frozen_fm` + `learned_rep` `.npz`: **workhorse GPU** (this contract is their spec).
- ⏳ J3 falsifier: laptop, once the two `.npz` land.
