# `dna_decode.forward` — the forward variant-effect cell

The **forward** direction of the decoder: given a genetic **edit** (a point mutation), predict the change in
**molecular phenotype**, with a calibrated magnitude + honest uncertainty, across the tree of life. This is
the "edit → effect" complement to the frozen determinant→R/S AMR decoder (which this package does NOT touch;
`verify_lock` stays green).

Validated per-variant against measured **Deep Mutational Scanning (DMS)** — the one place the project's label
wall doesn't bind (DMS gives free, independent, per-variant magnitude labels).

## Quickstart

```python
from dna_decode.forward import predict_effect
# score one edit on one protein (BLOSUM62 = deterministic, no deps)
p = predict_effect(protein_seq, "M69L", protein="TEM-1", method="blosum62")
print(p.predicted_effect, p.raw_score, p.confidence)   # preserved/damaging/uncertain + continuous score
```

## The regime router (`predict_edit`) — the capstone

Auto-routes an edit to the RIGHT predictor by G2P regime (`feedback_g2p_decoder_regime_boundary`):

| regime | edit hits… | predictor | note |
|---|---|---|---|
| **A** determinant | a curated AMR determinant | the determinant catalogue → R/S | resistance NEVER goes to a likelihood model (antagonistic-selection blindness) |
| **B** molecular | a protein's fitness/stability | BLOSUM / ESM2 / AlphaMissense / ESM-IF | the DMS-validated variant-effect methods |
| **C** organismal | a polygenic organism trait | — | **ABSTAIN** (closed negative) |

## Methods (Regime B) — validated numbers

Per-protein Spearman(prediction, measured DMS). Learned methods beat deterministic BLOSUM everywhere they run.

| method | scope | how | PTEN | TEM-1 (E.coli) |
|---|---|---|---:|---:|
| `blosum62` | universal, instant, no deps | substitution severity | 0.182 | 0.346 |
| `esm2` | universal (bacteria+eukaryote) | ESM2-650M masked-marginal | 0.518 | **0.732** |
| `alphamissense` | **human only** | AM pathogenicity (1−AM) | **0.539** | — |
| `esm_if` | structure-based | inverse-folding conditional-LL (Kaggle T4) | 0.479 | — |

Tree-of-life: ESM2 lifts BLOSUM by +0.28–0.39 on E.coli / human / yeast / Arabidopsis. AlphaMissense is
human-proteome-only. **Structure (ESM-IF 0.479) does NOT beat sequence (ESM2 0.518) on PTEN** — the
"+structure" margin doesn't materialize here (consistent with ProteinGym; the structure tier is led by
ProSST/SaProt, not ESM-IF). Full leaderboard: `wiki/forward_method_leaderboard_*.md`.

## Genome-level edit (`predict_genome_edit`)

Lifts the input from a protein mutation to a real **nucleotide edit in a CDS**: base substitution → codon →
AA change → the Regime-B predictor, classified silent / missense / nonsense. REF-base coordinate gate +
translated-WT double-check fail loudly. Validated end-to-end on a real blaTEM CDS (Spearman 0.76 vs measured
DMS over 1,715 real single-nt-accessible variants; `wiki/blatem_genome_demo_*.md`).

## Dosage head (`evaluate_dosage`) — calibrated MAGNITUDE, not just rank

Turns a rank-score into a calibrated magnitude prediction: isotonic score→effect calibrator + **split-
conformal prediction interval** at a target coverage. Reuses J2's conformal definition (`conformal_q` ==
`hiv_quantitative_calibration._conformal_q`). Cross-organism sweep (10 proteins): **calibrated 10/10**
(coverage ≈ target everywhere — the conformal guarantee is organism-agnostic), **informative 7/10**.
Key finding: **interval narrowing (dosage-informativeness) is DISTINCT from rank quality** — a method can
rank moderately (CcdB-ESM2 Spearman 0.49) yet NOT narrow the magnitude interval. "Ranks well" ≠ "pins the
dose." Memo: `wiki/forward_dosage_sweep_*.md`.

## The INVERSE (`dna-decode inverse` / `dna_decode.forward.inverse`) — effect → edit

The forward cell answers *edit → effect*. The inverse answers *effect → edit*, which is what a design loop
wants — using this cell as **label-free ground truth** (no phenotype label is ever consulted, which is how it
dodges the project's label wall).

```python
from dna_decode.forward.inverse import propose_edits
r = propose_edits(protein_seq, target_percentile=0.05, top_k=5, cds=cds_seq)  # the 5% most damaging
```

**It RANKS. It does not dose.** The narrowness is measured, not modest:

| question | verdict | why |
|---|---|---|
| hit a target **effect** (dose) | **not deployable** | needs a score→effect calibrator fit on the **target protein's own DMS** — and if you have that you already know every effect. Calibrators cannot transfer: the assays share no scale (CcdB's whole range [−9.00,−2.00] sits below TEM-1's minimum −3.56 — impossible *by construction*). The conformal interval is informative **0/6** splits: it brackets while proving nothing (coverage holds even for a useless model). |
| hit a target **percentile** (rank) | **ships** | needs no calibrator, no DMS, no label. **ESM** beats an exact no-oracle null 4/4 hand-picked proteins (~2–5 pct pts at top-5). **But the shipped `blosum62` default does NOT generalize** — N=200 ProteinGym: material on only **13.5%**, often worse than guessing (`wiki/proteingym_inverse_sweep_2026-07-17.md`). So: cheap first pass by default; ESM (GPU) for reliability; gate per protein. |

Three rails that ship inside every call (`does_not_support`, `evidence`, `notes`):

- **`propose k, assay k, keep the best`** — top-1 is ~4× worse than best-of-5.
- **`blosum62` is often the right answer, not a fallback** — the learned oracle earns its keep on only
  **3/4** proteins, and **utility does not track forward rank** (PTEN 0.5185 earns it; RL40A 0.5190 does
  not) → check per protein.
- **one edit per residue by default** — BLOSUM62 takes only **7 distinct scores** over 1,874 real blaTEM
  candidates (largest tie group 383), so a plain window returns *k shots at the same residue*. Diversity was
  measured free for ESM and **better** for BLOSUM (up to −0.07 pct pts).

Harnesses: `scripts/forward_inverse_{roundtrip,sweep,deployable}.py` →
`wiki/forward_inverse_*_2026-07-1{6,7}.md`.

## Honest scope (hard-won rails)

- **Regime B molecular only.** This predicts enzyme fitness/stability — NOT organism-level polygenic traits
  (Regime C, a closed negative → abstain), and NOT the antagonistic **clinical-resistance** direction (raw
  likelihood/severity scorers fail there — use the Regime-A catalogue; the resistance-conservativeness
  finding).
- The validated quantity is the **rank correlation** (per protein) + the **dosage coverage/narrowing**; the
  per-variant tier is a coarse read of a continuous score.
- Every coordinate is checked: a WT-vs-reference or REF-base mismatch fails loudly, never a silent wrong call.

## CLI (`dna-decode forward` / `dna-forward`)

The forward cell is a first-class command in the shipped tool (v0 = BLOSUM62, deterministic + wheel-only —
no network, no GPU; the same offline-safe posture as the blastn decoders):

```bash
dna-decode forward --mutation M69L --protein-seq MSIQHFRVALIPFFAAFCLPVFA...   # instant BLOSUM62
dna-decode forward --mutation S83L --protein-fasta gyrA.faa --protein gyrA --json
dna-forward       --mutation A2L  --protein-seq <seq>                          # same, direct entry
dna-decode list   # forward now appears with its DMS validation numbers
```

The WT-coordinate gate fails loudly (exit 2) on a residue/frame mismatch — never a silent wrong call. The
learned methods (ESM2/AlphaMissense/ESM-IF) beat BLOSUM everywhere but need a precomputed score table (they
run the model ONCE per protein), so they stay in the Python API (`predict_effect(..., method="esm2",
esm_table=...)`). **Scope: molecular fitness RANK, not clinical resistance — use `dna-decode amr` for R/S.**

## Run (scripts — validation harnesses)

```bash
uv run python scripts/tem1_forward_cell.py --dms-id <assay> [--method blosum62|esm2|alphamissense]
uv run python scripts/forward_leaderboard.py         # per-protein 4-method table
uv run python scripts/blatem_genome_demo.py          # real nucleotide-edit end-to-end
uv run python scripts/forward_router_demo.py         # regime routing on real predictors
uv run python scripts/forward_dosage_cell.py         # calibrated magnitude on PTEN
uv run python scripts/forward_dosage_sweep.py        # cross-organism dosage generalization
# ESM-IF needs the GNN stack (torch_scatter) — run on Linux/GPU (Kaggle kernel); see
# wiki/forward_esm_if_kaggle_result_2026-07-15.md for the 5 install traps.
```

Tests: `tests/test_forward_{variant_effect,genome_edit,router,alphamissense,structure,dosage}.py` (35).
Frozen decoder surface byte-unchanged throughout; `dna_decode/forward` is NON-frozen.
