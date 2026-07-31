# Dog "world model" masked reconstruction — F1′ CLEAN re-run (2026-07-31)

Supersedes the deprecated F1 smoke (`wiki/dog_masked_reconstruct_smoke_2026-07-31.md`). An adversarial
review found that smoke's `−0.15` headline biased against NT on two axes: (1) the Markov baseline was
fit on the SAME slice it scored (transductive label leakage → inflated Markov accuracy), and (2) NT was
scored by its single argmax 6-mer, discarding the per-base marginal distribution. This clean re-run
fixes both and re-runs on real NT weights.

## Method (leakage-free, fair)

- **NT per-base MARGINALS**, not argmax token: `P(base@offset j) = Σ` softmax mass over all vocab tokens
  whose j-th base is that base. Primary endpoint = **per-base negative log-likelihood (NLL)**.
- **Markov leakage control, two ways:** (a) **disjoint-fit** — trained on canFam4 chr1 bp 1200–6001,
  evaluated on bp 0–600 (no overlap); (b) **leave-one-out** — fit on the eval window but each target's
  own count excluded before predicting it.
- **Baseline given its best shot:** k-sweep 1–8; headline uses the k with the LOWEST Markov NLL (the
  HARDEST baseline for NT to beat) = k=3, disjoint.
- Both sides scored on the identical masked base set; strict guard raises on any dropped position.
- Substrate: canFam4 chr1 `NC_051804.1` slice; 600 bp eval window; 40 masked 6-mer tokens (240 bases).

## Result

| metric | deprecated F1 smoke | **F1′ clean** |
|---|---|---|
| Markov per-base accuracy | 0.4875 *(leaky same-slice, teacher-forced)* | **0.3417** *(disjoint/LOO)* |
| NT per-base accuracy | 0.3375 *(argmax-token)* | **0.3500** *(marginal)* |
| **accuracy delta (NT − Markov)** | **−0.150** | **+0.008 (a wash)** |
| NT per-base NLL | — | 1.4474 |
| Markov per-base NLL (disjoint, best k=3) | — | 1.3354 |
| **PRIMARY: NLL delta (Markov − NT), disjoint** | — | **−0.112 (NT slightly loses)** |
| NLL delta vs leave-one-out Markov | — | **+0.003 (tie)** |

disjoint NLL delta by k: `{1:−0.10, 2:−0.10, 3:−0.11, 4:−0.10, 5:−0.05, 6:−0.03, 7:−0.03, 8:−0.02}`.

## Interpretation (honest, corrected)

- **The −0.15 was almost entirely baseline label leakage.** Removing it drops Markov accuracy from
  0.49 → 0.34; NT was essentially unchanged (0.34 argmax → 0.35 marginal). The argmax-vs-marginal axis
  mattered far less than the leakage.
- **Corrected verdict: NT-v2-100M is at approximate PARITY with a low-order Markov chain** on canFam4
  masked reconstruction — a per-base accuracy wash (+0.008), and only a small NLL loss (−0.11 nats) to
  the strongest (data-rich disjoint) baseline; a tie (+0.003) against the leave-one-out baseline.
- This is **"parity", not "loses"** — a materially milder claim than the deprecated smoke. It is still
  an **unimpressive capability signal**: a 100M-parameter genomic LM that only ties a 3-mer Markov chain
  at reconstructing real dog sequence is not demonstrating deep sequence "understanding". But the honest
  statement is parity, and it is scoped to NT-v2-100M at this smoke scale.
- **Why disjoint Markov beats LOO Markov as the baseline:** the disjoint train region is ~4.8 kb vs the
  600 bp eval window, so its distributions are better-estimated (lower NLL) — the harder, fairer baseline.

## Scope + next (F2)

Smoke scale (40 tokens / 240 bases, one window, NT-v2-100M). A trustworthy general verdict needs F2:
region-stratified full-chromosome sweep (coding/intergenic/conserved), a compression/bidirectional
baseline alongside Markov, and — if the parity holds — a tightly-scoped bank, OR a bigger NT / single-nt
model (Caduceus) before any generalization. Do NOT generalize this to "genomic LMs don't work".

## Engine changes (this run)
- `markov_baseline.py`: `base_distribution` (Laplace-smoothed, LOO via `exclude_base`), `nll_on_masked`,
  `accuracy_on_masked(leave_one_out=…)`.
- `foundation.py`: NT + mock emit per-base `base_marginals`; `masked_token_predictions(strict=…)` raises
  on dropped positions.
- `masked_reconstruct.py`: per-base NLL primary; `score_from_predictions` (compute NT forward once, sweep
  Markov cheaply); disjoint/LOO modes.
- 10 offline tests (`tests/test_masked_reconstruct.py`); frozen AMR/forward surfaces byte-unchanged.
