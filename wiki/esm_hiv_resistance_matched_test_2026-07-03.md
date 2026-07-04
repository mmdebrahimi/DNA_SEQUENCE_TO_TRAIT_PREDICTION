# FAIR test — does a masked protein-LM carry AMR-resistance signal? (HIV, R-vs-S, 2026-07-03)

**The properly-powered replacement for the underpowered SARS-Mpro probe** (whose n=2 cross-position benign
set + within-R-only Spearman could not answer the question). Design: **single-mutant HIV isolates** (exactly
one non-consensus residue vs HXB2 → the PhenoSense fold is ATTRIBUTABLE to that one variant), balanced R and
S, AUC of ESM `-LLR` (WT-marginal deleteriousness) separating R (fold≥3) from S. No LD confound, no catalog
circularity, matched per-variant. Script: `scripts/esm_hiv_resistance_test.py`.

## Result — zero-shot ESM carries PARTIAL, POSITION-MECHANISM-DEPENDENT resistance signal

| class | gene/drug | R-vs-S AUC | n (R/S) | powered? | reading |
|---|---|---|---|---|---|
| **NRTI** | RT/3TC | **0.821** | 14/14 | yes | M184V at the conserved YMDD active site → deleterious AND resistant → ESM CAPTURES |
| **INSTI** | IN/RAL | **0.724** | 54/121 | yes | Q148K (LLR −9.3) / N155H at the integrase core → ESM CAPTURES |
| **PI** | PR/NFV | 1.0 | 1/32 | no (n_R=1) | I84V active-site (LLR −7.2) → deleterious; underpowered |
| **NNRTI** | RT/EFV | **0.244** | 33/24 | yes | K103N/Y181C in the TOLERANT NNRTI-binding pocket → evolutionarily fine → ESM ANTI-predictive |
| pooled | — | 0.587 | 102/191 | — | weakly positive, dragged down by NNRTI |

Verdict flag: `ESM_CARRIES_RESISTANCE_SIGNAL` (best powered-class AUC 0.821 ≥ 0.70). But the HONEST headline
is the **position-mechanism split**, not the binary.

## The corrected mechanism (supersedes the earlier blanket claim)

The earlier Mpro probe claimed a blanket law "resistance ⊥ evolutionary conservation." The fair test shows
that is **half right**:

- **Resistance ∥ conservation** when the resistance mechanism is a change at a **structurally-constrained /
  active-site** residue (NRTI M184V at YMDD; INSTI Q148/N155 at the integrase core; PI I84V). These ARE
  deleterious to fold-fitness, so a fitness-trained masked-LM SEES them. AUC 0.72–0.82.
- **Resistance ⊥ conservation** when the mechanism is **drug-pocket evasion** at an evolutionarily tolerant
  surface residue (NNRTI K103N/Y181C). These are fitness-neutral, so ESM is BLIND — indeed ANTI-predictive
  (AUC 0.24: ESM rates them MORE tolerated than susceptible variants). The Mpro all-blind result fits here:
  nirmatrelvir resistance is largely pocket-tolerant.

**So a masked-protein-LM is a partial, mechanism-gated resistance signal — strong where resistance rides
structural constraint, blind/anti where it rides pocket evasion.**

## What this greenlights (and what it does NOT)

- **Greenlit (bounded):** a learned scorer adds real value for CONSERVED-position resistance (extend the
  catalog to novel active-site variants it misses). It must NOT be trusted for pocket-evasion classes
  (NNRTI) — there the deterministic curated catalog remains essential and the learned score is anti-signal.
- **Motivates the SUPERVISED head:** zero-shot ESM = fitness; a head trained on fold LABELS could recover
  the NNRTI pocket signal that zero-shot misses (learn "resistance", not "fitness"). This is the decisive
  next test — does supervised beat BOTH the deterministic catalog and zero-shot ESM on held-out variants?
- Does NOT by itself justify the genotype-level world-model (the LD-confound argument is unchanged for pure
  self-supervised genotype models).

## Honest scope
- Single-mutant restriction shrinks N (the price of an attributable label); fold≥3 uniform cutoff; in-
  distribution HIVDB fold (independent-cohort fold is the eventual scale-up).
- PI/NNRTI single-mutant R sets are thin/absent for high-barrier drugs (LPV/DTG) — used lower-barrier
  drugs (NFV/RAL/3TC/EFV) to get R; PI still underpowered on R.
- Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9); READ-only over the frozen
  catalog + committed HXB2 references. ESM weights on D:.

## Reproduce
```bash
HF_HOME=D:/hf_cache TORCH_HOME=D:/dna_decode_cache/torch \
  uv run python scripts/esm_hiv_resistance_test.py --model esm2_t33_650M_UR50D
uv run pytest tests/test_esm_hiv_resistance.py -q   # offline helper tests (no model)
```
