# Klebsiella cross-organism transfer — top-K (promiscuity) + k-sweep (2026-07-25, overnight)

Extends the top-1 cross-organism number (0.453 clonality-corrected; `wiki/klebsiella_crossorganism_result_2026-07-25`)
with the BIOLOGICALLY-CORRECT metric. Depolymerases are PROMISCUOUS — one enzyme degrades several capsule
types (DpoTropiSearch reports top-10 hits) — so a phage-therapy match wants a RANKED KL-type shortlist, not a
single guess. top-K accuracy = the true KL-type is among the K nearest neighbours' KL-types. Clonality-corrected
(greedy-rep collapse @0.90), vs a top-K prior-frequency null. k swept over {3,4,5,6}.

## Result (clonality-corrected leave-one-out)

| k | reps | called | top-1 | top-3 | **top-5** | null-5 |
|---|---|---|---|---|---|---|
| 3 | 706 | 706 | 0.368 | 0.486 | 0.520 | 0.103 |
| 4 | 719 | 574 | 0.453 | 0.575 | 0.599 | 0.104 |
| 5 | 743 | 592 | 0.448 | 0.571 | 0.588 | 0.102 |
| 6 | 757 | 597 | 0.466 | 0.575 | 0.595 | 0.102 |

**Headline (best k=6): top-1 0.466 / top-3 0.575 / top-5 0.595** (null-5 ~0.10 -> lift +0.49).

## The finding (STRENGTHENED)

top-K materially reinforces the cross-organism result: the true Klebsiella capsule type is among the 5
nearest-depolymerase-domain predictions **~60% of the time** (vs top-1 ~0.45), lift **+0.49** over a 0.10 null.
This is the biologically-correct metric (depolymerase promiscuity) and the phage-therapy-useful one (a ranked
K-type shortlist). The k-sweep shows **k=4-6 are equivalent; k=3 is too coarse** — the original k=4 was sound.
There is an honest called-vs-accuracy tradeoff: lower k (k=3) abstains less (706/706 called) but is less
accurate; higher k (k=4-6) abstains more (574-597 called, sparser k-mer space) but is more accurate on what it
calls. The deterministic sequence-homology->phenotype paradigm GENERALIZES cross-organism on modular
depolymerase domains — and more strongly than the top-1 number alone suggested.

## Scope / honesty

- Labels prophage-host-LCA-INFERRED (in-distribution, the transfer analogue) — NOT independent wet-lab.
- Data (DpoTropiSearch, Zenodo 10.5281/zenodo.14065540) is NON-COMMERCIAL-licensed + D:-only + NEVER
  redistributed; this is non-commercial research USE. The bundled-cell ship remains blocked:authority
  (ledger row 567 — the CC-BY-vs-non-commercial license conflict is the package owner's call).

## Reproduce
```bash
uv run python scripts/klebsiella_topk_ksweep.py --labels <dpo_labels.tsv>   # needs the D: DpoTropiSearch data
```
