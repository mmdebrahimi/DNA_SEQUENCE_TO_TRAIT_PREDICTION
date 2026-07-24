# Phage RBP caller — CROSS-LAB independent number (2026-07-24)

**Corrects the same-day `/research` "data-blocked" verdict** (it missed phageReceptor). The user relayed a
DeepSeek pointer to the **phageReceptor** database (Zhang et al., *Bioinformatics* 2020 36(10):2975; Peng lab,
Hunan U) — a literature-curated phage→host-receptor DB **independent of LBNL/Arkin-Mutalik and BASEL/Maffei**.
It has 37 E. coli phages with outer-membrane-protein receptors in the caller's vocabulary; excluding classic
model phages (T4/T7/T5/λ/N4) + LBNL-name overlaps (M1,T1) leaves **24** independent test phages.

- **Reference (the caller):** committed LBNL RBP reference (data/phage_ref/rbp_reference.faa)
- **Test labels (independent):** phageReceptor DB (Zhang et al. Bioinformatics 2020; Peng lab, Hunan U) — INDEPENDENT of LBNL + BASEL
- **Protocol:** fetch each phage genome from GenBank → extract every 'tail fiber' CDS protein → best
  nearest-neighbour vs the LBNL RBP reference (`call_rbp_from_protein`) → compare to phageReceptor's measured receptor.

## Result

- Test 24 · called 11 · correct 4 → **cross-lab accuracy 0.364**
- Unscored: 2 no-genome, 2 no-tail-fiber-annotation, 9 INDETERMINATE (no LBNL RBP homolog ≥ 0.05)
- Per-receptor [correct/called]: { OmpC: 1/4, OmpA: 2/2, OmpF: 0/2, Tsx: 1/1, FadL: 0/1, FhuA: 0/1 }

| ok | phage | measured | predicted | sim |
|---|---|---|---|---|
| ✓ | Enterobacteria phage Bp7 | OmpC | OmpC | 1.0 |
| ✓ | Enterobacteria phage Mi | OmpA | OmpA | 0.817 |
| ✗ | Enterobacteria phage TuIa | OmpF | OmpC | 0.555 |
| ✗ | Enterobacteria phage vB_EcoM_IME281 | OmpF | Tsx | 0.867 |
| ✓ | Enterobacteria phage vB_EcoM_IME339 | Tsx | Tsx | 0.867 |
| ✓ | Enterobacteria phage vB_EcoM_IME340 | OmpA | OmpA | 0.747 |
| ✗ | Enterobacteria phage vB_EcoM_IME341 | FadL | OmpF | 0.63 |
| ✗ | Escherichia phage AR1 | OmpC | Tsx | 0.867 |
| ✗ | Escherichia phage PP01 | OmpC | Tsx | 0.867 |
| ✗ | Escherichia virus N15 | FhuA | LamB | 0.097 |
| ✗ | Phage 434 | OmpC | LamB | 0.676 |

## Honest reading (TWO corrections)

1. **My /research "data-blocked" verdict was WRONG.** phageReceptor IS a usable independent E. coli OMP-receptor
   source; the cross-lab number is achievable. The DeepSeek pointer was right; my searches missed it.
2. **The within-LBNL LOO (0.975) does NOT generalize cross-lab — the honest cross-lab number is 0.364.** The
   within-study number was optimistic. Failure modes: T4-like phages (AR1/PP01, true OmpC) mis-transfer to Tsx
   relatives at high k-mer similarity (RBP backbone similar, receptor different); 9/24 abstain (LBNL reference
   doesn't span these phages' RBP space).

**Caveat (do not over-read the exact 0.364):** the RBP is extracted by 'tail fiber' product annotation — a
phage carries several tail fibers and the best-match heuristic can pick the wrong one, so part of the drop is
extraction noise, not pure caller failure. The load-bearing finding is the DIRECTION + MAGNITUDE: within-study
0.975 → cross-lab ~0.36 is a large generalization gap (the "validate cross-lab, not within-study" lesson).

## Reproduce
```bash
uv run --with biopython python scripts/phagereceptor_crosslab_validate.py   --api http://www.computationalbiology.cn:18887/viralRecepetor
```
