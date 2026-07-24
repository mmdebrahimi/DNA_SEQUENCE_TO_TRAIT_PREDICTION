# Phage receptor-class cell v0 (2026-07-24)

**The first non-AMR, non-host-organism cell** — a bacteriophage-genome -> host-receptor-class decoder.

- **Axis:** bacteriophage host-receptor class (first non-AMR, non-host-organism cell)
- **Substrate:** BASEL E. coli phage collection (Maffei 2021 PLOS Biology 3001424; GenBank MZ501046-MZ501113)
- **Label source (FREE MEASURED):** experimentally-determined receptors (>50 single-gene K-12 mutants + EOP host-range), CC-BY
- **Method:** genome-homology receptor TRANSFER (nearest-BLAST-neighbour inherits its receptor)
- **Scope:** RECEPTOR-CLASS only (NOT the full phage x strain host-range matrix, which is polygenic/intractable from genome alone); clade-conserved clades only
- **Tier:** `IN_DISTRIBUTION` — closed for v0 - labels are clade-derived from the same BASEL Results the catalog is curated from; an INDEPENDENT number needs a held-out phage set with measured receptors

## Result (leave-one-out, native blastn)

- Genomes fetched: **68**; clean-labelled (clade-conserved): **29**;
  excluded RBP-variable: **39**.
- **Overall LOO accuracy: 1.000** (27/29 called; 2
  INDETERMINATE abstentions; **0 mis-calls**).

| receptor | LOO correct/called |
|---|---|
| BtuB | 9/9 |
| ECA | 12/12 |
| LPS_core | 4/4 |
| LptD | 2/2 |

## Honest reading

receptor-class TRANSFERS reliably along genome homology within clade-conserved clades (0 mis-calls); the caller ABSTAINS (INDETERMINATE) rather than mis-transfer when a phage has no reference homolog (e.g. the lone NfrA phage). The 100% is on the clade-conserved subset BY CONSTRUCTION - it validates the pipeline + abstention + catalog self-consistency, not a solved RBP->receptor map. RBP-variable clades are the documented tractability boundary.

**Excluded clades (the tractability boundary — receptor is receptor-binding-protein-determined, not clade-clean):**
- Tequatrovirus/Straboviridae (T-even: OmpC/FadL/Tsx vary by RBP)
- Drexlerviridae (FhuA/BtuB/YncD/TolC vary)
- Siphoviridae Dhillonvirus/Nonagvirus/Seuratvirus (LptD/FhuA/LamB vary)

## Reproduce

```bash
uv run python scripts/build_phage_receptor_report.py
# or just the number:
uv run python -c "from scripts.phage_receptor_caller import _load_manifest, leave_one_out; \
r=leave_one_out(*_load_manifest('data/phage_ref/basel_manifest.tsv','data/phage_ref/basel')); \
print(r.accuracy, r.per_receptor)"
```
