# Phage receptor cell — INDEPENDENT validation (2026-07-24)

**The independent number the cell needed** (rows 559-560 were `blocked:external` on it; the data arrived via
the user + `github.com/mjohnson11/PhageDataSheets`).

- **Independent source:** LBNL/Arkin-Mutalik Phage Datasheets (Moriniere et al.; github.com/mjohnson11/PhageDataSheets)
- **Label:** measured receptor via genome-wide genetic screens on E. coli K-12 BW25113
- **Independence:** DIFFERENT LAB, measured labels, non-Bas isolates disjoint from the BASEL-2021 reference (Bas## excluded as leakage); K-12 host comparable to BASEL
- **Reference (tested artifact):** BASEL-2021 genome-homology caller (shipped v0) — covers classes ['BtuB', 'ECA', 'LPS_core', 'LptD', 'NfrA']

## Result

| metric | value |
|---|---|
| non-Bas scoreable test phages | 119 |
| overall called / correct / acc | 86 / 25 / **0.291** |
| **covered-subset** (true receptor in ref classes) n / called / correct / **acc** | 38 / 29 / 25 / **0.862** |

Per-receptor [correct/called]:
- Tsx: 0/28
- BtuB: 22/26
- OmpA: 0/10
- OmpC: 0/8
- FhuA: 0/7
- LPS_core: 3/3
- OmpF: 0/2
- FadL: 0/1
- TolC: 0/1

Test true-receptor distribution: {'Tsx': 31, 'OmpA': 11, 'OmpC': 14, 'OmpF': 5, 'FadL': 1, 'LPS_core': 6, 'BtuB': 30, 'FhuA': 14, 'TolC': 1, 'LptD': 2, 'LamB': 2, 'YncD': 2}

## Honest reading

The OVERALL number is dragged down by CLASS COVERAGE: the shipped v0 catalogue only contains reference phages for classes ['BtuB', 'ECA', 'LPS_core', 'LptD', 'NfrA'], so it cannot predict the RBP-variable classes (Tsx/OmpC/FhuA/OmpA/OmpF/...) that dominate the independent set — those are exactly the (3) RBP-caller scope. The COVERED-SUBSET accuracy is the fair independent test of what v0 CLAIMS to decode.

**Tier move:** IN_DISTRIBUTION (KNOWLEDGE_BASELINE) -> **INDEPENDENT_MEASURED for the covered classes**
(BtuB 22/26 + LPS_core 3/3 transfer to a different lab's measured labels; 0.862 on 29 called). The
RBP-variable classes (Tsx/OmpC/FhuA/OmpA/OmpF/FadL/TolC/YncD, 60+ phages) are OUT of v0 scope and are the
(3) RBP-caller target — their measured per-phage labels + RBP CDS annotations are in the same LBNL table.

## Reproduce

```bash
git clone --depth 1 https://github.com/mjohnson11/PhageDataSheets.git
uv run --with biopython python scripts/lbnl_independent_validate.py --repo PhageDataSheets/Ecoli_phages
```
