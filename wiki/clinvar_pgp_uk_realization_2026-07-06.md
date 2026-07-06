# ClinVar/Mendelian decoder on REAL PGP-UK humans (2026-07-06)

**Extension (a):** the deterministic **ClinVar/Mendelian** decoder now runs **end-to-end on real, named,
open-consent human individuals** (PGP-UK) — the Mendelian sibling of the PGx real-people milestone
(`wiki/pgx_pgp_uk_realization_2026-07-05.md`). Non-duplication checked: no ClinVar-on-PGP-UK work exists on
the `mosfaer` branch (workhorse/Soraya-V3-5) — this is genuinely new.

## What ran
- **Decoder:** `dna_decode.data.clinvar.ClinVarDecoder` (curated-catalog deterministic Mendelian decoder;
  P/LP + B/LB across 10 canonical genes: CFTR, HBB, LDLR, BRCA1, BRCA2, TP53, PAH, G6PD, HFE, F8).
- **New per-VCF wrapper:** `scripts/clinvar_decode_vcf.py` — the decoder is a per-variant `call()`; this
  iterates a whole individual VCF, decodes every ALT the person actually **carries**, and collects the
  curated ClinVar classifications. Fail-closed: not-in-panel → INDETERMINATE (absence ≠ benign).
- **Build bridge (zero liftover):** the decoder is a pure dict lookup, so panel and VCF just need the SAME
  build. PGP-UK is GRCh37, so I built a **GRCh37 panel** (`data/clinvar/clinvar_panel_grch37.tsv`, 31,615
  variants) via the existing `scripts/capture_clinvar.py --vcf clinvar_GRCh37.vcf.gz` and intersect the
  GRCh37 individual VCF **natively** — no liftover of the 100 MB genome.

## Result — real curated Mendelian classifications on real people (N=3)

| sample | pathogenic (P/LP) | benign (B/LB) | not-in-panel (abstained) |
|---|---|---|---|
| FR07961000 | 0 | **338** | 4,954,697 |
| FR07961006 | 0 | **206** | 4,152,554 |
| FR07961007 | 0 | **364** | 4,182,670 |

The benign hits are real, curated ClinVar carrier variants — e.g. FR07961000: HFE (Hemochromatosis type 1,
2★), CFTR (Cystic fibrosis, 2★), HBB, PAH; FR07961006 also carries benign BRCA2 variants. **0 pathogenic**
across all 3 is the **honest, expected** result: a random individual rarely carries a P/LP variant in a
10-gene canonical panel, and the decoder abstains honestly on the ~4–5 M variants not in the panel (never
calling absence "benign").

## Honest tier (load-bearing)
- **Research demonstration**, NOT a clinical diagnosis. A ClinVar "Pathogenic/Benign" annotation is a curated
  database classification of the *allele*, not a clinical interpretation of the *person*. **NOT a clinical
  tool.**
- Regime-1 (curated-catalog → deterministic wins), like the AMR determinant catalog — NOT a learned predictor.
- The panel is a committed 10-gene subset; broader coverage = extend the gene list in `capture_clinvar.py`.

## Artifacts
- Wrapper: `scripts/clinvar_decode_vcf.py` (`--vcf <ind.vcf.gz> --panel <panel.tsv>`).
- GRCh37 panel: `data/clinvar/clinvar_panel_grch37.tsv` (committed; sibling of the GRCh38 `clinvar_panel.tsv`).
- Per-individual provenance: `wiki/pgp_uk_clinvar_results/FR0796100{0,6,7}_clinvar.json`.
- Full ClinVar GRCh37 VCF (192 MB) on `D:/dna_decode_cache/clinvar/clinvar_GRCh37.vcf.gz` (gitignored-class).
- Tests: `tests/test_clinvar_decode_vcf.py` (4 offline; the pathogenic + benign paths both fire).

**Frozen bacterial/viral/fungal/TB AMR surface byte-unchanged.**
