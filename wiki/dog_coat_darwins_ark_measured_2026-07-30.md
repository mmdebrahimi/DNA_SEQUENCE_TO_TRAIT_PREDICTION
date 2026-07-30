# Dog coat-colour cell — Darwin's Ark MEASURED concordance (2026-07-30)

**Verdict: BLACK validated (0.994, N=161); other base colours SUBSTRATE-LIMITED.** The deterministic
epistasis cell was scored per-individual on the free Darwin's Ark cohort (Dryad doi:10.5061/dryad.83bk3jb4r;
`canfam4_gp-0.70_biallelic` PLINK set, 3,277 dogs × 29,089,701 **biallelic SNVs**, imputed; N=1,930
owner-reported coat colours). The measured result **confirms the /probe's prediction**: a SNP-only imputed
panel does not carry the indel/structural/low-frequency causal variants needed to call most base colours.

## What was scorable (verified against the REAL .bim, not from memory)

Causal-variant coords were pinned from OMIA (canFam3.1) → lifted to canFam4 (UCSC canFam3ToCanFam4 chain,
pyliftover) → matched in the `.bim` (IDs are `chr:pos:ref:alt`) → **functionally validated by phenotype**:

| locus | variant | canFam4 | in panel? |
|---|---|---|---|
| **B** TYRP1 bc | c.121T>A | chr11:33376317 T:A | ✅ exact liftover match |
| **D** MLPH d1 | c.-22G>A | chr25:48403161 G:A | ✅ exact match |
| **D** MLPH d2 | c.705G>C | chr25:48431759 G:C | ✅ exact match |
| **E** MC1R e | c.916C>T (R306*) | ~chr5:64186728 | ❌ imputation gap (no clean SNP at/near the lifted pos; a +15 bp SNP was the WRONG variant — caught by the identity check: its "e/e" dogs were NOT the red dogs) |
| **B** TYRP1 bs | c.991C>T (Q331*) | ~chr11:33385200 | ❌ imputation gap (7 variants in the 4 kb window; the common brown allele) |
| **K** CBD103 K^B | ΔG23 (3 bp del) | — | ❌ indel — absent from a biallelic-SNV panel |
| **A** ASIP A^y/a^t | SINE insertion + coding | — | ❌ structural/multi — absent |
| **D** MLPH d3 | c.667_668insC | — | ❌ frameshift indel — absent |

## Measured concordance (435 single-colour dogs)

- **black 160/161 = 0.994** — the eumelanin default call (E-/no-brown/no-dilute), correctly made from the
  loci present. The one clean measured win.
- **blue/grey 11/31 = 0.355** — MLPH d1+d2 (verified present) partially recover dilution; misses reflect
  imputation noise + non-MLPH greys.
- **red/yellow — NOT SCORABLE.** MC1R e is imputation-gapped AND "red/yellow" is dominated by **A^y sable/fawn
  (ASIP)**, not e/e — a SNP panel structurally cannot carry ASIP's SINE, so the reds are unreachable by
  construction (a functional MC1R-region scan found NO SNP separating the 167 reds from black: best 55/167
  red-hom vs 3/161 black-hom).
- **brown/liver — NOT SCORABLE.** TYRP1 bs (the common brown allele) is imputation-gapped; only rare bc is
  present.

## Honest conclusion (confirms the /probe)

The deterministic cell is CORRECT in principle — black is called at 0.994 — but the free Darwin's Ark
substrate is a **biallelic-SNV imputed panel**, which by construction lacks the indel (CBD103 ΔG23, MLPH d3),
structural (ASIP SINE), and imputation-gapped (MC1R e, TYRP1 bs) causal variants that the *other* base
colours require. So the dog cell's measured tier on THIS substrate is **"black-validated + dilution-partial;
other colours substrate-limited"**, NOT a clean full-colour pass. This is a genuine finding, not a cell
defect: it is the exact regime the /probe flagged (indel/structural/low-freq loci missing from a SNP panel).

A full-colour measured validation would need a substrate that genotypes the causal indels/SVs directly
(e.g. Embark/VGL panels or WGS variant calls that include CBD103/ASIP/TYRP1-bs), not an imputed biallelic-SNV
set. The cell's ABSTAIN-on-missing-loci behavior handled the absence gracefully (no fabricated distribution
calls). Reproducibility: `scripts/dog_coat_darwins_ark_validate.py` (phenotype ingest, verified schema) +
`dna_decode/pigment/plink_io.py` (extractor) + the OMIA→liftover→.bim pinning above.
Frozen AMR/forward surfaces byte-unchanged (read-only scoring).
