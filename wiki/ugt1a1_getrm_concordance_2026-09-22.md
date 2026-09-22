# UGT1A1 tag-SNP surface vs GeT-RM CDC consensus (2026-09-22)

Scored **65** 1000G samples with CDC multi-lab consensus UGT1A1 diplotypes (32 carry `*28`). Bar frozen beforehand: `wiki/ugt1a1_tag_acceptance_bar.json`.

## The shipped tag arm (rs887829 / `*80` as a proxy for `*28`)

- sensitivity **1.0** · specificity **0.9394** · PPV **0.9412** · NPV **1.0** (n=65)
- exact `*28` **zygosity** agreement: 61/65 (0.9385)
- **verdict: TAG_CONFIRMED**

## The direct TA-repeat arm (the catalog's asserted structural wall)

- sensitivity **1.0** · specificity **1.0** · PPV **1.0** · NPV **1.0** (n=65)
- exact `*28` **zygosity** agreement: 64/65 (0.9846)
- **verdict: DIRECT_STRICTLY_DOMINATES_TAG**

**Scope limit.** A panel-specific result. It shows the repeat is callable FROM THE 1000G 30x PHASED PANEL, NOT that it is callable from an arbitrary short-read SNP VCF. Do NOT generalise to 'the structural wall is false'; the honest correction is 'the wall is panel-dependent and this panel clears it'.

## By superpopulation (the catalog names non-EUR LD breakdown as an unmeasured residual)

| pop | n | tag sens | tag PPV | direct sens |
|---|---|---|---|---|
| AFR | 21 | 1.0 | 0.9333 | 1.0 |
| AMR | 2 | None | 0.0 | None |
| EAS | 16 | 1.0 | 1.0 | 1.0 |
| EUR | 26 | 1.0 | 1.0 | 1.0 |

## POST-HOC: the clinical axis (reduced FUNCTION, not `*28` specifically)

*Not in the frozen bar. Added after diagnosing the frozen arm's errors; reported beside the frozen primary, never in place of it.*

**Diagnosed mechanism.** rs887829-T tags the reduced-function promoter haplotype class -- BOTH *28 (TA7) and *37 (TA8) -- not *28 alone. Confirmed directly: the three *28-specific tag errors (HG01190, NA19920, NA19239) are haplotypes where the C>CATAT (+2 TA = *37) allele is present. Against a *28-specific target these read as false positives; against reduced FUNCTION they are correct calls.

- tag vs reduced-function carriership (`*28` or `*37`): sensitivity **1.0** · specificity **1.0** · PPV **1.0** (n=65), zygosity-exact 64/65

## Unresolved panel-vs-consensus discrepancy

- **NA20509** — truth `(*28)/(*28 + *60) -- GeT-RM consensus = 2 copies of *28`; observed 1 copy, by BOTH independent assays (rs887829 tag GT 0|1 AND the C>CAT repeat indel GT 0|1, phase-consistent on the same haplotype); the other haplotype is reference TA6 at every repeat allele. Two independent assays of different variant types agree with each other and against the reference consensus. Reported as an UNRESOLVED panel-vs-consensus discrepancy -- this run does NOT adjudicate which is right, and it is the single error in the direct arm.

## Control: `*6` (rs4148323, a direct SNP -- not a tag)

- n=65 · sensitivity 1.0 · PPV 1.0 · zygosity-exact 1.0

## Honest limits

- GeT-RM consensus is caller-derived, not wet-lab -- bounds agreement, never correctness
- n=65 overlap; per-superpopulation strata will be small and several will be underpowered
- *60 is a separate promoter SNP in LD with *28 and is NOT in our catalog; it is not scored as a miss
- truth diplotypes carrying parenthesised/compound haplotypes are parsed per-haplotype; parse failures are counted and excluded, never silently defaulted
