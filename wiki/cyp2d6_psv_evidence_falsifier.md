# CYP2D6-CYP2D7 PSV evidence table — Phase A falsifier (read-level pileup)

_DIAGNOSTIC read-only evidence table; regional D6-fraction profile; NOT an identity caller. Flag contract: samtools mpileup -B -q 0 -Q 0 (permissive; matches the depth path). 117 PSVs from the real Cyrius CYP2D6_SNP_38.txt (GRCh38, paired-coord). Both pos_CYP2D6 AND pos_CYP2D7 counted per PSV._

**Falsifier: GO** — 3/3 hybrids signalled; 3/3 non-hybrids flat.

_The paired-coordinate PSV D6-fraction PROFILE reproduces the Cyrius hybrid signatures: *68 (5' high / 3' low), *13 (opposite), *36 (exon-9-tip dip); non-hybrids (normal/dup/deletion) stay flat. Directional 5'-3' cleanly catches *68 + *13; the *36 exon-9 conversion needs the dedicated exon9-tip feature (the subtlest case). PROOF-OF-SIGNAL on n=1/type -> GO to draft the Phase-B abstaining identity classifier; NOT a powered per-allele validation._

## Regional D6-fraction profile (3'->5'; region x sample)

| region | NA19239 (normal) | HG00436 (puredup) | NA12873 (del) | NA12878 (hyb68) | NA18563 (hyb36) | NA19785 (hyb13) |
|---|---|---|---|---|---|---|
| downstream_exon9 | 0.48 | 0.61 | 0.41 | 0.42 | 0.40 | 0.68 |
| exon9 | -- | -- | -- | -- | -- | -- |
| intron6 | 0.55 | 0.60 | 0.34 | 0.38 | 0.64 | 0.74 |
| exon7 | 0.55 | 0.66 | 0.23 | 0.40 | 0.56 | 0.73 |
| exon6 | 0.50 | 0.59 | 0.36 | 0.37 | 0.65 | 0.72 |
| intron5 | 0.49 | 0.60 | 0.40 | 0.34 | 0.56 | 0.69 |
| intron4 | 0.51 | 0.60 | 0.36 | 0.37 | 0.60 | 0.75 |
| exon3 | 0.47 | 0.54 | 0.36 | 0.38 | 0.57 | 0.75 |
| intron2 | 0.50 | 0.61 | 0.33 | 0.39 | 0.56 | 0.72 |
| exon2 | 0.59 | 0.78 | 0.25 | 0.41 | 0.59 | 0.75 |
| intron1 | 0.47 | 0.59 | 0.28 | 0.67 | 0.53 | 0.46 |
| exon1 | 0.49 | 0.59 | 0.29 | 0.68 | 0.62 | 0.41 |
| upstream_exon1 | 0.51 | 0.58 | 0.30 | 0.63 | 0.50 | 0.57 |
| **5'-3' shift** | -0.02 | -0.01 | -0.08 | +0.27 | -0.01 | -0.23 |
| **exon9-tip dip** | +0.02 | -0.01 | -0.08 | -0.02 | +0.17 | +0.04 |
| **signal** | flat_nonhybrid | flat_nonhybrid | flat_nonhybrid | directional_5p_high_3p_low | exon9_tip_dip | directional_5p_low_3p_high |
