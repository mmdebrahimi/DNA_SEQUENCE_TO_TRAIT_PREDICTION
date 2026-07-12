# Determinant co-occurrence / linkage world model + phenotype->genotype inverter (2026-07-11)

**Verdict: PASS_LINKAGE_STRUCTURE** — 143/144 within-organism testable determinants are LINKED (fraction 0.993; PASS bar 0.5). (raw)

Is there within-organism co-resistance LINKAGE (do a genome's other determinants predict a held-out determinant beyond base rates)? + the phenotype->genotype inversion.

Substrate: cached AMRFinder determinant calls (self-distillation). LINKED = a determinant's held-out presence is predicted from the OTHER determinants, WITHIN organism, with OOF-AUC 95% CI lower > 0.5.

| organism | genomes | testable | **linked** | frac | note |
|---|---|---|---|---|---|
| klebsiella | 307 | 60 | **60** | 1.0 |  |
| escherichia_coli_shigella | 240 | 45 | **44** | 0.978 |  |
| campylobacter | 100 | 4 | **4** | 1.0 |  |
| acinetobacter | 60 | 22 | **22** | 1.0 |  |
| salmonella | 60 | 13 | **13** | 1.0 |  |

## C3a — 'vice-versa' co-occurrence LIFT (what determinants travel together; verify-in-batch)
### klebsiella
- **fosA** → aac(3)-IVa(lift 1.07, n24), aac(6')-Ib(lift 1.07, n58), aac(6')-Ib'(lift 1.07, n8), ampC_Kaer-1(lift 1.07, n9)
- **oqxA** → aac(3)-IVa(lift 1.19, n24), ampC-Kaer(lift 1.19, n30), ampC_Kaer-1(lift 1.19, n9), aph(4)-Ia(lift 1.19, n24)
- **oqxB** → aac(3)-IVa(lift 1.67, n24), ampC-Kaer(lift 1.67, n30), ampC_Kaer-1(lift 1.67, n9), aph(4)-Ia(lift 1.67, n24)
- **parC_S80I** → aac(3)-IVa(lift 2.72, n24), aph(4)-Ia(lift 2.72, n24), ompK35_E42RfsTer47(lift 2.72, n72), gyrA_S83I(lift 2.66, n95)
- **sul1** → arr-3(lift 2.84, n10), cmlA5(lift 2.84, n8), dfrA12(lift 2.84, n61), aadA2(lift 2.69, n71)
- **gyrA_S83I** → aac(3)-IVa(lift 3.16, n24), aph(4)-Ia(lift 3.16, n24), ompK35_E42RfsTer47(lift 3.16, n72), sul3(lift 3.05, n27)
### escherichia_coli_shigella
- **gyrA_S83L** → gyrA_D87N(lift 1.51, n140), parC_A56T(lift 1.51, n8), parC_E84V(lift 1.51, n77), parE_L416F(lift 1.51, n11)
- **parC_S80I** → aac(6')-Ib-cr5(lift 1.7, n69), parC_A56T(lift 1.7, n8), parC_E84V(lift 1.7, n77), parE_S458A(lift 1.7, n33)
- **gyrA_D87N** → parC_A56T(lift 1.71, n8), parC_E84V(lift 1.71, n77), parE_L416F(lift 1.71, n11), parE_S458A(lift 1.71, n33)
- **uhpT_E350Q** → parC_E84V(lift 2.16, n77), parE_I529L(lift 2.16, n88), ptsI_V25I(lift 2.16, n94), blaCTX-M-27(lift 2.06, n21)
- **tet(A)** → bleO(lift 2.2, n8), cmlA1(lift 1.86, n11), blaCTX-M-27(lift 1.7, n17), dfrA12(lift 1.63, n20)
- **sul1** → aadA5(lift 2.14, n86), dfrA17(lift 2.08, n87), mph(A)(lift 1.99, n87), mrx(A)(lift 1.99, n87)
### campylobacter
- **blaOXA-193** → blaOXA-61_G-57T(lift 1.33, n17), aph(3')-IIIa(lift 1.25, n8), gyrA_T86I(lift 1.13, n40), tet(O)(lift 1.03, n38)
- **tet(O)** → aph(3')-IIIa(lift 1.92, n9), blaOXA-61_G-57T(lift 1.28, n12), gyrA_T86I(lift 1.23, n32), blaOXA-193(lift 1.03, n38)
- **gyrA_T86I** → blaOXA-61_G-57T(lift 1.67, n15), tet(O)(lift 1.23, n32), blaOXA-193(lift 1.13, n40)
- **blaOXA-61_G-57T** → gyrA_T86I(lift 1.67, n15), blaOXA-193(lift 1.33, n17), tet(O)(lift 1.28, n12)
### acinetobacter
- **parC_S84L** → aac(3)-IId(lift 1.5, n8), aph(3')-VIb(lift 1.5, n9), blaADC(lift 1.5, n8), blaADC-30(lift 1.5, n9)
- **sul2** → blaADC(lift 1.67, n8), blaOXA-100(lift 1.52, n10), ant(2'')-Ia(lift 1.34, n25), blaOXA-69(lift 1.27, n13)
- **ant(2'')-Ia** → blaADC(lift 1.94, n8), aph(3')-VIa(lift 1.83, n17), blaADC-79(lift 1.81, n14), blaOXA-100(lift 1.76, n10)
- **aph(3'')-Ib** → aac(3)-IId(lift 2.14, n8), aph(3')-VIb(lift 2.14, n9), blaADC-30(lift 2.14, n9), blaADC-76(lift 2.14, n12)
- **sul1** → aac(3)-IId(lift 2.4, n8), aadA1(lift 2.4, n12), blaCARB-2(lift 2.4, n8), blaPSE(lift 2.4, n8)
- **aph(3')-Ia** → blaTEM-1(lift 2.42, n8), aadA1(lift 2.27, n10), aac(3)-Ia(lift 2.1, n10), blaOXA-69(lift 1.93, n12)
### salmonella
- **tet(A)** → blaCTX-M-65(lift 2.31, n11), sul1(lift 2.19, n19), aadA1(lift 2.17, n16), aph(3')-Ia(lift 2.13, n12)
- **sul1** → blaCTX-M-65(lift 3.0, n11), aadA1(lift 2.82, n16), dfrA14(lift 2.73, n10), gyrA_D87Y(lift 2.68, n17)
- **gyrA_D87Y** → blaCTX-M-65(lift 3.16, n11), dfrA14(lift 3.16, n11), aadA1(lift 2.97, n16), sul1(lift 2.68, n17)
- **aadA1** → blaCTX-M-65(lift 3.53, n11), dfrA14(lift 3.21, n10), gyrA_D87Y(lift 2.97, n16), sul1(lift 2.82, n16)
- **aph(3')-Ia** → dfrA14(lift 4.62, n11), floR(lift 3.55, n10), aac(3)-IVa(lift 3.46, n9), aph(4)-Ia(lift 3.46, n9)
- **floR** → aac(3)-IVa(lift 4.23, n11), aph(4)-Ia(lift 4.23, n11), aph(3')-Ia(lift 3.55, n10), blaCTX-M-65(lift 3.36, n8)

## C3b — phenotype→genotype inversion (given a drug Class, the determinants that confer it)
### klebsiella
- **AMINOGLYCOSIDE** ← aph(6)-Id(78), aph(3'')-Ib(76), aadA2(75), aac(6')-Ib(58), aadA1(44)
- **AMINOGLYCOSIDE/QUINOLONE** ← aac(6')-Ib-cr5(48), aac(6')-Ib-cr(1)
- **BETA-LACTAM** ← blaTEM-1(94), blaSHV-11(91), blaCTX-M-15(87), ompK35_E42RfsTer47(72), blaSHV-12(71)
- **BLEOMYCIN** ← ble(9), bleO(1)
- **COLISTIN** ← crrB_L94M(3), mgrB_I12HfsTer13(3), crrB_S195N(2), crrB_Q10L(1), phoQ_D438N(1)
- **FOSFOMYCIN** ← fosA(287), fosA5(12), fosA7(9), fosA10(5)
- **LINCOSAMIDE/MACROLIDE** ← erm(42)(1)
- **LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN** ← erm(B)(3)
### escherichia_coli_shigella
- **AMINOGLYCOSIDE** ← aph(3'')-Ib(100), aadA5(90), aph(6)-Id(86), aac(3)-IIe(54), aadA1(40)
- **AMINOGLYCOSIDE/QUINOLONE** ← aac(6')-Ib-cr5(69), aac(6')-Ib-cr(1)
- **BETA-LACTAM** ← blaCTX-M-15(98), blaTEM-1(80), blaOXA-1(69), blaCTX-M-27(22), blaCTX-M-55(15)
- **BLEOMYCIN** ← bleO(8), ble(7)
- **COLISTIN** ← mcr-1.1(4), pmrB_C84Y(1), mcr-3.2(1), pmrB_T156M(1)
- **FOSFOMYCIN** ← uhpT_E350Q(111), ptsI_V25I(94), fosA(1), fosA3(1)
- **FOSMIDOMYCIN** ← cyaA_S352T(31)
- **LINCOSAMIDE** ← lnu(F)(2)
### campylobacter
- **AMINOGLYCOSIDE** ← aph(3')-IIIa(9), aadE-Cc(3), aad9(2), ant(6)-Ia(2), aph(2'')-Ig(1)
- **BETA-LACTAM** ← blaOXA-193(71), blaOXA-61_G-57T(18), blaOXA(4), blaOXA-461(3), blaOXA-591(3)
- **LINCOSAMIDE** ← lnu(C)(1)
- **MACROLIDE** ← rplV_A103V(5), 23S_A2075G(1), 23S_A2074T(1)
- **QUINOLONE** ← gyrA_T86I(50)
- **STREPTOTHRICIN** ← sat4(2)
- **TETRACYCLINE** ← tet(O)(52), tet(O/M/O)(1)
### acinetobacter
- **AMINOGLYCOSIDE** ← ant(3'')-IIa(60), ant(2'')-Ia(31), aph(3'')-Ib(28), aph(6)-Id(22), aph(3')-Ia(22)
- **BETA-LACTAM** ← blaOXA-23(19), blaOXA-58(19), blaOXA-69(17), blaADC-79(15), blaADC-76(12)
- **COLISTIN** ← pmrB_T232I(1), pmrB_S17R(1), pmrA_E8D(1)
- **EFFLUX** ← adeS_R152K(1), adeR_P116L(1), adeS_H189Y(1)
- **MACROLIDE** ← mph(E)(19)
- **MACROLIDE/STREPTOGRAMIN** ← msr(E)(19)
- **PHENICOL** ← floR(3), catA1(2), catB8(1), cmlA1(1)
- **QUINOLONE** ← gyrA_S81L(56), parC_S84L(40), parC_S84F(1), parC_E88K(1)
### salmonella
- **AMINOGLYCOSIDE** ← aadA1(17), aph(3')-Ia(13), aph(4)-Ia(12), aac(3)-IVa(12), aph(3'')-Ib(10)
- **BETA-LACTAM** ← blaTEM-1(11), blaCTX-M-65(11), blaCTX-M(2), blaCMY-2(1), blaNDM-1(1)
- **BLEOMYCIN** ← ble(3), bleO(2)
- **FOSFOMYCIN** ← fosA7(5), fosA3(5)
- **MACROLIDE** ← mph(E)(1)
- **MACROLIDE/STREPTOGRAMIN** ← msr(E)(1)
- **MULTIDRUG** ← ramR_Q11Ter(2), ramR_I21TerfsTer0(1), ramR_T18P(1), ramR_Q19Ter(1)
- **PHENICOL** ← floR(13), cmlA1(1), cmlA5(1)

## Honest caveats
- Cohorts are DRUG-R/S-SELECTED (resistance-enriched) -> co-occurrence reflects the curated cohorts + selection, NOT a random population sample. Linkage is 'within these cohorts'.
- Within-organism de-confounds SPECIES; --dedup-profiles is only a CRUDE clonality proxy (identical determinant sets). Full Mash-clonality correction is the follow-on (needs the assemblies + Docker).
- AMRFinder point-mutation determinants (gyrA_S83L) are chromosomal/organism-specific; acquired genes (sul1, tet(A)) are mobile/plasmid. Linkage mixes both mechanisms.
- A PASS (linkage exists) is the expected + valid finding — it demonstrates the world model captures real joint structure; a FAIL would mean within-organism conditional independence.
- Self-distillation from our own AMRFinder caller: associational, NOT a causal claim.