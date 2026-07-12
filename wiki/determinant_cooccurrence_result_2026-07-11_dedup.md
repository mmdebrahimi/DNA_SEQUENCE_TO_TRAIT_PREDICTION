# Determinant co-occurrence / linkage world model + phenotype->genotype inverter (2026-07-11)

**Verdict: PASS_LINKAGE_STRUCTURE** — 114/121 within-organism testable determinants are LINKED (fraction 0.942; PASS bar 0.5). (profile-deduped clonality proxy)

Is there within-organism co-resistance LINKAGE (do a genome's other determinants predict a held-out determinant beyond base rates)? + the phenotype->genotype inversion.

Substrate: cached AMRFinder determinant calls (self-distillation). LINKED = a determinant's held-out presence is predicted from the OTHER determinants, WITHIN organism, with OOF-AUC 95% CI lower > 0.5.

| organism | genomes | testable | **linked** | frac | note |
|---|---|---|---|---|---|
| klebsiella | 235 | 58 | **57** | 0.983 | profile-deduped 307->235 genomes (clonality proxy) |
| escherichia_coli_shigella | 163 | 42 | **39** | 0.929 | profile-deduped 240->163 genomes (clonality proxy) |
| campylobacter | 37 | 3 | **0** | 0.0 | profile-deduped 100->37 genomes (clonality proxy) |
| acinetobacter | 36 | 13 | **13** | 1.0 | profile-deduped 60->36 genomes (clonality proxy) |
| salmonella | 28 | 5 | **5** | 1.0 | profile-deduped 60->28 genomes (clonality proxy) |

## C3a — 'vice-versa' co-occurrence LIFT (what determinants travel together; verify-in-batch)
### klebsiella
- **fosA** → aac(3)-IVa(lift 1.08, n17), aac(6')-Ib(lift 1.08, n45), aac(6')-Ib'(lift 1.08, n8), aph(4)-Ia(lift 1.08, n17)
- **oqxA** → aac(3)-IVa(lift 1.21, n17), ampC-Kaer(lift 1.21, n10), aph(4)-Ia(lift 1.21, n17), blaOXA-48(lift 1.21, n8)
- **oqxB** → aac(3)-IVa(lift 1.73, n17), ampC-Kaer(lift 1.73, n10), aph(4)-Ia(lift 1.73, n17), ble(lift 1.73, n9)
- **parC_S80I** → aac(3)-IVa(lift 2.53, n17), aph(4)-Ia(lift 2.53, n17), ompK35_E42RfsTer47(lift 2.53, n54), gyrA_S83I(lift 2.46, n76)
- **sul1** → arr-3(lift 2.61, n10), dfrA12(lift 2.61, n48), aadA2(lift 2.48, n56), ant(2'')-Ia(lift 2.32, n8)
- **blaTEM-1** → arr-3(lift 2.8, n10), qnrB1(lift 2.5, n25), aac(3)-IIe(lift 2.35, n21), tet(A)(lift 2.07, n34)
### escherichia_coli_shigella
- **gyrA_S83L** → gyrA_D87N(lift 1.5, n92), parC_E84V(lift 1.5, n36), parE_L416F(lift 1.5, n10), parE_S458A(lift 1.5, n31)
- **parC_S80I** → aac(6')-Ib-cr5(lift 1.73, n43), parC_E84V(lift 1.73, n36), parE_S458A(lift 1.73, n31), gyrA_D87N(lift 1.72, n91)
- **gyrA_D87N** → parC_E84V(lift 1.77, n36), parE_L416F(lift 1.77, n10), parE_S458A(lift 1.77, n31), aac(6')-Ib-cr5(lift 1.73, n42)
- **blaTEM-1** → aac(3)-IId(lift 2.16, n16), dfrA14(lift 2.0, n20), aph(3')-Ia(lift 1.99, n13), floR(lift 1.68, n11)
- **sul1** → aadA5(lift 2.23, n50), dfrA17(lift 2.11, n51), aac(3)-IId(lift 2.08, n15), mph(A)(lift 1.94, n51)
- **sul2** → aph(3'')-Ib(lift 2.08, n60), aph(6)-Id(lift 2.04, n58), dfrA14(lift 1.95, n19), sat2(lift 1.89, n8)
### campylobacter
- **blaOXA-193** → blaOXA-61_G-57T(lift 1.64, n8), gyrA_T86I(lift 1.09, n10), tet(O)(lift 1.07, n11)
- **tet(O)** → gyrA_T86I(lift 1.15, n10), blaOXA-193(lift 1.07, n11)
- **gyrA_T86I** → tet(O)(lift 1.15, n10), blaOXA-193(lift 1.09, n10)
### acinetobacter
- **parC_S84L** → blaOXA-66(lift 1.44, n8), aph(3')-VIa(lift 1.33, n12), mph(E)(lift 1.3, n9), msr(E)(lift 1.3, n9)
- **sul2** → ant(2'')-Ia(lift 1.17, n14), aph(6)-Id(lift 1.15, n13), aph(3')-Ia(lift 1.1, n11), tet(B)(lift 1.09, n8)
- **ant(2'')-Ia** → aph(3')-VIa(lift 1.85, n12), blaOXA-69(lift 1.54, n10), sul2(lift 1.17, n14), aph(3')-Ia(lift 1.07, n8)
- **aph(3'')-Ib** → blaOXA-66(lift 2.12, n8), aph(6)-Id(lift 1.99, n16), tet(B)(lift 1.73, n9), parC_S84L(lift 1.1, n13)
- **aph(6)-Id** → blaOXA-66(lift 2.12, n8), aph(3'')-Ib(lift 1.99, n16), tet(B)(lift 1.93, n10), sul2(lift 1.15, n13)
- **aph(3')-Ia** → blaOXA-69(lift 1.66, n9), blaOXA-23(lift 1.6, n10), parC_S84L(lift 1.15, n12), sul2(lift 1.1, n11)
### salmonella
- **tet(A)** → sul1(lift 1.71, n11), aadA1(lift 1.68, n9), aph(3')-Ia(lift 1.66, n8), gyrA_D87Y(lift 1.56, n10)
- **gyrA_D87Y** → aadA1(lift 2.1, n9), sul1(lift 1.94, n10), tet(A)(lift 1.56, n10)
- **sul1** → aadA1(lift 2.1, n9), gyrA_D87Y(lift 1.94, n10), tet(A)(lift 1.71, n11)
- **aadA1** → gyrA_D87Y(lift 2.1, n9), sul1(lift 2.1, n9), tet(A)(lift 1.68, n9)

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
- **BETA-LACTAM** ← blaTEM-1(11), blaCTX-M-65(11), blaCTX-M(2), blaCMY-2(1), blaDHA(1)
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