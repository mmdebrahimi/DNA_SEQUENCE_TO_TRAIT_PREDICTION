# HIV epistasis vs additive — does interaction structure beat the additive rule? (2026-07-11)

**Verdict: FAIL_ADDITIVE_SUFFICES** — 2/24 powered drug-cells show a bootstrap-CI-positive interaction gain (fraction 0.083; PASS bar = 0.5).

Label = Stanford HIVDB PhenoSense fold-change (free, independent wet-lab; NOT Sierra interpretation). Does a constrained pairwise-interaction model beat the additive mutation-presence baseline out-of-sample on Stanford HIVDB PhenoSense fold-change?

Both models = ElasticNetCV in a nested 5-fold OOF harness on the SAME folds (fair paired comparison). `beat` = paired bootstrap delta_rho (interaction-additive) 95%-CI lower bound > 0.

| class:drug | n | add rho | int rho | **d_rho [95% CI]** | CI+ | add/int R2 |
|---|---|---|---|---|---|---|
| NNRTI:EFV | 2168 | 0.9108 | 0.9094 | **-0.0014 [-0.0038, 0.001]** | no | 0.8443/0.8426 |
| NNRTI:NVP | 2052 | 0.873 | 0.862 | **-0.011 [-0.0157, -0.0066]** | no | 0.7746/0.7899 |
| NNRTI:ETR | 998 | 0.8338 | 0.8303 | **-0.0035 [-0.0105, 0.0026]** | no | 0.7273/0.7216 |
| NNRTI:RPV | 311 | 0.7486 | 0.7415 | **-0.0071 [-0.0312, 0.0172]** | no | 0.642/0.6093 |
| NNRTI:DOR | 130 | 0.7841 | 0.7831 | **-0.001 [-0.0156, 0.011]** | no | 0.8131/0.8128 |
| NNRTI:CompMutList | 0 | — | — | too few isolates | — | — |
| NRTI:3TC | 1839 | 0.8957 | 0.8962 | **0.0005 [-0.0017, 0.0026]** | no | 0.9098/0.9212 |
| NRTI:ABC | 1731 | 0.8988 | 0.9055 | **0.0067 [0.0031, 0.0104]** | YES | 0.8181/0.8396 |
| NRTI:AZT | 1853 | 0.8822 | 0.8844 | **0.0021 [-0.0016, 0.0058]** | no | 0.7837/0.7987 |
| NRTI:D4T | 1846 | 0.8532 | 0.8508 | **-0.0024 [-0.0074, 0.0025]** | no | 0.7296/0.7276 |
| NRTI:DDI | 1849 | 0.8075 | 0.808 | **0.0005 [-0.0025, 0.0037]** | no | 0.7136/0.7114 |
| NRTI:TDF | 1548 | 0.7981 | 0.7991 | **0.0009 [-0.0073, 0.0096]** | no | 0.6177/0.6275 |
| NRTI:CompMutList | 0 | — | — | too few isolates | — | — |
| PI:FPV | 2052 | 0.9001 | 0.9006 | **0.0005 [-0.0023, 0.0034]** | no | 0.8628/0.8665 |
| PI:ATV | 1505 | 0.9125 | 0.9115 | **-0.0011 [-0.0039, 0.0015]** | no | 0.8825/0.8804 |
| PI:IDV | 2098 | 0.9123 | 0.9137 | **0.0015 [-0.0009, 0.0039]** | no | 0.8814/0.8833 |
| PI:LPV | 1807 | 0.9106 | 0.9107 | **0.0002 [-0.0021, 0.0022]** | no | 0.9086/0.9093 |
| PI:NFV | 2133 | 0.9087 | 0.9112 | **0.0025 [0.0002, 0.005]** | YES | 0.8676/0.8709 |
| PI:SQV | 2084 | 0.8929 | 0.8937 | **0.0008 [-0.0028, 0.0042]** | no | 0.873/0.8698 |
| PI:TPV | 1226 | 0.7702 | 0.7578 | **-0.0124 [-0.0269, 0.0028]** | no | 0.7203/0.7058 |
| PI:DRV | 993 | 0.8205 | 0.8052 | **-0.0153 [-0.0237, -0.0063]** | no | 0.8684/0.8616 |
| PI:CompMutList | 0 | — | — | too few isolates | — | — |
| INI:RAL | 753 | 0.8426 | 0.8516 | **0.009 [-0.0014, 0.0205]** | no | 0.896/0.9064 |
| INI:EVG | 754 | 0.8801 | 0.8825 | **0.0024 [-0.0066, 0.0112]** | no | 0.8223/0.8316 |
| INI:DTG | 370 | 0.5764 | 0.5703 | **-0.0062 [-0.0511, 0.0386]** | no | 0.4594/0.4156 |
| INI:BIC | 287 | 0.743 | 0.7454 | **0.0024 [-0.017, 0.0222]** | no | 0.4655/0.4796 |
| INI:CAB | 64 | 0.7396 | 0.6682 | **-0.0713 [-0.1417, -0.0155]** | no | 0.6137/0.5543 |
| INI:CompMutList | 0 | — | — | too few isolates | — | — |

## Verify-in-batch — top interaction terms (biology sanity check)
- **NNRTI:EFV** — 41L:215Y(-0.0277), 135T:70R(+0.0181), 200A:179I(-0.0175), 215Y:118I(-0.0172)
- **NNRTI:NVP** — 103N:181C(-0.15), 293V:228H(-0.0226), 184V:60I(+0.0222), 122E:245E(+0.0216)
- **NNRTI:ETR** — 184V:215Y(-0.0227), 184V:202V(-0.0187), 200A:228H(+0.0175), 272P:35I(-0.0173)
- **NNRTI:RPV** — 219Q:228H(+0.0527), 286A:118I(-0.0486), 286A:70R(+0.0316), 215Y:181C(-0.0257)
- **NRTI:3TC** — 184V:41L(-0.0402), 184V:67N(-0.0342), 184V:210W(-0.0338), 184V:215Y(-0.0314)
- **NRTI:ABC** — 184V:215Y(-0.0827), 70R:219Q(+0.0182), 215Y:74V(-0.0171), 41L:215Y(+0.0168)
- **NRTI:AZT** — 215Y:210W(+0.1125), 41L:215Y(+0.0855), 70R:219Q(+0.0829), 184V:215Y(-0.0687)
- **NRTI:D4T** — 196E:208Y(+0.0156), 70R:219Q(+0.0147), 41L:215Y(+0.0136), 219Q:215F(-0.012)
- **NRTI:DDI** — 215Y:210W(+0.01), 41L:215Y(+0.0079), 215Y:44D(+0.0064), 215Y:74V(-0.0064)
- **NRTI:TDF** — 70R:219Q(+0.0359), 184V:215Y(-0.0343), 135T:215F(+0.0238), 219Q:215F(-0.0233)
- **PI:FPV** — 46I:33F(-0.0311), 54V:15V(+0.0267), 36I:46L(+0.0216), 84V:32I(-0.0179)
- **PI:ATV** — 84V:10F(+0.0408), 36I:54V(+0.0357), 36I:46I(+0.0216), 63P:10I(+0.0199)
- **PI:IDV** — 82A:84V(-0.0436), 46I:54V(-0.0311), 84V:37D(+0.0204), 90M:20I(+0.0202)
- **PI:LPV** — 82A:84V(-0.0287), 46I:54V(-0.0256), 36I:35D(+0.0164), 71V:73S(+0.0154)
- **PI:NFV** — 36I:35D(+0.0353), 46I:33F(-0.035), 82A:84V(-0.0254), 90M:82A(-0.0225)
- **PI:SQV** — 33F:32I(-0.0323), 36I:35D(+0.0288), 54V:20R(+0.0251), 46I:82A(-0.0194)
- **PI:TPV** — 54V:84V(+0.0377), 35D:84V(+0.0355), 63P:84V(+0.0296), 84V:10F(-0.0239)
- **PI:DRV** — 71V:32I(+0.0271), 33F:46L(+0.0253), 54V:15V(+0.0246), 63P:32I(+0.0227)
- **INI:RAL** — 140S:148H(+0.1995), 140S:148R(+0.0493), 72V:97A(+0.0338), 148R:138K(+0.0305)
- **INI:EVG** — 140S:148H(+0.2014), 101I:155H(+0.0469), 72V:138K(-0.0373), 148R:138K(+0.0361)
- **INI:DTG** — 140S:148H(+0.1704), 124A:140S(+0.0462), 72V:138K(-0.04), 201I:140S(-0.0344)
- **INI:BIC** — 140S:148H(+0.0439), 148H:138K(-0.0336), 140S:138K(-0.0336), 124A:101I(-0.0224)
- **INI:CAB** — 140S:148H(+0.261), 148R:138K(+0.0798)

## Interpretation (how to read a FAIL)
A FAIL means the ADDITIVE mutation-presence model is out-of-sample-competitive with the explicit interaction model for RANK-prediction of fold-change — NOT that epistasis is absent. Read the top-interaction terms above: if they are known synergy pairs recovered with LARGE coefficients (e.g. INI G140S+Q148H `140S:148H`, NRTI TAMs `41L:215Y`/`215Y:210W`, PI `82A:84V`/`46I:54V`), the epistasis is REAL and correctly localized — it just does not improve the population rank-metric, because the additive main effects already rank double-mutants high (Spearman is rank-based) and the genuine primary-weak+accessory-weak synergy cases are rare in the population. The world-model implication: the curated ADDITIVE catalog wins in the quantitative/continuous regime too, extending the deterministic-decoder thesis beyond binary R/S. Interactions can even HURT on small cells (variance/overfit) — see any negative delta_rho.

## Honest caveats
- A mechanism-FEATURE interaction model, NOT a sequence embedding (the closed 0-for-5 arm).
- Censored folds ('>'/'<') kept at the numeric bound: rank-metric-safe, slightly biases R2.
- Paired bootstrap on FIXED OOF predictions is the standard held-out paired test (mildly optimistic; disclosed).
- Interactions constrained to top-K prevalent co-occurring (>=MIN_COOC) pairs + L1 — not the full O(p^2) space.
- A FAIL ('additive suffices for HIV DR') is a valid world-model finding, not a bug.

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset per HIVDB Terms of Use.