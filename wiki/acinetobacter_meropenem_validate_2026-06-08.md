# Acinetobacter meropenem — cross-organism validation — 2026-06-08

> Deployed dna-amr meropenem rule applied UNCHANGED. AMRFinder `-O Acinetobacter_baumannii`.
- NCBI group `Acinetobacter`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`
- First organism BEYOND the original Phase-1 four (E. coli, Klebsiella, Enterobacter, Salmonella).

## VERDICT: FAILS_BAR (generic rule does NOT transfer to Acinetobacter)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (meropenem, generic CARBAPENEM-class rule)** | 30 | **0.500** | 1.000 | **0.000** |
| naive AMRFinder | 30 | 0.500 | 1.000 | 0.000 |

The generic rule (threshold 1 + any CARBAPENEM-subclass determinant) calls **all 30 strains R** → no
better than calling everything R (naive = same). This is the **first organism where the deployed rule
fails to transfer out-of-the-box** (it held on the 4 Enterobacterales of Phase 1).

## Root cause — intrinsic OXA-51-family over-call

A. baumannii carries a **chromosomal, intrinsic blaOXA-51-like carbapenemase in EVERY isolate**
(OXA-66/68/69/100/312/113 are all OXA-51-family alleles). The broad CARBAPENEM-subclass count picks it
up in the **susceptible** strains too → spec collapses to 0. This is the same class of gotcha as the
Klebsiella OqxAB intrinsic-efflux over-call: **a broad drug-class bag counts intrinsic genes that don't
confer the phenotype.**

Per-strain carbapenemase content (all 30):
- **OXA-23** (strong, acquired): 10 strains — **all R** (true carbapenemase).
- **OXA-72** (OXA-24/40-like, strong, acquired): 1 strain — R.
- **OXA-58-like** (conditional/promoter-dependent): carried by **9 S strains AND 2 R strains** —
  presence alone does NOT predict phenotype (needs an ISAba2/ISAba3 upstream promoter to overexpress).
- **Intrinsic OXA-51-family**: in all 30 incl. every S strain — R only under ISAba1-driven
  overexpression, which gene-presence cannot see.

## Candidate refinement — literature-grounded carbapenemase-strength tiers (in-cohort)

Restrict the meropenem call to **strong acquired carbapenemases** (OXA-23-like, OXA-24/40-like incl.
OXA-72, OXA-143/235-like, NDM/IMP/VIM MBLs, KPC) — **excluding** intrinsic OXA-51-family AND the
conditional OXA-58-like. The strength hierarchy is established Acinetobacter AMR biology, **not fit to
these 30 labels**.

| rule | acc | sens | spec | FP | FN |
|---|---:|---:|---:|---:|---:|
| generic CARBAPENEM-class (deployed) | 0.500 | 1.000 | 0.000 | 15 | 0 |
| require any acquired carbapenemase (excl. OXA-51-fam) | 0.600 | 0.800 | 0.400 | 9 | 3 |
| **require STRONG acquired (excl. OXA-51-fam + OXA-58)** | **0.833** | 0.667 | **1.000** | **0** | 5 |

The strong-acquired rule: **0 false positives** (every strong-carbapenemase call is a true R), 5 FN.

### The 5 residual FN are a fundamental gene-presence blind spot, not a curation gap
- 3 R strains carry **only intrinsic OXA-51-family** (OXA-113 / OXA-312×2) → R via ISAba1-driven
  overexpression (an **expression-level**, not gene-presence, mechanism).
- 2 R strains carry **only OXA-58-like** → R only with an upstream promoter insertion (same class).

Gene-presence callers (AMRFinder, ResFinder, this decoder) **cannot** see IS-element-driven
overexpression of an intrinsic/weak gene. This is the honest, named ceiling for presence-based AMR on
Acinetobacter carbapenems — not a fixable rule gap.

## Honest scope / caveats
- 1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).
- The strength-tier refinement is evaluated **in-cohort** (same 30 strains). The tiers come from
  literature, but the acc/spec numbers need an **independent Acinetobacter cohort** before the
  refinement is productionized (per the validate-wrapper-vs-underlying-tool discipline).
- Refinement is a **candidate**, NOT wired into the deployed rule. Productionizing it = an
  organism-specific carbapenemase-strength catalog + external validation (a scoped follow-on).

## Takeaway
The transferability map now reads: the generic AMR-class rules transfer across **Enterobacterales**
(E. coli/Klebsiella/Enterobacter/Salmonella) but **NOT to Acinetobacter** without organism-specific
carbapenemase-strength curation — and even curated, presence-based calling has a principled ~33% FN
ceiling on Acinetobacter carbapenems (IS-element-driven expression mechanisms).
