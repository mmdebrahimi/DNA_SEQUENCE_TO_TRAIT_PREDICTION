# EXPRESSION-frontier falsifier: ISAba1→OXA-51 junction — 2026-06-10

> Tests the project's documented #1 limitation: the bacterial AMR decoder's **EXPRESSION abstain floor**
> (Acinetobacter carbapenem — gene-PRESENCE cannot see ISAba1-driven blaOXA-51 overexpression, the dominant
> false-negative cause). Hypothesis (/hypothesise strand #3): detecting an **ISAba1 insertion immediately
> upstream of blaOXA-51-family** in the assembly recovers the intrinsic-only-R strains gene-presence misses.
> Pure blastn over the 30 cached Acinetobacter genome assemblies; no money, no Docker. Falsifier script +
> raw JSON: `soraya_runs/2026-06-10-7qjq-expression-frontier-isaba1/`. Refs: `data/isaba1_ref/` (ISAba1 +
> OXA-51-family, both extracted from one GenBank junction record EU029998.1).

## VERDICT: SURVIVES — specific, but a PARTIAL floor-crossing

Method: blastn ISAba1 ref + OXA-51-family ref vs each genome; an OXA-51 hit is "junction-positive" if an
ISAba1 hit lies within 400 bp upstream (strand-aware) on the same contig. Cross-tabbed vs R/S + vs
strong-acquired-carbapenemase presence (the FN set = R strains with NO strong acquired carbapenemase).

| metric | value |
|---|---|
| junction-positive in **R** | **1 / 15** |
| junction-positive in **S** | **0 / 15** (zero false positives) |
| intrinsic-only-R (the FN ceiling) | 5 |
| of those, junction-positive (recovered) | **1 / 5** (`GCA_000692095.1`, 22 ISAba1 copies, one upstream of OXA-51) |
| false-rescue risk (S, junction-positive) | **0** |

## Reading — the floor is PARTIALLY crossable

1. **The expression-context signal is REAL and perfectly specific here.** When ISAba1 sits upstream of
   OXA-51, the strain is R — never S in this cohort (0/15 S). So reading *regulatory context* (IS position),
   not gene presence, genuinely crosses the abstain floor: it correctly upgrades 1 ABSTAIN→R with **zero**
   new false positives. This is the first demonstration that the floor is not absolute.
2. **But it only explains 1 of 5 false-negatives.** The other 4 intrinsic-only-R strains each carry a single
   ISAba1 copy NOT near OXA-51 → their R is via a *different* expression mechanism (efflux up-regulation /
   other promoter), still invisible to any sequence-presence/position signal. So the floor is *partially*
   crossable: IS-position rescues a minority; most expression-driven R remains uncrossable.
3. `GCA_000692095.1` (the rescue) has **22 ISAba1 copies** — a strain with massive IS expansion, consistent
   with the known ISAba1-amplification → OXA-51-overexpression mechanism.

## What this establishes / honest scope

- **Establishes:** a deterministic, presence-independent **expression-context feature** (ISAba1-upstream-of-
  target) is a valid, specific signal that can convert a known ABSTAIN into a correct R — quantifying that
  the EXPRESSION floor is real but ~20% crossable (1/5) by IS-position alone on this cohort.
- **Does NOT establish:** a production feature. N=30 (15R/15S), single cohort; ISAba1 ref is a 570 bp
  partial (the conserved transposase + IR — sufficient here since the 4 non-rescued FN strains have ISAba1
  *not near* OXA, so a fuller ref wouldn't manufacture absent junctions); 400 bp window + identity
  thresholds are tunable. Wiring ISAba1-junction → ABSTAIN-override into the deployed caller needs a fuller
  ISAba1 reference + an independent Acinetobacter cohort first (the validate-wrapper discipline).
- **Generalizes:** the same IS-element-upstream-of-target pattern is the mechanism for other intrinsic-gene
  overexpression (e.g. ISAba1/ISAba125 upstream of ampC/ADC; IS upstream of efflux regulators). This
  falsifier is the template for an `expression_context` signal class layered on the presence-based decoder.

## Bottom line

The bacterial decoder's hardest limitation (the carbapenem abstain floor) is **partially** addressable by
reading IS-element regulatory context — specifically (0 FP) but with low sensitivity (1/5 FN). Worth a
future `expression_context` feature with a fuller IS reference + independent validation; not worth wiring
into the deployed rule on N=30. The honest contribution: we now *quantify* how much of the floor is
sequence-visible (a minority) vs genuinely expression-only (the majority).
