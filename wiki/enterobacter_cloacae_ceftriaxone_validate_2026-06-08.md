# Enterobacter_cloacae ceftriaxone — cross-organism validation — 2026-06-08

> Deployed dna-amr ceftriaxone rule applied UNCHANGED. AMRFinder `-O Enterobacter_cloacae`.
- NCBI group `Enterobacter_cloacae`; cohort 11 (8R/3S), 11 runs; `ncbi/amr:4.2.7-2026-03-24.1`
- **SMALL-N / LABEL-LIMITED:** ceftriaxone-S *E. cloacae* with downloadable assemblies are RARE (the
  organism is intrinsically AmpC, so most clinical isolates are cef-R) — only 3 S found. The spec estimate
  rests on 2-3 strains; treat as DIRECTIONAL, not a graded verdict. The S-scarcity is itself a finding.

## VERDICT: FAILS_BAR (FN-dominated — an EXPRESSION blind-spot, the 3rd boundary flavor)

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (ceftriaxone, EXTENDED-SPECTRUM subclass rule)** | 11 | **0.455** | **0.375** | 0.667 |
| naive (≥1 any β-lactamase) | 11 | 0.727 | 1.000 | 0.000 |

Unlike Acinetobacter (spec→0 over-call) and Campylobacter (threshold under-count), this fails on
**FN: 5 of 8 R strains carry NO qualifying ceftriaxone determinant** (n=0 under the rule).

## Root cause — AmpC DEREPRESSION (expression-level), invisible to gene-presence

*E. cloacae* ceftriaxone resistance is dominated by **derepression / overexpression of the intrinsic
chromosomal AmpC (ampC / blaACT)** via ampD/ampR regulatory mutations — an **expression-level** mechanism.
Gene-presence callers see the *gene*, not its *expression level*, so a derepressed-AmpC R strain looks
determinant-free. Per-strain:

| label | n (rule) | determinants | note |
|---|---:|---|---|
| R | 3 | blaCTX-M-15, blaNDM-1, blaOXA-1 | acquired ESBL/carbapenemase — **caught** |
| R | 2 | blaACT (×2) | AmpC reported as a determinant — **caught** |
| R | 1 | blaCTX-M-3 | acquired ESBL — **caught** |
| **R × 5** | **0** | **(none)** | **derepressed intrinsic AmpC — MISSED (FN)** |
| S | 0 | (none) | correctly susceptible |

So the rule catches acquired determinants (CTX-M, NDM, reported blaACT) but misses the 5 R strains whose
resistance is pure derepression — exactly the "what gene-presence can't see" signature, here on a β-lactam.

## The subclass-refinement tension (same shape as the Campylobacter threshold)

naive (≥1 *any* β-lactamase) gets sens 1.000 / spec 0.000 — it over-calls (counts every β-lactamase incl.
ampicillin-only blaTEM) but catches the AmpC-carriers the EXTENDED-SPECTRUM subclass filter drops. The
filter that is **correct for E. coli** (excludes blaTEM-1 = ampicillin-R, not cef-R) is **mis-tuned for
E. cloacae**, where the intrinsic AmpC is the main cef driver. This is the organism-specific-tuning pattern
again — the right cef rule for an AmpC-intrinsic organism differs from the Enterobacterales default.

## Boundary-type map (updated — now THREE flavors)

| organism | drug | failure | boundary type |
|---|---|---|---|
| Acinetobacter | meropenem | spec→0 over-call (intrinsic OXA-51) | **CONTENT** |
| Campylobacter | cipro | sens→0 under-count (single gyrA T86I vs thr 2) | **TUNING** |
| **Enterobacter cloacae** | **ceftriaxone** | **sens→0.375 FN (derepressed AmpC)** | **EXPRESSION** |

CONTENT = counts the wrong genes; TUNING = right genes, wrong integer; EXPRESSION = right genes, can't see
their regulation. EXPRESSION is the deepest blind-spot (shared with the Acinetobacter ISAba1/OXA-51 FN
ceiling) — gene-presence fundamentally cannot read derepression; closing it needs promoter/regulatory-region
analysis from the assembly (the `/hypothesise` IS-element-upstream strand generalizes to ampD/ampR).

## Honest scope / caveats
- 1 organism, 1 drug, **N=11 (8R/3S)** — small, S-scarce; spec rests on 2-3 strains. DIRECTIONAL only.
- The FN signal (5/8 R determinant-free) is the robust part; the exact acc/spec are N-fragile.
- NCBI labels (different source/curation, not a different-lab study).
- No refinement wired — an AmpC-aware cef rule for *E. cloacae* would need to (a) count intrinsic AmpC and
  (b) ideally read derepression, which presence can't do. Documented, not built.
