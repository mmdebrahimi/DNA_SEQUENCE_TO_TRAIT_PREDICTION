# HIV NRTI within-subtype transfer check (2026-06-21)

Catalog = position-based v0 (consensus-B numbering). Label = Stanford HIVDB PhenoSense fold (independent wet-lab; NOT Sierra). Filter = Method=PhenoSense AND Type=Clinical; N = 3164.

**Subtype mix:** B = 3099, non-B (pooled) = 49. The data is B-dominated -> non-B is under-powered.

| Drug | B n (R) sens/spec/**balacc** | non-B n (R) sens/spec/**balacc** |
|---|---|---|
| lamivudine | 3051 (2001) 1.0/0.562/**0.781** | 49 (15) 0.933/0.912/**0.923** |
| abacavir | 2948 (1983) 0.999/0.601/**0.8** | 42 (10) 0.9/0.969/**0.934** |
| zidovudine | 3075 (1488) 0.999/0.371/**0.685** | 49 (11) 1.0/0.842/**0.921** |
| stavudine | 3086 (1345) 0.991/0.331/**0.661** | 49 (12) 0.917/0.838/**0.877** |
| didanosine | 3087 (1578) 0.984/0.373/**0.678** | 49 (15) 0.933/0.912/**0.923** |
| tenofovir | 2736 (718) 0.997/0.264/**0.63** | 31 (8) 1.0/0.609/**0.804** |

## Honest caveats
- the free HIVDB gp data is ~96% subtype B -> non-B transfer is UNDER-POWERED (a free-data limit, reported)
- per-clade non-B N is tiny; only B-vs-pooled-non-B is adequately powered, and even non-B is small
- a similar B vs non-B sens/spec is consistent with transfer; it does NOT prove non-B generalisation at scale

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from DRMcv.R.