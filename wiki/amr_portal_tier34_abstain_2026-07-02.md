# AMR-Portal Tier 3/4 — ABSTAIN / non-endorsed cells (Phase D, 2026-07-02)

Honest documentation of the Tier-3/4 cells that a determinant caller should NOT score — either empirically
demonstrated (the spec-guard rejected them) or class-C/D by mechanism. Per the plan, **an honest ABSTAIN is
a valid outcome, not a gap.** The frozen surface + deployed card are untouched; these are overlay cells whose
headline is INDETERMINATE/UNDERPOWERED, or class-D cells not attempted.

## Empirically demonstrated over-call / non-decodable (spec-guard rejected — Phase B)

| cell | result | why (mechanism) |
|---|---|---|
| **Serratia marcescens cef** | spec **0.0**, sens 1.0 → INDETERMINATE | intrinsic **AmpC** (SME) — every isolate flagged R; the textbook class-C over-call the guardrail exists for |
| **Enterobacter cloacae cef** | 117R/6S → UNDERPOWERED (S=6) | intrinsic AmpC; near-all-R, too few S to even test spec |
| **Enterobacter cloacae mero** | spec 0.61 → INDETERMINATE | carbapenem R = porin-loss + AmpC/ESBL **expression** (class D) — not gene-presence-decodable |
| **Serratia marcescens cipro** | sens 0.0 → INDETERMINATE | the default rule's determinants don't capture Serratia FQ resistance (chromosomal) |
| **Enterobacter cloacae gent** | spec 0.81 (< 0.85) → INDETERMINATE | borderline over-call; needs acquired-only aminoglycoside curation |
| Klebsiella TMP-SMX (experimental) | strata not reproduced → INDETERMINATE | the sul-AND-dfr pattern doesn't hold in Klebsiella |

## Class-C/D by mechanism — NOT attempted (documented ABSTAIN)

- **Streptococcus pneumoniae β-lactams** (cef, mero): resistance is **PBP mutation** (penicillin-binding
  proteins), not gene-presence — a determinant caller cannot decode it. (`organism_rules/pneumo_betalactam.py`
  exists as prior scaffolding but PBP typing is out of the gene-presence paradigm.)
- **Pseudomonas aeruginosa** non-meropenem (cipro/gent): efflux (MexAB) + porin **expression**-driven — not
  gene-presence-decodable. (Meropenem is already `EXPRESSION_FLOOR`/ABSTAIN in the deployed surface.)
- **Acinetobacter, Pseudomonas meropenem**: intrinsic OXA / oprD-loss carbapenem — already ABSTAINS_BY_DESIGN
  (Tier 2 of the triage).
- **Clostridioides difficile cef** / other anaerobes: sparse determinant curation; not attempted.

## Disposition

Tier-3/4 SCORED cells (Phase A+B+C + Tier-3 Campylobacter): **17 new** (NG cipro/tet, SA cipro/rif,
E. cloacae cipro/tet/TMP-SMX, Proteus gent/TMP-SMX, Serratia mero, E. faecium cipro/tet/gent,
C. jejuni tet/gent, C. coli tet/gent) + the 4 TMP-SMX experimental = **21 overlay SCORED**. The remaining
backlog is the ABSTAIN set above — a correct, honest terminal (the determinant paradigm genuinely can't
decode PBP/efflux/intrinsic-expression resistance). Any future attempt needs a different signal
(expression/regulatory), not another gene-presence rule.

Note on Campylobacter gent: the initial probe showed the AMINOGLYCOSIDE subclass dominated by aad9/spw
(spectinomycin/streptomycin, non-gent) — but the gent-R isolates DO carry true aph(2'')/aac(3) enzymes,
which the `_GENT_MARKERS`-only rule separates cleanly (spec 1.0 both species). The non-gent-marker exclusion
is the load-bearing piece (`tests/test_campylobacter_amr.py::test_gent_true_marker_only`).
