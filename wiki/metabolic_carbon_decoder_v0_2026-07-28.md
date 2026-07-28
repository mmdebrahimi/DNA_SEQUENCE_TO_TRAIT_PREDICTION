# Metabolic carbon-utilization decoder — v0 (2026-07-28)

`dna-decode metabolic` — E. coli carbon-source utilization from gene/operon presence. The AMR
determinant→phenotype paradigm applied to metabolism, plus the one twist metabolism adds that resistance
does not: **uptake-gating**.

## The rule

```
utilizes  iff  (all catabolic enzymes present) AND (a transporter present) AND (transporter expressed under the O2 condition)
```

A cell can catabolize a sugar iff it can (1) IMPORT it, (2) break it down, AND (3) the importer is expressed
under the queried aerobic/anaerobic condition. The naive AMR-style rule ("has the pathway genes → can use
it") is right for most sugars but wrong for the case metabolism is famous for.

## The differentiator — the citrate anchor (why this is not just the Nth R1 confirmation)

E. coli K-12 is **Cit⁻ aerobically**. It carries the full TCA cycle (it metabolizes citrate as an internal
intermediate constantly) **and** it carries the `citT` citrate/succinate antiporter — yet on an aerobic
citrate plate it cannot grow, because `citT` (the `citCDEFXGT` operon) is expressed **only anaerobically**
(CitAB two-component + anaerobiosis). A "has citrate genes → Cit+" rule says **+**; the measured aerobic
phenotype is **−**. The famous Lenski LTEE Cit+ mutants evolved *aerobic* `citT` expression via the
`rnk-citG` regulatory duplication (Blount et al. 2012 Nature) — i.e. they fixed exactly the expression gate
this rule models, not a new enzyme.

Real CLI:

```
$ dna-decode metabolic --source citrate --genes citD,citE,citF,citT --condition aerobic
citrate (aerobic): CANNOT UTILIZE  (confidence medium)
  - ...transporter is not expressed aerobically (expressed: anaerobic) -> uptake gate CLOSED...
  - a naive 'has the pathway genes -> can use it' rule MIS-CALLS this positive

$ dna-decode metabolic --source citrate --genes citD,citE,citF,citT --condition anaerobic
citrate (anaerobic): UTILIZES  (confidence high)
```

This is the metabolic analog of the flowering cell's Da(1)-12 anchor — the literature-sourced case a naive
rule gets wrong, pinned by `reference_integrity_ok()` as the anti-fabrication guard.

## Validation (measured E. coli K-12 MG1655 phenotypes; EcoCyc / Neidhardt)

| substrate | genes present | condition | call | measured K-12 | ✓ |
|---|---|---|---|---|---|
| lactose | lacZ,lacY | aerobic | utilizes | Lac+ | ✓ |
| lactose | lacY (lacZ⁻) | aerobic | cannot | Lac⁻ (enzyme KO) | ✓ |
| lactose | lacZ (lacY⁻) | aerobic | cannot | Lac⁻ (uptake KO) | ✓ |
| L-arabinose | araABD,araE | aerobic | utilizes | Ara+ | ✓ |
| D-glucose | pgi,ptsG | aerobic | utilizes | Glc+ | ✓ |
| **citrate** | **citDEF,citT** | **aerobic** | **cannot** | **Cit⁻** (naive rule says +) | ✓ |
| **citrate** | **citDEF,citT** | **anaerobic** | **utilizes** | **Cit+** | ✓ |

Catalog: lactose · L-arabinose · maltose · D-xylose · L-rhamnose · D-glucose · citrate (each entry sourced
in `dna_decode/metabolic/carbon_catalog.py`). Tier: **KNOWLEDGE_BASELINE**.

## Honest scope (load-bearing)

- **v0 = E. coli carbon catabolism only** (the cleanest curated determinant→phenotype metabolic map). N/S
  sources, auxotrophies, cross-organism transfer are out.
- Calls the **can/cannot DIRECTION** (aerobic/anaerobic), **not** growth rate / yield / lag.
- **Reads gene PRESENCE, not sequence integrity** — it cannot see a point mutation that silently inactivates
  a present gene (a v0.1 genome-mode + sequence follow-on, deliberately not fabricated here).
- Faithful-to-literature: applies published operon/transporter assignments; not a new model.

## Meta-note (the honest plateau signal)

This is the **8th** deterministic curated-catalog cell (bacterial/viral/fungal AMR → PGx → visible traits →
plant flowering → phage receptor → essentiality → metabolic). The R1 "curated-catalog → deterministic wins"
paradigm is now overwhelmingly established across phenotype domains and kingdoms. This cell earns its place
by encoding a **new rule shape** (uptake-gating, which the frozen count/OR AMR engine cannot express and a
naive has-the-genes rule gets wrong), not by re-confirming the paradigm. Further R1 cells are re-confirmation
motion; the genuine frontier is elsewhere (independent-label acquisition = authority/money fork; sequence-
integrity genome mode = a v0.1 build).

## Files

- `dna_decode/metabolic/carbon_catalog.py` — curated catalog + `call_carbon_utilization` + `reference_integrity_ok`
- `dna_decode/metabolic/cli.py` — `dna-decode metabolic` / `dna-metabolic`
- `tests/test_metabolic_carbon.py` (14) + `tests/test_metabolic_cli.py` (6)
- Wired: `cli.py` TRAITS+dispatch · `pyproject` `dna-metabolic` · `cell_registry` CellContract · `test_cli_dispatch` pin

Frozen AMR + forward surfaces byte-unchanged throughout. Offline, $0, deterministic.
