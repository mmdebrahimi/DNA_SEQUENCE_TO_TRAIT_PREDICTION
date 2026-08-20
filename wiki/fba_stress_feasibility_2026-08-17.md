# The stress axis is a structural NO-GO for iML1515 — and the reason is worth more than the panel

**Date:** 2026-08-17 · **Verdict:** `STRESS_AXIS_NOT_REPRESENTABLE_IN_iML1515`
**Artifact:** `wiki/fba_stress_feasibility_2026-08-17.json` · **Probe:** `scripts/fba_stress_feasibility_probe.py`

The nitrogen axis replicated all three carbon findings, so stress was the named next axis — 55 experiments
over 35 conditions, and a different perturbation *class* rather than just a different nutrient. The
checkpoint predicted exchange-mappability would be the blocker. It is, but the interesting part is *which
layer* fails.

## Three layers, each stricter

| layer | question | result |
|---|---|---|
| **L1** | does an exchange exist at all? | **8 of 35** (23 %) |
| **L2** | does adding it **reduce** growth? | **0 of 8** |
| **L3** | is the molecular target in the model? | **2 of 13** antibiotics |

**L1** is the shallow blocker. 27 of 35 conditions are antibiotics, ionic liquids, detergents or
cytotoxics — `1-ethyl-3-methylimidazolium chloride`, `MreB Perturbing Compound A22`, `Cisplatin`. These
are not metabolites, so iML1515 has no exchange for them and never could.

**L2 is the real finding, and it is not a counting problem.** The 8 that *do* map fail for a deeper
reason: an exchange models **supplementation**, which is the opposite of a stress.

| compound | exchange | Δ growth | |
|---|---|---|---|
| Sodium acetate | `EX_ac_e` | **+0.2445** | it is an extra carbon source |
| L-Lysine | `EX_lys__L_e` | **+0.0347** | it is an amino acid |
| Sodium Chloride | `EX_na1_e` | 0.0000 | inert |
| Sodium nitrite | `EX_no2_e` | 0.0000 | inert |
| Dimethyl Sulfoxide | `EX_dmso_e` | 0.0000 | inert |
| Cobalt / Nickel / Copper chloride | `EX_{cobalt2,ni2,cu2}_e` | 0.0000 | trace nutrients, already unlimited |

Zero reduce growth. Acetate — a *stressor* in the assay at high concentration — makes the model grow
**28 % faster**. So the medium-swap contract that carries the carbon and nitrogen axes is not merely
underpowered for stress; it has the **wrong sign**.

The cause is structural, not a modelling oversight: **FBA is stoichiometric.** It represents what a
network can produce from what it consumes. It has no representation of enzyme inhibition kinetics,
ribosome binding, membrane disruption or DNA damage — which is what every compound in this panel actually
does.

**L3 is the one narrow path, and it is honest about its size.** Exactly 2 of 13 panel antibiotics have a
*metabolic* target present in iML1515:

- **Phosphomycin** → MurA (`gene:murA`, `rxn:UAGCVT`)
- **D-Cycloserine** → alanine racemase / D-Ala-D-Ala ligase (`gene:alr`, `gene:ddlA`, `gene:ddlB`, `rxn:ALAALAr`)

Both act on peptidoglycan synthesis, which iML1515 does carry. They could in principle be modelled as
**target-directed reaction constraints** — a different experimental contract from the medium swap. The
other 11 target the ribosome, gyrase, PBP transpeptidation or the cytoskeleton: outside a metabolic model
**by construction**, not by omission.

**n = 2 is not a panel.** Carbon ran 25 conditions, nitrogen 13. Two conditions cannot support the
conditional-essentiality contract (a gene must be essential in ≥1 and dispensable in ≥1), and any
"replication" claim from it would be noise. Recorded as a real but unpursued path rather than a next step.

## What this changes

The three carbon findings replicate across **substrate** axes (carbon → nitrogen) but the method has a
hard boundary at the **perturbation class**. That boundary is a property of the model formalism, not of
the panel or the pipeline:

> **A stoichiometric model can vary what the cell is fed. It cannot represent something that poisons the
> machinery.**

That is a cleaner statement of the method's scope than another replication would have produced, and it
comes from an eight-minute probe rather than a full panel run.

## A wrong mapping I caught in my own draft

My first candidate map paired **`sodium fluoride` → `EX_fe2_e`**. `EX_fe2_e` is **ferrous iron**, not
fluoride. Had it survived, the probe would have "measured" an unrelated metabolite and manufactured a data
point. It is dropped rather than guessed at (iML1515 has no fluoride exchange), and
`test_fluoride_is_not_mapped_to_iron` pins it. `test_candidate_exchanges_are_all_metabolites_not_antibiotics`
pins the general form — no antibiotic may acquire a medium mapping.

## Next

Stress is closed as a **structural** wall — not code-closable; it would need a different model formalism
(kinetic / ME-model / whole-cell), which is a far larger undertaking than an axis.

The remaining named axis is **cross-organism** conditional essentiality: 47 organisms plus the Ortholog
table in `feba.db`. Its own feasibility question is the mirror of this one — not "can the condition be
represented" but "does a genome-scale model exist for these organisms at all", which is the next thing to
probe rather than assume.
