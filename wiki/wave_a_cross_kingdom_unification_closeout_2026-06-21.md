# Wave A closeout — cross-kingdom determinant unification (2026-06-21)

**Verdict: Wave A is effectively ALREADY SHIPPED. Do NOT build a genome-map cross-kingdom overlay tier.**

This closes the Wave A candidate from `wiki/phenotype_trait_tool_completion_assessment_2026-06-21.md`
("unify the shipped fungal/viral/TB/Plasmodium determinant decoders under the genome-map as overlay
tiers, the analog of the AMR + virulence tiers"). Produced via the planning chain: executor-run
`/idea-anchor` (flagged the literal framing as likely wrong) → multi-round `/brainstorm` (Codex, 2 rounds,
all findings grounded + repo-verified). Frozen AMR surface untouched throughout.

## The grounded findings (all repo-verified)

1. **The literal framing (Bakta cross-kingdom tiers) is REJECTED.** The genome-map re-tiers **Bakta**
   annotation, and Bakta annotates **bacteria only**. The AMR + virulence overlays work *because*
   AMRFinder/blastn return SUBJECT COORDINATES that join to a Bakta CDS via the length-reconciled
   contig-name map + coord-overlap. The four cross-kingdom engines return **substitution/mutation calls**
   (`ERG11:Y132F`, `NA:H275Y`, `K13:C580Y`, `rpoB:S450L`), NOT genome-space coordinate features — there is
   no Bakta CDS to join to and no coord-join to perform. Forcing them into Bakta tiers would manufacture
   fake join confidence AND inherit bacterial-only assumptions (`unknown_under_bakta_db_light`). [grounded]

2. **Unification is ~3/4 ALREADY SHIPPED.** `dna_decode/amr/cli.py::_target_site_record()` (line 65)
   already maps fungal **and** antiviral **and** antimalarial target-site calls onto the **uniform
   `amr-mechanism-call-v1` record — identical shape to the bacterial path** — and `dna-amr --drug` already
   ROUTES all of fungal/antiviral/antimalarial/bacterial. `dna_decode/data/shipped_decoder_surface.py`
   already unifies them as validation cells, and the report card renders them. So the four decoders are
   already unified at BOTH the record-schema and CLI-routing levels for three of four kingdoms. The
   net-new value of a "unification wave" is near zero for those three. [grounded, verified]

3. **TB is the genuine outlier.** TB is **not** routed through `dna-amr` (no TB import in `amr/cli.py`),
   is **absent** from `shipped_decoder_surface.py`, returns a `DrugCall` dataclass branded
   `KNOWLEDGE_BASELINE` / `NON-FROZEN`, and consumes **H37Rv VCF** input — a categorically different
   contract from FASTA→target-site-BLAST. The only honest build left in "Wave A" is a small non-frozen
   `DrugCall → amr-mechanism-call-v1` TB adapter — NOT a genome-map tier — and even that must NOT silently
   promote TB into the report-card "shipped/validated" surface (it is a knowledge-baseline). [grounded, verified]

## Decisions

- **Do not build a genome-map cross-kingdom overlay tier.** The coordinate-join host abstraction does not
  apply cross-kingdom. The genome-map stays a bacterial Bakta-feature host (AMR + virulence overlays).
- **The `/probe genome-map-cross-kingdom-determinant-unification` step is moot** — its central question
  ("do the engines return coordinates / is the Bakta host applicable") is already answered NO. Skip it.
- **A per-sample multi-drug/multi-engine report is DEFERRED** — it is a UX wrapper over the existing
  `amr-mechanism-call-v1` record, not a new decoder capability, and risks implying unsupported
  organism×drug cross-products are endorsed unless driven strictly from `shipped_decoder_surface.py`.
  Build only if a concrete "one sample, many drugs" operator workflow appears.

## Open fork (the only honest Wave-A build left) — AUTHORITY decision

**TB `DrugCall → amr-mechanism-call-v1` adapter — yes/no?**
- *Case for:* adds TB to the unified record/CLI surface (the one kingdom that's missing), small + non-frozen.
- *Case against / guardrail:* TB is a `KNOWLEDGE_BASELINE` (the WHO catalogue was built partly from
  CRyPTIC) and consumes VCF, not FASTA — it must stay branded non-frozen and stay OUT of the report-card
  "shipped/validated" cells. Lower value than Wave B.
- *Drafted recommendation:* **defer**. It's cosmetic unification of an honestly-different-contract,
  honestly-unvalidated cell; the higher-value move is Wave B. Revisit if a unified "decode any pathogen
  genome" CLI surface becomes a product goal.

## Redirect: Wave B is the higher-value net-new move (with a corrected frame)

The evidence says Wave A's capability is largely already shipped, so the higher-value net-new move is
**Wave B — stand up + validate a NEW HIV viral cell** against the free Stanford HIVDB genotype-phenotype
dataset (GO per `research_outputs/viral-antiviral-resistance-gp-datasets-2026-06-21.md`).

**Frame correction (grounded):** the currently-shipped *viral* cell is **influenza NA**
(`NO_FREE_PHENOTYPE_SOURCE`). A free HIV dataset validates a **NEW HIV DRM viral cell**, NOT the existing
influenza cell — the organism/drug/phenotype contracts don't match. Describe Wave B accurately as "a new,
validated HIV viral cell," never "validating the shipped influenza cell."

Wave B remains gated on: (1) a manual HIVDB dataset terms-of-use/license confirmation (the one
load-bearing verify item), and (2) the project's standard within-subtype de-confound precondition before
any AUROC claim. Both are real waves, not finishable in one autonomous batch.
