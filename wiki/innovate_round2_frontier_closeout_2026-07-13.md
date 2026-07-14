# /innovate round-2 frontier closeout — g8r2 + g5r2 are closed negatives (2026-07-13)

The `/innovate` round-2 pass (`wiki/innovate_genome_world_model_round2_2026-07-13.md`) left two survivors as "actionable but unbuilt": **g8r2 cross-organism catalog transfer** (the label-wall bypass) and **g5r2 self-supervised catalog** (flagged "must be scored, not built on"). Both are now resolved as **closed negatives** — honestly, before any build was committed on top of a refuted premise. This closes the entire round-2 forward frontier.

---

## g8r2 cross-organism catalog transfer — CLOSED (label-wall bypass refuted + already built)

A deep review (`/brainstorm`, 2 rounds, grounded-verified against the repo) established:

- **The core "zero target labels" claim is falsified by committed evidence.** `wiki/wider_amr_transferability_synthesis_2026-06-08.md` already ran the naive unchanged-transfer test and tagged the failures: **Campylobacter cipro → FAILS (TUNING)**, **Salmonella cipro → FAILS (TUNING+CONTENT)**, Acinetobacter/Pseudomonas meropenem → FAILS (CONTENT), Enterobacter ceftriaxone → FAILS (EXPRESSION). Its own conclusion: *"there is no family-wide rule; calibrate per organism."* You need the target's labels to KNOW the correct threshold (Campylobacter needs threshold 1, not E. coli's 2) — so transfer **relocates the label wall, it does not bypass it**.
- **The "4 QRDR cipro cells = one family-level catalog" premise is false.** `calibrated_amr_rules.json` gives them *different* rules: Campylobacter `qrdr_point` threshold 1; Klebsiella threshold 2 + `oqxA/oqxB` exclusions; Salmonella `broad` threshold 1. They were each calibrated per-organism.
- **The work is already largely built (non-duplication).** `scripts/klebsiella_cipro_transfer.py` + `wiki/klebsiella_cipro_transfer_2026-06-07.md` (first E. coli→Klebsiella transfer), `wiki/wider_amr_transferability_synthesis_2026-06-08.md` (the failure taxonomy), `scripts/organism_drug_validate.py` (the harness), and `scripts/amr_portal_explore_newcells.py` + `wiki/amr_portal_newcell_exploratory_result_2026-06-28.md` (the triage, already namespace-separated as EXPLORATORY). g8r2 re-treads all of it.
- **Organism-level catalog provenance-disjointness is ill-defined** for an AMRFinder-derived catalog (AMRFinder's reference catalog is pan-organism curated; `-O <target>` invokes a target-organism-specific detector — not "the E. coli catalog applied unchanged").

**Verdict:** g8r2 as pitched ("transfer a catalog UNCHANGED, zero target labels") is an already-closed negative. The only optional residual is a `TRANSFER_EXPERIMENTAL` reproducibility matrix that *systematizes* the 2026-06-08 synthesis — polish, not a frontier move, and ~80% covered by `organism_drug_validate.py`. Do NOT build the bypass.

---

## g5r2 self-supervised catalog — CLOSED (no negative class / circular where one exists)

The g5r2 claim: "train a generative effect-predictor SELF-SUPERVISED on the catalog's own edit→effect entries — no new wet-lab label required." Its round-2 kill-test was **precondition-only** (it proved the catalog is large enough, ≥100 entries — NOT that self-supervision works). Scored deeper here on committed data:

- **3 of 4 deployed catalog cells are ALL-POSITIVE — no negative class.** HIV NNRTI (16 major DRMs), SARS-CoV-2 Mpro (84 DRMs), fungal ERG11 (14/7) contain ONLY R-conferring mutations. A self-supervised discriminator has nothing to discriminate against — you cannot learn "this edit does NOT confer resistance" from a set that contains no non-resistance examples. Providing non-resistance examples requires external labels = **the label wall returns**.
- **The one graded catalog (WHO TB) makes the circularity explicit.** The WHO v2 catalogue *does* carry a negative class: **457 grade-1/2 (assoc-w-R) + 550 grade-4/5 (not-assoc) + 47,139 grade-3 "uncertain significance"** (88% of the catalogue is "uncertain" — the catalogue itself abstains on the bulk). BUT the WHO grades **are themselves distilled phenotype labels** (built from CRyPTIC + measured DST across labs). Training a predictor on WHO grades and validating on WHO grades is **circular** — it memorizes the deployed catalog, which IS the phenotype. The only non-circular version — generalize to propose NEW entries at un-catalogued positions — is exactly the **position-novelty self-awareness flag already built** (ledger row 423, `dna_decode/eval/position_novelty.py`), and its VALIDATION still requires independent phenotype labels for held-out edits (the label wall).

**Verdict:** g5r2 does not bypass the label wall. Self-supervision on catalog entries either has no negative class (3/4 cells) or is circular where one exists (WHO grades = distilled phenotype). It adds nothing beyond the shipped position-novelty flag. Do NOT build it.

---

## Net

The `/innovate` round-2 reframe (the catalog is a partial interventional world model) **stands** — but its two forward-build survivors both reduce to the project's one binding constraint: **LABELS, not models**. The honest highest-leverage frontier is unchanged from the reproducibility freeze: (1) **ACQUISITION** of a non-public wet-lab/clinical label source (a USER authority decision, clears the label gates by construction), or (2) **prospective-lock accrual** (already shipped + accruing; free, time-gated). Neither is an executor build task; both are already documented as the two non-foreclosed paths. Read-only over committed artifacts; frozen decoder surface byte-unchanged (verify_lock OK).
