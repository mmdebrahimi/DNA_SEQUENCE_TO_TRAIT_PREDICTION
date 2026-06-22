# Wave B grounding — a new HIV drug-resistance viral cell (2026-06-21)

> ## ✅ v0 SHIPPED + VALIDATED (2026-06-21) — the FIRST validated viral cell
> Built `dna_decode/data/hiv_amr.py` (v0 class-level NNRTI major-DRM catalog, sourced verbatim from the
> Stanford dataset page) + validated against the Stanford HIVDB PhenoSense fold-change (independent wet-lab,
> N=2272 isolates). **AUC (call separates fold): efavirenz 0.962 · nevirapine 0.985** (excellent) ·
> etravirine 0.75 · rilpivirine 0.70 · doravirine 0.56 (the honest class-level degradation on 2nd-gen
> NNRTIs — quantifies the v0.1 per-drug need). Result: `wiki/hiv_nnrti_v0_validation_2026-06-21.{md,json}`;
> harness `scripts/hiv_nnrti_validate.py` (+ tests). Circularity-safe (label = PhenoSense, NOT HIVDB Sierra).
> Frozen AMR surface byte-unchanged.
>
> ## ✅ EXTENDED — validate-vs-underlying-tool + a 2nd drug class (NRTI)
> **OLS baseline** (`scripts/hiv_nnrti_baseline.py`, a faithful Python reimplementation of Stanford's
> `DRMcv.R` — R not installed): the deterministic catalog **ties/beats** the full cross-validated OLS
> regression on EFV/NVP (Δbalacc +0.015 / −0.023) — the wrapper adds interpretability, not error; the
> regression's edge on 2nd-gen NNRTIs bounds the per-drug v0.1 headroom. **NRTI** (`hiv_amr` NRTI
> position-based catalog + `scripts/hiv_nrti_validate.py`, N=1867): sens ~0.98–1.0 but spec degraded
> (AZT/D4T ~0.49, TDF 0.41 — the *deliberate* T215-revertant/TAM over-call of position-based v0); OLS beats
> it +0.04…+0.16 balacc (the mutant-specific v0.1 headroom). Clinical cutoffs sourced from `DRMcv.R`.
> **v0.1 follow-ups (named):** mutant-specific catalogs (per-drug NNRTI + NRTI, data-derivable from the OLS
> coefficients) · within-subtype transfer check (unfiltered dataset's Subtype column) · report-card
> SCORED-cell integration · `dna-amr --drug` CLI route + HXB2-RT genome-mode caller for novel-FASTA input ·
> the remaining HIV classes (PI / INSTI / CAI — same datasets + cutoffs available).

Soraya `/probe`-by-hand + license-verification attempt, advancing Wave B (the highest-value net-new move
per `wiki/phenotype_trait_tool_completion_assessment_2026-06-21.md` + the Wave-A closeout) as far as
autonomy allows. **Result: the BUILD design is grounded + de-risked; the build itself is externally
walled.** Frozen AMR surface untouched (read-only grounding).

## What Wave B actually is (frame, corrected per the brainstorm)

Stand up a **NEW HIV drug-resistance viral cell** (RT / protease / integrase target-site DRMs → R/S),
validated against the FREE Stanford HIVDB genotype-phenotype dataset (PhenoSense IC50 fold-change). This
is the project's FIRST *validated* viral cell. It does **NOT** validate the shipped influenza NA cell
(different organism/drug/phenotype) — that cell stays `NO_FREE_PHENOTYPE_SOURCE`.

## Architecture: SUPPORTED — the existing viral template extends directly [grounded]

The HIV cell mirrors the shipped FOURTH-kingdom (viral) influenza pattern almost 1:1:

| Influenza NA cell (shipped) | New HIV cell (proposed) |
|---|---|
| `dna_decode/data/antiviral_amr.py` — `ANTIVIRAL_RESISTANCE_MUTATIONS` (drug→gene→subs), `AntiviralCall`, `call_from_observed_substitutions` | new `dna_decode/data/hiv_amr.py` — drug-class→gene→DRM set; same `Call` shape + `call_from_observed_substitutions` |
| `scripts/flu_na_caller.py` — BLAST NA-CDS-vs-genome → `observed_substitutions` (reuses the fungal caller) | new `scripts/hiv_dr_caller.py` — BLAST {RT, protease, integrase} CDS-vs-genome → `observed_substitutions` (SAME gene-generic caller) |
| routed via `dna-amr --drug oseltamivir` + `_target_site_record` → `amr-mechanism-call-v1` | routed via `dna-amr --drug <hiv-drug>` + the SAME `_target_site_record` |

The gene-generic `observed_substitutions` (in `scripts/fungal_erg11_caller.py`) already handles
multi-gene, intron-aware BLAST codon mapping — HIV RT/protease/integrase are intronless, so the colinear
single-HSP path is valid (like NA / K13). **No new engine architecture is needed.**

## The decisive design subtlety — the circularity trap [grounded, load-bearing]

**The HIVDB has TWO distinct products, and the validation MUST use the right one:**
1. **HIVDB GRT-IS / Sierra** — a rule-based genotype→susceptibility INTERPRETATION (penalty scores per
   DRM). This is a curated RULE system — the analog of the tool our caller competes against.
2. **The HIVDB genotype-phenotype DATASET** — PhenoSense IC50 fold-change, an INDEPENDENT wet-lab
   enzyme/replication assay measurement.

**Validate the new caller against (2), the PhenoSense fold-change — NEVER against (1), HIVDB's own
interpretation.** Scoring a DRM-catalog caller against HIVDB's GRT-IS output would be circular (rule vs
rule). The PhenoSense fold-change is independent of the DRM rules, so it is a legitimate
genotype→independent-phenotype test — this is exactly why HIV clears the project's circular-label gate.

Corollary (the project's standing "validate the wrapper vs the underlying tool" lesson): ship a NAIVE
HIVDB-Sierra-interpretation baseline ALONGSIDE the deterministic caller and **headline the delta** — does
our catalog add anything over just running Sierra? If not, the honest product is "use Sierra."

## The de-confound precondition [grounded — project standing discipline]

HIV **subtype** (B / C / CRF01_AE / …) is the lineage-confound analog that killed the bacterial embedding
arm. Any AUROC/accuracy claim must be reported **within-subtype** (concordance inside the same subtype),
not just overall — overall metrics conflate subtype structure with the resistance mechanism. Reuse the
`cohort_deconfound.py` discipline.

## Catalog scope [grounded]

Larger than influenza's 5-marker NA catalog: 4 drug classes (NRTI / NNRTI / PI / INSTI) × {RT, protease,
integrase}, hundreds of DRMs. But the lists are PUBLISHED + curatable offline (IAS-USA 2022 DRM list; the
WHO surveillance SDRM list; the Stanford mutation comments). A real but bounded curation effort — the
v0 slice could start with ONE drug class (e.g. NNRTIs: K103N / Y181C / G190A on RT — high-prevalence,
unambiguous) to prove the pipeline before the full catalog.

## External walls (what blocks the BUILD)

1. **Dataset license — RESOLVED 2026-06-21 (GO).** The user opened the authoritative Terms-of-Use page
   (`https://hivdb.stanford.edu/pages/FAQ/FAQ.html`, JS-rendered → read in-browser). Section I is a pure
   warranty disclaimer (no non-commercial / no no-redistribution / no research-use restriction); the FAQ
   states the data "belongs in the public domain", with "unfettered access", "publicly available" — "a
   curated public database". **Free for research/education download + validation use.** Required citation:
   **Rhee et al. 2003, Nucleic Acids Res 31:298-303**. ONE nuance: the footer is "© All Rights Reserved"
   with no explicit open-license grant → do NOT RE-HOST the curated dataset in the repo. The deliverable
   never needs that: follow the project's existing external-data pattern (refseq on D:, WHO TB catalogue
   gitignored w/ only checksums) — **download + validate + cite, gitignore the data, commit a fetch
   pointer + the cohort manifest**. So "redistributability" was never actually required; downloadability +
   research-use IS granted.
2. **Dataset download — EXTERNAL (user).** The genotype-phenotype dataset (`GENOTYPE-PHENO` nav →
   `https://hivdb.stanford.edu/pages/genopheno.dataset.html`) — confirm the download format (TSV/Excel) +
   that it pairs sequences with PhenoSense fold-change per drug.
3. **The DRM catalog must be SOURCED, not fabricated (provenance discipline).** The v0 NNRTI catalog
   (RT major DRMs) must come from a citable list — the HIVDB "NNRTI-associated DRMs" comment page (linked
   in the FAQ) and/or the IAS-USA mutations list — NOT curated from memory. This is a `/research`-with-
   sourcing step, the FIRST build step. (The genome-mode caller also needs an HXB2/consensus-B RT CDS
   reference FASTA, committed like the influenza `N1_NA` ref; the pure `--observed` wheel-only path needs
   neither the ref nor the dataset.)
4. **The build is a NEW wave** — sourced catalog (`data/hiv_amr.py`) + `hiv_dr_caller.py` + `dna-amr`
   routing + the de-confounded validation against PhenoSense + the Sierra baseline. Warrants the full
   planning chain (`/idea-anchor` → `/technical-plan` → pre-exec `/brainstorm` → `/save-plan` →
   `/execute-plan`) + the dataset + the sourced DRM list in hand. The pure-logic catalog + `call_from_
   observed_substitutions` slice is unit-testable WITHOUT the dataset once the DRM list is sourced.

## Drafted `/idea-anchor` (PARKED for user confirmation — the one skill Soraya does not self-anchor)

> **Formal rephrase:** Build a deterministic HIV-1 drug-resistance viral cell — a BLAST-based RT/protease/
> integrase target-site DRM caller (mirroring the shipped influenza NA cell) that emits the uniform
> `amr-mechanism-call-v1` record via `dna-amr --drug <hiv-drug>` — and validate it against the FREE
> Stanford HIVDB genotype-phenotype dataset's INDEPENDENT PhenoSense IC50 fold-change (NOT HIVDB's own
> Sierra interpretation), reported within-HIV-subtype, with a naive-Sierra baseline headlined as the delta.
>
> **The one foundational clarification:** start with a single high-prevalence drug class (NNRTI: K103N /
> Y181C / G190A) to prove the pipeline + the validation contract, OR curate the full 4-class catalog up
> front? (Drafted recommendation: single-class v0 slice first — proves the circularity-safe validation
> contract cheaply before the large catalog curation.)
>
> **Blunt opinion:** architecture is a solved, copy-the-influenza-template problem; the ONLY real risks
> are (a) the circularity trap (validate vs PhenoSense, not Sierra) and (b) the within-subtype de-confound
> — both named above. Value is real (first validated viral cell, on a free label that clears the gates
> that bound the learned arm). Gated only on the human license check + the dataset download.

## Recommended next step (for the user)

1. **Confirm the HIVDB dataset license** (open `https://hivdb.stanford.edu/pages/FAQ/FAQ.html` in a
   browser; the load-bearing go/no-go).
2. If usable: **ratify the drafted idea-anchor** (single-class v0 slice) → Soraya drives the planning
   chain + `/execute-plan` once the dataset is downloaded.
3. The actual HIV-cell build is a NEW wave — not finishable in this autonomous batch (externally walled on
   the license + the data).
