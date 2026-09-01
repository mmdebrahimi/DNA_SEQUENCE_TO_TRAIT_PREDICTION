# Open threads — what a fresh session should pick up

Short by design. **Transient** state only: what is in flight, what is waiting on the user, what was just
learned that isn't durable yet. Durable findings belong in `CLAUDE.md` / `wiki/`; scope facts belong in
`scripts/project_status.py` (derived, never written).

Prune aggressively. A stale entry here is worse than an empty file.

_Last updated: 2026-08-29._

---

## The system design (drafted 2026-08-31, awaiting ratification)

`plans/Hybrid_Decoder_Architecture_Plan.md`. **The hybrid is NOT "catalog + ML predictor"** — that framing
scored 0 survivors in the framing sweep. Measured shape: **CALL / DOUBT / EVIDENCE**.

| layer | status |
|---|---|
| **L1 CALL** — deterministic curated rules | shipped, 110 cells |
| **L2 DOUBT** — "this call may be incomplete, and why"; **never a competing call** | **mostly missing**, one unwired prototype |
| **L3 EVIDENCE** — de-confounding, nulls, denominators, leakage, provenance | built, under-exposed |
| **L4 LEARNED** — forward/inverse, orthogonal modalities; molecular + constructed ONLY | shipped, bounded |

**L2 is the innovation and it is cheap.** The catalog's failure mode is COMPLETENESS, not accuracy — found
twice the same way (gentamicin `rmt`, HIV NNRTI), both invisible until independent labels arrived. Nobody
in the field ships "my catalog might be wrong here."

**Critical path F-A → F-B. F-A IS COMPLETE 2026-08-31** — all four steps.
`wiki/doubt_layer_2026-08-31.md` is the memo; `dna_decode/eval/doubt.py` + `scripts/doubt_layer_per_cell.py`
+ 30 tests are the artifact. Headline: across 1,818 genomes / 1,279 uncounted determinant families /
six drugs, **exactly one family survives the family-wise correction — `rmtE1`, p = 4.11e-12 — and it is
the confirmed gap.** The raw signature flagged 5; the correction drops 4 and keeps the true one.

**Two things it is easy to get wrong here, both measured:**
- **Enrichment is the WRONG null.** A lower-tail binomial on the observed S count calls `aph(6)-Id`
  (62R/28S) STRONG at p≈5e-5 — and that is a CORRECT exclusion (a streptomycin gene travelling with
  gentamicin resistance by linkage). Every co-occurring determinant is R-enriched. The signature is
  **purity**; one S carrier ENDS the signal. Pinned by test.
- **A position-BASED catalog can never fire the position-novelty flag** (every substitution at a
  catalogued position is already called), so those cells report `not-applicable`, never "no doubt".

Registered augment-only on **both** trust surfaces — the inline badge and the standing report card
(a 4th disclosure layer beside lineage / source-concentration / prospective). Augment-only was
**verified by diff**: 27 cells before and after, zero non-doubt fields changed, all state counts
identical. Card rows are **drug-level**, which is why `rmtE1` renders against three organisms — it is
a property of the rule, not of those cohorts.

**Now unblocked: F-B (curation) has its measured baseline** — but both of its terminal moves are user
authority calls (below).

## The five project families exist now (seeded 2026-08-31, user-authorized)

Ledgers at `project_state/{doubt-layer,catalog-curation,evidence-surface,learned-narrow,label-acquisition}-2026-08-31.md`;
the frontier lives in `project_state/dna-decode-2026-05-11.md` under `## Project Families` +
`## Requirements Flow-down` and is **machine-readable** — `advance_ranker.rank()` parses it, ranks F-A
first, and correctly reports F-B blocked on F-A. Three families (F-C evidence-surface, F-D
learned-narrow, F-E label-acquisition) are eligible and untouched. Self-init account 6 / ceiling 25.

**F-D is a RESTRAINT family** — its deliverable is a boundary that stays enforced, not a build. Its
ledger records the corrected regime map (population design, not organism complexity) so the compression
error that has bitten three times has somewhere durable to live.

## F-C and F-D advanced 2026-08-31 (same run)

**F-C (evidence surface).** Four disclosure layers rendered on the report card; only **two** reached a
decoder call. `lineage` + `source_concentration` were card-only; `prospective` surfaced only when it
CONTRADICTED. All four now reach the record **and** print. The case that matters: e.coli x gentamicin
reports sens 0.893 from a cohort 95% one BioProject with zero `rmt` carriers (source-diverse: 0.523),
and that caveat is now on the call. Memo `wiki/evidence_surface_reachability_2026-08-31.md`.
*Not audited:* the HIV / TB / pgx cards — named follow-on.

**F-D (learned-narrow, the RESTRAINT family).** The regime boundary is now a function, not prose:
`dna_decode/eval/regime.py::screen_proposal(population, endpoint, method)`. **Exactly one regime
refuses** (natural x organism x ZERO-SHOT); a SUPERVISED natural-population proposal returns
`REQUIRES_DECONFOUNDING` **with conditions, never a refusal** — compressing that scope IS the error,
made three times. `uv run python scripts/regime_map.py` refuses to certify a regime whose cited
artifact is missing. Screen before proposing a learned decoder.

## PEAR is reclassified — F-E's premise was wrong (2026-08-31)

`wiki/pear_substrate_screen_2026-08-31.md`. Screened after a review flagged the ranking. Two corrections:

- **Not an L1 label source.** ~200k constructed MG1655 strains, single-copy `blaCTX-M-14` at a fixed
  chromosomal site, **relative growth** for ~23k under cefotaxime/ceftazidime, >90% of single mutations.
  That is `constructed_molecular` → regime **`WORKS`** (TEM-1 path, Spearman 0.761). **It cannot address
  the AMR label wall**, which lives in the natural×organism regime.
- **Not an acquisition.** Public and free — SRA + GitHub, no DUA, no money. The authority gate was a
  consequence of the misclassification.

**What it is:** the natural external replication of the one working learned regime — same shape as TEM-1,
different β-lactamase, different drugs. It **clears every applicable rejection gate**.

**The blocker is ARTIFACT FORMAT, not availability.** `PRJNA687219` resolves correctly but is 478 Gbases
of raw reads; the GitHub repo has no CSV/TSV at all — its two `.RData` files are serialized **ggplot2
plot objects** (`Figure.2A`), which is why `pyreadr` and `rdata` both fail. Needs R (not installed;
**C: at 99%**). **Cheapest untried route: the journal supplementary tables** (unread — PMC cookie-gated).
**Run `assay_degeneracy()` before believing any PEAR score** — a selection-growth assay has a floor, and
CcdB (79.3% tied at ceiling) posted the sweep's best number for exactly that reason.

## Doubt-layer firing rate measured (2026-08-31)

**0 STRONG across 747 candidate families on 4 drugs with no confirmed gap**, against **1 STRONG across
131 families** on the one drug that has one. The raw signature fired 4 times on those same 747; the
family-wise correction removed all 4. **Deliberately NOT a false-positive rate** — a hit on an
unconfirmed drug is ambiguous between a false positive and an undiscovered gap. Categories are
predeclared `confirmed` / `unconfirmed` / `unassessable`, never "clean".

## Waiting on the user (authority calls — not executor tasks)

1. ~~**v2 gentamicin lock**~~ — **DECIDED AND SHIPPED 2026-08-31** (user-authorized).
   `wiki/gentamicin_v2_lock_2026-08-31.md`. Rule is now `subclass_any={GENTAMICIN}` +
   `symbol_rescue=^(rmt[A-H]\d*|npmA\d*)$`; needed an ENGINE change (the first WIDENING refinement)
   because both prior refinements narrow and compose as AND. **E. coli N=131: sens 0.523 -> 0.892
   (+0.369), spec 0.985 unchanged.** Freeze + lock RETIRED, not broken —
   `prospective_lock_manifest_2026-08-31.json` supersedes; v1 manifest preserved.
   *Cost, now enforced in code:* the prospective evidence is SPENT for both E. coli cells; the card
   returns `superseded_by_surface_change` with numbers withheld and the clock restarts at 2026-08-31.
   *Honest limit, unchanged:* zero S-labelled `rmt` carriers anywhere, so the unchanged specificity is
   an absence, not a bound.

2. **Whether a single-source cell warrants more than disclosure.** 3 of 10 SCORED AMR cells rest on one
   BioProject (`wiki/provdisjoint_source_concentration_2026-08-28.md`). Current answer is *disclose*, in a
   namespace-separate layer. Demoting them is a scope decision.

3. **Whether to compress the other 17 long CLAUDE.md bullets (~9,100 words).** The file costs ~36,800
   tokens EVERY session; 18 bullets are long AND cite a resolvable memo, so their derivations could become
   pointers without losing anything. I compressed only my own (1,022 → 309). The rest is other sessions'
   institutional memory and rewriting it wholesale is a call about what belongs in the always-loaded
   surface. Measure first: `uv run python scripts/claude_md_weight.py`.
   *Two long bullets have NO external store and must stay whole — the tool already protects them.*

4. **The 7 unscreenable colour cells** — no existing evidence tier fits (`NO_FREE_SOURCE` is about labels;
   `NOT_CENSUSED` means never-scored). And whether curating the 40 unrecorded colour loci is worth doing
   (fabrication hazard unless every locus is OMIA/literature-sourced).

## The gene-LLM idea (raised 2026-08-30, prior-art-checked 2026-08-31)

Draft anchor `wiki/idea_anchor_genomic_language_model_2026-08-30.md` (**NOT ANCHORED** — user-confirmed
skill). Prior-art check `wiki/prior_art_genomic_language_models_2026-08-31.md`: **the field is crowded.**
gLM / gLM2 / GenSyntax already build gene-token models, and four of this repo's findings (650M peak,
fluency≠function, population-structure confounding, curated-rules-beat-ML) are already published.

**Cheapest decisive experiment, no training run:** score off-the-shelf `tattabio/gLM2_650M` on our
de-confounded benchmarks against the curated-catalog baseline. If a published gene-token model can't beat
a hand-written determinant catalog on constructed variation, that answers the idea for an inference pass
— and it is exactly the comparison the critique literature says nobody runs.

**Open user questions:** token level · training objective · natural vs constructed regime (drafted answers
in the anchor).

## Borrowable methods for the resistance blind spot (searched 2026-08-31, none tested)

`wiki/borrowable_methods_resistance_scoring_2026-08-31.md`. The mechanism (conservative substitutions at
average-conservation sites) says we need a signal that is NOT sequence plausibility. Ranked candidates:
**ΔΔG of BINDING** (Rosetta `flex_ddG` / FoldX — physics, maximally orthogonal, RMSE 1.2 kcal/mol
published, and it scores variants with no catalog entry) · **positive-selection scans** (HyPhy
`MEME`/`FUBAR`, free, no GPU — but codon scans MISS HGT-mediated resistance, so target-site cells only) ·
**drug-conditioned** ConPLex (PNAS 2023, code public) · Potts/DCA (fixes independence, not the
plausibility framing) · inverse folding (a *stability* specialist; resistance is *binding*).

~~**Proposed first move — score ΔΔG_bind on the catalog-negative subset**~~ — **KILLED 2026-08-31 by an
executed kill-test** (`wiki/innovate_blindspot_framing_sweep_2026-08-31.md`). A **zero-tool deterministic
position-novelty flag already recovers 60.4%** of the EFV blind spot (lift 4.69). ΔΔG's *premise* is sound
— the blind spot IS pocket-mediated, 3.05× burden-adjusted enrichment, `VERDICT: GO` — but it is **not the
cheapest move**, and the comparator is the free flag's **0.604**, not the catalog's 0.962.

**WINNING framing instead: the blind spot is a CURATION gap.** The drivers are named and counted
(V179D ×12, A98G ×10, H221Y ×7, F227C ×5, V108I ×4, V179E ×3) and sit at positions **absent from** the
deployed 8-position `NNRTI_RT_MAJOR_DRMS`. And `hiv_amr.py` is **NOT** pinned by the prospective lock, so
curating it does **not** invalidate the lock or the freeze — unlike the gentamicin `rmt` fix. **Whether to
edit a shipped catalog is a scope decision → user call.**

**Second, independent move (F3):** the position-novelty flag lives in `dna_decode/eval/` and is **never
surfaced** in `hiv_amr.py` or `cli.py`; `AbstentionVocab` already exists to carry it. Wiring only.

**Sixth family added after re-running a safeguard-blocked search (2026-08-31):** temporal /
frequency-trajectory selection inference (Wright-Fisher HMM, `WFABC`, Beta-with-Spikes; multinomial-logistic
clade models). **Regime split: fits our VIRAL cells** (direct precedent — applied to influenza drug
resistance) **and inherits our known-fatal confound on BACTERIAL cells** (a frequency rise in clonal
bacteria conflates "fitter variant" with "clone spread"). Only RELATIVE growth is ever identifiable.
Most plausible payoff: the prospective-lock arm (NCBI-PD carries collection dates). Also `Phylowave`
(*Nature* 2024) finds fitness-increased lineages **without predefined clades** — relevant to our reliance
on Mash/MLST/Napier as hand-chosen partitions.

**Acquisition target found:** PEAR — ~23,000 E. coli strains, each a unique single-copy `blaCTX-M-14`
variant, growth measured under cefotaxime/ceftazidime, with prospective/retrospective model split.
Constructed variation at scale on an AMR target with measured phenotype — the regime this repo says works.

## Cheap untried levers (executor work, no authority needed)

- ~~FBA conditional switch — continuous ratio as a ranking~~ **DONE 2026-08-29, bounded PASS.**
  Within-gene AUROC 0.7308 (non-flat, n=26, p=0.001); all-genes 0.5896 because 61% are flat. Oracle
  ceiling 11/67 exact-set vs deployed 3/67, and it ranks rather than calls. The failure is silence, not
  error. **REPLICATED on the 25-source Keio carbon axis: AUROC 0.8133 (n=69, p=0.0005), flatness 68.2%.**
  `wiki/fba_within_gene_ranking_2026-08-29.md`. **The follow-on is now CLOSED, not open:** on carbon the
  deployed rule already gets 23/217 and the oracle ceiling is 27/217 (+4 genes, +1.8pp), so estimating k
  would buy almost nothing. The 4-media "3->11" was a small-axis artifact. Do not build it.
- ~~FBA axis choice is a free lever~~ **SPENT + QUANTIFIED 2026-08-29.** All three axes measured; the
  lever is now a rule: **pick an axis whose WILDTYPE growth spreads** (distinct-growth fraction / CV,
  42 LP solves, seconds — `scripts/fba_axis_dynamic_range.py`). Flatness 61.2/68.2/75.5% tracks it
  monotonically on both summaries. Nitrogen was the worst available choice. n=3, so a direction that
  survived a common yardstick, not an established relationship.
- ~~Abstaining conditional-essentiality CLI~~ — checked 2026-08-30: there is **no shipped surface** to fix
  (`dna-fba`'s `conditional` is the variant-LOF-uncertain case; `carbon` is utilization). Shipping it would
  be a NEW phenotype-claiming surface on a published package = authority, listed above.
- **Staleness auditor** — one clean 110/110 corpus run at `TOTAL_TOKEN_BUDGET=5500` to verify the OOM
  mitigation. Still labelled unverified in `scripts/kaggle/staleness_corpus_kernel.py`.

## The FBA switch cell, as it now stands (2026-08-29)

Open, and its terms are now separated. **Direction is fine** — where the model varies it points the right
way (AUROC 0.71–0.81 on three substrates, all p≤0.001). **Accuracy when it commits is fine too** — 23/33
= 70% exact-set on carbon against a chance expectation of 0.78, ~30× chance. **Coverage is the whole
problem** — the model's call is constant for 85–94% of genes, so it is silent rather than wrong, and that
silence is structural (flat stoichiometry), unreachable by any readout change, and predicted in seconds by
the axis's own dynamic range. The readout lever is closed (+1.8 pp on the best-measured axis).

**Do not quote "10.6% exact-set" for this cell** — that scores the model against a target it cannot hit
for most genes. Quote the anatomy. The remaining bottleneck is the one already measured: the conditioning
signal is not measured in the conditions the phenotype data uses (PRECISE-1K ∩ Keio carbon = 11 of 28,
621 of 1,035 samples glucose).

## Known-stale / do not trust without re-deriving

- Any **cell count or trait count written in prose**. `scripts/project_status.py` is the authority — it
  caught two of my own figures wrong within an hour of writing them (46 traits → **44**; "~3x" → **4.1x**).
- `wiki/project_distillation_2026-08-29.md` says **46 CLI traits**. It is **44**. Left uncorrected as a
  worked example of exactly the drift this file exists to route around.
