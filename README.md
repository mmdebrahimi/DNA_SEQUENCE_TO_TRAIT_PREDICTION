# DNA_SEQUENCE_TO_TRAIT_PREDICTION

<!-- amr engine spans bacteria + M. tuberculosis + fungi + viruses (HIV / SARS-CoV-2 / influenza), routed by --drug; each cell carries its own honesty tier. -->
**`dna-decode` — a deterministic, interpretable genome→trait decoder.** Give it a genome (bacterial,
fungal, mycobacterial, or viral); it returns a phenotype call (antibiotic / antiviral / antifungal
resistance R/S, or E. coli pathotype) **plus the exact genes/mutations that drove the call** + its own
blind spots + provenance. Mechanism-feature based, not an embedding black-box. **Not a clinical tool.**

> **2026-06 — the independent-label breakthrough.** The project's binding constraint was a *free,
> independent, measured* phenotype label (everything else risks circularity — scoring a rule against
> another tool's predictions). It is now broken across **bacteria** (EBI AMR Portal measured AST:
> E. coli / Salmonella / Klebsiella / Shigella, acc 0.83–0.995), **M. tuberculosis** (WHO-2023-catalogue
> rule on N≈2,845 measured-AST isolates: rifampicin acc **0.937**, isoniazid **0.914**), and **HIV-1**
> (Stanford HIVDB PhenoSense wet-lab fold-change). One legible view across every validation surface —
> each surface's distinct independence tier preserved, never averaged into a misleading aggregate:
> [`wiki/cross_kingdom_validation_summary.md`](wiki/cross_kingdom_validation_summary.md).

## What it decodes (v0.5.0)

| Tool | Trait | Validation |
|---|---|---|
| `dna-decode amr` (bacterial) | antibiotic R/S — **cipro / cef / tet / gent / meropenem** across **E. coli, K. pneumoniae, P. aeruginosa, S. aureus** | 6 drugs × 4 organisms, in-cohort + held-out + cross-source (NCBI) + cross-organism; every per-drug rule beats naive AMRFinder. Capstone: `wiki/amr_multiorganism_capstone_2026-06-07.md` |
| `dna-decode amr` (**fungal**, v0.5.0) | azole / echinocandin R/S — **fluconazole / voriconazole / caspofungin / micafungin** for **Candida auris** (BLAST-ERG11/FKS1 target-site engine) | **kingdom-jump** — same determinant-scan method, validated on a de-confounded C. auris WGS+MIC cohort (Gate G1): sens 1.0 across clades (ERG11 Y132F/F126L), label-limited specificity. `wiki/fungal_ep7_g1_closeout_2026-06-08.md` |
| `dna-decode pathotype` | E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/…) compatibility + abstention | VirulenceFinder-marker resolver; ExPEC recall 0.917; rest documented scope-limit |
| `dna-decode plasmid` (**v0.5.0**) | plasmid Inc-replicon typing (IncF/IncH/IncI/IncX/IncN/…) — *is the resistance plasmid-borne?* | deterministic PlasmidFinder-blastn caller (identity 95 / coverage 60); faithful-to-tool (not an independent baseline); offline-safe degrade |
| `dna-decode serotype` (**new**) | E. coli **O:H serotype** (wzx/wzy/wzm/wzt O-antigen + fliC H-antigen) | deterministic SerotypeFinder-blastn caller (identity 85 / coverage 60); `O?/H?` when a locus is unresolved; offline-safe |
| `dna-decode resfinder` (**new**) | acquired **AMR genes** (ResFinder DB) — an **independent** cross-tool check vs `amr` | deterministic ResFinder-blastn caller (identity 90 / coverage 60); `caller_is_independent_baseline: true` (acquired genes only — no point-mutations/efflux); offline-safe |
| `dna-decode pointfinder` (**new**) | **chromosomal AMR point mutations** (PointFinder; v0 E. coli FQ QRDR gyrA/parC/gyrB/parE) | deterministic blastn + codon-position lookup vs `resistens-overview`; `caller_is_independent_baseline: true` (independent of `amr`'s AMRFinder POINT); offline-safe |
| `dna-decode disinfinder` (**new**) | **biocide/disinfectant resistance** genes (DisinFinder; qac/form quaternary-ammonium + formaldehyde) | deterministic DisinFinder-blastn caller (identity 90 / coverage 60); often plasmid-borne (pair with `coloc`); offline-safe |
| `dna-decode mlst` (**new**) | **MLST sequence type** (PubMLST; v0 E. coli Achtman 7-gene) — exact-allele → profile → ST | deterministic blastn 100/100 + PubMLST profile lookup; **validated: K-12 MG1655 → ST10**; `dna-mlst --fetch-db` installs the scheme; novel/incomplete → ST not guessed; offline-safe |

915+ tests green. **8 decoders** (shared curated-DB blastn engine `dna_decode/typing/blast_caller.py`
+ codon-mapping `dna_decode/typing/codon_map.py`) **+ 3 cross-decoder analyses** that compose them:

| Analysis | What | |
|---|---|---|
| `dna-decode concordance` | AMR cross-tool check — **AMRFinder (`amr`) vs ResFinder (`resfinder`)** acquired-gene calls, gene-family level + Jaccard agreement | the independent second-opinion `resfinder` was built for |
| `dna-decode profile` | **run-all** — every assembly-FASTA decoder (pathotype+serotype+plasmid+resfinder) on one genome → one unified report | the "tell me everything" UX; each section degrades independently |
| `dna-decode coloc` | **AMR×plasmid co-localization** — is *this* acquired resistance gene on the same contig as a plasmid replicon (likely plasmid-borne)? | turns "both present" into "the gene sits on the plasmid"; same-contig is suggestive, not proof | The deterministic rules live in `dna_decode/eval/amr_rules.py::DRUG_RULE` (per-drug
threshold + AMRFinder-Subclass / QRDR-point / gene-prefix refinement). Engineering principle that held
across every organism: **count the drug's specific resistance determinants, not the broad drug-class bag.**

### Validated resistance cells — by kingdom (each with its own honesty tier)

Beyond the E. coli-family bacterial decoders above, the same determinant-scan method now spans four
kingdoms. Each cell is validated on the strongest *free* label available and carries its **honest
independence tier** — an `in-distribution` knowledge baseline is never relabelled as `independent`, and
each kingdom keeps a **namespace-separate** standing report card so the tiers can't be conflated:

| Kingdom | Cell | Independent validation | Tier |
|---|---|---|---|
| **Bacteria** | cipro / cef / tet / gent / meropenem × E. coli · Klebsiella · Salmonella · Shigella | EBI AMR Portal **measured AST** (free; BioSample/GCA-disjoint), acc **0.83–0.995** | provenance-disjoint, measured — non-circular (`wiki/amr_portal_independent_report_card.md`) |
| **M. tuberculosis** | rifampicin (`rpoB`) + isoniazid (`katG`/`inhA`) — WHO-2023 catalogue rule | EBI AMR Portal measured AST, **N≈2,845**: RIF acc **0.937**, INH **0.914** | independent, measured (`wiki/tb_report_card.md`) |
| **Virus — HIV-1** | NNRTI / NRTI / PI / INSTI / CAI (RT · protease · integrase · capsid) | Stanford HIVDB **PhenoSense** wet-lab fold-change (NNRTI EFV AUC **0.962**) | free, independent, isolate-level wet-lab label (`wiki/hiv_decoder_report_card.md`) |
| **Virus — SARS-CoV-2** | nirmatrelvir / ensitrelvir (Mpro / 3CLpro) | Stanford CoV-RDB fold-change — **in-distribution, underpowered** (37R/5S) | knowledge baseline, honestly labelled (`wiki/sarscov2_mpro_validation_result_2026-06-23.md`) |
| **Fungus — C. auris** | fluconazole / voriconazole / caspofungin / micafungin (ERG11 / FKS1) | de-confounded WGS+MIC cohort, sens 1.0 across clades; spec label-limited | kingdom-jump, G1-validated (`wiki/fungal_ep7_g1_closeout_2026-06-08.md`) |

The bacterial 6-drug deployed surface is under a **reproducibility freeze**
(`wiki/reproducibility_freeze_2026-06-13.md`) — frozen, sha-pinned, one-command reproducible. A
learned-embedding alternative to the deterministic rules was tested to a decisive verdict and is a
**closed 0-for-4 negative** (it learned population structure, not mechanism); the deterministic decoder
suite is the validated, shippable artifact (`wiki/negative_results_map_2026-06-13.md`).

**Every prediction carries its own trust badge inline.** As of the productization pass, each `dna-amr`
call emits a `validation:` line (and a `validation` block in the JSON record) reporting that exact cell's
honest tier + headline metric + the standing report card it came from — e.g. `validation:
INDEPENDENT_MEASURED -- acc 0.919 (N=8778)` for E. coli ceftriaxone, `INDEPENDENT_WETLAB -- AUC 0.962` for
HIV efavirenz, `IN_DISTRIBUTION` for SARS-CoV-2, `NO_FREE_PHENOTYPE_SOURCE` for the fungal cells,
`ABSTAINS_BY_DESIGN` for the carbapenem abstainers. The honesty discipline is now *user-facing* at the
CLI, not buried in the wiki (`dna_decode/data/trust_surface.py`; tiers never averaged, metrics never
fabricated).

### Organism-aware AMR calling (`dna-amr --organism`)

The per-drug `DRUG_RULE` is E. coli-tuned. Cross-organism validation (6 organisms × cipro/meropenem,
N≈30 each, NCBI AST) found it fails to transfer in three distinct ways — a **boundary taxonomy**:
**CONTENT** (counts intrinsic genes that don't confer R: Acinetobacter OXA-51, Pseudomonas nalC/oprD →
over-call), **TUNING** (threshold wrong where a single mutation suffices: Campylobacter cipro), and
**EXPRESSION** (regulation/derepression-driven R that gene-presence can't see: Enterobacter AmpC). Map +
evidence: `wiki/wider_amr_transferability_synthesis_2026-06-08.md`.

`dna_decode/eval/calibrate_organism.py` auto-selects the per-organism config (determinant **counter** ×
**threshold** + intrinsic gene-family exclusions) from a ≥15R/15S labeled cohort by leave-one-out balanced
accuracy, and **abstains** (`EXPRESSION_FLOOR`) when no presence-based config clears the floor (one-class/
under-powered cohorts → `INSUFFICIENT_EVIDENCE`). Validated configs ship in the committed registry
`dna_decode/data/calibrated_amr_rules.json` (independent-cohort out-of-sample validated:
`wiki/calibrated_registry_independent_validation_2026-06-09.md`). Pass **`dna-amr --organism <name>`** to
use a calibrated config when one exists (Campylobacter / Klebsiella / Salmonella cipro); an
`EXPRESSION_FLOOR` organism (Acinetobacter / Pseudomonas carbapenem) prints **`CALL: ABSTAIN`** rather than
over-calling. The registry is **opt-in** — the default (no `--organism`, or an organism with no entry)
uses the unchanged `DRUG_RULE`; calibrated configs are NCBI-AST in-sample-derived and stay opt-in pending
a different-lab cohort.

## Install

```bash
uv sync          # or: pip install -e .
# AMR genome mode also needs Docker + an AMRFinderPlus DB (see Gotchas); cached-run mode is pure-Python.
```

## Quickstart (verified output)

```text
$ uv run dna-decode list
dna-decode 0.5.0 - deterministic genotype->phenotype decoders
  amr        antibiotic resistance R/S - bacterial (cipro/cef/tet/gent/meropenem) + FUNGAL azole/echinocandin (C. auris) + ANTIMALARIAL artemisinin/K13 + chloroquine/pfcrt-K76T (P. falciparum)
  pathotype  E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/...) compatibility call + abstention

$ uv run dna-decode amr --drug ceftriaxone --amrfinder-run data/amrfinder_runs/GCA_008727135.1
sample: GCA_008727135.1  drug: ceftriaxone
CALL: R  [MODERATE | 1 determinant(s)]
  driven by: blaCMY-2  (CEPHALOSPORIN, 100.00% id)

$ uv run dna-amr --drug fluconazole --observed ERG11:Y132F --sample-id isolate1   # fungal, pure (no BLAST)
sample: isolate1  drug: fluconazole  organism: Candida_auris
CALL: R  [high | 1 determinant(s)]
  driven by: ERG11:Y132F
```

```bash
# Plasmid replicon typing on a genome assembly (blastn + PlasmidFinder DB; composes with amr):
uv run dna-decode plasmid path/to/assembly.fna --sample-id MY_STRAIN
# (downloads the DB once: curl -sSL https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/HEAD/enterobacteriales.fsa -o data/plasmidfinder_db/enterobacteriales.fsa)

# Pathotype on a genome assembly (pure-stdlib, no Docker):
uv run dna-decode pathotype path/to/assembly.fna --sample-id MY_STRAIN

# AMR on a novel genome (genome mode — runs AMRFinder via Docker; --organism selects the AMRFinder -O):
uv run dna-decode amr --drug ciprofloxacin --genome-fasta X.fna --organism Klebsiella_pneumoniae
```

Full capability table + validation provenance: **[Shipped decoders](#shipped-decoders-v040--two-interpretable-e-coli-genometrait-tools)** below. The rest of this README is project history (how the tool was arrived at).

---

## Project history — Phase 1 → v0.4.0 (how we got here)

> The sections below are the chronological research record (embedding-thesis exploration, Evidence
> Packets, the deterministic pivot). For *using the tool*, the section above is all you need.

## Status: Phase 1 — CLOSED 2026-05-17 (infrastructure + cross-drug architectural finding)

Phase 1 evidence collection closed 2026-05-17. Cross-drug architectural finding synthesis at `wiki/ep1_ep2_cross_drug_architectural_finding_2026-05-17.md`:

> At 12-strain smoke fidelity, frozen-NT-whole-genome-pooling PASSES on concentrated-signal AMR mechanisms (cipro QRDR point mutations: AUROC 0.750; cef plasmid acquired-gene β-lactamases: AUROC 0.833) AND FAILS on distributed mobile-element mechanisms (tet tet-family efflux + ribosomal protection: AUROC 0.400, anti-predictive). The architecture's failure mode appears mechanism-class-bounded, largely independent of drug identity at smoke fidelity.

EP1 cipro closed internally (`wiki/cipro_ep1_closeout_2026-05-17.md`) with a 4-tier adversarial audit infrastructure (mechanism × MIC × opacity merge with structurally-enforced SUSPEND gate). EP2 cef + tet smoke fired (cef PASS, tet FAIL, H17 falsified). No Databricks burst spent. External publication deferred per PC1=`internal_closeout`.

Phase 1 code: all 18 implementation steps shipped Wave 0-7 (2026-05-11 → 2026-05-12) + 3 hardening waves; cross-drug Evidence Packet evidence collection completed 2026-05-17 per the Evidence Packets framing reset 2026-05-15. **Phase 2 entry fired 2026-05-18**: BV-BRC strict-MIC 4-drug feasibility census ran (`scripts/bvbrc_strict_mic_4drug_census.py` + `wiki/bvbrc_strict_mic_4drug_census_2026-05-18.{md,json}`) — NO drug clears N=150 per-class at either strict-MIC or relaxed-MIC bars; structural bottleneck is `assembly_accession`. North star clarified: AI DNA decoder tool, not papers. v0 UX + success criteria LOCKED at `wiki/decoder_v0_ux_and_success_criterion.md` (CLI via `pipeline.py predict`, LOSO AUROC ≥ 0.70, cipro v0 / cef v0.1, JSON + markdown sidecar). 3 of 5 v0 criteria green via 24 new tests; 2 gated on Databricks N=147 cipro cache landing.

**Phase 2 in-flight (2026-05-22 → 2026-05-24)**: cipro interpretability audit completed on Precision 7780 (RTX 3500 Ada) by parallel Codex CLI session. Bounded-falsifier coordination plan + post-falsifier ship-path technical plan covered all 4 verdict branches × 3 gate states pre-committed. Codex on Precision 7780 ran the falsifier 2026-05-23 — **verdict = FAIL** (ranking-only rescue did not improve the ELX-family failure case on 12-strain Bucket B). Per the FAIL branch + north star, **v0 shipped 2026-05-24** as a cached-strain cipro predictor (`scripts/pipeline.py predict --strain-id ...`) with a documented scope-limit. v0 spec RELOCKED at `wiki/decoder_v0_ux_and_success_criterion.md` to match the implemented cached-strain surface (not the original genome-input decoder concept). Leakage-safe retrain on `leave_one_accession_out` CV yielded **AUROC 0.8697**. v0 closeout handoff: `wiki/dna_decoder_v0_closeout_handoff_2026-05-24.md`.

**v0.1 cipro genome-input slice 1 LANDED (2026-05-25)**: Codex shipped `pipeline.py predict --genome-fasta X.fna --annotations Y.gff3` end-to-end. Cross-path concordance validated on 4-strain mixed panel + same-strain parity (max prob delta 0.011599). Live embedding in batched chunks (OOM fix). Audit fallback to cohort-level framing when sample missing from per_strain.

**v0.1 cef cached-strain LANDED + validated (2026-05-25 → 2026-05-26)**: Codex shipped a dedicated 67-strain NT cache + trained cef classifier (CV AUROC 0.895 / AUPRC 0.838 on N=49 usable, 25R/24S). Duplicate-accession audit PASS (no LOSO leakage). Full-panel cached-vs-genome-input validation: **49/49 prediction concordance, 47/49 label alignment, max prob delta 0.063**. 2 shared model misses (562.28389 FP, 562.7695 FN) at decision-boundary probabilities. Currently debug-mode (no audit sidecar yet).

**v0.1 cef audit-aware packet (IN FLIGHT 2026-05-26)**: closeout slice per `plans/Cef_Audit_Aware_Packet_Design.md` — 5 artifacts (AMRFinder mechanism audit + cef MIC tier audit + new `scripts/drug_mechanism_phenotype_merge.py` + 4 canonical predict examples in canonical_audit_aware mode + release packet update with pre-committed verdict-branch wording). ~3.5-4 hr Precision 7780 compute. Closes the last debug-mode gap.

**Long-horizon roadmap drafted (2026-05-26)**: `plans/Trait_Decoding_Roadmap.md` maps Phase 0 (v0 cipro) → Phase 6 (eukaryotic organisms) with per-phase terminal claims + dataset prerequisites + falsifier triggers. EP-4 first non-AMR phenotype scoping: pathotype prediction (EnteroBase substrate; multiclass EPEC/EHEC/ETEC/UPEC/EAEC/commensal) per `plans/EP_4_Non_AMR_Phenotype_Candidates.md`.

**Test count:** 847 green (as of v0.4.0, 2026-06-07).

See `plans/EP1_EP2_Cross_Drug_Synthesis_Plan.md` for the synthesis plan; `plans/Cipro_Decision_Bundle_Plan.md` + `plans/Cipro_Decision_Bundle_Technical_Plan.md` for the EP1 closeout planning chain; `plans/EP2_Cef_Tet_Smoke_Design_Plan.md` for the EP2 design. See `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` for the original Phase 1 contracted ship-path. See `plans/Ecoli_G2P_Platform_Technical_Plan.md` for the full Phase 1 plan with Tier 1-5 attribution-success framework. See `docs/ARCHITECTURE.md` for the module map.

What runs end-to-end today:

| Surface | Entry point | Notes |
|---|---|---|
| **Pilot gate (HARD)** | `python -m scripts.pilot_gate --ast-tsv <path>` | Validates per-drug strain counts before ingestion fires. Exit 0=GO, 1=NO-GO, 2=PilotGateError, 3=no source. |
| **Full pipeline** | `python -m scripts.pipeline {ingest, train, predict, attribute}` | Single CLI with 4 subcommands; shared config-driven path resolution. |
| **Smoke regression** | `python scripts/smoke_pipeline.py` | <60s synthetic-fixture end-to-end via MockFoundationModel; asserts AUROC ≥0.85 + top-1 attribution = seeded gene. |
| **Leaderboard fan-out** | `python scripts/leaderboard.py --drugs ... --models evo,dnabert2` | Loops pipeline.py train per (model × drug); writes `data/processed/leaderboard.md`. |
| **Quant-fidelity check** | `python scripts/quantize_fidelity_check.py --full-precision-attributions <manifest.json> --quantized-attributions <manifest.json>` | One-time 4-bit vs full-precision ISM concordance check; gates whether Phase 1 attribution numbers are quantization-conditional. |
| **Viz** | `dna_decode.viz.browser.render_attribution_plot` + `export_attribution_tsv` | matplotlib PNG + TSV export; pygenometracks deferred to Phase 2. |
| **BV-BRC strict-MIC 4-drug feasibility census** | `python -m scripts.bvbrc_strict_mic_4drug_census` | Phase 2 entry (2026-05-18). Per-drug feasibility at strict + relaxed bars for cipro/cef/tet/gent. Writes `wiki/bvbrc_strict_mic_4drug_census_<date>.{md,json}`. Imports from `dna_decode/data/mic_tiers.py` (shared per-drug catalogs). |
| **v0 decoder predict** | `python -m scripts.pipeline predict --strain-id X --model-path M.pkl --cache C.h5 --annotations G.gff3 --audit-merge-json A.json --output Y.json` | v0 schema per `wiki/decoder_v0_ux_and_success_criterion.md` (2026-05-18 LOCKED). Emits JSON + markdown sidecar with prediction + calibrated_probability + confidence_tier + top_k_attribution + audit_verdict (SUSPEND propagation) + provenance. |
| **Provenance-disjoint validation** | `python scripts/provenance_disjoint_validate.py --drug X --organism Y ...` | Scores a deployed decoder on a submitter/lab-disjoint NCBI-PD cohort. Leakage exclusion via the data-driven `dna_decode/eval/cohort_manifest.py` (EXACT-self identity over all raw+parquet cohorts); FAILS CLOSED on an incomplete manifest (`INCOMPLETE_MANIFEST`, exit 2) unless `--allow-incomplete-manifest`. Powering census `scripts/ncbi_pd_provenance_census.py` self-persists to `wiki/provdisjoint_census_results.json`. |
| **Decoder-suite validation report card** | `python scripts/build_validation_report_card.py` | Read-only roll-up (exit 0 always — a report, not a gate). Rows = deployed-claim surface (`dna_decode/data/shipped_decoder_surface.py`) ∪ observed cells; honest per-cell tier, no aggregate headline. Renders a **Lineage disclosure (clonality-corrected)** table from the lineage sidecar when present. Writes `wiki/decoder_validation_report_card.{md,json}`. |
| **Lineage-disclosure metrics** | `python scripts/compute_lineage_metrics.py` | Recomputes per-cohort clonality-corrected sens/spec — the report card's raw sens/spec counts one vote per *isolate*, so an over-sampled clone inflates it. Greedy-representative Mash clustering (`dna_decode/eval/clonality.py`; chaining-resistant, NOT single-linkage), cluster-weighted confusion + Wilson CI + effective-lineage-N, graded lineage bucket. M4 reconciles raw sens/spec vs the committed artifact before trusting any weighted number. Needs Docker (Mash). Writes `wiki/provdisjoint_lineage_metrics.json` (schema `provdisjoint-lineage-metrics-v1`). |
| **External re-validation — Gate-0 preflight** | `python scripts/external_cohort_preflight.py --project PRJNA604975 --cohort-name oxford [--mic-open\|--mic-gated]` | Wave-0 go/no-go for re-validating the FROZEN decoder on an INDEPENDENT measured-MIC cohort. Bidirectional Entrez-primary/ENA-fallback BioSample resolver (`dna_decode/eval/biosample_resolver.py`) closes the accession-string leakage blind spot; emits assembly-availability (FREE vs ASSEMBLY-REQUIRED) + a BioSample-level leakage verdict (FAIL-CLOSED if any tuning overlap, >5% unresolved, or Entrez/ENA disagreement) + MIC-openness. Writes `wiki/external_preflight_<cohort>_<date>.json`. |
| **External re-validation — scorer** | `python scripts/external_cohort_revalidate.py --cohort oxford --drug ciprofloxacin --labels-dir D --preflight-json P.json` | Mirrors `provenance_disjoint_validate` (ensure_run → `call_resistance` → `_conf`) on the external cohort, into a SEPARATE `external-validation-v1` namespace (`evidence_tier=external_clinical`). Organism triple pinned VERBATIM from the frozen E. coli cells (AMRFinder `-O Escherichia`, registry `Escherichia_coli_Shigella`). MIC→R/S via `dna_decode/data/external_mic_labels.py` (`classify_tier` strict=HIGH_R/HIGH_S primary, relaxed=+DECISIVE secondary, bucket counts — NOT naive thresholding). Fail-closed unless preflight PASS. Writes `wiki/external_validation_<cohort>_<drug>_<date>.json`. |
| **External re-validation — roll-up** | `python scripts/build_external_validation_report.py --run-id R [--allow-degraded]` | Run-scoped (refuses glob-all without `--run-id`/`--artifacts`/`--allow-unscoped-glob`); globs the separate `external_validation_*` namespace → `wiki/external_validation_report_card.{md,json}`; per cell raw strict/relaxed + cluster-weighted STRICT sens/spec + Wilson CI + effective-lineage-N (reuses `clonality.py` math inline; Docker Mash, degrades to raw on non-Docker hosts). Skips `powering.hard_fail` cells; degraded only with `--allow-degraded`. FROZEN decoder report card + `compute_lineage_metrics` NOT touched (Fix C). |
| **Oxford ingestion — W0 probe** | `python scripts/oxford_w0_probe.py --project PRJNA604975 --mic-table T --key-col K --drug-col COL=alias` | Pins the MIC-table schema + crosswalk feasibility BEFORE the ingester codes against it: ENA candidate-field cardinality, row/key/dup summary, operator/censoring distribution, MIC-key→BioSample resolution rate → `wiki/oxford_w0_probe_<date>.json`. |
| **Oxford ingestion — labels + manifest** | `python scripts/build_oxford_labels.py --project PRJNA604975 --mic-table T --key-col K --drug-col COL=alias --run-id R` | Ingest → alias→BioSample crosswalk (ABORTS on hard conflict) → per-drug `selected_{strict,relaxed}.tsv` + `buckets` + the single `cohort_manifest_external_<run_id>.json` (the exact scored-cohort definition). |
| **Oxford ingestion — one-command driver** | `python scripts/run_oxford_revalidation.py --project PRJNA604975 --mic-table T --key-col K --drug-col COL=alias --drugs ciprofloxacin` | Chains W0 probe → labels → exact-set preflight (abort != PASS) → per-drug scorer → roll-up ONLY IF every drug run is acceptable (driver gating). Exit = worst child (3 hard-fail > 1 degraded > 2 gate > 0). |

Module map: `dna_decode/data/` (ingestion) + `dna_decode/models/` (foundation wrappers + cache + classifiers + classical baselines; `cache.verify_complete` integrity gate added 2026-05-15) + `dna_decode/interp/` (ISM + Tier 1-5 attribution) + `dna_decode/eval/` (CV + metrics + batched-call phylogeny + clade-only baseline) + `dna_decode/viz/` (browser) + `tools/` (Stage 2 bioinformatics-tool runner via Docker Desktop — Mash + AMRFinderPlus + Bakta).

## Phase 1 scope

| Aspect | Value |
|---|---|
| Organism | E. coli |
| Phenotypes | Ciprofloxacin, ceftriaxone, tetracycline binary resistance |
| Foundation models | Evo (primary), DNABERT-2, Nucleotide Transformer, GENA-LM (leaderboard) |
| Classical baselines | AMRFinder gene calls, k-mer logreg + XGBoost, gene-presence XGBoost (Step 18) |
| Baseline ML | Frozen foundation-model embeddings + XGBoost per drug |
| Attribution | In-silico mutagenesis (gene-level + nucleotide-level saturation) |
| CV | Leave-one-Mash-clade-out + clade-only baseline + per-clade reporting |
| Target | AUROC ≥0.80 SLO / ≥0.85 stretch; clade-baseline-gap ≥0.10 on ≥75% of held-out clades; ≥3pp gap vs best classical baseline on ≥2 of 3 drugs |
| Horizon | 3 months Phase 1; 12 months Phase 1+2+3 |
| Compute | Local GTX 860M (4 GiB Maxwell, NT v2 only — verified 2026-05-13) + Databricks burst for larger cohorts. 4-bit Evo unavailable (bitsandbytes requires CC ≥ 7.0). Original target was RTX 4090 + 4-bit Evo; never materialized. |

## Long-term vision

Multimodal genotype-phenotype platform — start with bacterial AMR (Phase 1), expand toward eukaryotes + image-paired phenotype data in later phases. NOT a direct stepping stone to "DNA → animal image" prediction; that would require a parallel multimodal track.

## Setup

```bash
# 1. Install uv (if not already on PATH)
#    Windows PowerShell:  irm https://astral.sh/uv/install.ps1 | iex
#    Linux/macOS:         curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Sync deps. pytest is in default deps (Wave 1.5 hardening fix):
uv sync

# 3. Run the test suite
uv run pytest tests/ -v

# 4. Optional: install dev tooling (ruff + pytest-cov)
uv sync --extra dev

# 5. (Advanced, gated on hardware) install bitsandbytes for 4-bit Evo quantization
#    Requires CC ≥ 7.0 GPU (Ampere / Ada / Hopper). NOT compatible with the project's
#    actual GTX 860M (CC=5.0). Skip unless running on A100+ or similar.
uv sync --extra quantize
```

## Phase 1 quickstart

End-to-end Phase 1 run (assumes `uv sync` + BV-BRC AST TSV downloaded + Mash CLI installed):

```bash
# 1. HARD gate: confirm you have enough labeled strains per drug
uv run python -m scripts.pilot_gate \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --target-per-drug 150 \
  --ast-tsv path/to/bvbrc_ast.tsv

# 2. Smoke test: <60s end-to-end on synthetic fixtures (sanity check before real run)
uv run python scripts/smoke_pipeline.py

# 3. Ingest: build cohort + download cohort genomes
#    For real-data runs, pass --assembly-metadata-csv pointing at the BV-BRC
#    Genomes-tab export (CSV). The adapter at dna_decode/data/bvbrc_genome.py
#    feeds contig_count + N50 + MLST + assembly_accession into
#    candidates_from_bvbrc_ast. Both --assembly-metadata (legacy YAML) and
#    --assembly-metadata-csv are mutually exclusive.
uv run python -m scripts.pipeline ingest \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --ast-tsv path/to/BVBRC_genome_amr.csv \
  --assembly-metadata-csv path/to/BVBRC_genome.csv \
  --download-genomes

# 4. Populate the embedding cache (deferred — see ARCHITECTURE.md for the wiring;
#    embedding cache populate is invoked from a Phase 2 helper script that hasn't
#    shipped yet; Phase 1 callers populate the cache externally via cache.populate()).

# 5. Train per-drug classifier + run CV + emit clade-only baseline + validation gate
uv run python -m scripts.pipeline train \
  --drug ciprofloxacin --model evo --include-clade-baseline

# 6. Run ISM attribution + Tier 1-5 classification for one strain
uv run python -m scripts.pipeline attribute \
  --strain-id <bvbrc-strain-id> \
  --drug ciprofloxacin \
  --card-path path/to/card.json \
  --amrfinder-path path/to/amrfinder.tsv \
  --output data/processed/attribution_report.json

# 7. Build leaderboard across foundation models + classical baselines
uv run python scripts/leaderboard.py \
  --drugs ciprofloxacin,ceftriaxone,tetracycline \
  --models evo,dnabert2

# 8. (Optional, gated on CC ≥ 7.0 GPU) Validate that 4-bit Evo attribution matches full-precision
uv run python scripts/quantize_fidelity_check.py \
  --full-precision-attributions full_manifest.json \
  --quantized-attributions quantized_manifest.json \
  --drug ciprofloxacin
```

## Decoder v0 quickstart (Phase 2 in-flight)

The v0 AI DNA decoder operates on **cached strains** — a strain whose NT embeddings already live in the HDF5 cache (built by `pipeline ingest` + the Databricks N=147 cipro populate). UX + success criteria locked in `wiki/decoder_v0_ux_and_success_criterion.md`.

```bash
# Predict cipro R/S for a cached strain, with top-K attribution + audit-verdict propagation.
uv run python -m scripts.pipeline predict \
  --model-path data/processed/models/ciprofloxacin_nucleotide_transformer.pkl \
  --strain-id 562.12345 \
  --cache D:/dna_decode_cache/embeddings/nt_n147_cipro.h5 \
  --annotations D:/dna_decode_cache/refseq/GCF_xxx.x/annotations.gff3 \
  --audit-merge-json wiki/cipro_mechanism_phenotype_merge_2026-05-17.json \
  --output result.json
```

Writes `result.json` + `result.md` (markdown sidecar) per the v0 schema:

- `prediction` (R/S) + `calibrated_probability` + `confidence_tier` (HIGH/MEDIUM/LOW)
- `top_k_attribution` — gene-level ISM hits with resistance-catalog tier labels (Tier 1–5)
- `audit_verdict` — propagated from the merge gate; explicit `suspend_gate_fired` flag + verdict explanation when training cohort had `SUSPEND_CONDITION_4`
- `provenance` — model, training cohort, LOSO AUROC, trained-on date

**Not a clinical decision support tool.** Audit verdict + provenance must accompany any downstream interpretation. See `wiki/decoder_v0_ux_and_success_criterion.md` for full v0 schema + success criteria.

## Shipped decoders (v0.4.0) — two interpretable E. coli genome→trait tools

The project's delivered value is **two deterministic, interpretable decoders** (installable console
commands after `uv sync` / `pip install -e .`). Both take a genome assembly and emit a call + the exact
genes/mutations that drove it + provenance — biologically interpretable, not embedding black-boxes.

| Command | Trait | What it reports | Validation |
|---|---|---|---|
| `dna-pathotype` | E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/…) | virulence-cluster compatibility call + abstention + canonical-VirulenceFinder diff | compatibility resolver; ExPEC/EPEC/ETEC supported, rest documented scope-limit |
| `dna-amr` | antibiotic resistance (R/S) — cipro / cef / tet / gent / **meropenem**; **E. coli + Klebsiella + Pseudomonas** (`--organism`) | R/S call + the curated AMRFinder resistance determinants driving it (e.g. `gyrA_S83L`, `blaCTX-M-15`, `aac(3)-IIa`) + `undetectable_mechanisms` blind-spots on S calls | **cipro** E. coli N=147 acc 0.925 (cross-organism QRDR-POINT rule); **cef** N=60 0.933; **gent** N=128 0.945; **tet** N=12 0.833. **Cross-source (NCBI, zero BV-BRC overlap): cipro 1.0 / cef 0.864 / gent 1.0 / tet 0.909**, beating naive AMRFinder. **Cross-ORGANISM: full Klebsiella 5-drug matrix** — cipro 1.0 / cef 0.80 / gent 0.867 / meropenem 0.867 (all ✅); tet 0.80 (sens-limited by efflux). **3rd organism Pseudomonas cipro N=30 acc 0.867; 1st Gram-positive S. aureus oxacillin/mecA sens 1.0** (genotype transfers; oxacillin-label-confounded spec — `wiki/staphylococcus_aureus_oxacillin_validate_2026-06-07.md`). Per-drug rules in `amr_rules.py::DRUG_RULE`; genome mode takes `--organism` |

```bash
# Unified entry (dispatches to the trait decoders):
uv run dna-decode list                                  # what it decodes + per-trait validation status
uv run dna-decode pathotype path/to/assembly.fna --sample-id MY_STRAIN
uv run dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x

# The per-decoder entries remain independently usable:
uv run dna-pathotype path/to/assembly.fna --sample-id MY_STRAIN
uv run dna-amr --drug ceftriaxone --genome-fasta X.fna   # genome mode needs Docker + AMRFinder DB
```

**Why deterministic, not embeddings:** the frozen-genome-embedding (NT-mean-pool) thesis was tested to a
decisive verdict and found to have **no E. coli AMR niche** — on the cleanest substrate (cipro) it lost
to the QRDR-POINT knowledge baseline and within-lineage scored at chance (it learned lineage, not
mechanism). See `plans/AMR_embedding_niche_decision_2026-06-05.md`. For AMR, mechanism features win and
are interpretable; the embedding architecture's remaining open frontier is non-AMR phenotypes lacking a
curated knowledge baseline (gated on a de-confounded labeled substrate; `wiki/HANDOFF_session_2026-06-05.md`).

## Pathotype resolver (E. coli) — v0 tool (SHIPPED 2026-06-04, tag `pathotype-v0`)

A self-contained, **pure-stdlib** CLI that takes an E. coli genome assembly (FASTA) and emits an
auditable **pathotype-compatibility** call with virulence-cluster provenance + a side-by-side diff
against canonical VirulenceFinder. Honest framing: it is a marker-based **compatibility resolver with
abstention, NOT a clinical predictor**. Supported (externally-valid) classes = ExPEC / EPEC / ETEC;
EAEC / commensal / clean-EHEC are a documented scope-limit (the resolver reports their modules but
flags low external validity).

```bash
# After `uv sync` (or `pip install -e .`), the `dna-pathotype` command is available:
uv run dna-pathotype path/to/assembly.fna --sample-id MY_STRAIN --out result.json
# or equivalently:  uv run python -m dna_decode.pathotype path/to/assembly.fna ...
```

Emits provenance JSON + a human summary:
- `derived_call` — 11-class honest surface (EHEC/STEC/tEPEC/aEPEC/ETEC/EAEC/UPEC/HYBRID/AMBIGUOUS/
  UNCLASSIFIED/COMMENSAL) with `confidence_tier` + `external_validity` + abstention rules. `--legacy-6class`
  preserves the original 6-class promise.
- `cluster_profile` + `marker_hits` — which virulence clusters drove the call (k=15 k-mer-seed coverage
  over the VirulenceFinder E. coli allele DB; ≥0.80 = confident).
- `vf_diff` — **canonical VirulenceFinder side-by-side** via real `blastn` over the SAME VF DB: per-gene +
  per-cluster concordance. **HONESTY:** both callers use the same DB, so `caller_is_independent_baseline:
  false` + a same-DB caveat ship in every diff — it is an AUDIT of the fast caller, not independent
  validation. Degrades to `status: unavailable` (never dropped) when `blastn` is absent. Use `--no-vf-diff`
  to skip.

**BLAST+ for the diff:** install NCBI BLAST+ (`blastn` + `makeblastdb`) and either put it on PATH or set
`$BLASTN_BIN`. Without it, the resolver still runs fully; only the canonical diff degrades to `unavailable`.

**Marker DB:** `data/virulencefinder_db/virulence_ecoli.fsa` (fetch from the VirulenceFinder Bitbucket DB;
the CLI prints the exact `curl` command if it's missing). DB checksum is pinned in every provenance record.

## Phase 1 success criteria

Phase 1 ships when:

- Smoke pipeline passes (`scripts/smoke_pipeline.py` returns exit 0)
- LOMO-clade-out CV AUROC ≥0.80 SLO / ≥0.85 target per drug
- Embedding model AUROC ≥0.10 above clade-only baseline on ≥75% of held-out clades
- Top-K=20 attribution-tier distribution: cipro ≥40% Tier 1-3 hits; ceftriaxone ≥25%; tet ≥30%; all ≤20% Fail
- Best foundation model beats best classical baseline by ≥3pp AUROC on ≥2 of 3 drugs
- Quantization-fidelity check returns GO (mean Spearman ≥0.7, intersection ≥0.6)

Phase 2 redesign trigger: classical baselines win on ≥2 drugs. The Step 18 classical-baselines control wires this empirically — see `plans/Ecoli_G2P_Platform_Technical_Plan.md` validation-gate section.

## Pilot gate alternate inputs

- `BVBRC_AST_TSV=path/to/ast.tsv` env var
- `bvbrc_ast.local_tsv_path: path/to/ast.tsv` in `config/datasources.yaml`

## Optional: route caches to a USB drive

Phase 1 runtime needs ~25GB (foundation models + strain genomes + embeddings). If your C: drive is tight, route caches to external storage:

```bash
# Replace E: with your drive letter (Windows) or /mnt/usb (Linux)
export HF_HOME=E:/hf_cache               # HuggingFace tokenizer + model cache
export DNA_DECODE_CACHE_ROOT=E:/dna_decode_cache
```

Then edit `config/datasources.yaml` to point `cache_dir` fields at the USB-backed path.

## Project workflow

Built using a personal Claude Code skill ladder for project planning:
- `/idea-anchor` → `/project-init` → `/brainstorm` ×3 → `/technical-plan` → `/probe` → `/execute-plan`
- Project ledger maintained via `/project-state` skill
- Execution state tracked in `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`
- All planning artifacts captured as audit trail

See `project_state/dna-decode-2026-05-11.md` for full decision history (17 hypotheses, 12+ decisions made, 54+ action-log entries as of 2026-05-17; "Phase 1 / 2 / 3" labels retrospective-only — new work tracked as Evidence Packets per the 2026-05-15 framing reset; Phase 1 evidence collection CLOSED 2026-05-17 with the cross-drug architectural finding synthesis).
