# Horesh 2021 File F1 — label-provenance audit (H1 frontier)

**Date:** 2026-06-04. **Substrate:** Horesh et al. 2021 curated E. coli collection, File F1
(`F1_genome_metadata.csv`, 10,146 records, Figshare DOI 10.6084/m9.figshare.13270073, downloaded to the
gitignored `data/horesh_2021/`). **Question (H1):** does a de-confounded, non-circular ExPEC+EPEC label
substrate exist so the v0.1 embedding-classifier upgrade can be honestly validated against the v0
deterministic resolver?

## TL;DR — the de-confounded cohort is NOT buildable; the labels frontier for pathotype discrimination is closed

Independent (isolation-source) labels DO exist robustly for ExPEC, but the ExPEC-vs-EPEC discrimination
task is **intrinsically confounded with sampling context** along every control axis. The handoff's
"within-study / lineage-matched ExPEC+EPEC set" cannot be constructed — not for lack of data quality, but
because **the pathotype distinction operationally IS the sampling context** (extraintestinal vs
intestinal). The prior 24-genome k-mer baseline of 0.514 (= chance) was a correct null, not a tooling bug.

## Label provenance (the per-record `(predicted)` suffix is the flag)

Per [[feedback_per_record_label_provenance_audit]] — provenance was hidden in a suffix:

| Class | N | Share |
|---|---|---|
| curated / non-`(predicted)` (independent-ish) | 2,077 | 20.5% |
| `(predicted)` — marker/lineage-derived (CIRCULAR with resolver) | 5,295 | 52.2% |
| Not determined | 2,774 | 27.3% |

The 52.2% `(predicted)` exactly matches the prior ledger decision (row 143). **On the full collection H1
FAILS** (independent fraction 20.5% < 70% threshold). The curated subset is independent by construction
but is the minority.

## Curated subset — isolation source is decisive

| Class | N (curated) | Blood+Urine | Feces | distinct STs | distinct PopPUNK lineages |
|---|---|---|---|---|---|
| **ExPEC** | 1,574 | **1,567 (99.6%)** | 7 | 231 | 242 |
| **EPEC** | 269 | 0 | 269 (100%) | 84 | 79 |
| **ETEC** | 183 | 0 | 182 | 82 | 83 |

ExPEC labels = blood/urine isolation (extraintestinal) → genuinely **independent** of the virulence-marker
rules the resolver uses. EPEC/ETEC = feces (intestinal); their curated labels still ultimately rest on
markers (eae / enterotoxins) → residual circularity on those axes.

## The de-confound check — confound persists on ALL three axes

| Control axis | ExPEC vs EPEC overlap | Verdict |
|---|---|---|
| **Source study** | shared pool = 7 ExPEC + 59 EPEC (ExPEC ← Kallonen/Salipante bacteremia-UTI; EPEC ← Ingle/Hazen diarrheal) | disjoint — batch≈class |
| **Sampling context** | ExPEC = blood/urine (100%); EPEC = feces (100%) | disjoint **by definition of the pathotype** |
| **Phylogeny (PopPUNK)** | 5 shared lineages, 7 ExPEC + 20 EPEC strains | too small to balance |

ETEC is effectively single-source (von Mentzer 2014, 181/183) → batch-confounded with any contrast.

**Root cause (the durable lesson):** ExPEC and EPEC are *sampled from different body sites because they
cause different diseases*. The label (pathotype) and the batch variable (which study / which isolation
site) are the same axis. You cannot matched-pair your way out of a confound that is the definition of the
classes. This is intrinsic, not a data-collection gap — no amount of additional public data fixes it.

## Implications

1. **The handoff §3 "de-confounded ExPEC+EPEC cohort" is infeasible.** Don't keep hunting for it.
2. **The v0.1 embedding-classifier upgrade for pathotype is a dead-end on available labels** — not because
   NT embeddings are weak, but because no honest validation substrate exists for the discrimination task.
   Any AUROC would measure sampling context, exactly as the 24-genome confound did.
3. **This STRENGTHENS the v0 deterministic-resolver decision.** The marker-based compatibility resolver +
   abstention is the honest tool for pathotype; there is no validatable classifier uplift to chase here.
4. **Cross-project (roadmap-level) lesson:** for any phenotype *operationally defined by sampling context*
   (clinical site, disease presentation), study==class confounding is INTRINSIC and irremediable by
   within-collection matching. The embedding-upgrade frontier should target phenotypes whose labels are
   **sampling-independent lab measurements** (AMR MIC; quantitative assays) — NOT sampling-defined
   categories. See [[feedback_sparse_auroc_can_be_assembly_batch]] + [[feedback_threshold_vs_null_baseline_sanity_check]].

## What Horesh F1 IS still good for

- A large, independent-isolation-source **ExPEC positive set** (1,574 blood/urine; 242 lineages) — usable
  as a clade-diverse cohort for *resolver conformance* + abstention-quality evaluation (H4), and as a real
  out-of-cohort specificity check for the cross-axis ExPEC rule (the K=1 in-sample risk noted in
  `markers.py`). NOT for a discrimination AUROC.
- Per-lineage pathotype distribution for any future clade-balanced work (H3 infrastructure already exists
  via `scripts/mash_cluster_n147.py`).

## Provenance
`data/horesh_2021/F1_genome_metadata.csv` (Figshare 13270073, file 25552514, 4.83 MB). Analysis inline
(pandas value-counts + source/lineage set overlap). No model run; no AUROC computed (intentionally — the
finding is that an honest AUROC is not constructible for this task).
