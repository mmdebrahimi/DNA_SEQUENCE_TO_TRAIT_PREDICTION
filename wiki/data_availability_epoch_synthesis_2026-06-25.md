# Data-availability epoch — synthesis + acquisition brief (2026-06-25)

Banks the "where/how do we get the data we need" arc and turns the standing frontier (non-public
acquisition) into an actionable spec. The acquisition RELATIONSHIP is the user's decision; this is the
gate-clearing data SPEC + candidate source types + the ask — ready for when that decision is made.

## What the epoch delivered (the "labels not models" wall cracked 3 ways)
| Win | Number | Tier |
|---|---|---|
| **Pneumo capsular serotype vs phenotypic Quellung** | serogroup **0.939** / exact 0.661 (n=230) | INDEPENDENT measured-label (the GPS Poland Quellung cohort) |
| **Pneumo gene-presence AMR** (macrolide/tet) | macrolide **0.961** / tet **0.932** | NEAR-INDEPENDENT (label independent; gene-presence ~AMRFinder-equiv) |
| **(earlier) HIV PhenoSense** | the prior epoch's win | INDEPENDENT measured-label |

Plus two new GREEN typing cells shipped (Salmonella serovar + pneumo serotype callers), the pneumo
β-lactam breakpoint foundation, and the 4-reservoir data-acquisition map. **Key reframe:** the labels wall
is NOT uniform — phenotypic serology + viral phenotype DBs crack it cleanly; bacterial-AMR-MIC is saturated
for public data; and even a "free measured label" AMR cell (pneumo) carries a real engine + breakpoint +
single-cohort cost.

## The 4-reservoir verdict (free independent measured labels)
1. **Phenotypic serology (Quellung) → WIN** (serogroup 0.939). The cleanest free crack this epoch.
2. **TB CRyPTIC MIC → gated** (independent number needs a hand-curated post-2023 gold set).
3. **HCV/HBV viral → no free measured dataset** (geno2pheno = rule tool, faithful-to-tool only).
4. **Bacterial AMR MIC → closed for public** (the 8-gate negative-results map; only acquisition reopens it).

## Deferred (named, scoped — NOT abandoned)
- **Pneumo β-lactam AMR cell** — the CDC PBP-type→MIC engine (`wiki/pneumo_betalactam_assessment_2026-06-25.md`).
  Target proven (penicillin@meningitis 0.973); a multi-session Docker-gated build. Foundation banked.
- **Salmonella serovar real DB** — SeqSero2 antigen-DB engineering (`wiki/salm_serovar_report_card.md`).
- **Independent-AMRFinder swap** for the pneumo gene-presence cell (near→fully independent; Docker, post-ktype).

---

## ACQUISITION BRIEF — the frontier, made actionable (user-owned relationship)

**Why acquisition is THE frontier:** the negative-results map proves public bacterial-AMR is label-bound by 8
gates. A non-public, isolate-level, wet-lab-measured source clears them BY CONSTRUCTION — it is the only move
that both (a) gives a fully-independent number and (b) reopens the learned-model bet that 0-for-4'd on public
data (which always learned population structure, never the causal signal, because public labels are
study-confounded).

**The data shape that clears all 8 gates** (screen any offered source against this):
- **Isolate-level records**, each = `{genome (assembly or reads)}` + `{a MEASURED phenotype}` + `{provenance fields}`.
- **Measured phenotype = wet-lab** (BMD/ETEST MIC, disc zone, or a functional assay) — NOT a gene-call/tool
  output (clears G1 circular), NOT a sampling descriptor like isolation-site (clears G3 sampling-defined).
- **Provenance fields populated** (submitter / lab / collection date / source) so a leakage-clean
  provenance-disjoint or temporal split is buildable (clears G7).
- **Lineage diversity** — ideally >1 BioProject/clone per class, else disclose clonality (clears G8).
- **Powered** — ≥ ~20–30 per class per drug after QC (clears G4 surveillance-domination + powering).

**Candidate source TYPES** (the user picks the actual relationship; not specific orgs):
- Hospital/clinical **microbiology labs** — routine genome + AST, the richest measured-label stream.
- **National surveillance programs** with a data-sharing arm — the RAW lab AST + assemblies (not the
  already-public dedup'd dump that trips G2/G4).
- **Culture biobanks** (e.g. type-collection repositories) that ship strains WITH measured AST.
- **Academic collections** / a collaborator's curated cohort with measured phenotype + metadata.
- **Diagnostics companies** with genome+AST panels.

**The ask (template):** "Isolate-level {genome assembly or reads} + {measured AST/MIC or assay result} +
{collection metadata: lab, date, source}, N≥30 per drug/class, with permission to publish AGGREGATE
validation metrics (no patient identifiers)." A modest first cohort (one drug, one organism, one lab) is
enough to prove the fully-independent cell.

**What it unlocks:** the first fully-independent, learned-model-eligible cell — and reopening the
bacterial-AMR expansion the public-data map closed.

## Status
- **Scope of the "acquisition is the frontier" claim (corrected 2026-06-25):** acquisition is the frontier for
  a **FULLY-INDEPENDENT bacterial-AMR learned-model** cell specifically — NOT a global "no more free data"
  claim. Several **public/code-closable** paths remain open: TB CRyPTIC gold-set curation, the pneumo
  AMRFinder swap (near→fully-independent), fuller GPS measured-AST exploitation, and the β-lactam PBP engine.
- **Banked.** The epoch's wins are committed; the deferred cells are scoped with named walls; the acquisition
  frontier is now a concrete spec + ask, awaiting the user's relationship decision.
- **No code/research task remains on the acquisition path** — it is, by design, a relationship decision. The
  next *code* work is any of the named deferred cells, when chosen.
