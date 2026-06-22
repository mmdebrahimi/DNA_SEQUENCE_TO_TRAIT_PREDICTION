<!-- memo-schema: 0.4 -->
# Non-circular phenotype label sources — ranked acquirability shortlist (2026-06-22)

> Supported memo. Validated from `noncircular-label-sources-2026-06-22.raw.md` (intake run by hand within the /research orchestrator). All rows medium-confidence (WebSearch-summary provenance — confirm by primary-source fetch before acquisition labor). Screened against the 8 rejection gates in `wiki/negative_results_map_2026-06-13.md` and the 4 acquisition criteria in `wiki/next_epoch_idea_anchor_prompt_2026-06-13.md`.

## Research Context (problem anchor)

The deterministic AMR decoder is the banked terminal product of the public-label track; the binding
constraint is **labels, not models**. This survey asked: is there an acquirable phenotype label source —
sampling-independent, non-circular (wet-lab/clinical measured, not genome-tool-derived), ≥~100-150
same-organism isolates with downloadable assemblies, and provenance-disjoint-feasible outside the US
surveillance ecosystem — that reopens the decoder beyond the saturated AMR grid?

**Headline finding (honest):** the survey largely **RE-CONFIRMS the wall** for E. coli. The free +
measured + genome-linked + non-circular pool is essentially the open INSDC / BV-BRC set the project already
exhausted. Two corrections + a small number of real leads emerged (below). **The bottleneck is acquisition,
not analysis — and most of the would-be "non-public clinical" unlock turns out not to exist** (bacterial
AST+genomes are open INSDC, not DUA-gated).

## Ranked shortlist (by acquirability; free first)

| # | Source | Measured & non-circular? | Genomes public? | N-after-assembly-filter (est.) | Outside US surveillance? | Gate verdict | Access path | Acquirability |
|---|---|---|---|---|---|---|---|---|
| 1 | **CRyPTIC TB MIC compendium** | YES — BMD MIC, 13 drugs, single assay | YES — EBI FTP + ENA | ~12,289 (June-2022) → 44,405 WGS / 56,405 pDST (2025) | YES — academic, 23–27 countries | clears G1/G3/G4/G5/G8; **in-distribution** for the project's WHO-catalogue TB cell (the catalogue was built partly from CRyPTIC) | **FREE download** (ftp.ebi.ac.uk/.../cryptic + Zenodo) | **HIGH — free, now** |
| 2 | **BV-BRC lab-derived AST** (project's existing substrate) | YES — filter `Lab Typing Method ≠ Computational Prediction` | YES — NCBI | E. coli already mined; provdisjoint **exhausted / surveillance-dominated** | mixed | clears G1; **trips G4** on the E. coli disjoint pool | FREE API/download | already used (saturated) |
| 3 | **ARESdb / 11,087 clinical isolates** (Ares Genetics) | YES — culture-based AST MIC, 21 drugs | **NO — VERIFIED proprietary** (GEAR-base: registered academic-only; "cannot make all data fully accessible for batch download"; raw seqs need a special request; no public accessions; Ares Genetics "sole owner") | ≤11,087 / 18 species; E. coli subset NOT stated | YES — "concerted", non-NARMS, bias-reduced | **trips G5 + acquisition wall** (no downloadable genome+label pairs) | **commercial / academic-MTA** | **LOW — closed for a solo non-academic project** |
| 4 | **von Mentzer ETEC / Whittam DECA** | NO — label is toxin-gene/clade-defined | YES — ENA PRJEB33365 (~362) | ~362 (ETEC) | YES | **trips G1/G3** (circular pathotype label) | FREE download (ENA) | genomes free, **label circular → no unlock** |
| 5 | **Pfizer ATLAS / SENTRY** | YES — BMD MIC, 6.5M MICs / ~1M isolates | NO — Pfizer owns isolates; US excluded | ~0 joinable (no public WGS) | surveillance-adjacent | **trips G5** | FREE (registration) | measured but **no genomes → unusable** |
| 6 | **EUCAST MIC distributions** | YES — but aggregate | NO | 0 (aggregate per-species×drug only) | mixed | **trips G5 / no linkage** | FREE | aggregate only → unusable |
| 7 | **BacDive** | YES — biochemical | partial (mostly type strains) | <100 / organism (type-strain depth) | YES | **trips organism-depth** | FREE API | depth-limited |

**Two corrections to the anchor's assumptions** (verified this pass):
- von Mentzer ETEC genomes are **public on ENA** (PRJEB33365 / ERP004228), **not** MTA-gated — but the pathotype label is circular, so it still doesn't unlock pathotype decoding.
- There is **no DUA-gated bacterial AST treasure trove**: EGA/dbGaP controlled access is for *human* data; bacterial WGS+AST live in *open* INSDC. The "non-public clinical label" lane is much narrower than assumed.

## Decisions for Human Confirmation

| # | Decision (authority / money — yours) | Candidate use | Verification needed before acting |
|---|---|---|---|
| 1 | **TB independent validation via CRyPTIC** — free + downloadable now. The decoder's TB cell exists; an INDEPENDENT score needs a determinant rule NOT derived from CRyPTIC OR a hand-curated post-2023 gold set (the project's existing TB deliverable-(b), `wiki/tb_independent_goldset_acquisition_2026-06-17.md`). This is **compute-bound (the ~1.6 TB regeno → D: wall), not acquisition-bound.** Highest-leverage FREE move. | independent-label validation of the shipped TB cell | confirm a non-CRyPTIC-derived rule source; provision the D: regeno compute |
| 2 | **ARESdb / 11,087-isolate dataset** — **VERIFIED proprietary** (GEAR-base academic-only, no batch download, no public genome accessions; Ares Genetics sole owner). Effectively **closed** for this solo non-academic project barring a commercial/academic-MTA arrangement. | would have been a new measured-MIC+WGS non-surveillance substrate, but the genome+label pairs are not downloadable | DONE — verified PMC6624217 data-availability; no further free path |
| 3 | **Accept public-label saturation for E. coli** → pivot to **prospective-lock validation** of the frozen decoder (needs no new label; `wiki/reproducibility_freeze_2026-06-13.md`). | validate the existing decoder on later-arriving independent data, no acquisition | none — executor-eligible whenever you choose |
| 4 | **Contact-gated reference collections** (Whittam DECA; any clinical lab AST export). I can DRAFT the request; **sending is yours** (send/authority gate). Only worth it if paired with a measured (non-definitional) phenotype. | a bespoke non-circular cohort | identify a collection whose label is a *measurement*, not a gene-definition |

## Honest gaps

- No new **free E. coli** measured-phenotype+genome source surfaced — re-confirms "labels not models" for the project's home organism.
- ARESdb access model (free vs MTA vs paid) + its E. coli depth are **unverified** — the one lead that could be new requires a primary-source fetch + a money/authority decision.
- Non-AMR quantitative-assay phenotypes (growth/fitness/biochemical) at ≥100 same-organism depth with public genomes: **none found** (open gap, not a located source).

Citation discipline: all source URLs are real https/FTP; scale figures cross-attested; medium-confidence pending primary-source confirmation.
