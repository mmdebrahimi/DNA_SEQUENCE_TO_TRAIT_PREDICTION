# Data-source census across the tree of life — new decoder-cell substrates (2026-07-13)

**Question.** What NEW *free, independent, isolate-level* genotype→phenotype label sources exist across **every family of life**, to seed new cells of the deterministic DNA decoder?

**Method.** 8 parallel kingdom scouts (WebSearch + WebFetch), each screening candidates against the decoder's winning **regime** (a curated target-site mutation catalog meeting a free/independent/isolate-level/wet-lab phenotype label) and the **8 rejection gates** from `wiki/negative_results_map_2026-06-13.md`. Machine-readable companion: `wiki/data_source_census_tree_of_life_2026-07-13.json` (22 scored candidates).

**Gates.** G1 circular-label (label from a genome tool) · G2 study==class · G3 sampling-defined label · G4 surveillance domination · G5 assembly/access attrition · G6 measurement censoring at breakpoint · G7 provenance not separable · G8 dedup/clonality collapses balance.

---

## The load-bearing meta-finding

**Across nearly every kingdom, the curated catalog (regime fit) already exists — but the FREE + INDEPENDENT + ISOLATE-LEVEL label is the binding constraint.** This is the project's north-star lesson ("labels not models") reconfirmed across the entire tree of life, not just AMR. Two label archetypes clear the gates:

- **(a) a measured MIC / EC50 / bioassay that is INDEPENDENT of the catalog** (TB CRyPTIC MIC, N. gonorrhoeae MIC, CMV recombinant-phenotyping fold-change, nematode egg-hatch assay, WWARN clinical outcome); or
- **(b) DMS / deep mutational scanning** (MaveDB/ProteinGym, baloxavir-PA) — but this is the *molecular-property* regime where a learned scorer (ESM) can win and a deterministic catalog is a less natural fit (see the project's G2P regime-boundary lesson).

The recurring **trap** is the circular label (G1): almost every public "resistance database" (geno2pheno-HCV/HBV, HCV-GLUE, CPIC, ClinVar, Pf7's own inferred resistance status, PBP-typed pneumococcal MIC, gene-presence VRE) is the *catalog's own interpretation*, and scoring against it proves nothing. The honest label must be a separate wet-lab/clinical measurement.

---

## Ranked shortlist — the highest-VOI NEW cells

| # | Kingdom | Organism × phenotype | Label (independent?) | Regime fit | Dominant gate risk | Verdict |
|---|---|---|---|---|---|---|
| **1** | Bacteria | **M. tuberculosis** × second-line + new drugs (bedaquiline, linezolid, moxifloxacin, amikacin, ethionamide, clofazimine, delamanid) | wet-lab MIC (CRyPTIC, ~12k–21k isolates) ✅ | **STRONG** — WHO catalogue v2 **already pinned in-repo**; same VCF pipeline as shipped RIF/INH | G6-minor (ECOFF), G8 (TB clonal; barcode collapse already in-repo) | **IMMEDIATE — lowest-friction new cells on the tree** |
| **2** | Virus | **CMV** × ganciclovir (UL97), foscarnet/cidofovir (UL54), letermovir (UL56) | wet-lab EC50 fold-change via **marker transfer** (Chou catalog) ✅ clears G1 | **STRONG** — drops into the existing HIV/SARS-CoV-2 "catalog scored vs fold-change" harness | per-mutation not per-isolate | **TOP VIRAL — build first among viral; exact match to the validated pattern** |
| **3** | Bacteria | **N. gonorrhoeae** × ceftriaxone/azithromycin/cipro | wet-lab MIC (Euro-GASP + Pathogenwatch, several thousand) ✅ | **STRONG** — mature curated penA/mtrR/gyrA/23S catalog + downloadable ENA assemblies | G4 (Euro bias — pool GISP), G8 (clonal complexes) | **HIGH — new bacterial genus, clean measured MIC** |
| **3.5** | Animalia (human) | **MaveDB / ProteinGym DMS** (BRCA1, PTEN, TP53, TPMT…) | DMS functional score (wet-lab MAVE) ✅ G1-clean by construction | cleanest free+independent+wet-lab of all — **but molecular-property regime** | (regime, not gate) — learned-scorer territory | **CLEANEST LABEL; partly in-repo (ProteinGym) — a learned-scorer track, not a catalog cell** |
| **4** | Protozoa | **P. falciparum** × sulfadoxine-pyrimethamine (dhfr/dhps) | WWARN clinical outcome + pyrimethamine IC50 ✅ (NOT Pf7's inferred label) | **EXCELLENT** — cleanest enzyme active-site catalog after Kelch13; MalariaGEN Pf7/Pf8 free genotype at scale | **G1 if you use Pf7's inferred status** (must pair with WWARN), G4/G8 | **TOP PROTOZOA — cost = building genotype↔outcome linkage** |
| **5** | Fungi | **Aspergillus fumigatus** × azole (CYP51A TR34/L98H) | wet-lab EUCAST MIC (poolable free WGS studies) ✅ | **EXCELLENT** — FunResDB + FungAMR + MARDy catalog; TR-mechanism is target-site/promoter | G4/G8 (TR34 clonal sweeps — dedup mandatory), G5 (some reads-only) | **TOP FUNGAL — closest to a FREE fungal MIC label** (existing fungal cells are `no_free_source`) |
| **6** | Virus | **Influenza A** × baloxavir (PA I38/L106R) | DMS + focus-reduction EC50 ✅ | HIGH via DMS — small target ideal for mutant-level catalog | engineered library (not isolate-level) | **HIGH via DMS — modern PhenoSense-analog for small viral targets** |
| **6.5** | Invertebrate (nematode) | **Haemonchus contortus** + human STH × benzimidazole (β-tubulin 167/198/200/134) | wet-lab egg-hatch / larval-development EC50 / FECRT ✅ clears G1/G3 | **BEST invertebrate** — mature curated catalog + independent bioassay + single target gene | G8 (mostly pooled larvae), G6 (LDA dose-censoring) | **TRUE AMR ANALOG in a NEW KINGDOM — needs data-engineering (single-worm subset / harvest supplements)** |

Full 22-candidate list with URLs, N estimates, and per-candidate gate screens: the JSON sidecar.

---

## Per-kingdom summary

- **Bacteria** — richest new-cell vein. **TB second-line (CRyPTIC)** is the single lowest-friction win (same pipeline + WHO catalogue already in-repo). **N. gonorrhoeae** is the strongest new-genus candidate. **S. pneumoniae (GPS)** is huge but conditional — must filter to the *measured*-MIC subset (default PBP→MIC is G1). H. pylori has the cleanest single-SNP catalog but N is fragmented. Enterococcus vanco is a gene-presence identity (low novelty).
- **Archaea** — **AMR/isolate side is a truthful CLOSED NEGATIVE** (no CLSI breakpoints, no catalog, no usable isolate label; substantial archaeal data — Haloferax/M. maripaludis transposon fitness, Sulfolobus forward genetics — is the wrong regime; Madin trait DBs are taxon-confounded G3/G7). **But the DMS side is NOT a total void** (correction verified 2026-07-13 against the vendored ProteinGym v1.3 CSV, not inferred): there are **exactly 2 archaeal DMS assays already in ProteinGym**, both *Saccharolobus (Sulfolobus) solfataricus* (`SACS2`) — `DN7A_SACS2_Tsuboyama_2023` (Sso7d stability, ESM2-650M ρ=0.337) and `TRPC_SACS2_Chan_2017` (indole-3-glycerol-P synthase, ρ=0.653). The thermophile look-alikes (Thermotoga/Thermus/Aquifex) are correctly excluded as bacteria. So archaea is a **low-priority molecular-property (DMS) extension with a tiny in-repo substrate**, NOT an acquisition-only gap — the regime-compatible path already exists, it just has N=2.
- **Fungi** — the catalog side is strong (FungAMR 35k entries, MARDy, FunResDB); the label is distributed across per-study supplements, not aggregated. **A. fumigatus azole (CYP51A)** is the closest to a free MIC label. **C. glabrata echinocandin (FKS1/FKS2)** has the best label density and a clean target-site catalog (azole is efflux-confounded). Cryptococcus (heteroresistance/aneuploidy) and DMS are regime-mismatched.
- **Protozoa** — **P. falciparum SP (dhfr/dhps)** is the clear winner: free MalariaGEN genotype at scale + independent WWARN clinical/IC50 phenotype (critically, NOT Pf7's own inferred resistance status = G1). Piperaquine (plasmepsin2/3 CNV) and mefloquine (pfmdr1 CNV) are CNV-driven + GMS-clonal. P. vivax / Leishmania / Trichomonas / Trypanosoma all reject (no catalog, or label-phenotype decoupling).
- **Viruses** — every public "database" (geno2pheno, HCV-GLUE, HBVdb, FluSurver) is a **G1-circular rules engine**. The independent phenotype lives in two non-DB forms: **recombinant/marker-transfer phenotyping** (CMV UL97/UL54/UL56 — the top viral pick, matches the HIV/SARS pattern exactly) and **DMS** (baloxavir-PA). MPXV F13L is the only isolate-level+measured viral find but small-N + G3 + access-uncertain. HCV/HBV have catalogs but no free independent phenotype set (unlike HIV's PhenoSense).
- **Invertebrates** — **nematode β-tubulin benzimidazole** is the true AMR analog (independent bioassay clears G1/G3; mature 167/198/200/134 catalog), blocked mainly by pooled-larvae design (G8) and data being un-aggregated. Insect insecticide resistance (kdr/ace-1/Rdl) has strong catalogs but a **structural G3/G8**: Ag1000G has individuals but no phenotype; VectorBase PopBio pairs genotype+phenotype only at *cohort* resolution (ecological regression, not an isolate cell).
- **Plants (weeds)** — regime mirrors AMR (Murphy & Tranel 2019 = the curated ALS/ACCase/EPSPS/PSII/PPO rule), but the binding constraint is again labels: the dominant gate is **G8 (pooled populations, not isolates)** plus **G1** (many TSR-diagnostic papers *call* resistance from the mutation). weedscience.org is a catalog/ontology (aggregate counts), not a paired table. Cleanest pairing = Lolium Chile whole-plant dose-response + tri-gene genotyping (small N).
- **Animalia / Human** — **MaveDB/ProteinGym DMS** is the cleanest free+independent+wet-lab label of the entire census, but it is the *molecular-property* regime (learned-scorer territory, already partly in-repo). **OMIA** livestock/companion single-gene traits are a real catalog + observed phenotype but hit **G8 breed-clonality** (the project's clonality-inflation trap) at small N. **PharmGKB/CPIC** = catalog-circular + independent label access-gated (dbGaP/biobank, G5). **ClinVar** = G1-severe (labels embed undisclosed computational evidence).

---

## Recommended next moves (for user ratification — this is a scouting memo, not a build)

1. **TB second-line (CRyPTIC)** — the immediate, lowest-risk expansion: same substrate + the WHO catalogue already pinned in-repo. New cells: bedaquiline, linezolid, moxifloxacin, amikacin.
2. **CMV UL97/UL54/UL56** — the top viral expansion; drops into the existing HIV/SARS fold-change harness with an independent (non-circular) label.
3. **N. gonorrhoeae** — the strongest new-bacterial-genus candidate with clean measured MIC.
4. **P. falciparum SP (dhfr/dhps)** — highest-VOI parasite cell, contingent on building the MalariaGEN↔WWARN linkage (and explicitly avoiding Pf7's circular label).

Each of the above is an ACQUISITION/BUILD decision (cohort-assembly labor), not a plug-in. The three closed/blocked findings (archaea = no label; HCV/HBV/CPIC/ClinVar = circular or access-gated; insect/weed = structural G3/G8) are recorded so the tree isn't re-walked.

*Discovery-tier: several dataset accessions + N estimates are flagged unverified in the JSON (per-candidate) and in the scout outputs — confirm data-availability before any cohort build. Content-filter false-positives blocked a few resistance-topic web queries this session (Vibrio, Mycoplasma genitalium, some mechanism queries); those cells are marked unverified.*
