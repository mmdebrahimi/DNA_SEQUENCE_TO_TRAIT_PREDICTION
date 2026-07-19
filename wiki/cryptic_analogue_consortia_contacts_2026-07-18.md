# CRyPTIC-analogue consortia + data-access contacts — the named, gate-screened shortlist

**Date:** 2026-07-18 · **Status:** DECISION AID (acquisition = user authority; but see the headline below —
the two top targets are largely FREE/public, so the "money-gated" framing softens) ·
**Supersedes-by-sharpening:** `wiki/data_acquisition_voi_memo_2026-07-18.md` (that memo named the *shape*;
this one names the *programs, access routes, and contacts* and screens each against the 8 gates) ·
**Grounding:** a 2026-07-18 targeted web-research pass (Vivli AMR Register, CDC/FDA AR Isolate Bank, NARMS,
Euro-GASP/PubMLST, GHRU) + `wiki/negative_results_map_2026-06-13.md` (the 8 rejection gates).

## Headline (the update that matters)

The prior memo framed acquisition as a money/partner call. Research sharpens that: **the two highest-VOI
targets are largely FREE and already reachable** — a CDC/FDA request (no cost for data + isolates) and a
convergent-pathogen WGS+MIC corpus (gonorrhoea) that is mostly *public* in ENA/PubMLST. And it surfaces one
**trap**: the biggest measured-MIC platform (Vivli/ATLAS, 6M isolates) is the **wrong shape** for us — it has
MIC but **no full genome sequences**, so our decoder pipeline can't run on it. The label wall's cheapest
opening is a *request/fetch*, not a licence.

## The join our decoder needs (the screen)

CRyPTIC's power was **measured BMD-MIC ⋈ full WGS, per isolate**. A candidate must supply BOTH — a real
phenotype measurement (G1: not a genome-derived call) AND a downloadable genome (not just resistance-gene
genotype). "MIC + β-lactamase-gene list" is NOT enough — the genome is what `call_resistance` +
the supervised complement consume. Convergent-evolution organism is the VOI multiplier (the supervised
complement is HIV-proven on convergent resistance; it transfers).

## Ranked, named shortlist

### 🥇 Tier 1a — CDC & FDA Antimicrobial Resistance (AR) Isolate Bank — *free, actionable now*
- **What:** a curated biorepository of bacterial + fungal isolates. Each isolate has **reference broth
  microdilution MIC** (CLSI method — real G1 phenotype) + interpretation + known mechanisms, and — for
  sequenced isolates — a **BioSample accession linking to the WGS** (so MIC ⋈ genome is provided). CDC/FDA
  commit to sequencing *all* bank isolates.
- **Gate screen:** G1 ✓ (BMD, not a genome call). G5 ✓ (genomes exist — no assembly-attrition). Physical
  isolates ship free → you can **re-measure MIC yourself**, which *fully* de-circularizes the label (the
  strongest possible G1 clearance). Caveats: the panel is **resistance-enriched/selected** (great for a
  determinant catalog, but not prevalence-representative → don't quote population sens/spec off it) and
  **clinical metadata is intentionally sparse** (IRB avoidance) → provenance partial (G7). Not dedup-fragile
  because it's already a de-duplicated curated panel.
- **Access:** register with **institutional** contact at `wwwn.cdc.gov/ARIsolateBank` → order up to 5 panels
  → data (MIC + mechanisms + BioSample/WGS links) is **free**; only FedEx/UPS shipping (dry ice, ~15–20 lb)
  is on the requester. Restricted to institutions with a biosafety officer; strict no-redistribution.
- **Contact:** **`ARbank@cdc.gov`**. Convergent targets in-bank: *Neisseria gonorrhoeae*, Enterobacterales,
  *Acinetobacter*, *Candida auris* (already a validated fungal cell).

### 🥇 Tier 1b — Euro-GASP / GISP gonorrhoea corpus (via ENA + PubMLST) — *convergent + mostly public*
- **What:** *N. gonorrhoeae* is a textbook **convergent-resistance** pathogen. The European Gonococcal AMS
  Programme (Euro-GASP) + US GISP have produced **WGS ⋈ MIC** collections: the 2013 pan-European genomic
  survey (~20 countries), a global 481-isolate WGS+MIC set (agar-dilution MICs, 5 countries, 2004–2019), and
  large aggregations (tens of thousands of genomes). Genomes in **ENA**; MICs + NG-STAR typing in article
  supplements and **PubMLST (`pubmlst.org/neisseria`)**.
- **Gate screen:** G1 ✓ (measured MIC). Convergent ✓ → supervised complement transfers. Less
  surveillance-single-clone-dominated than Enterobacterales bloodstream (G4 softer). Main friction is
  **stitching MIC (supplements/PubMLST) to reads (ENA)** per study, and agar-dilution vs BMD method mix
  (note the method per source). MIC censoring at breakpoints applies (G6) — model operator-aware, as we do.
- **Access:** ENA (open) + PubMLST *Neisseria* (open, hosts genome+MIC+NG-STAR) + per-paper supplements.
  Mostly a **fetch**, minimal acquisition. **Contacts:** PubMLST *Neisseria* is Oxford-curated (Jolley/Maiden
  group, `pubmlst.org` contact form); Euro-GASP is ECDC/UKHSA-coordinated; GISP is US CDC.

### 🥈 Tier 2 — NARMS BMD subset (FDA/USDA/CDC) — *a free re-cut, partly already-mined*
- **What:** enteric convergent pathogens (*Salmonella*, *E. coli*, *Campylobacter*) with **Sensititre BMD
  MIC** (15 drugs) + WGS in SRA; MIC↔accession tables in the NARMS ML/genotype papers and in **NCBI Pathogen
  Detection**. This is a *real* CRyPTIC-analogue that is already public.
- **Gate screen:** the honest caveat — this overlaps the **NCBI Pathogen Detection structured-antibiogram**
  source our 10 SCORED cells + prospective-lock already use → **surveillance-domination (G4)** for the drugs
  we tested. BUT the **BMD-method-labelled NARMS subset** is a *sharper, cleaner cut* than the broad PD pull
  (method-known, exact MIC). Low-cost lever: re-mine NCBI PD filtered to NARMS BMD records for a drug/organism
  the deterministic decoder doesn't yet cover — not a new acquisition, a better filter.
- **Access:** free — NCBI SRA + Pathogen Detection + FDA "NARMS Now" + paper supplements.

### ⛔ Deprioritize — Vivli AMR Register (Pfizer ATLAS, Merck SMART, Venatorx GEARS, Shionogi SIDERO-WT, Paratek KEYSTONE, GSK SOAR, …)
- **Why it's the trap:** enormous **measured MIC** at scale (ATLAS alone >900k–6M isolates, free via request),
  BUT the programs provide **MIC + at most a resistance-mechanism genotype (β-lactamase genes only for ATLAS;
  none for most)** — **no full WGS**. Our determinant/embedding pipeline needs the genome; a gene list can't
  feed it. So Vivli is **the wrong data shape for the decoder** (excellent for MIC-trend / spatiotemporal
  modeling — the Vivli Data Challenge space — not for us).
- **Access (for the record):** free request via the **Vivli AMR Register** (`amr.vivli.org` /
  `searchamr.vivli.org`) — sign up, 300-word EOI, team form; Merck SMART redirects to `globalsmartsite.com`;
  ATLAS also at `atlas-surveillance.com`.

### ⛔ Deprioritize — GHRU / Pathogenwatch (Wellcome Sanger / CGPS)
- WGS-rich LMIC *Klebsiella*/*E. coli* surveillance, but the AST is **submitted-metadata-dependent** and the
  headline output is **Kleborate/Pathogenwatch = genotype-derived** → **circular (G1)** where it's dense.
  Same verdict as the prior memo. (`pathogen.watch`; CGPS at the Big Data Institute, Oxford.)

## What this changes for the decision

1. **The cheapest opening isn't money — it's a request/fetch.** Tier-1a (CDC data free; email `ARbank@cdc.gov`)
   and Tier-1b (ENA/PubMLST, mostly open) are reachable without a licence or a paid partnership. That lowers
   the bar the prior memo set.
2. **A named trap is now ruled out:** don't chase Vivli/ATLAS for the decoder — no WGS. (Revisit only if a
   contributor ever adds genomes.)
3. **Convergent-first is concrete:** gonorrhoea (Euro-GASP) is the convergent pathogen where the
   supervised-complement bet has the best transfer odds and the data is closest to free.

## The concrete next step (still user authority, but cheaper than framed)

Two low-cost moves, both a *request/fetch* rather than a purchase:
- **(a)** email **`ARbank@cdc.gov`** from an institutional address for a gonorrhoea/Enterobacterales panel's
  MIC + BioSample/WGS links (free data), OR
- **(b)** fetch the **Euro-GASP global 481-isolate WGS+MIC** set (ENA + supplement) + the PubMLST *Neisseria*
  MIC/NG-STAR metadata and run the existing preflight → `call_resistance` → supervised-complement → report
  card, exactly the validated pipeline.
Everything downstream of a clean measured-AST ⋈ WGS set is already built; the moment such a set is in hand the
decoder extends by construction.

## Honest scope

This memo names programs + access routes + contacts and screens them against the project's own gates; it does
NOT perform a request or acquisition (that + any institutional/biosafety attestation = user authority). The
free-vs-gated read is a genuine update to the prior memo, not a reopening of the closed free-public-MIC
negative — the difference is **WGS-linked** measured MIC from a *specific named* program vs the broad
surveillance pull that tripped G4. Frozen decoder surface untouched (no code changed).

## Sources

- Vivli AMR Register / programs: https://amr.vivli.org/resources/research-programs/ ·
  https://amr.vivli.org/faq/atlas/ · https://searchamr.vivli.org/
- CDC & FDA AR Isolate Bank: https://wwwn.cdc.gov/ARIsolateBank · https://wwwn.cdc.gov/ARIsolateBank/QA ·
  FDA overview https://www.fda.gov/medical-devices/in-vitro-diagnostics/cdc-fda-antibiotic-resistance-isolate-bank ·
  JCM 2017 https://journals.asm.org/doi/10.1128/jcm.01415-17 · contact `ARbank@cdc.gov`
- NARMS WGS+BMD-MIC: https://pmc.ncbi.nlm.nih.gov/articles/PMC6355527/ (ML MIC, per-isolate MIC↔SRA) ·
  https://journals.asm.org/doi/10.1128/aac.01030-16 (McDermott 2016) · NCBI Pathogen Detection
  https://www.ncbi.nlm.nih.gov/pathogens
- Euro-GASP / gonorrhoea WGS+MIC: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6010626/ (2013 survey) ·
  https://www.omicsdi.org/dataset/biostudies-literature/S-EPMC10319951 (global 481-isolate WGS+MIC) ·
  PubMLST https://pubmlst.org/organisms/neisseria-spp
- GHRU / Pathogenwatch: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8634497/ · https://pathogen.watch/
