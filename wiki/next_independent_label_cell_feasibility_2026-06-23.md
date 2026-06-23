# "Hunt the next HIV" — free-independent-label cell feasibility (2026-06-23)

Direction-1 feasibility check (user-chosen). HIV is the project's one cell validated against a FREE,
INDEPENDENT, isolate-level, non-circular wet-lab phenotype (Stanford HIVDB PhenoSense). This memo scores
candidate next cells against the SAME gate **before any build** (the TB lesson: verify the label exists +
is free + independent + non-circular FIRST). Label-gate evidence is verified; powering/build details are
build-phase work.

## The "next-HIV" rubric (a candidate must clear ALL 5)
1. **Free** — publicly downloadable, no DUA, no author-contact.
2. **Independent / non-circular** — a wet-lab MEASUREMENT (fold-change / IC50 / selection), NOT a
   genotype→prediction tool's own output (the circularity that killed the public TB label).
3. **Isolate/variant-level** — per-sample (or per-mutation) genotype ↔ measured phenotype, not aggregate.
4. **Determinant-tractable** — a target-site / catalog mechanism we can encode deterministically (the
   HIV RT/PR/IN, fungal ERG11, TB rpoB pattern).
5. **Genotype-downloadable** — sequence/variant callable from a reference (we already do BLAST→codon-map).

## Candidate scorecard
| Candidate | 1 Free | 2 Indep | 3 Isolate | 4 Determinant | 5 Genotype | Verdict |
|---|---|---|---|---|---|---|
| **SARS-CoV-2** nirmatrelvir (Mpro/3CLpro) + remdesivir (RdRp) — Stanford **CoV-RDB** | ✅ | ✅ | ✅ | ✅ | ✅ | **GREENLIT — the next NEW cell** |
| **HIV PI + INSTI** deconfounded v0.1 (deepen existing) | ✅ | ✅ | ✅ | ✅ | ✅ | **BUILD-READY NOW (zero label risk)** |
| Influenza NA-inhibitor (deepen existing NA cell) | ✅ | ✅ | ✅ | ✅ | ✅ | viable validation substrate (medium) |
| HCV DAA (NS3/NS5A/NS5B) | ⚠️ | ✅ | ⚠️ | ✅ | ✅ | DEFER — no single free Stanford-style isolate-level DB (scattered per-paper replicon data = the TB scattered-source risk) |
| HBV (pol/RT) | ⚠️ | ⚠️ | ⚠️ | ✅ | ✅ | DEFER — resistance is known-mutation; little free isolate-level measured-phenotype |
| DMS / ProteinGym / MaveDB | ✅ | ✅ | ✅ | ❌ | ✅ | OUT of this direction — free+measured but CONTINUOUS variant→fitness = the LEARNED/embedding shape (a closed 0-for-4 negative), not a deterministic R/S decoder |

## Verified evidence (this scan)
- **CoV-RDB is free + open + downloadable** (GitHub `hivdb/covid-drdb` SQLite daily + raw CSV
  `hivdb/covid-drdb-payload`, MIT). Same Stanford lab/data-shape as HIVDB.
- **It carries small-molecule target-site resistance data, lab-measured (not predicted):**
  `tables/invitro_selection_results/` has e.g. `agostini18` = remdesivir/GS-441524 → **RdRP L480, L557**;
  `jochmans23`/`iketani22c`/`fenwick22`/`heilmann22` = 3CLpro inhibitors → **Mpro 50F/166/167**; plus FDA
  EUA in-vitro tables. Schema `(ref, rx_name, gene, position, amino_acid)` = exactly the determinant
  catalog substrate. A separate susceptibility (fold-change) layer exists for validation.
- **HIV PI/INSTI datasets already on disk** (`data/raw/hiv/PI_DataSet.txt`, `INI_DataSet.txt`) — the
  deconfounded-v0.1 pattern already shipped for NRTI; PI/INSTI is pure engineering, no label risk.
- HCV/HBV: EuropePMC shows only scattered per-paper replicon datasets — no single free isolate-level
  genotype↔measured-phenotype DB surfaced (discovery-tier; would need a deeper dig before committing).

## Recommendation (ranked)
1. **Build the SARS-CoV-2 Mpro/nirmatrelvir (+ RdRp/remdesivir) cell — the next HIV.** Label gate CLEARED:
   free, independent, downloadable, determinant-tractable (single target enzyme, clean resistance positions).
   Same engine pattern as HIV/fungal (BLAST a committed Mpro/RdRp reference → codon-map → catalog). It puts a
   **3rd validated viral pathogen** on the board via the same free-independent-label mechanism that made HIV
   the crown jewel — and a 4th gene-target family.
   - **Make-or-break for the BUILD (not the label):** powering. Clinical nirmatrelvir/remdesivir resistance
     is RARE, so the independent per-isolate fold-change validation cohort is smaller than HIV's. Expect an
     honestly-UNDERPOWERED first cell (like our CAB/2nd-gen-INSTI cells) — the catalog will be well-founded
     from measured data; the validation N is the limit, reported as such. NOT a circularity problem.
2. **In parallel (zero-risk quick win): HIV PI + INSTI deconfounded v0.1.** Datasets in hand, validated
   pattern, finishes the HIV class set (the CLAUDE.md "deferred: mutant-specific deconfounded v0.1 for
   PI/INSTI" item). Certain, fast, no external dependency.
3. Defer HCV/HBV (feasibility-risky); skip DMS for this direction (wrong shape).

## Honesty tiering
- VERIFIED: CoV-RDB free/open/downloadable + small-molecule target-site measured data present (fetched the
  repo tree + two real CSVs). HIV PI/INSTI datasets present on disk.
- DISCOVERY-TIER (not yet verified, flagged): exact size/clinical-fraction of the CoV-RDB antiviral
  fold-change validation cohort; HCV/HBV free-DB existence. These are build-phase / deeper-dig items, not
  blockers to the recommendation.

## Provenance
Scan 2026-06-23 (Soraya, Direction-1). Sources: GitHub `hivdb/covid-drdb` + `covid-drdb-payload` (tree +
`agostini18-invitro.csv`, `jochmans23-invitro.csv`); EuropePMC (HCV); local `data/raw/hiv/`. WebSearch was
usage-filter-blocked on antiviral-resistance phrasing; routed via GitHub API + EuropePMC + repo CSVs.
