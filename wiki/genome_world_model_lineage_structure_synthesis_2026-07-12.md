# Genome world model — lineage-structure of the multi-axis co-occurrence network (2026-07-12)

Autonomous Soraya day-run synthesis. Consolidates today's lineage de-confound work into one organizing
finding, grounds it in the literature, and states the decoder implications. **Frozen decoder surface
byte-unchanged throughout** (`verify_lock OK`); all analyses local/CPU on cached assemblies, no new
labels / GPU / embeddings.

## The question

The overnight multi-axis-C run (ledger rows 399–400) found that an E. coli genome's AMR + plasmid + virulence
features predict each other at high AUC (0.79–0.98) — resistance, plasmids, and virulence are NOT independent.
But the honest caveat was: **is this a real biological co-occurrence, or just clonal co-inheritance?** (ExPEC
ST131-type clones carry everything together, so any within-cohort model can "predict" one axis from another
by implicitly recognizing the clone.) The proper control is to hold out whole lineages.

## Method

Mash-cluster the 240 cached E. coli/Shigella genomes (→ 85 clades at Mash 0.005, largest clade 34.2%), then
re-run each axis's cross-prediction under **leave-one-clade-out (GroupKFold) CV** — a genome's entire clade is
absent from its training fold, so the model cannot exploit clonal identity. A clade-surviving AUC means the
co-occurrence generalizes *across* lineages (real beyond clonal co-inheritance); a collapse means it was
lineage-memorization. `scripts/crossaxis_lineage_deconfound.py` (`--target-axis {determinant,virulence,plasmid}`),
9 CI-safe tests.

## The finding — a coherent 3-axis ordering

| axis (predict → from) | verdict | median clade-grouped AUC | reading |
|---|---|---|---|
| **determinant ↔ determinant** | GENERALIZES | **0.908** (39/45, E. coli n=232) | resistance-gene co-occurrence is real beyond lineage |
| **AMR+plasmid → virulence** | SPLIT | **0.676** | core functions generalize; accessory genes clonal |
| **AMR+virulence → plasmid** | LINEAGE-MEDIATED | **0.615** (8/23) | plasmid backbone identity is clade-fixed |

(The determinant number is 0.908 on the harvest-consistent E. coli cohort of 232 genomes; 0.922 on the
240-genome virulence-cache cohort — same GENERALIZES verdict either way.)

The three axes fall in a clean, mechanistically-sensible order — **0.922 > 0.676 > 0.615** — and each end of
the ordering tells a distinct, non-obvious story:

- **Resistance-gene CONTENT transfers across lineages (0.922).** AMR determinants co-occur in cassette/integron
  blocks (sul1→aadA/dfrA, QRDR clusters) that predict each other even when the target's whole clade is held
  out. The cassette is the mobile unit; its internal linkage is intrinsic, not clonal. This **upgrades the
  C-core co-occurrence finding (row 397) from a determinant-profile-dedup *proxy* to a real *phylogenetic*
  control — and it passes strongly.** The only exceptions are chromosomal QRDR point mutations (`parC_E84V`,
  `gyrA_*`) whose positives all live in one clade → flagged `clade_concentrated` (they are lineage signatures,
  the determinant analog of hlyA-in-C2 below).

- **The plasmid VEHICLE is lineage-locked (0.615).** This **falsified the pre-registered prediction** that
  mobile conjugative plasmids would generalize *more* than chromosomal islands. Instead the specific IncFII
  sub-variants (pRSB107 0.867→0.369, pHN7A8 0.843→0.353, pCoo 0.793→0.24, pAMA1167-NDM-5 0.897→0.404) collapse
  hardest; only broad-host-range backbones (IncN, IncX1, IncFIA/B) survive. E. coli plasmid↔host pairings are
  clade-fixed by stable co-inheritance (post-segregation killing / addiction systems keep an IncF backbone
  vertically maintained within a lineage). *Which* resistance rides *which* plasmid backbone is a clade
  signature, not free horizontal mixing.

- **Virulence splits by prevalence (0.676).** High-prevalence core functions (P_FIMBRIAE 0.901, SIDEROPHORES
  0.867, CAPSULE_SERUM 0.763) generalize across clades; low-prevalence accessory genes collapse (HEMOLYSIN
  0.796→0.286, CNF1, AFA_DRA). A function carried across many clades has cross-clade signal to learn; a
  clade-restricted accessory gene does not.

**One-line world-model insight:** *the ACQUIRED resistance gene content is transferable across Enterobacterales
lineages (confirmed in E. coli AND Klebsiella), but the plasmid vehicle that carries it, the chromosomal
virulence context around it, and the intrinsic/chromosomal determinants are all lineage-locked.*

## Cross-organism replication (E. coli + Klebsiella)

The determinant↔determinant result is not an E. coli quirk. Re-running it per-organism (each organism
Mash-clustered separately, `--organism`):

| organism | n | clades | verdict | median clade AUC | generalize / concentrated |
|---|---|---|---|---|---|
| E. coli/Shigella | 232 | 80 | GENERALIZES | 0.908 | 39/45, 1 conc |
| Klebsiella | 307 | 118 | GENERALIZES | **0.913** | 46/60, 3 conc |
| Salmonella | 42 | — | under-powered (< MIN_GENOMES 60) | — | skipped (honest) |

Klebsiella lands at essentially the same median (0.913 vs 0.908) — resistance-gene co-occurrence transfers
across lineages in a second Enterobacterales genus. And it sharpens the mechanism: in Klebsiella the split is
cleanly **ACQUIRED-vs-INTRINSIC**. Acquired/mobile determinants generalize (blaTEM-1 0.815, qnrS1 0.835,
aac(3)-IIe 0.835, blaOXA 0.79, dfrA1 0.702 — plasmid/integron-borne), while the genes that COLLAPSE are exactly
Klebsiella's intrinsic/chromosomal ones — **fosA** (intrinsic fosfomycinase, present 287/307, → 0.674),
**blaSHV-11** (intrinsic SHV, 0.579), and chromosomal gyrA/ompK point mutations. This independently recovers
the documented "intrinsic genes are core-genome/lineage-structured, acquired genes are mobile" distinction
([[feedback_intrinsic_genes_break_broad_amr_class_rules]]) — the de-confound sees it with no gene annotation,
purely from the clade-generalization signal. So the headline tightens: **it is specifically the ACQUIRED
resistance content that transfers across lineages; intrinsic/chromosomal determinants are lineage-locked.**

## Literature grounding

- **ST131 reviews** (Johnson / Nicolas-Chanoine; PMC3916147, PMC4135879, PMC8487868) independently support the
  virulence + plasmid halves: resistance↔virulence co-occurrence is driven by **clonal expansion** ("now the
  dominant mechanism"), hlyA is documented **C2-clade-restricted** co-occurring with aac6Ib/blaCTX-M-15
  *within* clade, and IncF plasmids are **stably maintained** in the lineage. Our leave-one-clade-out control
  recovered all three independently — HEMOLYSIN collapses hardest, and the IncFII backbones are the most
  clade-locked.
- **MIC-ML literature** (PMC9280306 N. gonorrhoeae, PMC10995625 K. pneumoniae) flags censoring (left + right)
  and two-fold-dilution uncertainty as *the* core MIC-modeling problem and notes distribution-free coverage
  "would be especially useful" — yet no published paper combines conformal prediction with MIC. The Family-B
  censoring-aware conformal MIC intervals (CRyPTIC RIF/INH, ~0.90 coverage) are methodologically apt and
  ahead of the published intersection.

## Decoder implications

1. **A resistance-gene co-occurrence predictor is defensible cross-lineage.** "If a genome carries sul1, expect
   aadA/dfrA" is a real cassette-structure inference, not a clone artifact — it survives the phylogenetic
   control (0.922). This is a usable, honest decoder capability.
2. **Do NOT infer plasmid backbone or accessory virulence from resistance content alone.** Those axes are
   lineage-locked (0.615 / accessory-half of 0.676); a cross-lineage prediction there is memorizing the clone.
   Any such output must be gated on explicit lineage context.
3. **Prevalence is the honest confidence axis for virulence context.** Core functions (fimbriae/siderophores/
   capsule) are safe to associate; accessory PAI genes (hemolysin/CNF1/afimbrial adhesin) are clade-specific
   and should abstain without lineage.

## Honest scope

- E. coli/Shigella only (the axis with a curated virulence caller). Cohorts are drug-R/S-selected — associational,
  not causal. Mash 0.005 resolves lineages but not sub-clade structure; a within-clade residual can still hide
  below Mash resolution (the surviving AUCs are "not PURELY clonal", not "clonality-free").
- Cross-organism replication done for Klebsiella (n=307); Salmonella (n=42) is below the MIN_GENOMES=60
  both-class/5-fold floor and honestly skipped. A broader multi-genus sweep would need larger per-organism
  cohorts.

## Artifacts

- `scripts/crossaxis_lineage_deconfound.py` (+ `tests/test_crossaxis_lineage_deconfound.py`, 9 tests)
- `wiki/crossaxis_lineage_deconfound_2026-07-12.{json,md}` (virulence, SPLIT)
- `wiki/crossaxis_lineage_deconfound_plasmid_2026-07-12.{json,md}` (plasmid, LINEAGE — prediction falsified)
- `wiki/crossaxis_lineage_deconfound_determinant_2026-07-12.{json,md}` (E. coli determinant, GENERALIZES 0.908)
- `wiki/crossaxis_lineage_deconfound_determinant_klebsiella_2026-07-12.{json,md}` (Klebsiella, GENERALIZES 0.913)
- commits 0eff2f3 (virulence) → 2e07995 (plasmid) → f0c0206 (determinant) → cross-organism; frozen surface byte-unchanged.
