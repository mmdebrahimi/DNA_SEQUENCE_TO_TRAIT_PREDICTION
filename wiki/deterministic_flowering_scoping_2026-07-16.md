# Deterministic flowering-time decoder scoping (B-3, 2026-07-16)

The open question from the path-B probe: the Arabidopsis flowering-time EMBEDDING test was a closed NEGATIVE
(2026-06-12, PlantCaduceus de-confounded within-group r² −0.13 — the embedding learned population structure,
not the causal signal; the 3rd de-confounded embedding failure across the kingdom boundary). **Can a
DETERMINISTIC curated-causal-locus decoder work where the embedding failed?** A `/research`-grounded scope, no
build.

## Finding: a curated causal-locus catalog for flowering DOES exist

Arabidopsis flowering-time natural variation is driven predominantly by DISCRETE, catalogued functional
variants at two loci — the textbook winter-annual vs rapid-cycling switch:
- **FRI (FRIGIDA) loss-of-function alleles** — the major determinant. Catalogued functional variants with a
  known effect DIRECTION (LoF → early flowering; functional → late/winter-annual unless vernalized):
  **FRI-Col** (premature stop), **FRI-Ler** (start-codon deletion), plus independent substitution alleles
  (Cvi, Wil-2 — first-intron + in-frame stop). FRI disruption accounts for **~40–70% of long-day
  flowering-time variation** (study-dependent).
- **FLC (FLOWERING LOCUS C) natural nulls** — for accessions with functional FRI: **Van-0** (nonsense),
  **Bur-0** (aberrant splice → null). Adds the FRI-independent early-flowering cases.

This is structurally the SAME shape as the AMR/pigmentation curated catalogs: a small set of curated causal
functional variants → a phenotype-direction call. The embedding failed by learning lineage; a catalog reads
the mechanism directly.

## Verdict: BUILDABLE candidate — the plant analog of the AMR/pigmentation cells

A deterministic **FRI/FLC flowering-habit decoder** (curated LoF-variant catalog → early / late [winter-annual]
call) is buildable from free data and would likely BEAT the embedding on the concentrated large-effect axis —
exactly the AMR pattern (curated catalog catches the concentrated-signal mechanism; the embedding couldn't).
It expands the deterministic decoder to a NEW kingdom (plant) in the winning paradigm.

**Honest caveats (load-bearing — the same disciplines that shaped AMR):**
1. **PARTIAL decoder.** FRI/FLC explains ~40–70%; the rest is polygenic + strongly environment-dependent
   (photoperiod, vernalization). It decodes flowering HABIT / DIRECTION (early vs late, winter-annual vs
   rapid-cycling), NOT the quantitative days-to-flower. Frame it as a habit call with an abstain/uncertain
   tier for the polygenic residue — NOT a continuous predictor (mirrors AMR R/S vs MIC-continuous).
2. **Multi-locus AND + epistasis.** FRI-LoF → early requires FLC to be up-regulatable (counterexample Lz-0:
   FRI deletion but LATE due to high FLC). So the rule is `FRI-LoF AND not(high-FLC-independent)` — a
   multi-locus interaction the frozen count/OR engine can't express (like the TMP-SMX `sul AND dfr` overlay
   or the TB `organism_rules` cell). It needs a NON-frozen `organism_rules`-style module.
3. **Build cost (v0.1+, not this session):** curate the FRI/FLC functional-variant catalog (Col/Ler/Cvi/
   Wil-2/Van-0/Bur-0 with genomic coords from the primary papers) + a variant caller (detect the LoF variants
   from an accession genome/VCF) + validate the habit call against the FREE AraPheno flowering data (1135
   accessions, 1001 Genomes; CSV/PLINK). All free; no money.

## Recommendation

**GO as a v0.1+ candidate** — the deterministic FRI/FLC flowering-habit decoder is a genuine curated-causal-
locus cell where the embedding failed, buildable from free data, expanding the winning paradigm to plants. It
is a real build (catalog curation + caller + AraPheno validation), not a this-session task. Sequence it within
the free curated-catalog frontier (alongside the pigmentation hair/skin extension) BEFORE any acquisition (A).
It also stands as the honest counterpoint to the closed embedding negative: **the flowering signal was
decodable all along — deterministically, not by embedding.**

Sources: [FRI/FLC role in flowering variation, PubMed 15908596](https://pubmed.ncbi.nlm.nih.gov/15908596/) ·
[Molecular analysis of FRIGIDA, PubMed 11030654](https://pubmed.ncbi.nlm.nih.gov/11030654/) ·
[FRIGIDA-independent variation, PMC1451178](https://pmc.ncbi.nlm.nih.gov/articles/PMC1451178/) ·
[FRI-Ler functional analysis, PMC4158083](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4158083/) ·
AraPheno free data: [arapheno.1001genomes.org](https://arapheno.1001genomes.org/study/12/).
