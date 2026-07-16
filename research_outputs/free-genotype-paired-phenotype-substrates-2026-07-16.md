<!-- memo-schema: 0.4 -->
# Free genotype-paired phenotype substrates for a deterministic decoder — supported memo (2026-07-16)

## Research Context

**Problem anchor (verbatim intent):** Path B of the VOI menu — the cheap, reversible probe run *before* the
expensive gated commitment (A = label acquisition). Question: does a FREE, independent, genotype-paired
phenotype substrate exist in an *unexploited* modality that (i) clears the dna_decode label wall by
construction AND (ii) fits the project's *winning* deterministic curated-catalog paradigm? If B finds one →
new capability, no money, A unnecessary. If B confirms the wall → A is provably the only ceiling-raiser.

## Supported findings

1. **Free genotype+phenotype data is ABUNDANT — "no free labels" was never the wall.**
   - **PGS Catalog:** an open database of **5,022 polygenic scores across 656 traits**, each shipped as a
     deterministic apply-from-VCF formula ("weights/effect sizes, effect allele, genome build") + a REST API
     + `pgsc_calc` pipeline. `[medium]` scale / `[high]` format+license.
   - **AraPheno / 1001 Genomes:** 1,135 Arabidopsis accessions with free common-garden phenotypes (flowering
     DTF1) in CSV/PLINK/JSON. `[medium]`
   - **Yeast 1011 / Y1000+:** 1,011 genomes × 223 lab-measured life-cycle/growth traits (under drugs, stress,
     carbon/nitrogen); Y1000+ adds 1,154 strains × 24 conditions. `[medium/low]`

2. **The abundant free space is exactly the polygenic / ancestry-confounded regime — the project's OWN wall,
   baked in.** Polygenic-score accuracy **decays with genetic distance**: Martin et al. report relative
   accuracy **1.6× lower (Hispanic/Latino), 2.5× (East Asian), 4.9× (African)**; a 2025 systematic analysis
   finds predictive power drops to **24% (African) / 37% (East Asian) / 51% (South Asian)** of European.
   `[medium]` This is *lineage/ancestry* confound — the identical "captures population structure, not
   mechanism" failure the project mapped **0-for-5** across the kingdom boundary. So the free-polygenic space
   carries the closed-negative regime by construction.

3. **The free polygenic space is already OWNED by an incumbent the project has no edge to beat.** PGS Catalog
   already provides free, deterministic, published-weight predictors for 656 traits. Building a dna_decode
   cell there means *beating PGS Catalog* — and the project's own disciplines (beat-the-domain-knowledge-
   baseline; validate-wrapper-vs-underlying-tool) say that requires a demonstrated edge, of which there is
   none for polygenic traits.

4. **The ONLY free substrate that fits the project's WINNING paradigm is a curated causal-locus trait.** The
   deterministic decoder wins where a curated catalog of *causal* loci exists (AMR determinants, PGx star
   alleles, viral/fungal drug-resistance mutations). Among unexploited traits, exactly one free such
   substrate surfaced — and it was already found in the D scoping: **visible pigmentation** (HIrisPlex-S: 41
   curated SNPs + free logistic coefficients, eye AUC ~0.9) validated on the free **openSNP** cohort (6,401
   genotypes × 668 phenotypes incl. eye/hair colour). B finds **no additional** free ceiling-raiser.

## Verdict (path B)

**B confirms the wall — it does NOT find a new free substrate that clears it.** Free genotype+phenotype data
is plentiful, but it lives in the polygenic/ancestry-confounded regime that is (a) the project's own
closed-negative (0-for-5), (b) already owned by PGS Catalog with no edge to beat, and (c) not
curated-causal-locus-decodable. The single free substrate that fits the winning paradigm — visible pigmentation
(HIrisPlex-S/openSNP) — was already identified in the D probe; B adds none.

**Therefore A (label acquisition) is provably the only ceiling-raiser** for genuinely new capability — with a
sharpened target: the useful acquisition is NOT "more phenotype data" (abundant, confounded, incumbent-owned)
but **mechanism-resolved / curated-causal-locus labels** for a trait class where the project's deterministic-
catalog edge applies. This is exactly the "if B confirms the wall → draft A on solid ground" branch.

## Decisions for Human Confirmation

| # | Decision | Candidate use | Verification needed |
|---|---|---|---|
| 1 | Build the visible-pigmentation cell (HIrisPlex-S + openSNP)? | The one reachable free curated-catalog cell (from D) | DNA-11 lane coordination + application-framing authority call |
| 2 | Proceed to draft the A acquisition target-list? | B confirmed A is the only new-capability lever; target = mechanism-resolved labels | User go-ahead (draft is reversible; acquisition is money+authority) |
| 3 | Deterministic flowering-time decoder (FLC/FRI loci, AraPheno)? — open, not a found substrate | Possible curated-locus eukaryotic cell where embedding failed | Confirm curated causal-locus catalog exists + is decodable (untested) |

## Honest gaps

See raw memo. Key: PGS live count not fetched (dated-pub value, medium); Arabidopsis deterministic-locus
decoding untested (embedding was a closed negative); yeast drug-resistance-locus subset not enumerated.

Sources: [PGS Catalog](https://www.pgscatalog.org/about/) · [PGS portability review PMC9391275](https://pmc.ncbi.nlm.nih.gov/articles/PMC9391275/) · [ancestry accuracy PMC12622163](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12622163/) · [AraPheno](https://arapheno.1001genomes.org/study/12/) · [yeast 1011 PMC12583546](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12583546/) · [openSNP](https://opensnp.org/) · [HIrisPlex-S FSI:G 2018](https://www.sciencedirect.com/science/article/abs/pii/S1872497318302205).
