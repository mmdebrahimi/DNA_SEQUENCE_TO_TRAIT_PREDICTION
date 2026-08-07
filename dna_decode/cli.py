"""Unified `dna-decode` tool entry — one command, the trait decoders underneath.

Turns the two validated deterministic decoders into a single coherent tool (per the project north star:
"AI DNA decoder tool, not papers"). Thin dispatcher — each subcommand delegates argv verbatim to the
existing decoder `main(argv)`, so the per-decoder CLIs (`dna-amr`, `dna-pathotype`) stay independently
usable and their behavior/output is unchanged.

    dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x
    dna-decode amr --drug ceftriaxone  --genome-fasta X.fna --sample-id X   # needs Docker + AMRFinder DB
    dna-decode pathotype assembly.fna --sample-id MY_STRAIN
    dna-decode list        # what this tool decodes + per-trait validation status
    dna-decode --version

Deterministic mechanism-feature decoders — NOT embeddings (the frozen-genome-embedding thesis was tested
to a decisive FAIL on every reachable de-confounded substrate; see plans/AMR_embedding_niche_decision_2026-06-05.md
+ CHANGELOG 0.3.0). NOT a clinical decision tool.
"""
from __future__ import annotations

import argparse
import sys

# Per-trait capability registry: `dna-decode` SUBCOMMAND -> (delegate dotted-path main, one-line
# capability + validation). This is INTENTIONALLY a different namespace from
# `dna_decode.data.cell_registry.cli_routable_manifest()`, which maps top-level CONSOLE SCRIPTS
# (dna-amr / dna-pgx / dna-hla / dna-clinvar / traits) to their routable cells. The two are
# orthogonal by design — do NOT "unify" TRAITS to be generated from the cell registry; they answer
# different questions (subcommand dispatch table vs console-script->cell manifest).
TRAITS = {
    "amr": {
        "summary": "antibiotic resistance R/S - bacterial (cipro/cef/tet/gent/meropenem; E.coli/Klebsiella/Pseudomonas/S.aureus) + M. tuberculosis (rif/inh) + FUNGAL azole/echinocandin (fluconazole/voriconazole/caspofungin/micafungin; C. auris) + VIRAL target-site (HIV NNRTI/NRTI/PI/INSTI/CAI, SARS-CoV-2 Mpro, influenza NA, HCMV herpesvirus ganciclovir/cidofovir/foscarnet/letermovir via --observed) via --drug",
        "validation": "bacterial: cipro 0.925 (held-out 0.862, cross-source 1.0) | cef 0.933 | gent 0.945 | tet 0.833 | mero 0.867; cross-organism (capstone). fungal C. auris fluconazole G1: sens 1.0 across clades, label-limited spec (wiki/fungal_ep7_g1_closeout_2026-06-08)",
    },
    "pathotype": {
        "summary": "E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/...) compatibility call + abstention",
        "validation": "VirulenceFinder-marker resolver; ExPEC recall 0.917; rest documented scope-limit",
    },
    "plasmid": {
        "summary": "plasmid Inc-replicon typing (IncF/IncH/IncI/IncX/IncN/... via PlasmidFinder allele DB) - composes with amr (is the resistance plasmid-borne?)",
        "validation": "deterministic PlasmidFinder-blastn caller (identity 95 / coverage 60); faithful-to-tool, not an independent baseline; offline-safe degrade",
    },
    "serotype": {
        "summary": "E. coli O:H serotype (wzx/wzy/wzm/wzt O-antigen + fliC H-antigen via SerotypeFinder allele DB)",
        "validation": "deterministic SerotypeFinder-blastn caller (identity 85 / coverage 60); faithful-to-tool; O?/H? when a locus is unresolved; offline-safe degrade",
    },
    "resfinder": {
        "summary": "acquired AMR genes (ResFinder allele DB) - an INDEPENDENT cross-tool check vs amr (AMRFinder DB)",
        "validation": "deterministic ResFinder-blastn caller (identity 90 / coverage 60); caller_is_independent_baseline=True (acquired genes only, no point-mutations/efflux); offline-safe degrade",
    },
    "pointfinder": {
        "summary": "chromosomal AMR point mutations (PointFinder; v0 E. coli FQ QRDR gyrA/parC/gyrB/parE) - INDEPENDENT vs amr's AMRFinder POINT",
        "validation": "deterministic blastn + codon-position lookup vs resistens-overview; caller_is_independent_baseline=True; epistasis recorded not enforced; offline-safe degrade",
    },
    "disinfinder": {
        "summary": "biocide/disinfectant resistance genes (DisinFinder; qac/form... quaternary-ammonium + formaldehyde) - often plasmid-borne (pair with coloc)",
        "validation": "deterministic DisinFinder-blastn caller (identity 90 / coverage 60); faithful-to-tool; offline-safe degrade",
    },
    "mlst": {
        "summary": "MLST sequence type (PubMLST; v0 E. coli Achtman adk/fumC/gyrB/icd/mdh/purA/recA) - exact-allele -> profile -> ST",
        "validation": "deterministic blastn 100/100 exact-allele + PubMLST profile lookup; novel/incomplete -> ST not guessed; `dna-mlst --fetch-db` to install; offline-safe degrade",
    },
    "ktype": {
        "summary": "Klebsiella K-antigen (capsule) type via the wzi allele scheme (BIGSdb Pasteur, Kleborate-bundled) - the serotype sibling",
        "validation": "deterministic wzi-blastn caller (identity 90 / coverage 80); faithful-to-tool; wzi->K ~94% NOT one-to-one; full-locus Kaptive more accurate; offline-safe degrade",
    },
    "salmserovar": {
        "summary": "Salmonella enterica serovar via the Kauffmann-White antigenic formula (O + H1=fliC + H2=fljB; SeqSero2-style antigen DB)",
        "validation": "deterministic antigen-blastn + formula lookup (identity 90 / coverage 80); faithful-to-tool (SeqSero2/Kauffmann-White); serovar only when formula resolves uniquely; free measured label = traditional serotyping (validate vs wet-lab, not the tool); offline-safe degrade",
    },
    "pneumoserotype": {
        "summary": "S. pneumoniae capsular serotype via the cps-locus reference scheme (PneumoCaT/SeroBA-style)",
        "validation": "INDEPENDENT vs phenotypic Quellung (GPS Poland n=230): serogroup 0.939 / exact 0.661 (QUELLUNG-subset n=42: serogroup 0.952). deterministic cps-reference-blastn (id 90/cov 70); serogroup-reliable v0, within-serogroup (6A/6B,19A/19F) needs allele logic (v0.1); offline-safe degrade",
    },
    "pgx": {
        "summary": "HUMAN pharmacogenomics (--gene): CYP2C19 / CYP2C9 diplotype + CPIC metabolizer phenotype, or VKORC1 warfarin sensitivity, from a phased VCF (GRCh38) -- the first human cells",
        "validation": "deterministic VCF->defining-SNP->star-allele->diplotype->CPIC phenotype. GeT-RM consensus concordance on real 1000G (caller independent of the consensus tools): CYP2C19 core 72/72, CYP2C9 core 73/73. CALLING independently validatable; PHENOTYPE faithful-to-CPIC (ref tool PharmCAT). v0 core SNP set; non-core star -> CYP2C19 withholds (sentinel), CYP2C9 mis-calls *1 (sentinel=v0.1). VKORC1 = single-SNP rs9923231 (minus-strand). NOT a clinical tool",
    },
    "forward": {
        "summary": "FORWARD variant-effect (--mutation M69L --protein-seq/--protein-fasta): a protein point mutation -> predicted MOLECULAR-phenotype change (Regime B, enzyme fitness/stability) - the edit->effect complement to `amr`. v0 CLI = BLOSUM62 (deterministic, offline); learned methods (ESM2/AlphaMissense/ESM-IF) via the Python API",
        "validation": "IN-DISTRIBUTION vs measured Deep Mutational Scanning (ProteinGym): ESM2-650M Spearman TEM-1 0.732 / PTEN 0.518, AlphaMissense(human) 0.539, BLOSUM62 weaker (0.35/0.18) but instant+offline; calibrated-magnitude dosage head coverage 10/10 organisms. PROSPECTIVE (leakage-free: MaveDB DMS whose genes are NOT in the ProteinGym benchmark; R2 has no population-structure confound): ESM2-650M median |Spearman| 0.478 over 2383 held-out assays (0.492 on 978 human proteins; pharmacogenes CYP2C19/2C9/G6PD/NUDT15/VKOR 0.547), and beats BLOSUM62 90% paired (p=5e-15) (wiki/mavedb_full_esm2_2026-07-22 + wiki/mavedb_esm_vs_blosum_paired_2026-07-21). LEAKAGE-FREE HYBRID AT SCALE (N=76 held-out, Kaggle T4): ESM2 0.538 / ProSST 0.596 / hybrid 0.602 median |Spearman| -- the hybrid BEATS BOTH components PAIRED (70/76 vs ESM2 +0.063; 52/76 vs ProSST +0.011, sign-test p=0.0009 -- significant, confirmed by doubling N from 38), and structure (ProSST) is the strongest single modality, above AM 0.502 (wiki/mavedb_holdout_hybrid_2026-07-23). 3-WAY (ESM2+GEMME+ProSST, GEMME=evolution via the finalized Docker toolchain, TEM-1 0.719): on the held-out GEMME-covered subset (N=25) the 3-way beats the 2-way 21/25 (sign-p=0.0005) + beats GEMME-alone 22/25 -- adding evolution LIFTS the hybrid; modest N, orthogonal signal (wiki/gemme_threeway_holdout_2026-07-23). CLINICAL (ClinVar path/benign AUROC, actionable human genes): fitness-alignment CEILING (DMS-itself) TP53 0.996 / MSH2 0.955 vs BLOSUM floor 0.707/0.832; the deployable LEARNED decoders fill the gap near the top -- AlphaMissense 0.986/0.936 (no-GPU, best on TP53), shipped ESM2+ProSST hybrid 0.918/0.937 (wins MSH2); winner is gene-dependent. in-distribution-clinical NOT held-out; single-class genes (BRCA1/PTEN) AUROC-inapplicable by design (wiki/clinical_variant_effect_validation_2026-07-22 + wiki/clinical_am_hybrid_auroc_2026-07-22). Regime B molecular fitness RANK, NOT clinical resistance (use `amr` for R/S)",
    },
    "inverse": {
        "summary": "INVERSE design (--protein-seq/--protein-fasta --target-percentile 0.05 [--cds-fasta]): effect -> EDIT. Proposes the edits at a target percentile of predicted molecular damage, using the DMS-validated forward oracle as LABEL-FREE ground truth (no phenotype label consulted). The effect->edit complement to `forward`",
        "validation": "graded NON-circularly against MEASURED wet-lab DMS (calibrate on held-out positions; grade on the proposed variant's measured value, never the model's re-score): beats an exact no-oracle null on 4/4 usable proteins across 4 kingdoms, ~2-5 percentile points at top-5. RANKS, NEVER DOSES -- the magnitude version needs a calibrator fit on the TARGET protein's own DMS (which would make the inverse unnecessary; and calibrators cannot transfer -- the assays share no scale), and its conformal interval is uninformative even where it brackets. The learned oracle beats plain BLOSUM62 on only 3/4, so the blosum62 default is often right, not a fallback; utility does NOT track forward rank (PTEN 0.5185 earns keep, RL40A 0.5190 does not) -> per-protein check required. Regime B molecular fitness only, NOT clinical resistance (use `amr`)",
    },
    "coatcolor": {
        "summary": "DOG coat colour (--loci E=e/e,K=KB/KB,A=at/at,B=B/b,D=D/d): eumelanin/phaeomelanin pigment type + eumelanin colour (black/brown/blue/isabella) + distribution (solid/sable/agouti/tan-points) via the FIVE classic OMIA loci (E/MC1R, K/CBD103, A/ASIP, B/TYRP1, D/MLPH) resolved in fixed epistatic order - the first PHYSICAL/visible-trait animal cell + the deterministic curated-catalog form of 'DNA->appearance'",
        "validation": "deterministic epistatic curated-catalog rule (literature-anchored Little 1957 / Schmutz & Berryere; OMIA-sourced causal loci); reference-integrity biology-checked incl. the E-locus epistasis anchor a naive has-the-allele rule mis-calls (e/e is red/yellow even when K^B + b/b). KNOWLEDGE_BASELINE; per-individual scoring vs the free Darwin's Ark/Dryad cohort (N=1930 owner-reported coat colour + N=3277 canFam4 genotypes, doi:10.5061/dryad.83bk3jb4r) = the v0.1 measured tier (scripts/dog_coat_darwins_ark_validate.py). Calls COLOUR not shade/length/spotting; pattern loci (merle/spotting) ABSTAIN. Companion-animal, NOT human/forensic",
    },
    "morphology": {
        "summary": "DOG body SIZE + EAR type (--dosages IGF1=2,HMGA2=2,STC2=1,GHR=1,EAR=2 | --vcf dog.vcf): relative size rank (toy/small..large/giant) from a 4-locus additive polygenic score (IGF1/HMGA2/STC2/GHR) + ear type (MSRB3 erect/drop) — the pinned + Darwin's-Ark-VALIDATED quantitative/visible-trait sibling of `coatcolor`. Input is per-locus big-allele DOSAGE (0/1/2) OR a canFam4 dog genome VCF (pinned SNPs called by coordinate); coat length/curl + leg length + the 4 rerun morph traits ABSTAIN",
        "validation": "MEASURED relative-signal on the free Darwin's Ark cohort (Dryad doi:10.5061/dryad.83bk3jb4r, canFam4 imputed, N=3277): body-size polygenic score r=+0.619 (R2=0.383) vs owner-reported height Q121; ear MSRB3 lead chr10:8612500 r=+0.543 vs Q125, cleanly resolved from the HMGA2 body-size SNP (the Morrill 2022 confound). Causal SNP coords OMIA/lit canFam3.1 -> canFam4 liftover -> .bim-verified -> functionally validated (unlike the coat indels, the body-size SNPs ARE in-panel). RELATIVE rank NOT absolute inches; ear erect/drop naming is MSRB3-literature-anchored (medium confidence). Companion-animal, NOT human/forensic. See wiki/dog_morphology_darwins_ark_validated_2026-07-30 + dog_body_size_darwins_ark_pinned_2026-07-30",
    },
    "rabbitcolor": {
        "summary": "RABBIT coat colour (--loci A=A/a,B=B/b,C=C/C,D=D/d,E=E/e): the textbook A-E mammalian series — agouti/tan/self (A/ASIP) x black/chocolate (B/TYRP1) x full-colour/chinchilla/Himalayan/albino (C/TYR) x dense/dilute (D/MLPH) x extension/steel/red (E/MC1R). Via the shared mammalian-colour engine",
        "validation": "deterministic curated OMIA epistatic rule (classic rabbit A-E series); reference-integrity biology-checked incl. anchors a naive rule mis-calls (albino c/c masks all; e/e red hides agouti; Ed self-black). KNOWLEDGE_BASELINE: no free per-individual validation substrate. Companion/lab-animal, NOT human/forensic",
    },
    "mousecolor": {
        "summary": "MOUSE coat colour (--loci A=a/a,B=b/b,C=C/C,D=D/d,P=p/p,E=E/E): the FOUNDATIONAL mammalian pigment loci — agouti (A/ASIP) x brown (B/Tyrp1) x albino (C/Tyr) x dilute (D/Myo5a) x pink-eyed dilution (P/Oca2) x extension (E/Mc1r). Via the shared engine",
        "validation": "deterministic curated OMIA epistatic rule (the century-old mouse pigment genetics); reference-integrity biology-checked (albino c/c masks; e/e yellow hides agouti). KNOWLEDGE_BASELINE. Lab-animal, NOT human/forensic",
    },
    "cattlecolor": {
        "summary": "CATTLE coat colour (--loci E=ED/e,PMEL=Dc/n): base black/red/wild via MC1R Extension (ED dominant-black > E+ > e recessive-red) + incompletely-dominant PMEL/SILV dilution (Charolais Dc / Highland Dh -> dun/silver). Via the shared engine",
        "validation": "deterministic curated OMIA epistatic rule (MC1R ED/E+/e; PMEL Dc/Dh dosage dilution); reference-integrity biology-checked (ED dominant black; e/e red; PMEL incomplete-dominant). KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "pigcolor": {
        "summary": "PIG coat colour (--loci KIT=I/i+,E=ED/e): KIT Dominant-White (I, epistatic — masks all colour) + MC1R Extension (ED dominant-black > E+ > e recessive-red). Via the shared engine (OMIA 001199-9823 extension)",
        "validation": "deterministic curated OMIA epistatic rule (KIT dominant white epistatic over MC1R E-series); reference-integrity biology-checked (KIT I masks; ED black; e/e red). KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "sheepcolor": {
        "summary": "SHEEP coat colour (--loci A=AWt/a,E=ED/E+): ASIP Agouti (A^Wt dominant white/tan from a 190kb duplication > a recessive black) + MC1R Extension (ED dominant-black overrides ASIP white). Via the shared engine",
        "validation": "deterministic curated OMIA epistatic rule (OMIA 000201-9940 agouti; ASIP dominant-white duplication vs recessive-black LOF; MC1R ED dominant black); reference-integrity biology-checked incl. ED-overrides-ASIP-white. KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "camelcolor": {
        "summary": "CAMEL (dromedary) coat colour (--loci MC1R=W/E+,A=A/a): MC1R c.901C>T DOMINANT white (dominant-negative, heterozygote sufficient) + ASIP recessive black (exon-2 frameshift); wild = light brown. Via the shared engine (Almathen 2018)",
        "validation": "deterministic curated OMIA rule (MC1R dominant-white / ASIP recessive-black, Almathen 2018 / Alshanbari 2019); reference-integrity biology-checked. KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "minkcolor": {
        "summary": "AMERICAN MINK coat colour (--loci C=C/c,B=B/b,D=D/d): TYR albino/Himalayan + TYRP1 American-Palomino brown + MLPH Silverblue dilute; wild = dark. Via the shared engine (OMIA 000202/000031-452646, Manakhov 2019)",
        "validation": "deterministic curated OMIA rule (TYR albino nonsense; TYRP1 Palomino intron-2 insertion; MLPH Silverblue splice c.901+1G>A); reference-integrity biology-checked. KNOWLEDGE_BASELINE. Fur-farm animal, NOT human/forensic",
    },
    "roedeercolor": {
        "summary": "ROE DEER coat colour (--loci A=A/a): ASIP c.33G>T p.Leu11Phe -- A chestnut (G, phaeomelanin) > a recessive black (T). Via the shared engine (OMIA 000201-9858, Reissmann 2020)",
        "validation": "deterministic curated OMIA rule (ASIP c.33G>T; TT black / GG-GT chestnut, Reissmann 2020); reference-integrity biology-checked. KNOWLEDGE_BASELINE. Wildlife, NOT human/forensic",
    },
    "guineapigcolor": {
        "summary": "GUINEA PIG coat colour (--loci A=A/a,B=B/b,C=C/C,D=D/d,E=E/e): the classic A/B/C/D/E series (ASIP agouti/non-agouti black c.181delTTCA / TYRP1 brown / TYR / MLPH dilute / MC1R e recessive-red). Via the shared engine (OMIA 000201-10141, 001199-10141)",
        "validation": "deterministic curated OMIA epistatic rule; reference-integrity biology-checked. KNOWLEDGE_BASELINE. Lab/companion animal, NOT human/forensic",
    },
    "foxcolor": {
        "summary": "FOX (red/silver) coat colour (--loci E=EA/E+,A=A/a): the non-epistatic silver-fox system -- MC1R EA Alaska-Silver dominant-black (Cys125Arg gain) + ASIP A wild-red / a Standard-Silver recessive-black (both routes -> dark). Via the shared engine (OMIA 000201-9627)",
        "validation": "deterministic curated OMIA rule (Vage 1997 non-epistatic ASIP x MC1R); reference-integrity biology-checked (Alaska + Standard silver both black via different loci; wild red). KNOWLEDGE_BASELINE. NOT human/forensic",
    },
    "donkeycolor": {
        "summary": "DONKEY coat colour (--loci E=E+/e,A=A/a,C=C/c): MC1R e recessive-red (c.629T>C, Abitbol 2014) + ASIP light-points/grey-dun vs solid-black-no-light-points (c.349T>C, Sun 2017) + TYR Asinara-white albinism. Via the shared engine (OMIA 000201/001199-9793)",
        "validation": "deterministic curated OMIA rule; reference-integrity biology-checked (grey-dun light-points / solid-black / red / Asinara albino). ASIP heterozygote variability -> a residual gene is unidentified (documented). KNOWLEDGE_BASELINE. NOT human/forensic",
    },
    "buffalocolor": {
        "summary": "WATER BUFFALO coat colour (--loci A=AW/a): ASIP A^W DOMINANT white via a 2809-bp LINE-1 insertion (10x ASIP overexpression, Liang 2020) > a black. MC1R is monomorphic in buffalo so ASIP is the sole driver. Via the shared engine (OMIA 000213-89462)",
        "validation": "deterministic curated OMIA rule (LINE-1 ASIP dominant-white; convergent with cattle); reference-integrity biology-checked. KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "goatcolor": {
        "summary": "GOAT coat colour (--loci A=AWt/a,B=B/b): ASIP Agouti (A^Wt dominant white/tan, the CNV-driven many-pattern hub > a recessive nonagouti black) + TYRP1 brown (Copperneck). Via the shared mammalian engine (OMIA 000201-9925)",
        "validation": "deterministic curated OMIA epistatic rule (ASIP dominant-white/tan vs recessive-black; TYRP1 brown); reference-integrity biology-checked. Goat MC1R association is incomplete in the literature so ASIP is the modeled driver. KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "alpacacolor": {
        "summary": "ALPACA/llama fleece colour (--loci E=E/e,A=A/a): MC1R Extension (E coloured; e/e = recessive WHITE regardless of ASIP -- the camelid twist) + ASIP Agouti (A functional -> fawn/agouti > a loss-of-function -> black). Via the shared engine",
        "validation": "deterministic curated OMIA epistatic rule (MC1R E/e camelid recessive-white; ASIP fawn-vs-black); reference-integrity biology-checked incl. the ee-recessive-white anchor (white regardless of ASIP). KNOWLEDGE_BASELINE. Livestock, NOT human/forensic",
    },
    "pigeoncolor": {
        "summary": "PIGEON plumage colour (--loci B=BA/B+,E=E+/e,D=D/d,C=C/+ [--sex female]): base ash-red/blue/brown (B/TYRP1, Z-linked) + recessive-red (E/SOX10, epistatic) + dilute dun/khaki/ash-yellow (D/SLC45A2, Z-linked) + wing pattern T-check/checker/bar/barless (C/NDP). One of the best-characterised colour systems in any organism (Shapiro lab). A 2nd BIRD cell; B/D Z-linked -> FEMALE (ZW) hemizygous",
        "validation": "deterministic curated epistatic rule; molecularly-confirmed causal genes (TYRP1 B-locus Domyan 2014; SOX10 recessive-red; SLC45A2 dilute; NDP wing-pattern Vickrey 2018 eLife). reference-integrity biology-checked incl. anchors a naive rule mis-calls: (1) SOX10 e/e is RED regardless of the TYRP1 base; (2) Z-linked reversed hemizygosity (FEMALE is ZW). KNOWLEDGE_BASELINE: no free per-individual validation substrate. Calls base+dilute+wing-pattern not modifiers/shade. Hobby/livestock, NOT human/forensic",
    },
    "plumage": {
        "summary": "CHICKEN plumage colour (--loci E=E/E,B=B/b+,S=S/s+,I=i+/i+,BL=bl+/bl+ [--sex male]): eumelanin canvas (extended-black/birchen/wheaten/partridge via E/MC1R) + Z-LINKED barring (B/CDKN2A) + Z-linked silver/gold (S/SLC45A2) + dominant white (I/PMEL17) + blue/splash (Bl) + lavender (MLPH) + recessive white (TYR). The B/S loci are Z-linked so a FEMALE (ZW) is HEMIZYGOUS (reversed from mammals). A 4th-organism (bird) visible-trait cell",
        "validation": "deterministic epistatic curated-catalog rule; OMIA-sourced causal variants (MC1R E-locus series OMIA 000374-9031; CDKN2A sex-linked barring OMIA 000102-9031 Hellstrom/Schwochow; SLC45A2 silver OMIA 000370-9031 Gunnarsson; PMEL17 dominant white OMIA 000373-9031 Kerje; Bl blue; MLPH lavender). reference-integrity biology-checked incl. the epistasis anchors a naive rule mis-calls: (1) EXTENSION is the canvas (barring/blue barely show on a wheaten bird); (2) Z-LINKED barring/silver with REVERSED hemizygosity (the FEMALE is ZW-hemizygous, mirror of cat's X-linked orange); (3) dominant/recessive white mask eumelanin. KNOWLEDGE_BASELINE: no free per-individual validation substrate. Calls canvas+major modifiers not fine pattern/lacing/pencilling. Livestock, NOT human/forensic",
    },
    "catcolor": {
        "summary": "CAT coat colour (--loci O=O/o,A=a/a,B=B/B,D=D/D,C=C/C,W=w/w [--sex female]): base (black/chocolate/cinnamon + dilute blue/lilac/fawn) + X-LINKED orange (red/cream) + TORTOISESHELL/CALICO mosaic + tabby + colorpoint (Siamese/Burmese) + white spotting + dominant white, via the OMIA loci (W/KIT, O/ARHGAP36 X-linked, A/ASIP, B/TYRP1, D/MLPH, C/TYR). The O locus is X-linked (1 allele=male, 2=female); a female O/o is a tortoiseshell mosaic. A 3rd-organism visible-trait cell",
        "validation": "deterministic epistatic curated-catalog rule; OMIA-sourced causal variants incl. the 2025-identified X-linked ORANGE gene (ARHGAP36 5.1-kb deletion, Toh/Kaelin Current Biology 2025), plus KIT dominant-white/spotting (FERV1), ASIP agouti, TYRP1 brown, MLPH dilute, TYR albino-series (Siamese cs). reference-integrity biology-checked incl. THREE epistasis anchors a naive rule mis-calls: (1) W dominant-white masks ALL colour; (2) a female O/o is a TORTOISESHELL mosaic (X-inactivation), not uniform — +white spotting=CALICO; (3) orange is EPISTATIC over brown (a b/b orange cat is red, not chocolate). KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate. Calls colour+major pattern not tabby sub-pattern/shade/spotting-extent. Companion-animal, NOT human/forensic",
    },
    "horsecolor": {
        "summary": "HORSE coat colour (--loci E=E/e,A=A/a,CR=Cr/N,D=nd1/nd1,G=n/n): base (chestnut/bay/black) + cream dilution (palomino/buckskin/cremello/perlino) + dun (red dun/grullo) + grey (progressive), via the five OMIA loci (E/MC1R, A/ASIP, CR/SLC45A2, D/TBX3, G/STX17) resolved in fixed epistatic order — the best-characterised animal coat system, a 2nd-organism visible-trait cell alongside dog `coatcolor`",
        "validation": "deterministic epistatic curated-catalog rule (OMIA-sourced causal variants: MC1R S83F chestnut / ASIP a black / SLC45A2 c.457G>A cream / TBX3 dun / STX17 dup grey); reference-integrity biology-checked incl. the TWO epistasis anchors a naive rule mis-calls (e/e is chestnut even when A/A bay; a G/n horse greys out regardless of base). KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate (unlike the dog cell's Darwin's Ark). Calls COLOUR not sooty/flaxen shade or spotting extent; champagne/silver/pearl/roan/tobiano/appaloosa ABSTAIN. v0 = allele-call input; genome-mode = v0.1. Livestock/companion, NOT human/forensic",
    },
    "flowering": {
        "summary": "PLANT trait — Arabidopsis thaliana flowering HABIT (--fri/--flc allele calls): summer-annual-early vs winter-annual-late (vernalization-requiring), from the curated FRI/FLC causal loci. The deterministic counterpart to the CLOSED-NEGATIVE flowering EMBEDDING test (which learned lineage, not mechanism)",
        "validation": "deterministic curated-causal-allele rule (late iff functional FRI AND strong FLC; FLC is downstream so a weak/null FLC calls early regardless of FRI). Literature-anchored (Johanson 2000 FRI / Michaels 2003 PNAS weak-FLC / Werner 2005 FRI-independent); reference-integrity biology-checked incl. the Da(1)-12 anchor a naive FRI-only rule mis-calls. PARTIAL: FRI/FLC ~40-70% of long-day variation -> HABIT/direction only, NOT days-to-flower; FRI-route confidence capped by the Lz-0 counterexample. v0 = allele-call input; genome-mode = v0.1",
    },
    "pigment": {
        "summary": "HUMAN visible-trait pigmentation (--trait eye/hair/skin --genotypes rsID=GT,...): eye colour (IrisPlex 6-SNP -> blue/intermediate/brown), hair colour (blond/brown/red/black), skin colour (very-pale..dark-black) -- the deterministic multinomial-logistic form of 'DNA->appearance'. Benign visible-trait genetics, NOT a forensic tool",
        "validation": "EYE = Walsh-2011 IrisPlex coefficients (irisplex.py), reference-integrity biology-checked (HERC2 GG->blue/AA->brown) + POPULATION-VALIDATED on real 1000G (EUR blue 0.468; AFR/EAS/SAS brown ~1.0). HAIR+SKIN = HIrisPlex-S deployed models RECOVERED from the erasmusmc webtool (papers publish betas but NOT the intercepts -> webtool-only) via a designed-genotype-basis query + LS-fit, VALIDATED on 20 held-out genotypes: max |dP| eye 6e-15 / hair 6e-16 / skin 9e-3 (reproduces the deployed webtool). Population geography also confirmed on 1000G. Population-level, NOT per-individual (openSNP deleted 2025-04-30)",
    },
    "phage": {
        "summary": "BACTERIOPHAGE genome/lineage/RBP -> host-RECEPTOR class: --lineage <genus> (wheel-only catalogue) | --genome-fasta X.fna (genome-homology transfer, needs blastn) | --rbp-fasta X.faa (tail-fiber RBP k-mer transfer, wheel-only — covers the RBP-variable mixed clades Tsx/OmpC/FhuA/OmpA that --lineage abstains on). The first non-AMR, non-host-organism cell (a virus-of-bacteria host-tropism axis)",
        "validation": "INDEPENDENT_MEASURED (covered classes): on the LBNL/Arkin-Mutalik Phage Datasheets (github.com/mjohnson11/PhageDataSheets; measured receptors on K-12 BW25113, non-Bas isolates disjoint from the BASEL reference) the BASEL-2021 catalogue = 25/29=0.862 on covered classes (BtuB 22/26, LPS_core 3/3); in-dist LOO on 29 clade-conserved BASEL phages = 27/27. RBP-level caller (--rbp-fasta) leave-one-out = 156/160=0.975 overall, 0.961 on the RBP-variable classes genome-homology got 0/N (Tsx 35/35, FhuA 22/22, OmpC 18/18). Receptor-CLASS only. See wiki/phage_{independent_result,rbp_caller_result}_2026-07-24",
    },
    "kleb": {
        "summary": "KLEBSIELLA phage depolymerase -> host CAPSULE (KL-type), ranked top-K (--depolymerase-fasta X.faa): the CROSS-ORGANISM cell (E. coli phage-receptor paradigm transferred to a different host + phenotype). FETCH-ONLY — bundles NO data; build the local reference via scripts/fetch_dpotropisearch.py (DpoTropiSearch/Zenodo; CC-BY record / repo non-commercial license — verify your use). Offline/no-reference -> actionable INDETERMINATE",
        "validation": "KNOWLEDGE_BASELINE / in-distribution (DpoTropiSearch prophage-LCA labels; Concha-Eloko/Nat Commun 2025, Zenodo 10.5281/zenodo.14065540 — NOT independent wet-lab). Clonality-corrected leave-one-out (greedy-rep @0.90): top-1 ~0.45 / top-5 ~0.60 over 147-165 KL-types, lift +0.49 over a 0.10 prior null -> the deterministic sequence-homology->phenotype paradigm GENERALIZES cross-organism on modular depolymerase domains. See wiki/klebsiella_{crossorganism_result,topk_ksweep}_2026-07-25",
    },
    "essentiality": {
        "summary": "single-gene KO -> ESSENTIAL / non-essential (--gene X --product '...' | --feature-table X.txt.gz): the deterministic conserved-core decoder. Predicts essentiality from gene FUNCTION (translation / replication / transcription / cell-envelope / division catalogue), label-independent + offline. High-precision, conservative-recall (universal core; the learned E3 complement lifts the tail).",
        "validation": "KNOWLEDGE_BASELINE / validated vs gold-standard. E. coli AUROC 0.695 genome-wide (Goodall 2018 mBio TraDIS Table S1); composition matches the known essentialome (208/4318, translation/envelope/replication-dominated). Cross-organism transfer to human (BAGEL CEG2/NEG) AUROC 0.580. The learned E3 complement (aa-composition+length+core, 5-fold CV) lifts it: E. coli 0.795 / human 0.911. NOT clinical. See wiki/essentiality_{decoder_v0,ecoli_v0_1_auroc,e4_transfer,e3_learned,e3_human}_2026-07-28",
    },
    "metabolic": {
        "summary": "E. coli carbon-source utilization (--source lactose/citrate/... --genes lacZ,lacY | --feature-table X.txt.gz): the deterministic UPTAKE-GATED catabolism decoder. utilizes iff (catabolic enzymes present) AND (a transporter present) AND (transporter expressed under the O2 condition). The uptake-gate is what a naive AMR-style has-the-genes rule misses.",
        "validation": "KNOWLEDGE_BASELINE / validated vs measured E. coli K-12 MG1655 phenotypes (EcoCyc/Neidhardt). Anchors: lac+ ara+ mal+ xyl+ rha+ glc+, and the CITRATE anchor (Blount 2012 Nature LTEE) -- Cit- aerobic / Cit+ anaerobic, the case a naive has-the-genes rule mis-calls positive. Reads gene presence not sequence integrity; calls can/cannot DIRECTION not growth rate. NOT clinical. See wiki/metabolic_carbon_decoder_v0_2026-07-28",
    },
    "motility": {
        "summary": "flagellar SWIMMING motility from gene presence (--genes flhD,flhC,fliC,motA,... | --feature-table X.txt.gz): the first NON-metabolic trait catalog. MOTILE iff all 5 flagellar modules present (master flhDC -> sigma-28 fliA -> flagellin fliC/fljB -> motor motAB -> basal-body/export fliF/fliG/flhA/fliI); chemotaxis (cheA/W/Y/Z) reported SEPARATELY (a che-mutant still swims). The determinant->phenotype paradigm (like amr/metabolic) applied to a physical behaviour.",
        "validation": "KNOWLEDGE_BASELINE / curated flagellar catalog vs literature anchors: E. coli K-12 MG1655 + Salmonella motile (all modules) vs Shigella flexneri non-motile (flagellar pseudogenes) + flhDC/fliC/motAB KO non-motile. Presence-based DIRECTION (swim/no-swim), NOT speed; cannot see a present-but-inactivated gene (the K-12 flhD IS-insertion) -> sequence-mode is v0.1. NOT clinical. See wiki/motility_catalog_v0_2026-08-03",
    },
    "fba": {
        "summary": "gene edit -> QUANTITATIVE cell-level trait via genome-scale flux-balance analysis (iML1515 E. coli default; --organism saureus|salmonella|pputida|yeast generalizes -- each alias loads a model of THAT organism; P. aeruginosa has no BiGG reconstruction and is refused, never substituted). KO mode: --gene b0720|gltA | --knockout b0720,b0721 | --wildtype -> growth rate (/h) + essential/non-essential over ANY model gene. SYNTHETIC-LETHALITY: --knockout A,B --synthetic-lethality -> is the PAIR lethal though neither single is? DESIGN mode (the INVERSE direction -- product -> edits): --design-target succ -> searches knockouts that make producing the target NECESSARY for growth (growth-coupled strain design; two-sided LP at a fixed growth floor). POINT-MUTATION mode (composes `forward`): --gene gltA --mutation D362A --protein-seq S -> forward scores the missense (LOF?) -> if damaging, model as KO -> cell trait; uncertain -> reported both-ways (never forced). The first GENERAL edit->quantitative-trait rung; computes from stoichiometry, sidesteps population-structure confounding.",
        "validation": "KNOWLEDGE_BASELINE (in-distribution) / genome-wide single-gene-deletion essentiality (glucose M9 aerobic) vs the free Keio-collection mutant-fitness gold standard (Bernstein 2023 method, fitness<-2): accuracy 0.954 / MCC 0.652. Point-mutation mode inherits forward's DMS validation (missense->LOF) + this Keio validation (LOF->trait); the LOF binarization is forward's own method-aware threshold (heuristic, not a calibrated LOF probability). METABOLIC traits only -- NOT virulence/regulation. cobrapy required (`uv pip install cobra`); iML1515 auto-fetched from BiGG + cached. See wiki/fba_keio_validation_2026-08-03 + wiki/fba_variant_compose_2026-08-03",
    },
}


def _version() -> str:
    try:
        from importlib.metadata import version
        return version("dna_decode")
    except Exception:
        return "unknown"


def _delegate(trait: str, rest: list[str]) -> int:
    if trait == "amr":
        from dna_decode.amr.cli import main as amr_main
        return amr_main(rest)
    if trait == "pathotype":
        from dna_decode.pathotype.cli import main as patho_main
        return patho_main(rest)
    if trait == "plasmid":
        from dna_decode.plasmid.cli import main as plasmid_main
        return plasmid_main(rest)
    if trait == "serotype":
        from dna_decode.serotype.cli import main as serotype_main
        return serotype_main(rest)
    if trait == "resfinder":
        from dna_decode.resfinder.cli import main as resfinder_main
        return resfinder_main(rest)
    if trait == "pointfinder":
        from dna_decode.pointfinder.cli import main as pointfinder_main
        return pointfinder_main(rest)
    if trait == "disinfinder":
        from dna_decode.disinfinder.cli import main as disinfinder_main
        return disinfinder_main(rest)
    if trait == "mlst":
        from dna_decode.mlst.cli import main as mlst_main
        return mlst_main(rest)
    if trait == "ktype":
        from dna_decode.ktype.cli import main as ktype_main
        return ktype_main(rest)
    if trait == "salmserovar":
        from dna_decode.salmserovar.cli import main as salmserovar_main
        return salmserovar_main(rest)
    if trait == "pneumoserotype":
        from dna_decode.pneumoserotype.cli import main as pneumoserotype_main
        return pneumoserotype_main(rest)
    if trait == "pgx":
        from dna_decode.pgx.cli import main as pgx_main
        return pgx_main(rest)
    if trait == "forward":
        from dna_decode.forward.cli import main as forward_main
        return forward_main(rest)
    if trait == "inverse":
        from dna_decode.forward.inverse_cli import main as inverse_main
        return inverse_main(rest)
    if trait == "pigment":
        from dna_decode.pigment.cli import main as pigment_main
        return pigment_main(rest)
    if trait == "coatcolor":
        from dna_decode.pigment.coat_cli import main as coat_main
        return coat_main(rest)
    if trait == "morphology":
        from dna_decode.pigment.morphology_cli import main as morphology_main
        return morphology_main(rest)
    if trait == "horsecolor":
        from dna_decode.pigment.horse_coat_cli import main as horse_main
        return horse_main(rest)
    if trait == "catcolor":
        from dna_decode.pigment.cat_coat_cli import main as cat_main
        return cat_main(rest)
    if trait == "plumage":
        from dna_decode.pigment.chicken_plumage_cli import main as plumage_main
        return plumage_main(rest)
    if trait == "pigeoncolor":
        from dna_decode.pigment.pigeon_plumage_cli import main as pigeon_main
        return pigeon_main(rest)
    if trait in ("rabbitcolor", "mousecolor", "cattlecolor", "pigcolor", "sheepcolor", "goatcolor", "alpacacolor",
                 "guineapigcolor", "foxcolor", "donkeycolor", "buffalocolor", "camelcolor", "minkcolor",
                 "roedeercolor"):
        from dna_decode.pigment import mammal_color_cli as mcc
        return {"rabbitcolor": mcc.rabbit_main, "mousecolor": mcc.mouse_main, "cattlecolor": mcc.cattle_main,
                "pigcolor": mcc.pig_main, "sheepcolor": mcc.sheep_main, "goatcolor": mcc.goat_main,
                "alpacacolor": mcc.alpaca_main, "guineapigcolor": mcc.guineapig_main, "foxcolor": mcc.fox_main,
                "donkeycolor": mcc.donkey_main, "buffalocolor": mcc.buffalo_main, "camelcolor": mcc.camel_main,
                "minkcolor": mcc.mink_main, "roedeercolor": mcc.roedeer_main}[trait](rest)
    if trait == "flowering":
        from dna_decode.organism_rules.flowering_cli import main as flowering_main
        return flowering_main(rest)
    if trait == "phage":
        from dna_decode.phage.cli import main as phage_main
        return phage_main(rest)
    if trait == "kleb":
        from dna_decode.kleb.cli import main as kleb_main
        return kleb_main(rest)
    if trait == "essentiality":
        from dna_decode.essentiality.cli import main as essentiality_main
        return essentiality_main(rest)
    if trait == "metabolic":
        from dna_decode.metabolic.cli import main as metabolic_main
        return metabolic_main(rest)
    if trait == "fba":
        from dna_decode.fba.cli import main as fba_main
        return fba_main(rest)
    if trait == "motility":
        from dna_decode.motility.cli import main as motility_main
        return motility_main(rest)
    if trait == "concordance":
        from dna_decode.concordance.cli import main as concordance_main
        return concordance_main(rest)
    if trait == "profile":
        from dna_decode.profile.cli import main as profile_main
        return profile_main(rest)
    if trait == "coloc":
        from dna_decode.colocalization.cli import main as coloc_main
        return coloc_main(rest)
    raise ValueError(f"unknown trait: {trait}")


# Cross-decoder ANALYSES (compose the decoders; NOT new traits/DBs — kept out of TRAITS so the
# decoder registry contract stays the 5-decoder set).
ANALYSES = {
    "decode": "ROUTER: point at any file (FASTA/VCF) -> which decoders apply + the exact command for each",
    "concordance": "AMR cross-tool concordance (AMRFinder vs ResFinder acquired-gene calls)",
    "profile": "unified genome profile - run all assembly-FASTA decoders in one report",
    "coloc": "resistance-gene x plasmid co-localization (is this acquired AMR gene plasmid-borne?)",
}


def _print_list() -> int:
    print(f"dna-decode {_version()} - deterministic genotype->phenotype decoders\n")
    for name, meta in TRAITS.items():
        print(f"  {name:11} {meta['summary']}")
        print(f"  {'':11} validation: {meta['validation']}")
    print("\nanalyses (compose the decoders):")
    for name, summary in ANALYSES.items():
        print(f"  {name:11} {summary}")
    print("\nrun `dna-decode <trait|analysis> --help` for options.")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode",
        description="Unified DNA trait decoder — deterministic, interpretable, mechanism-feature based.",
        epilog=(f"{len(TRAITS)} traits: " + ", ".join(TRAITS) + f".  {len(ANALYSES)} analyses: "
                + ", ".join(ANALYSES) + ".  Run `dna-decode list` for every command + its validation "
                "status. Zero-setup decodes (no Docker/BLAST/downloads): `forward`, `inverse`, "
                "`flowering`, and `amr --drug <hiv/fungal drug> --observed ...`."),
    )
    ap.add_argument("--version", action="version", version=f"dna-decode {_version()}")
    # metavar was hardcoded "{amr,pathotype,list}" -- a lie: it hid 16 of the 19 commands from the usage
    # line, the first thing `dna-decode --help` shows. Honest placeholder; the full set is in the body + epilog.
    sub = ap.add_subparsers(dest="trait", metavar="<command>")
    # Register thin pass-through subparsers; real arg parsing happens in each decoder's main().
    for name, meta in TRAITS.items():
        sub.add_parser(name, add_help=False, help=meta["summary"])
    for name, summary in ANALYSES.items():
        sub.add_parser(name, add_help=False, help=summary)
    sub.add_parser("list", help="show what this tool decodes + per-trait validation status")

    # Split argv at the subcommand so the rest passes through verbatim (incl. --help) to the decoder.
    if not argv:
        ap.print_help()
        return 0
    trait = argv[0]
    if trait in ("-h", "--help"):
        ap.print_help()
        return 0
    if trait == "--version":
        print(f"dna-decode {_version()}")
        return 0
    if trait == "list":
        return _print_list()
    if trait == "decode":
        # input-aware router: `dna-decode decode <file>` -> which decoders apply + the exact commands;
        # `dna-decode decode <file> --run` -> actually RUN the auto-runnable ones + report the rest.
        rest = argv[1:]
        if not rest or rest[0] in ("-h", "--help"):
            print("usage: dna-decode decode <input.fasta|input.vcf> [--run]\n"
                  "  Detects the input kind (nucleotide/protein FASTA or VCF) and lists every applicable\n"
                  "  decoder with its claim, honest tier, and the exact command to run.\n"
                  "  --run: actually run the auto-runnable decoders (genome -> profile; protein -> inverse)\n"
                  "         and report the ones that need a specific parameter (--mutation / --gene).")
            return 0
        do_run = "--run" in rest
        files = [a for a in rest if not a.startswith("-")]
        if not files:
            print("error: decode needs an input file", file=sys.stderr)
            return 2
        if do_run:
            from dna_decode.decode_router import run_decode_plan
            return run_decode_plan(files[0])
        from dna_decode.decode_router import render_decode_plan
        print(render_decode_plan(files[0]))
        return 0
    if trait not in TRAITS and trait not in ANALYSES:
        ap.error(f"unknown subcommand {trait!r}; traits: {', '.join(TRAITS)}; "
                 f"analyses: {', '.join(ANALYSES)} (or `list`)")
    return _delegate(trait, argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
