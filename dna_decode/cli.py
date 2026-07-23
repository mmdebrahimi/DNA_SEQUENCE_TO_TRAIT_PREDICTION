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
        "summary": "antibiotic resistance R/S - bacterial (cipro/cef/tet/gent/meropenem; E.coli/Klebsiella/Pseudomonas/S.aureus) + FUNGAL azole/echinocandin (fluconazole/voriconazole/caspofungin/micafungin; C. auris) via --drug",
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
        "validation": "IN-DISTRIBUTION vs measured Deep Mutational Scanning (ProteinGym): ESM2-650M Spearman TEM-1 0.732 / PTEN 0.518, AlphaMissense(human) 0.539, BLOSUM62 weaker (0.35/0.18) but instant+offline; calibrated-magnitude dosage head coverage 10/10 organisms. PROSPECTIVE (leakage-free: MaveDB DMS whose genes are NOT in the ProteinGym benchmark; R2 has no population-structure confound): ESM2-650M median |Spearman| 0.478 over 2383 held-out assays (0.492 on 978 human proteins; pharmacogenes CYP2C19/2C9/G6PD/NUDT15/VKOR 0.547), and beats BLOSUM62 90% paired (p=5e-15) (wiki/mavedb_full_esm2_2026-07-22 + wiki/mavedb_esm_vs_blosum_paired_2026-07-21). LEAKAGE-FREE HYBRID AT SCALE (N=38 held-out, Kaggle T4): ESM2 0.519 / ProSST 0.601 / hybrid 0.586 median |Spearman| -- the hybrid BEATS BOTH components PAIRED (34/38 vs ESM2 +0.060; 26/38 vs ProSST +0.006; read paired, NOT the medians), and structure (ProSST) is the strongest single modality, above AM 0.502 (wiki/mavedb_holdout_hybrid_2026-07-23). CLINICAL (ClinVar path/benign AUROC, actionable human genes): fitness-alignment CEILING (DMS-itself) TP53 0.996 / MSH2 0.955 vs BLOSUM floor 0.707/0.832; the deployable LEARNED decoders fill the gap near the top -- AlphaMissense 0.986/0.936 (no-GPU, best on TP53), shipped ESM2+ProSST hybrid 0.918/0.937 (wins MSH2); winner is gene-dependent. in-distribution-clinical NOT held-out; single-class genes (BRCA1/PTEN) AUROC-inapplicable by design (wiki/clinical_variant_effect_validation_2026-07-22 + wiki/clinical_am_hybrid_auroc_2026-07-22). Regime B molecular fitness RANK, NOT clinical resistance (use `amr` for R/S)",
    },
    "inverse": {
        "summary": "INVERSE design (--protein-seq/--protein-fasta --target-percentile 0.05 [--cds-fasta]): effect -> EDIT. Proposes the edits at a target percentile of predicted molecular damage, using the DMS-validated forward oracle as LABEL-FREE ground truth (no phenotype label consulted). The effect->edit complement to `forward`",
        "validation": "graded NON-circularly against MEASURED wet-lab DMS (calibrate on held-out positions; grade on the proposed variant's measured value, never the model's re-score): beats an exact no-oracle null on 4/4 usable proteins across 4 kingdoms, ~2-5 percentile points at top-5. RANKS, NEVER DOSES -- the magnitude version needs a calibrator fit on the TARGET protein's own DMS (which would make the inverse unnecessary; and calibrators cannot transfer -- the assays share no scale), and its conformal interval is uninformative even where it brackets. The learned oracle beats plain BLOSUM62 on only 3/4, so the blosum62 default is often right, not a fallback; utility does NOT track forward rank (PTEN 0.5185 earns keep, RL40A 0.5190 does not) -> per-protein check required. Regime B molecular fitness only, NOT clinical resistance (use `amr`)",
    },
    "flowering": {
        "summary": "PLANT trait — Arabidopsis thaliana flowering HABIT (--fri/--flc allele calls): summer-annual-early vs winter-annual-late (vernalization-requiring), from the curated FRI/FLC causal loci. The deterministic counterpart to the CLOSED-NEGATIVE flowering EMBEDDING test (which learned lineage, not mechanism)",
        "validation": "deterministic curated-causal-allele rule (late iff functional FRI AND strong FLC; FLC is downstream so a weak/null FLC calls early regardless of FRI). Literature-anchored (Johanson 2000 FRI / Michaels 2003 PNAS weak-FLC / Werner 2005 FRI-independent); reference-integrity biology-checked incl. the Da(1)-12 anchor a naive FRI-only rule mis-calls. PARTIAL: FRI/FLC ~40-70% of long-day variation -> HABIT/direction only, NOT days-to-flower; FRI-route confidence capped by the Lz-0 counterexample. v0 = allele-call input; genome-mode = v0.1",
    },
    "pigment": {
        "summary": "HUMAN visible-trait pigmentation (--genotypes rsID=GT,...): v0 = IrisPlex EYE colour (6 curated SNPs -> P(blue)/P(intermediate)/P(brown) + call) - the deterministic curated-catalog form of 'DNA->appearance'. Benign visible-trait genetics, NOT a forensic tool",
        "validation": "deterministic multinomial-logistic Walsh-2011 IrisPlex coefficients (curated, provenance brianbhsu/eye-color; re-verify vs Walsh Table = v0.1); reference-integrity biology-checked (HERC2 GG->blue, AA->brown). Eye pigmentation AUC ~0.9 (HIrisPlex-S lit). Hair/skin (full 41-SNP) + VCF input + openSNP scoring = v0.1 follow-ons",
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
    if trait == "flowering":
        from dna_decode.organism_rules.flowering_cli import main as flowering_main
        return flowering_main(rest)
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
    if trait not in TRAITS and trait not in ANALYSES:
        ap.error(f"unknown subcommand {trait!r}; traits: {', '.join(TRAITS)}; "
                 f"analyses: {', '.join(ANALYSES)} (or `list`)")
    return _delegate(trait, argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
