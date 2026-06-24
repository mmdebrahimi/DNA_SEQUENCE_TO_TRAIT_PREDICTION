"""Deterministic AMR mechanism decoder — in-package CLI (console entry `dna-amr`).

Genome FASTA (or a cached AMRFinder run) -> R/S call per drug + the curated resistance determinants that
drove it + provenance. Mechanism-feature decoding, NOT embeddings (per
`plans/AMR_embedding_niche_decision_2026-06-05.md`). Sibling of `dna_decode.pathotype.cli` (dna-pathotype).

In-package so it ships in the wheel. Cached-run mode is pure (reads main.tsv via amr_rules — no Docker).
Genome mode lazily imports the AMRFinder Docker runner from `scripts/` (repo-only; needs Docker + a
Docker-readable DB) and errors cleanly if unavailable — so the console entry installs + imports without
the scripts/ dir.

    dna-amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x
    dna-amr --drug ciprofloxacin --genome-fasta X.fna --sample-id X      # needs Docker + data/amrfinder_db

NOT a clinical decision tool. cipro N=147 op-chars (threshold=2): acc 0.939 / sens 0.931 / spec 0.947.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.data.antimalarial_amr import (
    call_from_observed_substitutions as antimalarial_call_from_observed,
    gene_for_drug,
    supported_antimalarial_drugs,
)
from dna_decode.data.antiviral_amr import (
    call_from_observed_substitutions as antiviral_call_from_observed,
    supported_antiviral_drugs,
)
from dna_decode.data.fungal_amr import (
    call_from_observed_substitutions,
    supported_fungal_drugs,
)
from dna_decode.data.sarscov2_amr import (
    all_supported_sarscov2_drugs,
    call_sarscov2_observed,
    gene_for_sarscov2_drug,
)
from dna_decode.data.hiv_amr import (
    all_supported_hiv_drugs,
    call_hiv_observed,
    gene_for_hiv_drug,
)
from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.data.trust_surface import one_line, trust_block
from dna_decode.eval.amr_rules import AMRFINDER_IMAGE_PINNED, call_resistance

# Fungal target-site path (BLAST ERG11/FKS1 -> catalog), NOT AMRFinder (no AMRFinder-for-fungi). Routed by
# drug: fluconazole/voriconazole/caspofungin/micafungin -> fungal engine. G1-validated on C. auris
# (2026-06-08, wiki/fungal_ep7_g1_closeout): method transfers, sens 1.0 across clades; label-limited spec.
_DEFAULT_ERG11_REF = Path(__file__).resolve().parent.parent.parent / "data" / "fungal_ref" / "Cauris_ERG11_cds.fna"
# Antimalarial target-site path (BLAST Pfkelch13 -> WHO-validated marker catalog), the 3rd kingdom
# (protozoan). Routed by drug: artemisinin/artesunate/dihydroartemisinin -> K13 engine.
_DEFAULT_K13_REF = Path(__file__).resolve().parent.parent.parent / "data" / "antimalarial_ref" / "Pf3D7_K13_cds.fna"
_DEFAULT_PFCRT_REF = Path(__file__).resolve().parent.parent.parent / "data" / "antimalarial_ref" / "Pf3D7_pfcrt_cds.fna"
# Antiviral target-site path (BLAST influenza NA -> CDC/WHO-recognized NAI marker catalog), the 4th kingdom
# (viral). Routed by drug: oseltamivir/peramivir/zanamivir -> NA engine. Reference is N1 (WT His275).
_DEFAULT_NA_REF = Path(__file__).resolve().parent.parent.parent / "data" / "antiviral_ref" / "N1_NA_NC026434_cds.fna"

# HIV-1 target-site path (BLAST the class's gene CDS -> Stanford/HIVDB-sourced major-DRM catalog). Routed
# by drug -> gene (RT NNRTI/NRTI, PR PI, IN INSTI, CA CAI). References are HXB2 (K03455.1) in-frame CDS,
# consensus-B WT at every catalogued DRM position.
_HIV_REF_DIR = Path(__file__).resolve().parent.parent.parent / "data" / "hiv_ref"
# SARS-CoV-2 Mpro (3CLpro) target-site path — the next free-independent-label viral cell (validated vs the
# Stanford CoV-RDB measured fold-change). Committed Wuhan-Hu-1 NC_045512.2 nsp5 reference, WT at every
# catalogued position (catalytic dyad H41/C145 + E166 verified).
_DEFAULT_SARSCOV2_MPRO_REF = (Path(__file__).resolve().parent.parent.parent / "data" / "sarscov2_ref"
                              / "SARSCoV2_Mpro_NC045512_cds.fna")
_DEFAULT_HIV_RT_REF = _HIV_REF_DIR / "HIV1_RT_HXB2_cds.fna"
_DEFAULT_HIV_PR_REF = _HIV_REF_DIR / "HIV1_PR_HXB2_cds.fna"
_DEFAULT_HIV_IN_REF = _HIV_REF_DIR / "HIV1_IN_HXB2_cds.fna"
_DEFAULT_HIV_CA_REF = _HIV_REF_DIR / "HIV1_CA_HXB2_cds.fna"


def _parse_observed(observed: str) -> dict[str, set[str]]:
    """'ERG11:Y132F,ERG11:K143R,FKS1:S639F' -> {'ERG11': {'Y132F','K143R'}, 'FKS1': {'S639F'}}."""
    out: dict[str, set[str]] = {}
    for tok in (t.strip() for t in observed.split(",") if t.strip()):
        gene, _, sub = tok.partition(":")
        if not sub:
            raise ValueError(f"bad --observed token {tok!r}; expected GENE:SUBSTITUTION (e.g. ERG11:Y132F)")
        out.setdefault(gene.strip(), set()).add(sub.strip())
    return out


def _target_site_record(call, sample_id: str, drug: str, organism: str, provenance: dict, *,
                        caller_name: str, source: str) -> dict:
    """Map a target-site Call (fungal OR antimalarial — same shape) onto the uniform
    amr-mechanism-call-v1 record (same shape as the bacterial path)."""
    dets = [{"symbol": d.split(":", 1)[0], "subclass": d.split(":", 1)[1] if ":" in d else "",
             "class": "TARGET_SITE_MUTATION", "pct_identity": None} for d in call.determinants]
    return {
        "sample_id": sample_id, "drug": drug,
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-mechanism-call-v1",
        "prediction": call.prediction,
        "confidence": "high" if call.prediction == "R" else ("n/a" if call.prediction == "INDETERMINATE" else "screen"),
        "n_determinants": len(call.determinants), "determinants": dets,
        "resistance_threshold": 1,
        "undetectable_mechanisms": list(call.undetectable_mechanisms),
        "caller": {"name": caller_name, "rule": call.rule, "source": source,
                   "caller_is_independent_baseline": False},
        "caveat": call.caveat,
        "validation": trust_block(drug, organism),
        "provenance": {**provenance, "organism": organism},
    }


def _fungal_main(args) -> int:
    """Fungal target-site decoder branch (routed when --drug is a fungal drug)."""
    # the --organism default 'Escherichia' is the bacterial AMRFinder default; on the fungal path, relabel
    # to the validated fungal organism unless the user explicitly set one.
    if args.organism == "Escherichia":
        args.organism = "Candida_auris"
    if args.observed is not None:
        obs = _parse_observed(args.observed)
        call = call_from_observed_substitutions(args.drug, obs)
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.fungal_erg11_caller import call_erg11  # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: fungal genome mode needs scripts/fungal_erg11_caller + BLAST+ ({e}). "
                  "Use --observed with known substitutions for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_erg11(str(args.genome_fasta), str(args.erg11_ref), args.drug)
        prov = {"mode": "blast-erg11", "erg11_ref": str(args.erg11_ref)}
    else:
        print("ERROR: fungal drug needs --genome-fasta OR --observed GENE:SUB[,...]", file=sys.stderr)
        return 2

    rec = _target_site_record(call, sample_id, args.drug, args.organism, prov,
                              caller_name="dna_decode-fungal-target-mutation-v0",
                              source="hand-curated fungal catalog (no AMRFinder-for-fungi)")
    return _emit_target_site(rec, call, sample_id, args)


def _emit_target_site(rec: dict, call, sample_id: str, args) -> int:
    """Shared output for the fungal + antimalarial target-site branches."""
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  drug: {args.drug}  organism: {args.organism}")
        print(f"CALL: {call.prediction}  [{rec['confidence']} | {len(call.determinants)} determinant(s)]")
        for d in call.determinants:
            print(f"  driven by: {d}")
        if not call.determinants:
            print("  driven by: (no catalogued target-site resistance mutation)")
        if call.undetectable_mechanisms:
            print(f"  blind spots (an S call can't rule out): {', '.join(call.undetectable_mechanisms)}")
        print(f"  {call.caveat}")
        print(f"  {one_line(rec['validation'])}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if call.prediction != "INDETERMINATE" else 4


def _antimalarial_main(args) -> int:
    """Antimalarial K13 target-site decoder branch (routed when --drug is an antimalarial drug)."""
    if args.organism == "Escherichia":            # relabel the bacterial default on this path
        args.organism = "Plasmodium_falciparum"
    if args.observed is not None:
        call = antimalarial_call_from_observed(args.drug, _parse_observed(args.observed))
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        # Genome mode uses the BLAST codon-mapper (intron-aware / multi-HSP as of 2026-06-10), so both the
        # intronless K13 and the 13-exon pfcrt work. Pick the committed CDS reference by target gene.
        gene = gene_for_drug(args.drug)
        ref_by_gene = {"K13": args.k13_ref, "pfcrt": args.pfcrt_ref}
        ref = ref_by_gene.get(gene)
        if ref is None or not Path(ref).exists():
            print(f"ERROR: genome mode for {args.drug} (gene {gene}) needs a committed {gene} CDS reference"
                  f"{f' at {ref}' if ref else ' (none configured)'}. Use --observed {gene}:SUB for a "
                  f"wheel-only call.", file=sys.stderr)
            return 3
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.pf_kelch13_caller import call_kelch13   # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: antimalarial genome mode needs scripts/pf_kelch13_caller + BLAST+ ({e}). "
                  "Use --observed GENE:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_kelch13(str(args.genome_fasta), str(ref), args.drug, gene=gene)
        prov = {"mode": f"blast-{gene.lower()}", f"{gene.lower()}_ref": str(ref)}
    else:
        print("ERROR: antimalarial drug needs --genome-fasta OR --observed K13:SUB[,...]", file=sys.stderr)
        return 2
    rec = _target_site_record(call, sample_id, args.drug, args.organism, prov,
                              caller_name="dna_decode-antimalarial-k13-target-mutation-v0",
                              source="hand-curated WHO-validated Pfkelch13 catalog (no AMRFinder-for-Plasmodium)")
    return _emit_target_site(rec, call, sample_id, args)


def _hiv_main(args) -> int:
    """HIV-1 target-site decoder branch (5 classes / 4 genes; validated vs HIVDB PhenoSense).

    Routed by drug -> gene: NNRTI/NRTI->RT, PI->PR (protease), INSTI->IN (integrase), CAI->CA (capsid).
    Wheel-only `--observed GENE:SUB[,...]` (e.g. RT:K103N, PR:V82A, IN:Q148H, CA:M66I), OR genome-FASTA mode
    (`--genome-fasta X.fna`) which BLASTs the committed HXB2 CDS reference for that gene vs the assembly and
    codon-maps the substitutions (scripts/hiv_rt_caller; needs BLAST+), mirroring the influenza NA path."""
    if args.organism == "Escherichia":            # relabel the bacterial default on this path
        args.organism = "HIV-1"
    gene = gene_for_hiv_drug(args.drug)            # RT / PR / IN / CA
    ref_by_gene = {"RT": args.hiv_rt_ref, "PR": args.hiv_pr_ref,
                   "IN": args.hiv_in_ref, "CA": args.hiv_ca_ref}
    if args.observed is not None:
        call = call_hiv_observed(args.drug, _parse_observed(args.observed))
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed, "gene": gene}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        ref = ref_by_gene.get(gene)
        if ref is None or not Path(ref).exists():
            print(f"ERROR: genome mode for {args.drug} (gene {gene}) needs a committed HIV-1 {gene} CDS "
                  f"reference at {ref}. Use --observed {gene}:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.hiv_rt_caller import call_hiv_target   # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: HIV genome mode needs scripts/hiv_rt_caller + BLAST+ ({e}). "
                  f"Use --observed {gene}:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_hiv_target(str(args.genome_fasta), str(ref), args.drug, gene)
        prov = {"mode": "blast-hiv-target", "gene": gene, "hiv_ref": str(ref)}
    else:
        print(f"ERROR: HIV drug needs --observed {gene}:SUB[,...] (e.g. {gene}:K103N) OR --genome-fasta X.fna.",
              file=sys.stderr)
        return 2
    rec = _target_site_record(call, sample_id, args.drug, args.organism, prov,
                              caller_name="dna_decode-" + call.rule.replace("_", "-"),
                              source="HIVDB-PhenoSense-validated (in-distribution; see hiv_decoder_report_card)")
    return _emit_target_site(rec, call, sample_id, args)


def _antiviral_main(args) -> int:
    """Antiviral NA target-site decoder branch (routed when --drug is an NA-inhibitor drug). 4th kingdom."""
    if args.organism == "Escherichia":            # relabel the bacterial default on this path
        args.organism = "Influenza_A_virus"
    if args.observed is not None:
        call = antiviral_call_from_observed(args.drug, _parse_observed(args.observed))
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        if not Path(args.na_ref).exists():
            print(f"ERROR: genome mode for {args.drug} (gene NA) needs a committed N1 NA CDS reference at "
                  f"{args.na_ref}. Use --observed NA:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.flu_na_caller import call_neuraminidase   # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: antiviral genome mode needs scripts/flu_na_caller + BLAST+ ({e}). "
                  "Use --observed NA:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_neuraminidase(str(args.genome_fasta), str(args.na_ref), args.drug, gene="NA")
        prov = {"mode": "blast-na", "na_ref": str(args.na_ref)}
    else:
        print("ERROR: antiviral drug needs --genome-fasta OR --observed NA:SUB[,...]", file=sys.stderr)
        return 2
    rec = _target_site_record(call, sample_id, args.drug, args.organism, prov,
                              caller_name="dna_decode-antiviral-na-target-mutation-v0",
                              source="hand-curated CDC/WHO-recognized influenza NA marker catalog (no AMRFinder-for-influenza)")
    return _emit_target_site(rec, call, sample_id, args)


def _sarscov2_main(args) -> int:
    """SARS-CoV-2 Mpro target-site decoder branch (nirmatrelvir/ensitrelvir/lufotrelvir; validated vs CoV-RDB).

    Wheel-only `--observed Mpro:E166V[,...]`, OR genome-FASTA mode (`--genome-fasta X.fna`) which BLASTs the
    committed Wuhan-Hu-1 Mpro CDS reference vs the assembly and codon-maps the substitutions
    (scripts/sarscov2_caller; needs BLAST+), mirroring the HIV / influenza-NA paths."""
    if args.organism == "Escherichia":            # relabel the bacterial default on this path
        args.organism = "SARS-CoV-2"
    gene = gene_for_sarscov2_drug(args.drug)       # Mpro
    if args.observed is not None:
        call = call_sarscov2_observed(args.drug, _parse_observed(args.observed))
        sample_id = args.sample_id or "observed"
        prov = {"mode": "observed-substitutions", "observed": args.observed}
    elif args.genome_fasta is not None:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        ref = args.sarscov2_mpro_ref
        if not Path(ref).exists():
            print(f"ERROR: genome mode for {args.drug} (gene {gene}) needs a committed SARS-CoV-2 Mpro CDS "
                  f"reference at {ref}. Use --observed {gene}:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            from scripts.sarscov2_caller import call_sarscov2_target   # repo-only; needs BLAST+
        except ImportError as e:
            print(f"ERROR: SARS-CoV-2 genome mode needs scripts/sarscov2_caller + BLAST+ ({e}). "
                  f"Use --observed {gene}:SUB for a wheel-only call.", file=sys.stderr)
            return 3
        call = call_sarscov2_target(str(args.genome_fasta), str(ref), args.drug, gene)
        prov = {"mode": "blast-sarscov2-mpro", "gene": gene, "sarscov2_mpro_ref": str(ref)}
    else:
        print(f"ERROR: SARS-CoV-2 drug needs --observed {gene}:SUB[,...] (e.g. {gene}:E166V) OR "
              "--genome-fasta X.fna.", file=sys.stderr)
        return 2
    rec = _target_site_record(call, sample_id, args.drug, args.organism, prov,
                              caller_name="dna_decode-sarscov2-mpro-target-mutation-v0",
                              source="Stanford CoV-RDB selection-derived Mpro catalog (validate vs measured fold-change)")
    return _emit_target_site(rec, call, sample_id, args)


def _run_amrfinder_for_genome(fasta: Path, sample_id: str, out_root: Path, db: Path,
                              organism: str = "Escherichia") -> Path:
    """Genome mode: lazily import the repo's AMRFinder Docker runner (not in the wheel).

    `organism` selects AMRFinder's `-O` (organism-specific point-mutation detection — gyrA/parC QRDR calls
    are organism-specific, so a Klebsiella genome MUST use 'Klebsiella_pneumoniae' or its QRDR is missed)."""
    try:
        import scripts.drug_mechanism_audit as dma  # repo-only; needs Docker + DB
        from scripts.drug_mechanism_audit import _run_amrfinder
    except ImportError as e:
        raise RuntimeError(
            "genome mode needs the repo's AMRFinder runner (scripts/drug_mechanism_audit) + Docker + a "
            "Docker-readable DB at --amrfinder-db; not available in a wheel install. Use --amrfinder-run "
            f"with a precomputed run instead. ({e})"
        ) from e
    out_dir = out_root / (sample_id or fasta.stem)
    out_dir.mkdir(parents=True, exist_ok=True)
    if db:
        dma.AMRFINDER_DB = str(db)
    _run_amrfinder(fasta, out_dir, organism=organism)
    return out_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-amr",
                                 description="Deterministic AMR R/S decoder from AMRFinder curated determinants")
    ap.add_argument("--drug", required=True,
                    choices=sorted(set(supported_drugs()) | set(supported_fungal_drugs())
                                   | set(supported_antimalarial_drugs()) | set(supported_antiviral_drugs())
                                   | set(all_supported_hiv_drugs()) | set(all_supported_sarscov2_drugs())),
                    metavar="DRUG", help="bacterial (AMRFinder engine), fungal (BLAST-ERG11 engine), "
                                         "antimalarial (BLAST-Pfkelch13 engine), or antiviral "
                                         "(BLAST-influenza-NA engine) drug")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--amrfinder-run", type=Path, help="[bacterial] existing AMRFinder run dir (main.tsv)")
    src.add_argument("--genome-fasta", type=Path, help="genome FASTA (bacterial: AMRFinder via Docker; "
                                                       "fungal: BLAST ERG11 via scripts/fungal_erg11_caller)")
    src.add_argument("--observed", default=None, help="[fungal] known substitutions GENE:SUB[,...] "
                                                      "(e.g. ERG11:Y132F) — pure, no BLAST")
    ap.add_argument("--erg11-ref", type=Path, default=_DEFAULT_ERG11_REF,
                    help="[fungal] in-frame ERG11 CDS reference FASTA (default: committed C. auris allele)")
    ap.add_argument("--k13-ref", type=Path, default=_DEFAULT_K13_REF,
                    help="[antimalarial] in-frame Pfkelch13 CDS reference FASTA (default: committed 3D7 allele)")
    ap.add_argument("--pfcrt-ref", type=Path, default=_DEFAULT_PFCRT_REF,
                    help="[antimalarial] in-frame pfcrt CDS reference FASTA (default: committed 3D7 allele; "
                         "genome mode is intron-aware so a genomic pfcrt allele works)")
    ap.add_argument("--na-ref", type=Path, default=_DEFAULT_NA_REF,
                    help="[antiviral] in-frame influenza N1 NA CDS reference FASTA (default: committed "
                         "NC_026434.1 A/California/07/2009 allele, WT His275)")
    ap.add_argument("--hiv-rt-ref", type=Path, default=_DEFAULT_HIV_RT_REF,
                    help="[HIV] in-frame HIV-1 RT CDS reference FASTA (default: committed HXB2 "
                         "K03455.1:2550-4229 allele, consensus-B WT at every DRM position)")
    ap.add_argument("--hiv-pr-ref", type=Path, default=_DEFAULT_HIV_PR_REF,
                    help="[HIV] in-frame protease CDS reference FASTA (default: committed HXB2 K03455.1:2253-2549)")
    ap.add_argument("--hiv-in-ref", type=Path, default=_DEFAULT_HIV_IN_REF,
                    help="[HIV] in-frame integrase CDS reference FASTA (default: committed HXB2 K03455.1:4230-5093)")
    ap.add_argument("--hiv-ca-ref", type=Path, default=_DEFAULT_HIV_CA_REF,
                    help="[HIV] in-frame capsid CDS reference FASTA (default: committed HXB2 K03455.1:1186-1878)")
    ap.add_argument("--sarscov2-mpro-ref", type=Path, default=_DEFAULT_SARSCOV2_MPRO_REF,
                    help="[SARS-CoV-2] in-frame Mpro (3CLpro/nsp5) CDS reference FASTA (default: committed "
                         "Wuhan-Hu-1 NC_045512.2:10055-10972 allele, WT at every catalogued position)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--organism", default="Escherichia",
                    help="AMRFinder -O organism for genome mode (organism-specific QRDR point-mutation "
                         "detection). E.g. Escherichia (default), Klebsiella_pneumoniae. Validated "
                         "cross-organism for E. coli + K. pneumoniae.")
    ap.add_argument("--amrfinder-db", type=Path, default=Path("data/amrfinder_db"),
                    help="AMRFinder DB root (Docker-readable; default data/amrfinder_db)")
    ap.add_argument("--out-root", type=Path, default=Path("data/amrfinder_runs"))
    ap.add_argument("--resistance-threshold", type=int, default=None,
                    help="min #curated determinants for an R call. Default: per-drug validated config "
                         "(cipro=2 QRDR; cef=1 + extended-spectrum refinement; tet/gent=1). Pass an int to override.")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    # Route by drug: fungal -> BLAST-ERG11; antimalarial -> BLAST-Pfkelch13; bacterial -> AMRFinder.
    if args.drug in supported_fungal_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; fungal drugs use --genome-fasta or --observed",
                  file=sys.stderr)
            return 2
        return _fungal_main(args)

    if args.drug in supported_antimalarial_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; antimalarial drugs use --genome-fasta or --observed",
                  file=sys.stderr)
            return 2
        return _antimalarial_main(args)

    if args.drug in supported_antiviral_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; antiviral drugs use --genome-fasta or --observed",
                  file=sys.stderr)
            return 2
        return _antiviral_main(args)

    if args.drug in all_supported_hiv_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; HIV drugs use --observed GENE:SUB[,...]",
                  file=sys.stderr)
            return 2
        return _hiv_main(args)

    if args.drug in all_supported_sarscov2_drugs():
        if args.amrfinder_run is not None:
            print("ERROR: --amrfinder-run is bacterial-only; SARS-CoV-2 drugs use --observed Mpro:SUB[,...]",
                  file=sys.stderr)
            return 2
        return _sarscov2_main(args)

    if args.observed is not None:
        print("ERROR: --observed is fungal-only; bacterial drugs use --amrfinder-run or --genome-fasta",
              file=sys.stderr)
        return 2

    if args.amrfinder_run:
        run_dir = args.amrfinder_run
        sample_id = args.sample_id or run_dir.name
    else:
        if not args.genome_fasta.exists():
            print(f"ERROR: genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            return 2
        sample_id = args.sample_id or args.genome_fasta.stem
        try:
            run_dir = _run_amrfinder_for_genome(args.genome_fasta, sample_id, args.out_root,
                                                args.amrfinder_db, organism=args.organism)
        except Exception as e:
            print(f"ERROR: AMRFinder run failed ({type(e).__name__}: {e}).", file=sys.stderr)
            return 3

    # Forward --organism so a calibrated registry entry (opt-in) is used when the organism is known. The
    # default 'Escherichia' has no registry entry -> DRUG_RULE default (unchanged); an explicit
    # Campylobacter/Klebsiella/Salmonella resolves its independent-cohort-validated config, and an
    # EXPRESSION_FLOOR organism (Acinetobacter/Pseudomonas carbapenem) returns prediction 'ABSTAIN'.
    # Pass the genome FASTA (genome mode only; None for --amrfinder-run) so the EXPRESSION_FLOOR
    # expression-context override can read the assembly when its registry block is enabled (opt-in).
    call = call_resistance(run_dir / "main.tsv", args.drug, args.resistance_threshold,
                           organism=args.organism, genome_fasta=args.genome_fasta)
    rec = {
        "sample_id": sample_id, "drug": args.drug,
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-mechanism-call-v1",
        "prediction": call["prediction"], "confidence": call["confidence"],
        "n_determinants": call["n_determinants"], "determinants": call["determinants"],
        "resistance_threshold": call.get("resistance_threshold"),
        "undetectable_mechanisms": call.get("undetectable_mechanisms", []),
        "caller": {"name": "dna_decode-amr-rules-v1", "rule": call["rule"],
                   "source": "AMRFinderPlus curated main.tsv", "caller_is_independent_baseline": False},
        "caveat": call["caveat"],
        "validation": trust_block(args.drug, args.organism),
        "provenance": {"amrfinder_run": str(run_dir), "amrfinder_image": AMRFINDER_IMAGE_PINNED,
                       "amrfinder_organism": args.organism},
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  drug: {args.drug}  organism: {args.organism}")
        if call["prediction"] == "ABSTAIN":
            print("CALL: ABSTAIN  [gene-presence cannot decode this organism×drug]")
            print(f"  {call['caveat']}")
        else:
            nd = call["n_determinants"]
            print(f"CALL: {call['prediction']}  [{call['confidence']} | {nd} determinant(s)]")
            for x in call["determinants"]:
                print(f"  driven by: {x['symbol']}  ({x['subclass'] or x['class']}, {x['pct_identity']}% id)")
            if not call["determinants"]:
                print("  driven by: (no curated resistance determinants for this drug)")
            print(f"  {call['caveat']}")
        print(f"  {one_line(rec['validation'])}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return {"INDETERMINATE": 4, "ABSTAIN": 5}.get(call["prediction"], 0)


if __name__ == "__main__":
    raise SystemExit(main())
