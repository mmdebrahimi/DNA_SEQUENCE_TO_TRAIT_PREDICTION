"""`dna-decode forward` — the forward variant-effect decoder as a first-class CLI command.

The forward direction of the decoder: given a genetic EDIT (a point mutation on a protein), predict the change
in MOLECULAR phenotype (Regime B — enzyme fitness/stability), validated per-variant against measured Deep
Mutational Scanning (the one place the project's label wall does not bind). This is the "edit -> effect"
complement to the deterministic determinant->R/S AMR decoder (`dna-decode amr`), which this does NOT touch.

    dna-decode forward --capabilities                                            # what can this host run?
    dna-decode forward --mutation M69L --protein-seq MSIQ...                      # BLOSUM62, instant, no deps
    dna-decode forward --mutation M69L --protein-seq MSIQ... --method esm2        # learned (needs torch)
    dna-decode forward --mutation M69L --uniprot P00552 --method hybrid           # ESM2+ProSST (validated best)
    dna-decode forward --mutation M69L --protein-seq MSIQ... --method auto        # strongest runnable here
    dna-decode forward --mutation c.205G>A --cds-fasta cds.fna                    # a GENOME edit (nt -> codon -> AA)
    dna-decode forward --mutation c.248C>T --genome-fasta g.fna --annotations a.gff3 --gene gyrA

The nucleotide forms lift the input from a protein mutation to a real genome edit: resolve the HGVS coding
substitution to its codon, classify it SILENT / MISSENSE / NONSENSE, and (for the two that change the protein)
fall through to exactly the same predictor path — so a genome edit inherits every method and the same
degradation behaviour. `--genome-fasta` + `--annotations` + `--gene` extracts that gene's CDS from a whole
genome (strand-aware: a minus-strand gene is reverse-complemented), matching `--gene` on `gene_symbol` first —
the CROSS-STRAIN identifier — then `locus_tag`, then `gene_id`.

Everything here fails CLOSED, because a wrong coordinate is the failure mode this path exists to catch: the
reference base must match the CDS, the CDS must be in frame, a separately supplied protein must agree with the
translation, only single-base `c.` substitutions parse, and a gene name matching several CDS features (a
multi-copy gene, alternative products of one locus, or a joined multi-segment CDS such as the dnaX programmed
frameshift) is REFUSED with the field that actually separates them — never resolved to an arbitrary copy.
A SILENT edit returns the coding-consequence call and makes no effect prediction.

The STRONG methods (esm2 / prosst / gemme / hybrid) are now first-class from the CLI — they compute their
per-protein score table ONCE (heavy: minutes on CPU for ESM2, needs a structure for ProSST, Docker+an MSA for
GEMME) and are validated to beat BLOSUM everywhere (wiki/mavedb_holdout_hybrid_2026-07-23.md). When a method's
dependency is missing, the CLI DEGRADES to the strongest runnable method and says so (--no-degrade to error
instead). `--capabilities` is the honest preflight. Precomputed tables (--esm-table / --prosst-table / ...)
skip the heavy compute.

HONEST SCOPE: this predicts MOLECULAR fitness rank (the DMS-validated quantity), NOT clinical antibiotic
resistance — the antagonistic-selection direction fails for likelihood/severity scorers; use `dna-decode amr`
(the frozen determinant catalogue) for R/S. NOT a clinical tool.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

STRONG = ("esm2", "prosst", "gemme", "hybrid", "auto")


def _read_fasta_seq(path: str) -> str:
    seq = "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))
    return "".join(c for c in seq.upper() if c.isalpha())


def _load_json(path: str | None):
    if not path:
        return None
    return json.loads(Path(path).read_text())


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode forward",
        description="Forward variant-effect decode: a point mutation -> predicted molecular-phenotype change "
                    "(Regime B, DMS-validated). BLOSUM62 (deterministic, no deps) + the strong learned methods "
                    "(esm2/prosst/gemme/hybrid) with a capability preflight + honest graceful degradation.",
        epilog="scope: molecular fitness rank, NOT clinical resistance (use `dna-decode amr` for R/S).",
    )
    ap.add_argument("--capabilities", action="store_true",
                    help="print the capability preflight (installed deps -> runnable methods) and exit")
    ap.add_argument("--mutation", help="point mutation, e.g. M69L (WT-pos-ALT, 1-based)")
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--protein-seq", help="protein amino-acid sequence (verifies the WT residue at the position)")
    src.add_argument("--protein-fasta", help="path to a protein FASTA whose sequence to use")
    cds_src = ap.add_mutually_exclusive_group()
    cds_src.add_argument("--cds-fasta", help="path to a CDS nucleotide FASTA — decode a genome edit given as "
                                             "HGVS `c.<pos><REF>><ALT>` (e.g. c.205G>A)")
    cds_src.add_argument("--cds-seq", help="CDS nucleotide sequence (alternative to --cds-fasta)")
    cds_src.add_argument("--genome-fasta", help="whole-genome nucleotide FASTA -- with --annotations and "
                                                "--gene, extracts that gene's CDS (strand-aware)")
    ap.add_argument("--annotations", help="GFF3 for --genome-fasta")
    ap.add_argument("--gene", help="gene to decode from --genome-fasta (matched on gene_symbol, e.g. gyrA, "
                                   "then locus_tag, then gene_id)")
    ap.add_argument("--protein", default="protein", help="protein name/label (default: protein)")
    ap.add_argument("--method", default="blosum62",
                    choices=["blosum62", "esm2", "prosst", "gemme", "hybrid", "auto"],
                    help="blosum62 (default, no deps) | esm2 | prosst | gemme | hybrid | auto (strongest "
                         "runnable). Learned methods beat blosum62 everywhere but need torch / a structure / "
                         "Docker; run `--capabilities` to see what this host can run.")
    ap.add_argument("--regime", default="B_molecular", choices=["B_molecular", "C_organismal"],
                    help="B_molecular = enzyme fitness/stability (default); C_organismal = abstain (closed negative)")
    # structure sources for ProSST / hybrid
    ap.add_argument("--uniprot", help="UniProt accession -> fetch an AlphaFold structure for ProSST")
    ap.add_argument("--pdb", help="path to a PDB/mmCIF structure for ProSST")
    ap.add_argument("--structure-tokens", help="path to a JSON list of ProSST structure tokens")
    ap.add_argument("--msa", help="path to an MSA (a3m/FASTA) for GEMME")
    # precomputed tables (skip the heavy per-protein compute)
    ap.add_argument("--esm-table", help="path to a precomputed ESM2 score table (JSON)")
    ap.add_argument("--prosst-table", help="path to a precomputed ProSST score table (JSON)")
    ap.add_argument("--gemme-table", help="path to a precomputed GEMME score table (JSON)")
    ap.add_argument("--am-table", help="path to a precomputed AlphaMissense score table (JSON)")
    ap.add_argument("--no-degrade", action="store_true",
                    help="error instead of degrading to the strongest runnable method when a dep is missing")
    ap.add_argument("--json", action="store_true", help="emit the prediction as JSON")
    args = ap.parse_args(argv)

    # --- capability preflight (no mutation required) ---
    if args.capabilities:
        from .capabilities import render_capabilities
        print(render_capabilities())
        return 0

    if not args.mutation:
        ap.error("--mutation is required (unless --capabilities)")

    seq = ""
    if args.protein_fasta:
        try:
            seq = _read_fasta_seq(args.protein_fasta)
        except OSError as e:
            print(f"error: cannot read --protein-fasta {args.protein_fasta}: {e}", file=sys.stderr)
            return 2
    elif args.protein_seq:
        seq = "".join(c for c in args.protein_seq.upper() if c.isalpha())

    # --- CDS / genome-edit path: `c.205G>A` on a coding sequence -> codon -> AA change ---------------
    # Resolves the nucleotide edit to its amino-acid consequence, then FALLS THROUGH to the protein path
    # below, so the genome-edit input inherits every method (esm2/prosst/gemme/hybrid + degradation)
    # rather than becoming a second-class branch.
    cds = ""
    if args.cds_fasta:
        try:
            cds = _read_fasta_seq(args.cds_fasta)
        except OSError as e:
            print(f"error: cannot read --cds-fasta {args.cds_fasta}: {e}", file=sys.stderr)
            return 2
    elif args.cds_seq:
        cds = "".join(c for c in args.cds_seq.upper() if c.isalpha())
    elif args.genome_fasta:
        if not (args.annotations and args.gene):
            ap.error("--genome-fasta requires --annotations <gff3> and --gene <name>")
        # heavy import (pandas) deferred to the branch that needs it
        from dna_decode.data.annotations import extract_cds_sequences, parse_gff3
        from dna_decode.forward import cds_record_key, select_gene_cds
        try:
            table = parse_gff3(args.annotations)
            rec = select_gene_cds(table.to_dict("records"), args.gene)
            cds = extract_cds_sequences(args.genome_fasta, table)[cds_record_key(rec)]
        except (OSError, ValueError, KeyError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        args.protein = args.protein if args.protein != "protein" else args.gene
        print(f"gene lookup: {args.gene} -> {rec.get('seqid')}:{rec.get('start')}-{rec.get('end')} "
              f"strand {rec.get('strand')}  (locus_tag {rec.get('locus_tag') or '-'})")

    if cds:
        from dna_decode.forward import cds_point_edit, parse_hgvs_c, translate_cds
        try:
            nt_pos, ref_base, alt_base = parse_hgvs_c(args.mutation)
            edit = cds_point_edit(cds, nt_pos, ref_base, alt_base)   # verifies REF against the CDS, loudly
            translated = translate_cds(cds)
        except (ValueError, KeyError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 2

        wt_aa, alt_aa, aa_pos = edit["wt_aa"], edit["alt_aa"], edit["aa_pos"]
        consequence = "silent" if wt_aa == alt_aa else ("nonsense" if alt_aa == "*" else "missense")
        print("genome edit -> codon consequence")
        print(f"  CDS: {len(cds)} nt ({len(translated)} aa)   edit: c.{nt_pos}{ref_base}>{alt_base}  "
              f"(REF verified against the CDS)")
        print(f"  codon {aa_pos}: {edit['wt_codon']} -> {edit['alt_codon']}   "
              f"{wt_aa} -> {alt_aa}   consequence: {consequence.upper()}")

        if consequence == "silent":
            print("  synonymous — no protein change, so no variant-effect prediction is made.")
            print("  [scope: this is a CODING-consequence call only. A synonymous edit can still matter via "
                  "codon usage / mRNA structure / splicing — none of which this tool models.]")
            return 0

        # An independently-supplied protein must AGREE with the translation (double coordinate check);
        # otherwise the translation IS the reference.
        if seq and seq[:len(translated)] != translated[:len(seq)]:
            print("error: --protein-* sequence disagrees with the translated --cds-* sequence "
                  "(coordinate/frame mismatch — refusing rather than picking one)", file=sys.stderr)
            return 2
        if not seq:
            seq = translated
        args.mutation = f"{wt_aa}{aa_pos}{alt_aa}"
        print(f"  -> protein mutation: {args.mutation}\n")

    want_strong = args.method in STRONG and args.regime != "C_organismal"

    # --- BLOSUM62 / organismal-abstain: the original wheel-only path (unchanged) ---
    if not want_strong:
        from dna_decode.forward import predict_effect
        try:
            pred = predict_effect(seq, args.mutation, protein=args.protein, method="blosum62",
                                  regime=args.regime)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        d = pred.as_dict()
        d.update({"method_requested": args.method, "method_used": d["method"], "degraded": False,
                  "degrade_reason": None})
        return _emit(d, seq, args.json)

    # --- strong path: resolve against host capabilities, compute tables, degrade honestly ---
    from .deploy import predict_effect_deployable
    try:
        tokens = _load_json(args.structure_tokens)
        d = predict_effect_deployable(
            seq, args.mutation, method=args.method, protein=args.protein,
            uniprot=args.uniprot, pdb_path=args.pdb, structure_tokens=tokens,
            esm_table=_load_json(args.esm_table), prosst_table=_load_json(args.prosst_table),
            gemme_table=_load_json(args.gemme_table), am_table=_load_json(args.am_table),
            msa=args.msa, degrade=not args.no_degrade)
    except (ValueError, RuntimeError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return _emit(d, seq, args.json)


def _emit(d: dict, seq: str, as_json: bool) -> int:
    if as_json:
        print(json.dumps(d, indent=2))
        return 0
    score = d.get("raw_score")
    print("forward variant-effect decode")
    print(f"  protein: {d['protein']}   mutation: {d['mutation']} ({d['wt']}->{d['alt']} at {d['pos']})")
    used = d.get("method_used", d["method"])
    if d.get("degraded"):
        print(f"  method: requested {d['method_requested']} -> USED {used}   (DEGRADED)")
        print(f"  degrade: {d['degrade_reason']}")
    else:
        print(f"  method: {used}   regime: {d['regime']}")
    print(f"  predicted_effect: {d['predicted_effect']}   "
          f"raw_score: {'nan' if score is None or score != score else f'{score:.3f}'}   "
          f"confidence: {d['confidence']}")
    print(f"  phenotype axis: {d['phenotype_axis']}")
    for n in d.get("notes", []):
        print(f"  note: {n}")
    if not seq:
        print("  note: no protein sequence supplied — WT residue NOT verified against a reference")
    if used == "blosum62":
        print("  [BLOSUM62 = deterministic substitution severity — the WEAKEST measured method; ~2.5x below "
              "ESM2 on held-out MaveDB. Run `dna-decode forward --capabilities` to see the stronger methods "
              "this host can run.]")
    else:
        print(f"  [{used} is a learned method — measured to beat BLOSUM62 "
              "(wiki/mavedb_holdout_hybrid_2026-07-23.md: hybrid ~ prosst > esm2 >> blosum62)]")
    print("  [scope: molecular fitness rank, NOT clinical resistance — use `dna-decode amr` for R/S]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
