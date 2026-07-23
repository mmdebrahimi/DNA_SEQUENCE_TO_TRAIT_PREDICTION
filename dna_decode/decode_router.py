"""Input-aware decode router — the 'point at your DNA -> what can this tool decode?' coherence layer.

The tool ships 22 console scripts (dna-amr / dna-pgx / dna-forward / dna-mlst / ...), but a user with a file
in hand has no single place that says WHICH decoders apply to it. This module closes that gap: detect the
input KIND (VCF / protein FASTA / nucleotide FASTA) from the file, then list every applicable decoder with its
one-line claim + honest evidence tier + the exact command to run. Pure + offline (file sniff + a grounded
table + a registry read); no network, no heavy deps.

Design: the applicable-decoder table is curated but GROUNDED -- every route is a real console script and every
example is a real invocation. Evidence tiers are read live from `cell_registry` where a route maps to a cell,
so the honest trust surface travels with the recommendation (never a bare 'run this').
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# nucleotide vs protein alphabets (for FASTA sniffing)
_NT = set("ACGTUN-")
_AA_ONLY = set("EFILPQZ*")  # letters that appear in protein but never in a DNA/RNA sequence


@dataclass(frozen=True)
class Decoder:
    route: str          # console script, e.g. "dna-amr"
    what: str           # one-line: what it decodes
    example: str        # a real example invocation
    track: str          # registry track (for the tier lookup), or "" if not registry-backed


# input_kind -> the decoders that apply. Grounded in the real console scripts (pyproject [project.scripts]).
DECODERS: dict[str, list[Decoder]] = {
    "vcf_human": [
        Decoder("dna-pgx", "pharmacogenomics: CYP2C19/2C9 diplotype + CPIC metabolizer phenotype, VKORC1 warfarin",
                "dna-pgx --gene CYP2C19 --vcf sample.vcf", "pgx"),
        Decoder("dna-clinvar", "Mendelian: curated ClinVar pathogenic/benign calls for the variants carried",
                "dna-clinvar --vcf sample.vcf", "mendelian"),
        Decoder("dna-hla", "HLA typing from the VCF region", "dna-hla --vcf sample.vcf", "hla"),
    ],
    "protein_fasta": [
        Decoder("dna-forward", "variant-effect: a protein point mutation -> molecular-phenotype change (fitness rank)",
                "dna-forward --mutation M69L --protein-fasta prot.fasta", "finder"),
        Decoder("dna-inverse", "inverse design: a target effect -> a ranked list of candidate edits (ranks, never doses)",
                "dna-inverse --protein-fasta prot.fasta --target-percentile 90", "finder"),
    ],
    "nucleotide_fasta": [
        Decoder("dna-amr", "antibiotic resistance R/S (bacterial cipro/cef/tet/gent/mero; also fungal/viral via --drug)",
                "dna-amr --drug ciprofloxacin --genome-fasta genome.fna --sample-id X", "amr"),
        Decoder("dna-pathotype", "E. coli pathotype (EPEC/EHEC/ETEC/UPEC/...) compatibility call",
                "dna-pathotype genome.fna --sample-id X", "finder"),
        Decoder("dna-mlst", "MLST sequence type (7-locus Achtman -> ST)", "dna-mlst genome.fna", "typing"),
        Decoder("dna-serotype", "E. coli O:H serotype (O-antigen + fliC)", "dna-serotype genome.fna", "typing"),
        Decoder("dna-salmserovar", "Salmonella serovar (Kauffmann-White antigenic formula)", "dna-salmserovar genome.fna", "typing"),
        Decoder("dna-ktype", "Klebsiella K (capsule) type via wzi", "dna-ktype genome.fna", "typing"),
        Decoder("dna-pneumo-serotype", "S. pneumoniae capsular serotype (cps locus)", "dna-pneumo-serotype genome.fna", "typing"),
        Decoder("dna-plasmid", "plasmid Inc-replicon typing (is resistance plasmid-borne?)", "dna-plasmid genome.fna", "finder"),
        Decoder("dna-resfinder", "acquired AMR genes (ResFinder DB; independent cross-check vs dna-amr)", "dna-resfinder genome.fna", "finder"),
        Decoder("dna-pointfinder", "chromosomal AMR point mutations (QRDR)", "dna-pointfinder genome.fna", "finder"),
        Decoder("dna-disinfinder", "biocide/disinfectant resistance genes", "dna-disinfinder genome.fna", "finder"),
        Decoder("dna-coloc", "co-localization: is a resistance + a virulence/plasmid marker on the same contig?",
                "dna-coloc genome.fna", ""),
    ],
}
# nucleotide FASTA also admits the forward cell on a CDS (edit->effect on the translated protein)
DECODERS["nucleotide_fasta"].append(
    Decoder("dna-forward", "variant-effect on a CDS via genome-nt edit (translate -> molecular phenotype)",
            "dna-forward --mutation c.205G>A --genome-fasta cds.fna", "finder"))


def detect_input_kind(path: str | Path, *, max_bytes: int = 65536) -> str:
    """Sniff a file -> 'vcf_human' | 'protein_fasta' | 'nucleotide_fasta' | 'unknown'. PURE (reads a prefix).

    VCF by the `##fileformat=VCF` header; FASTA split into protein vs nucleotide by the residue alphabet of the
    first record (a protein-only letter like E/F/I/L/P/Q => protein; otherwise nucleotide)."""
    p = Path(path)
    try:
        head = p.read_text(encoding="utf-8", errors="replace")[:max_bytes]
    except (OSError, UnicodeError):
        return "unknown"
    stripped = head.lstrip()
    if stripped.startswith("##fileformat=VCF") or "\n#CHROM\t" in head:
        return "vcf_human"
    if stripped.startswith(">"):
        # collect residues from the sequence lines of the first record
        seq = []
        started = False
        for ln in head.splitlines():
            if ln.startswith(">"):
                if started:
                    break
                started = True
                continue
            if started:
                seq.append(ln.strip().upper())
        residues = set("".join(seq))
        residues.discard("")
        if residues & _AA_ONLY:
            return "protein_fasta"
        if residues and residues <= _NT | {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}:  # + IUPAC nt ambiguity
            return "nucleotide_fasta"
        # a FASTA whose alphabet is ambiguous (only A/C/G/T/N which are also valid AAs) -> nucleotide by default
        return "nucleotide_fasta" if residues else "unknown"
    return "unknown"


def applicable_decoders(kind: str) -> list[Decoder]:
    """The decoders that apply to an input kind (empty for 'unknown')."""
    return DECODERS.get(kind, [])


def _tier_for(track: str, route: str) -> str:
    """Best-effort honest evidence tier for a route, read LIVE from the cell registry (or '' if not backed)."""
    if not track:
        return ""
    try:
        from dna_decode.data import cell_registry as reg
        # collect all contracts across the registry's builders
        builders = ["_amr_contracts", "_viral_contracts", "_typing_finder_contracts", "_hla_contracts"]
        cells = []
        for b in builders:
            fn = getattr(reg, b, None)
            if callable(fn):
                try:
                    cells += list(fn())
                except Exception:
                    pass
        cells += list(getattr(reg, "_PGX_CONTRACTS", []) or [])
        cells += list(getattr(reg, "_MENDELIAN_CONTRACTS", []) or [])
        cells += list(getattr(reg, "_TRAIT_CONTRACTS", []) or [])
        tiers = {str(getattr(c, "evidence_tier", "")).split(".")[-1] for c in cells
                 if getattr(c, "route", None) == route}
        tiers.discard("")
        return "/".join(sorted(tiers)) if tiers else ""
    except Exception:
        return ""


def render_decode_plan(path: str | Path) -> str:
    """Human-readable 'what can I decode from this file?' report + the exact commands. PURE-ish (registry read)."""
    kind = detect_input_kind(path)
    labels = {"vcf_human": "human VCF", "protein_fasta": "protein FASTA",
              "nucleotide_fasta": "nucleotide FASTA (genome / CDS)", "unknown": "unrecognized"}
    lines = [f"dna-decode: input '{Path(path).name}' detected as {labels.get(kind, kind)}", ""]
    decs = applicable_decoders(kind)
    if not decs:
        lines.append("  No decoder recognizes this input kind. Supported: a nucleotide/protein FASTA, or a VCF.")
        lines.append("  (Sniff is content-based: VCF header, or the FASTA residue alphabet.)")
        return "\n".join(lines)
    lines.append(f"  {len(decs)} applicable decoder(s):")
    for d in decs:
        tier = _tier_for(d.track, d.route)
        tier_s = f"  [tier: {tier}]" if tier else ""
        lines.append(f"  - {d.route:20s} {d.what}{tier_s}")
        lines.append(f"      run: {d.example}")
    lines.append("")
    lines.append("  Honest scope: each decoder emits its own trust surface + abstains when out-of-scope; run "
                 "`dna-decode list` for per-trait validation status.")
    return "\n".join(lines)
