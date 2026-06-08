"""PointFinder chromosomal point-mutation caller — blastn ref gene vs genome + codon-position lookup.

For each reference gene CDS (.fsa), blastn vs the assembly, recover the subject AA at each catalogued codon
(typing.codon_map, gap-aware), and emit a resistance mutation when the subject AA is one of the Res_codon
values for that (gene, codon_pos) in resistens-overview.txt. Reuses pathotype/vf_runner's BLAST resolvers.
Offline-safe: no blastn / no DB -> status 'unavailable' (never raises).

resistens-overview.txt columns: Gene_ID  Gene_name  Codon_pos  Ref_nuc  Ref_codon  Res_codon  Resistance  ...
(Res_codon + Resistance may be comma-lists). 'Required_mut' (epistasis) is recorded but NOT enforced in v0.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from dna_decode.pathotype.vf_runner import _find_makeblastdb, find_blastn
from dna_decode.typing.codon_map import subject_aa_by_codon, translate

POINT_IDENTITY_THRESHOLD = 90.0   # min HSP identity to trust the gene alignment


def parse_overview(path: str | Path) -> dict[tuple[str, int], dict]:
    """resistens-overview.txt -> {(gene, codon_pos): {ref_aa, res_aas:set, resistances:set, required:set}}."""
    out: dict[tuple[str, int], dict] = {}
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 7:
            continue
        gene = f[1].strip() or f[0].strip()
        try:
            pos = int(f[2])
        except ValueError:
            continue
        ref_aa = f[4].strip()
        res_aas = {a.strip() for a in f[5].split(",") if a.strip()}
        resist = {r.strip() for r in f[6].split(",") if r.strip()}
        required = {x.strip() for x in (f[10].split(",") if len(f) > 10 else []) if x.strip()}
        key = (gene, pos)
        if key in out:
            out[key]["res_aas"] |= res_aas
            out[key]["resistances"] |= resist
            out[key]["required"] |= required
        else:
            out[key] = {"ref_aa": ref_aa, "res_aas": res_aas, "resistances": resist, "required": required}
    return out


def _read_single_fasta(path: str | Path) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def _best_hsp(genome_fasta, ref_fasta, blastn, makeblastdb, timeout):
    """Return (qstart, qseq, sseq) of the best (highest-bitscore) identity-passing HSP, or None."""
    with tempfile.TemporaryDirectory(prefix="pf_") as td:
        dbp = str(Path(td) / "g")
        subprocess.run([makeblastdb, "-in", str(genome_fasta), "-dbtype", "nucl", "-out", dbp],
                       check=True, capture_output=True, text=True, timeout=timeout)
        out = subprocess.run(
            [blastn, "-query", str(ref_fasta), "-db", dbp,
             "-outfmt", "6 qstart qend sstart send qseq sseq bitscore",
             "-perc_identity", str(POINT_IDENTITY_THRESHOLD), "-max_target_seqs", "5", "-evalue", "1e-20"],
            check=True, capture_output=True, text=True, timeout=timeout).stdout
    rows = [ln.split("\t") for ln in out.strip().splitlines() if ln.strip()]
    rows = [r for r in rows if len(r) >= 7]
    if not rows:
        return None
    best = max(rows, key=lambda r: float(r[6]))
    qstart = int(best[0])
    if qstart > int(best[1]):    # query should be plus-strand; skip exotic orientation in v0
        return None
    return qstart, best[4], best[5]


def call_point_mutations(genome_fasta: str | Path, gene_refs: dict[str, str | Path],
                         overview: dict[tuple[str, int], dict], *, blastn_bin: str | None = None,
                         timeout: int = 600) -> dict:
    """Call catalogued point mutations across `gene_refs` (gene -> ref CDS .fsa). Returns mutations +
    the union of conferred resistances. Offline-safe."""
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return {"status": "unavailable", "tool": "blastn", "mutations": [], "resistances": [],
                "reason": "blastn not found (set $BLASTN_BIN or install NCBI BLAST+)"}
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return {"status": "unavailable", "tool": "makeblastdb", "mutations": [], "resistances": [],
                "reason": "makeblastdb not found alongside blastn"}

    mutations = []
    genes_aligned = []
    for gene, ref in gene_refs.items():
        if not Path(ref).exists():
            continue
        ref_prot = translate(_read_single_fasta(ref))
        try:
            hsp = _best_hsp(genome_fasta, ref, blastn, makeblastdb, timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
            return {"status": "unavailable", "tool": "blastn", "mutations": [], "resistances": [],
                    "reason": f"blastn invocation failed: {type(e).__name__}"}
        if hsp is None:
            continue
        genes_aligned.append(gene)
        qstart, qseq, sseq = hsp
        aa_by_codon = subject_aa_by_codon(qseq, sseq, qstart, len(ref_prot))
        for (g, pos), spec in overview.items():
            if g != gene:
                continue
            subj_aa = aa_by_codon.get(pos)
            if subj_aa and subj_aa in spec["res_aas"]:
                ref_aa = spec["ref_aa"] or (ref_prot[pos - 1] if pos - 1 < len(ref_prot) else "?")
                mutations.append({"gene": gene, "mutation": f"{ref_aa}{pos}{subj_aa}", "codon_pos": pos,
                                  "resistances": sorted(spec["resistances"]),
                                  "required_mut": sorted(spec["required"])})
    resistances = sorted({r for m in mutations for r in m["resistances"]})
    return {"status": "ok", "tool": "blastn", "method": "pointfinder_blastn_codonmap_v0",
            "genes_aligned": sorted(genes_aligned),
            "mutations": sorted(mutations, key=lambda m: (m["gene"], m["codon_pos"])),
            "resistances": resistances}
