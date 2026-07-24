"""Phage genome -> host-receptor-class caller (genome-homology receptor TRANSFER).

The v0 caller for the phage receptor cell (`dna_decode/data/phage_receptor.py`): a query phage
inherits the receptor of its nearest genome-BLAST neighbour among a reference set of phages with
known receptors. This measures how well receptor usage transfers along genome similarity - the
honest scientific question - rather than claiming a solved receptor-binding-protein -> receptor map.

Reuses the project's `blastn` resolver (`dna_decode.pathotype.vf_runner.find_blastn`; native BLAST+
at C:/Users/Farshad/ncbi-blast/bin on this host). Offline-safe: no blastn -> status INDETERMINATE,
never a fabricated call.

Leave-one-out validation (`leave_one_out`) is the in-distribution number: each labelled phage is
predicted from the OTHERS and checked against its measured receptor. Because receptor is genus-
conserved but VARIES within the T-even Tevenvirinae by receptor-binding protein, transfer accuracy
is reported per-receptor as well as overall (honest, not a single flattering headline).
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from dna_decode.pathotype.vf_runner import _find_makeblastdb, find_blastn


@dataclass(frozen=True)
class ReceptorCall:
    status: str                 # "CALLED" | "INDETERMINATE"
    predicted_receptor: str | None
    nearest_label: str | None
    percent_identity: float | None
    aln_bitscore: float | None
    reason: str = ""
    method: str = "genome_homology_transfer_v0"


def _read_fasta_ids(path: str | Path) -> list[str]:
    ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def _write_labeled_db_fasta(refs: dict[str, str], out_fa: Path, exclude: str | None = None) -> dict[str, str]:
    """Concatenate reference genomes into ONE fasta, renaming each record's id to `label`.
    Returns {db_seq_id -> label}. A multi-contig reference gets label, label__1, ... all mapped back.
    """
    id_to_label: dict[str, str] = {}
    with open(out_fa, "w", encoding="utf-8") as out:
        for label, fna in refs.items():
            if exclude is not None and label == exclude:
                continue
            n = 0
            with open(fna, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        seq_id = label if n == 0 else f"{label}__{n}"
                        id_to_label[seq_id] = label
                        n += 1
                        out.write(f">{seq_id}\n")
                    else:
                        out.write(line)
    return id_to_label


def _run_blastn(blastn: str, query: str, db: str) -> list[tuple[str, float, float]]:
    """Return [(subject_id, pident, bitscore), ...] sorted by bitscore desc (outfmt 6)."""
    cmd = [blastn, "-query", query, "-db", db, "-outfmt", "6 sseqid pident bitscore",
           "-max_target_seqs", "5", "-evalue", "1e-10"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    rows: list[tuple[str, float, float]] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                rows.append((parts[0], float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    rows.sort(key=lambda r: r[2], reverse=True)
    return rows


def call_receptor(query_fna: str | Path, refs: dict[str, str], receptors: dict[str, str],
                  *, exclude: str | None = None, blastn_bin: str | None = None) -> ReceptorCall:
    """Predict a phage's receptor by nearest-genome-homology transfer.

    refs: {label -> reference genome fasta path}. receptors: {label -> receptor class}.
    exclude: a label to drop from the reference DB (leave-one-out). Returns INDETERMINATE when
    blastn is absent or no homolog clears the threshold.
    """
    blastn = blastn_bin or find_blastn()
    if not blastn:
        return ReceptorCall("INDETERMINATE", None, None, None, None,
                            reason="blastn not found (set $BLASTN_BIN or install NCBI BLAST+)")
    makeblastdb = _find_makeblastdb(blastn)
    if not makeblastdb:
        return ReceptorCall("INDETERMINATE", None, None, None, None, reason="makeblastdb not found")

    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        db_fa = tdp / "refdb.fna"
        id_to_label = _write_labeled_db_fasta(refs, db_fa, exclude=exclude)
        if not id_to_label:
            return ReceptorCall("INDETERMINATE", None, None, None, None,
                                reason="empty reference set after exclusion")
        db_path = tdp / "refdb"
        mk = subprocess.run([makeblastdb, "-in", str(db_fa), "-dbtype", "nucl", "-out", str(db_path)],
                            capture_output=True, text=True, timeout=120)
        if mk.returncode != 0:
            return ReceptorCall("INDETERMINATE", None, None, None, None,
                                reason=f"makeblastdb failed: {mk.stderr.strip()[:120]}")
        hits = _run_blastn(blastn, str(query_fna), str(db_path))
        if not hits:
            return ReceptorCall("INDETERMINATE", None, None, None, None,
                                reason="no genome homolog cleared e-value threshold")
        sseqid, pident, bitscore = hits[0]
        label = id_to_label.get(sseqid, sseqid)
        receptor = receptors.get(label)
        if receptor is None:
            return ReceptorCall("INDETERMINATE", None, label, pident, bitscore,
                                reason=f"nearest phage {label} has no receptor label")
        return ReceptorCall("CALLED", receptor, label, pident, bitscore)


@dataclass
class LOOResult:
    n_total: int = 0
    n_called: int = 0
    n_correct: int = 0
    per_receptor: dict[str, list[int]] = field(default_factory=dict)  # receptor -> [correct, called]
    predictions: list[dict] = field(default_factory=list)

    @property
    def accuracy(self) -> float | None:
        return (self.n_correct / self.n_called) if self.n_called else None


def leave_one_out(refs: dict[str, str], receptors: dict[str, str],
                  *, blastn_bin: str | None = None) -> LOOResult:
    """Leave-one-out receptor-transfer accuracy over the labelled genome set."""
    res = LOOResult()
    for label, true_receptor in receptors.items():
        res.n_total += 1
        call = call_receptor(refs[label], refs, receptors, exclude=label, blastn_bin=blastn_bin)
        correct = call.status == "CALLED" and call.predicted_receptor == true_receptor
        if call.status == "CALLED":
            res.n_called += 1
            if correct:
                res.n_correct += 1
            bucket = res.per_receptor.setdefault(true_receptor, [0, 0])
            bucket[1] += 1
            if correct:
                bucket[0] += 1
        res.predictions.append({
            "phage": label, "true": true_receptor, "predicted": call.predicted_receptor,
            "nearest": call.nearest_label, "status": call.status,
            "pident": call.percent_identity, "correct": bool(correct),
        })
    return res


def _load_manifest(manifest_tsv: str | Path, genome_dir: str | Path) -> tuple[dict[str, str], dict[str, str]]:
    """Load a BASEL-style manifest (accession, phage_name, genus, family) + attach receptors via the
    catalog's genus/family lineage lookup. Returns (refs{label->fna}, receptors{label->receptor}),
    keeping only phages whose lineage resolves to a catalogued receptor.
    """
    from dna_decode.data.phage_receptor import label_receptor_for_lineage
    refs: dict[str, str] = {}
    receptors: dict[str, str] = {}
    gdir = Path(genome_dir)
    with open(manifest_tsv, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            acc = parts[idx["accession"]]
            genus = parts[idx.get("genus", -1)] if "genus" in idx else ""
            family = parts[idx.get("family", -1)] if "family" in idx else ""
            fna = gdir / f"{acc}.fna"
            if not fna.exists():
                continue
            # only clade-conserved taxa yield a label (RBP-variable clades excluded, not mislabelled)
            receptor = label_receptor_for_lineage([genus, family])
            if receptor is None:
                continue
            refs[acc] = str(fna)
            receptors[acc] = receptor
    return refs, receptors


if __name__ == "__main__":  # pragma: no cover - thin CLI
    import argparse
    import json

    ap = argparse.ArgumentParser(description="Phage genome -> receptor-class caller + LOO validation")
    ap.add_argument("--manifest", help="BASEL manifest TSV (accession,phage_name,genus,family)")
    ap.add_argument("--genome-dir", help="dir of <accession>.fna genomes")
    ap.add_argument("--query", help="single query phage genome fasta (predict its receptor)")
    ap.add_argument("--loo", action="store_true", help="run leave-one-out over the manifest")
    args = ap.parse_args()

    if args.manifest and args.genome_dir:
        refs, receptors = _load_manifest(args.manifest, args.genome_dir)
        print(f"loaded {len(refs)} labelled phages across {len(set(receptors.values()))} receptors")
        if args.query:
            print(json.dumps(call_receptor(args.query, refs, receptors).__dict__, indent=2))
        if args.loo:
            r = leave_one_out(refs, receptors)
            print(json.dumps({"accuracy": r.accuracy, "n_called": r.n_called, "n_total": r.n_total,
                              "per_receptor": r.per_receptor}, indent=2))
