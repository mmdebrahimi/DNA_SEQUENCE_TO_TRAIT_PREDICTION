"""EP-7 G1 — assemble C. auris SRA runs into genomes for the deterministic ERG11 caller.

For each row of a cohort TSV (with a `run` SRA accession column), this:
  1. prefetch + fasterq-dump the run (sra-tools Docker image) -> paired FASTQs,
  2. de-novo assemble (SPAdes Docker image, --isolate) -> contigs,
  3. copy the contigs to `<assemblies>/<isolate_id>.fna`.
Then writes an AUGMENTED cohort TSV adding a `genome_fasta` column pointing at each assembly, ready for
`scripts/build_fungal_cohort.py` (which BLASTs the ERG11 CDS reference vs each assembly -> R/S call).

Restartable: a run whose `<assemblies>/<isolate>.fna` already exists is skipped. Reads/intermediate SPAdes
dirs live under `--workdir` (default D:, which has the free space) and can be pruned after.

Docker images are pinned + configurable (--sra-image / --spades-image); the project routes all
bioinformatics tools through Docker on this Windows host (see tools/docker_runner). Small MiSeq C. auris
genomes (~12.5 Mb) assemble in ~10-20 min each on a few cores.

This is COMPUTE (CPU), the user-greenlit "minor compute on this laptop" path (no GPU, no money).
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import time  # noqa: E402

from tools import docker_runner  # noqa: E402
from tools.docker_runner import DockerRunnerError  # noqa: E402

SRA_IMAGE = "quay.io/biocontainers/sra-tools:3.1.1--h4304569_2"
SPADES_IMAGE = "quay.io/biocontainers/spades:3.15.5--h95f258a_1"
SKESA_IMAGE = "quay.io/biocontainers/skesa:2.5.1--hdcf5f25_1"


def _read_cohort(tsv: Path) -> list[dict]:
    return list(csv.DictReader(tsv.read_text(encoding="utf-8").splitlines(), delimiter="\t"))


def fetch_reads(run: str, reads_dir: Path, sra_image: str, timeout: float) -> tuple[Path, Path]:
    """prefetch + fasterq-dump <run> -> reads_dir/<run>_1.fastq, _2.fastq. Returns the pair."""
    reads_dir.mkdir(parents=True, exist_ok=True)
    r1 = reads_dir / f"{run}_1.fastq"
    r2 = reads_dir / f"{run}_2.fastq"
    if r1.exists() and r2.exists():
        return r1, r2
    mounts = {str(reads_dir): "/work"}
    # prefetch into /work/<run>/<run>.sra, then fasterq-dump it
    docker_runner.run(sra_image, ["prefetch", "-O", "/work", run],
                      mounts=mounts, check=True, timeout=timeout)
    docker_runner.run(sra_image,
                      ["fasterq-dump", "--split-files", "-e", "4", "-O", "/work", f"/work/{run}/{run}.sra"],
                      mounts=mounts, check=True, timeout=timeout)
    if not (r1.exists() and r2.exists()):
        raise FileNotFoundError(f"fasterq-dump did not produce paired FASTQs for {run} in {reads_dir}")
    return r1, r2


def assemble(run: str, r1: Path, r2: Path, asm_dir: Path, out_fna: Path,
             spades_image: str, threads: int, mem_gb: int, timeout: float) -> Path:
    """SPAdes --isolate -> contigs.fasta copied to out_fna."""
    asm_dir.mkdir(parents=True, exist_ok=True)
    mounts = {str(r1.parent): "/reads", str(asm_dir): "/asm"}
    docker_runner.run(
        spades_image,
        ["spades.py", "--isolate", "-1", f"/reads/{r1.name}", "-2", f"/reads/{r2.name}",
         "-o", "/asm", "-t", str(threads), "-m", str(mem_gb)],
        mounts=mounts, check=True, timeout=timeout)
    contigs = asm_dir / "contigs.fasta"
    if not contigs.exists():
        raise FileNotFoundError(f"SPAdes produced no contigs.fasta for {run} in {asm_dir}")
    out_fna.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(contigs, out_fna)
    return out_fna


MINIMAP2_IMAGE = "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"


def map_erg11(run: str, r1: Path, r2: Path, work_dir: Path, out_fna: Path, erg11_ref: Path,
              minimap2_image: str, samtools_image: str, threads: int, min_depth: int,
              timeout: float) -> Path:
    """TARGETED: map reads to the ERG11 CDS reference + call consensus -> out_fna (an ERG11-locus FASTA).

    Far faster than whole-genome assembly when only ERG11 is needed: at high coverage only the ERG11 reads
    align to the tiny reference, so this is seconds-minutes/isolate regardless of total depth. The consensus
    is over the reference's coordinates (intronless C. auris ERG11 -> colinear), and `call_erg11` BLASTs the
    same ref CDS against it to extract substitutions — identical downstream contract to an assembly.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = erg11_ref.parent
    mounts = {str(ref_dir): "/ref:ro", str(r1.parent): "/reads:ro", str(work_dir): "/work"}
    docker_runner.run(
        minimap2_image,
        ["minimap2", "-ax", "sr", "-t", str(threads), f"/ref/{erg11_ref.name}",
         f"/reads/{r1.name}", f"/reads/{r2.name}", "-o", "/work/aln.sam"],
        mounts=mounts, check=True, timeout=timeout)
    docker_runner.run(
        samtools_image, ["samtools", "sort", "-@", str(threads), "/work/aln.sam", "-o", "/work/aln.bam"],
        mounts=mounts, check=True, timeout=timeout)
    docker_runner.run(
        samtools_image,
        ["samtools", "consensus", "-f", "fasta", "--min-depth", str(min_depth),
         "/work/aln.bam", "-o", "/work/erg11_cons.fa"],
        mounts=mounts, check=True, timeout=timeout)
    cons = work_dir / "erg11_cons.fa"
    if not cons.exists() or cons.stat().st_size == 0:
        raise FileNotFoundError(f"samtools consensus produced nothing for {run} (no ERG11 coverage?)")
    out_fna.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(cons, out_fna)
    return out_fna


def assemble_skesa(run: str, r1: Path, r2: Path, asm_dir: Path, out_fna: Path,
                   skesa_image: str, threads: int, mem_gb: int, timeout: float) -> Path:
    """SKESA (fast isolate assembler) -> contigs copied to out_fna. ~5-10x faster than SPAdes multi-K."""
    asm_dir.mkdir(parents=True, exist_ok=True)
    mounts = {str(r1.parent): "/reads", str(asm_dir): "/asm"}
    docker_runner.run(
        skesa_image,
        ["skesa", "--reads", f"/reads/{r1.name},/reads/{r2.name}",
         "--cores", str(threads), "--memory", str(mem_gb), "--contigs_out", "/asm/contigs.fasta"],
        mounts=mounts, check=True, timeout=timeout)
    contigs = asm_dir / "contigs.fasta"
    if not contigs.exists() or contigs.stat().st_size == 0:
        raise FileNotFoundError(f"SKESA produced no contigs for {run} in {asm_dir}")
    out_fna.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(contigs, out_fna)
    return out_fna


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True, help="cohort TSV with isolate_id + run columns")
    ap.add_argument("--workdir", default="D:/dna_decode_cache/fungal_g1",
                    help="root for reads/ + asm/ + assemblies/ (needs free space)")
    ap.add_argument("--out-tsv", default=None,
                    help="augmented cohort TSV with genome_fasta column (default <cohort>.assembled.tsv)")
    ap.add_argument("--method", choices=["map", "assemble"], default="map",
                    help="map = targeted ERG11 read-mapping+consensus (fast, default); assemble = full WGS")
    ap.add_argument("--erg11-ref",
                    default="data/fungal_ref/Cauris_ERG11_cds.fna",
                    help="ERG11 CDS reference for --method map")
    ap.add_argument("--min-depth", type=int, default=3, help="min read depth for consensus base (map mode)")
    ap.add_argument("--assembler", choices=["skesa", "spades"], default="skesa",
                    help="(assemble mode) skesa = fast isolate assembler (default); spades = thorough multi-K")
    ap.add_argument("--sra-image", default=SRA_IMAGE)
    ap.add_argument("--spades-image", default=SPADES_IMAGE)
    ap.add_argument("--skesa-image", default=SKESA_IMAGE)
    ap.add_argument("--minimap2-image", default=MINIMAP2_IMAGE)
    ap.add_argument("--samtools-image", default=SAMTOOLS_IMAGE)
    ap.add_argument("--threads", type=int, default=4, help="threads PER SPAdes job")
    ap.add_argument("--jobs", type=int, default=1, help="concurrent isolates (RAM ~3GB/job; mem-gb is per job)")
    ap.add_argument("--mem-gb", type=int, default=12, help="SPAdes RAM cap PER job (GB)")
    ap.add_argument("--fetch-timeout", type=float, default=2400.0)
    ap.add_argument("--assemble-timeout", type=float, default=3600.0)
    ap.add_argument("--limit", type=int, default=0, help="assemble only the first N rows (0 = all)")
    ap.add_argument("--keep-reads", action="store_true", help="don't delete FASTQs after assembly")
    ap.add_argument("--retries", type=int, default=2,
                    help="per-isolate retries on transient DockerRunnerError (e.g. daemon 125)")
    a = ap.parse_args(argv)

    work = Path(a.workdir)
    reads_root = work / "reads"
    asm_root = work / "asm"
    asmb_dir = work / "assemblies"
    rows = _read_cohort(Path(a.cohort))
    if a.limit:
        rows = rows[:a.limit]
    out_tsv = Path(a.out_tsv) if a.out_tsv else Path(a.cohort).with_suffix(".assembled.tsv")

    def process_one(i: int, row: dict):
        iso, run = row["isolate_id"], row["run"]
        out_fna = asmb_dir / f"{iso}.fna"
        tag = f"[{i}/{len(rows)}] {iso} ({run})"
        if out_fna.exists():
            print(f"{tag}: assembly exists, skip", flush=True)
            row["genome_fasta"] = str(out_fna)
            return ("done", row)
        last_err = None
        for attempt in range(1, a.retries + 2):  # retries+1 total attempts
          try:
            if attempt > 1:
                print(f"{tag}: retry {attempt-1}/{a.retries} after transient error", flush=True)
                time.sleep(10 * attempt)  # back off (lets a hiccuping daemon recover)
            print(f"{tag}: fetching reads...", flush=True)
            r1, r2 = fetch_reads(run, reads_root / run, a.sra_image, a.fetch_timeout)
            print(f"{tag}: {a.method} ({a.assembler if a.method=='assemble' else 'minimap2+samtools'})...",
                  flush=True)
            if a.method == "map":
                map_erg11(run, r1, r2, asm_root / iso, out_fna, Path(a.erg11_ref),
                          a.minimap2_image, a.samtools_image, a.threads, a.min_depth, a.assemble_timeout)
            elif a.assembler == "skesa":
                assemble_skesa(run, r1, r2, asm_root / iso, out_fna, a.skesa_image, a.threads,
                               a.mem_gb, a.assemble_timeout)
            else:
                assemble(run, r1, r2, asm_root / iso, out_fna, a.spades_image, a.threads, a.mem_gb,
                         a.assemble_timeout)
            row["genome_fasta"] = str(out_fna)
            print(f"{tag}: OK -> {out_fna}", flush=True)
            if not a.keep_reads:
                shutil.rmtree(reads_root / run, ignore_errors=True)
                shutil.rmtree(asm_root / iso, ignore_errors=True)
            return ("done", row)
          except DockerRunnerError as e:
            last_err = e  # transient (e.g. daemon 125) -> retry
            print(f"{tag}: DockerRunnerError (attempt {attempt}): {str(e)[:120]}", flush=True)
            continue
          except Exception as e:
            print(f"{tag}: FAILED — {type(e).__name__}: {str(e)[:200]}", flush=True)
            return ("failed", (iso, run, str(e)[:120]))
        print(f"{tag}: FAILED after {a.retries} retries — {str(last_err)[:160]}", flush=True)
        return ("failed", (iso, run, f"retries exhausted: {str(last_err)[:100]}"))

    done, failed = [], []
    if a.jobs <= 1:
        results = [process_one(i, row) for i, row in enumerate(rows, 1)]
    else:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            results = list(ex.map(lambda ir: process_one(ir[0], ir[1]),
                                  [(i, row) for i, row in enumerate(rows, 1)]))
    for status, payload in results:
        (done if status == "done" else failed).append(payload)

    # write augmented TSV (only successfully-assembled rows; build_fungal_cohort reads genome_fasta)
    if done:
        cols = ["isolate_id", "run", "clade", "fluconazole_mic", "genome_fasta"]
        with open(out_tsv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in done:
                w.writerow(row)
        print(f"\nwrote {out_tsv}  ({len(done)} assembled, {len(failed)} failed)")
    if failed:
        print("FAILED:", failed)
    return 0 if done and not failed else (1 if failed else 2)


if __name__ == "__main__":
    raise SystemExit(main())
