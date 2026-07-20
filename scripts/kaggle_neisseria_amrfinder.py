#!/usr/bin/env python
"""SELF-CONTAINED Kaggle CPU notebook: SRA-assemble + AMRFinder the AR-Bank N. gonorrhoeae isolates.

The AR Bank deposits many N. gonorrhoeae isolates as SRA reads only -> the neisseria_amr cell
cannot be validated on the assembly-download path used for the Gram-negatives. Locally, the SRA-assembly +
AMRFinder path wedges Docker Desktop (WSL2 mount corruption on concurrent bioconda-in-Docker load). Kaggle
runs bioconda NATIVELY (no Docker nesting) with a 12 h CPU budget -> the right host for this.

Per isolate (checkpointed -- each {biosample}.amrfinder.tsv written immediately, so a partial/timed-out run
is still usable): Entrez BioSample->SRR -> prefetch+fasterq-dump -> skesa assemble -> amrfinder
-O Neisseria_gonorrhoeae -> save the full main.tsv to /kaggle/working. Scoring (parse Element symbols ->
call_ng_cefixime / call_ng_ciprofloxacin -> vs CDC labels) is done LOCALLY in the repo (verifiable),
not here -- this notebook is only the AMRFinder factory.

HOW TO RUN ON KAGGLE
--------------------
Settings: Accelerator = None (CPU); Internet = ON. Paste this file into ONE cell. Run.
Outputs: /kaggle/working/<biosample>.amrfinder.tsv (21 files) + assembly_manifest.json.
~15 min/isolate (SRA fetch + skesa + amrfinder) x 21 ~= 5 h; fits the 12 h CPU limit.
"""
import json
import os
import subprocess
import sys
import time
import urllib.request

# 21-isolate powering-optimal AR-Bank N. gonorrhoeae subset: cefixime 13R/8S + ciprofloxacin 12R/9S.
# AMRFinder yields ALL determinants per assembly, so one run scores both drugs downstream.
BIOSAMPLES = [
    "SAMEA3165649",
    "SAMEA3165668",
    "SAMEA3165653",
    "SAMEA3165655",
    "SAMEA3165656",
    "SAMEA3165648",
    "SAMEA3165247",
    "SAMEA3165657",
    "SAMEA3165659",
    "SAMN35332250",
    "SAMN35332251",
    "SAMEA3165293",
    "SAMEA3165270",
    "SAMEA3165275",
    "SAMEA3165272",
    "SAMEA3165236",
    "SAMEA3165292",
    "SAMEA3165249",
    "SAMEA3165250",
    "SAMEA3165241",
    "SAMEA3165296",
]

WORK = "/kaggle/working"
ENVP = "/kaggle/working/bioenv"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def sh(cmd, **kw):
    print(f"$ {cmd}", flush=True)
    return subprocess.run(cmd, shell=True, check=kw.pop("check", True), **kw)


def setup_bioconda():
    """Install sra-tools + skesa + AMRFinderPlus via micromamba (fast solver), fetch the AMRFinder DB."""
    if os.path.exists(f"{ENVP}/bin/amrfinder"):
        print("bioenv already present"); return
    sh("curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba")
    sh(f"./bin/micromamba create -y -p {ENVP} -c conda-forge -c bioconda "
       f"sra-tools skesa ncbi-amrfinderplus", timeout=3600)
    sh(f"{ENVP}/bin/amrfinder -u", timeout=1200)   # download the AMRFinder DB


def entrez_srr(biosample):
    """BioSample -> first SRR run accession (Entrez esearch db=sra + efetch runinfo)."""
    u = f"{NCBI}/esearch.fcgi?db=sra&term={biosample}&retmode=json"
    ids = json.load(urllib.request.urlopen(u, timeout=60))["esearchresult"]["idlist"]
    time.sleep(0.5)
    if not ids:
        return None
    u2 = f"{NCBI}/efetch.fcgi?db=sra&id={ids[0]}&rettype=runinfo&retmode=text"
    txt = urllib.request.urlopen(u2, timeout=60).read().decode()
    time.sleep(0.5)
    for ln in txt.splitlines():
        if ln[:3] in ("SRR", "ERR", "DRR"):
            return ln.split(",")[0]
    return None


def process(biosample):
    out_tsv = f"{WORK}/{biosample}.amrfinder.tsv"
    if os.path.exists(out_tsv):
        print(f"  {biosample}: already done"); return "cached"
    srr = entrez_srr(biosample)
    if not srr:
        print(f"  {biosample}: NO SRR"); return "no_srr"
    d = f"/kaggle/tmp/{biosample}"
    os.makedirs(d, exist_ok=True)
    env = f"{ENVP}/bin"
    try:
        sh(f"{env}/prefetch -O {d} {srr}", timeout=1800)
        # --split-3: paired -> {srr}_1/_2.fastq; SINGLE-end -> {srr}.fastq (NO _1 suffix).
        # The old --split-files + hardcoded {srr}_1.fastq silently missed single-end runs (fasterq
        # writes {srr}.fastq for those) -> skesa got a nonexistent file -> CalledProcessError. The
        # AR-Bank E. faecium S-side block (SAMN15040xxx) is largely single-end, so this gutted the S class.
        sh(f"{env}/fasterq-dump --split-3 -O {d} {d}/{srr}/{srr}.sra", timeout=1800)
        r1, r2, single = f"{d}/{srr}_1.fastq", f"{d}/{srr}_2.fastq", f"{d}/{srr}.fastq"
        if os.path.exists(r1) and os.path.exists(r2):
            reads = f"{r1},{r2}"           # paired
        elif os.path.exists(single):
            reads = single                 # single-end ({srr}.fastq)
        elif os.path.exists(r1):
            reads = r1                      # unpaired _1 fallback
        else:
            raise FileNotFoundError(f"fasterq-dump produced no FASTQ for {srr} in {d}")
        sh(f"{env}/skesa --reads {reads} --cores 4 --memory 12 --contigs_out {d}/contigs.fasta",
           timeout=2400)
        sh(f"{env}/amrfinder -n {d}/contigs.fasta -O Neisseria_gonorrhoeae -o {out_tsv}", timeout=900)
        # free disk: drop reads/SRA immediately (Kaggle /kaggle/tmp is small)
        sh(f"rm -rf {d}", check=False)
        return "ok"
    except Exception as e:  # noqa: BLE001
        print(f"  {biosample}: FAILED {type(e).__name__}: {str(e)[:120]}")
        sh(f"rm -rf {d}", check=False)
        return f"fail:{type(e).__name__}"


def main():
    setup_bioconda()
    manifest = {}
    for i, bs in enumerate(BIOSAMPLES, 1):
        print(f"\n[{i}/{len(BIOSAMPLES)}] {bs}", flush=True)
        manifest[bs] = process(bs)
        with open(f"{WORK}/assembly_manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2)
    print("\nDONE:", json.dumps(manifest, indent=2))
    ok = sum(1 for v in manifest.values() if v in ("ok", "cached"))
    print(f"{ok}/{len(BIOSAMPLES)} isolates assembled + AMRFinder-scanned")


main()
