"""Quick pre-flight: check how many of the N=40 cipro cohort strains have
their genome.fna + annotations.gff3 already downloaded under the RefSeq cache.
Reports the count + a list of missing accessions so we know whether to download
before populate.
"""
from __future__ import annotations

import sys
from pathlib import Path

from dna_decode.data.cohort import load_cohort
from dna_decode.data.refseq import fasta_path, gff_path

REFSEQ_CACHE = Path("D:/dna_decode_cache/refseq")
COHORT = Path("data/processed/gate_b_n40_cipro_cohort.parquet")


def main() -> int:
    cohort = load_cohort(COHORT)
    print(f"[check] cohort: {len(cohort.strains)} strains")
    missing_fna: list[str] = []
    missing_gff: list[str] = []
    ok_count = 0
    for s in cohort.strains:
        acc = getattr(s, "assembly_accession", None)
        if not acc:
            missing_fna.append(f"{s.strain_id} (no assembly_accession)")
            continue
        fna_ok = fasta_path(acc, REFSEQ_CACHE).exists()
        gff_ok = gff_path(acc, REFSEQ_CACHE).exists()
        if fna_ok and gff_ok:
            ok_count += 1
        else:
            if not fna_ok:
                missing_fna.append(f"{s.strain_id} ({acc})")
            if not gff_ok:
                missing_gff.append(f"{s.strain_id} ({acc})")
    print(f"[check] ready: {ok_count} / {len(cohort.strains)}")
    print(f"[check] missing genome.fna: {len(missing_fna)}")
    for m in missing_fna[:20]:
        print(f"  - {m}")
    print(f"[check] missing annotations.gff3: {len(missing_gff)}")
    for m in missing_gff[:20]:
        print(f"  - {m}")
    return 0 if not (missing_fna or missing_gff) else 1


if __name__ == "__main__":
    sys.exit(main())
