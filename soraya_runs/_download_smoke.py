import sys, csv
sys.path.insert(0, r"C:/Users/Farshad/PythonProjects/dna_decode")
from dna_decode.data.refseq import download_genome
DNA=r"C:/Users/Farshad/PythonProjects/dna_decode"
cache="D:/dna_decode_cache/refseq"
rows=list(csv.DictReader(open(f"{DNA}/research_outputs/pathotype_smoke_lt_vs_st_cohort_2026-05-31.csv",encoding='utf-8')))
for r in rows:
    acc=r['gca_accession']
    try:
        p=download_genome(acc,cache)
        import os
        ok=(os.path.exists(os.path.join(p,'genome.fna')))
        print(f"OK {acc} ({r['binary_label']}) fna={ok}")
    except Exception as e:
        print(f"FAIL {acc}: {type(e).__name__} {str(e)[:120]}")
print("downloads done")
