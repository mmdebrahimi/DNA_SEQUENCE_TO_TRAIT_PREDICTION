import sys
sys.path.insert(0, r"C:/Users/Farshad/PythonProjects/dna_decode")
from tools.docker_runner import run, DockerRunnerError
acc = "GCA_003984445.1"
try:
    proc = run(
        image="oschwengers/bakta:v1.11.4",
        args=["--db","/db/db-light","--output","/out","--prefix",acc,
              "--skip-plot","--force","--threads","4","/data/genome.fna"],
        mounts={
            f"D:/dna_decode_cache/refseq/{acc}": "/data:ro",
            "C:/Users/Farshad/dna_decode_stage2/bakta_db": "/db:ro",
            "D:/dna_decode_cache/bakta_out": "/out",
        },
        capture_output=True, check=False, timeout=1800,
    )
    print("EXIT:", proc.returncode)
    print("STDERR tail:", (proc.stderr or "")[-1200:])
except DockerRunnerError as e:
    print("DockerRunnerError:", str(e)[:500])
