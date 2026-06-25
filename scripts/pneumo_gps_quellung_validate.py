"""GREEN-VALIDATED number for dna-pneumo-serotype vs PHENOTYPIC (Quellung) serotype — the GPS Poland cohort.

This is the real independent-measured-label validation (the HIV-style win, for a typing trait): per isolate,
fetch the GPS-deposited ENA ASSEMBLY (contig.fa.gz), call the serotype with the deterministic cps-reference
caller (native blastn), and score concordance vs the WET-LAB phenotypic serotype from the GPS pipeline-paper
supplement (Nat Commun 2025, Supplementary Data 1; column `Phenotypic_serotype`, method `QUELLUNG`/unspecified).

Clears the circularity rail: the label is a wet-lab serology MEASUREMENT, independent of the genome (NOT the
in-silico Monocle field). Assemblies are GPS-deposited (independent of OUR caller). Native blastn (no Docker).

Cohort TSV: `ERR<tab>ERS<tab>phenotypic_serotype` (built from Supplementary Data 1). Assemblies resolved via
ENA: ERS -> analysis ERZ -> contig.fa.gz. Checkpointed jsonl; contigs deleted per-isolate (bound disk).

    BLASTN_BIN=C:/Users/Farshad/ncbi-blast/bin/blastn.exe uv run python scripts/pneumo_gps_quellung_validate.py \
        --cohort-tsv D:/dna_decode_cache/pneumo_gps/poland_quellung_cohort.tsv \
        --db-dir data/pneumoserotype_db --work D:/dna_decode_cache/pneumo_gps --limit 10
"""
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from dna_decode.pneumoserotype.runner import call_pneumo_serotype  # noqa: E402


def _norm(st: str) -> str:
    """'04'->'4', '06B'->'6B', '19F'->'19F' (strip leading zeros from the numeric prefix)."""
    st = st.strip()
    i = 0
    while i < len(st) and st[i].isdigit():
        i += 1
    return (str(int(st[:i])) if i else "") + st[i:]


def _serogroup(st: str) -> str:
    i = 0
    while i < len(st) and st[i].isdigit():
        i += 1
    return st[:i]


def _ena_assembly_ftp(ers: str, timeout: int = 30) -> str | None:
    """ERS -> the contig.fa.gz https URL via the ENA analysis filereport (None if no assembly)."""
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ers}"
           "&result=analysis&fields=analysis_accession,generated_ftp&format=tsv")
    try:
        out = subprocess.run(["curl", "-s", "--max-time", str(timeout), url],
                             capture_output=True, text=True, timeout=timeout + 5).stdout
    except (subprocess.TimeoutExpired, OSError):
        return None
    for line in out.splitlines()[1:]:
        parts = line.split("\t")
        ftp = parts[-1] if parts else ""
        for token in ftp.split(";"):
            if token.endswith((".fa.gz", ".fasta.gz", "contig.fa.gz")):
                return "https://" + token if not token.startswith("http") else token
    return None


def _fetch_assembly(ers: str, dest: Path, timeout: int = 120) -> bool:
    ftp = _ena_assembly_ftp(ers)
    if not ftp:
        return False
    gz = dest.with_suffix(".fa.gz")
    try:
        subprocess.run(["curl", "-sL", "--max-time", str(timeout), ftp, "-o", str(gz)],
                       check=True, capture_output=True, timeout=timeout + 10)
        with gzip.open(gz, "rb") as fi:
            dest.write_bytes(fi.read())
        gz.unlink(missing_ok=True)
        return dest.exists() and dest.stat().st_size > 1000
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError, gzip.BadGzipFile):
        return False


def load_cohort(tsv: Path) -> list[dict]:
    out = []
    for line in tsv.read_text(encoding="utf-8").splitlines():
        p = [x.strip() for x in line.split("\t")]
        if len(p) >= 3 and p[0].startswith(("ERR", "SRR", "DRR")):
            out.append({"err": p[0], "ers": p[1], "measured": p[2]})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-tsv", type=Path, required=True)
    ap.add_argument("--db-dir", type=Path, default=REPO / "data" / "pneumoserotype_db")
    ap.add_argument("--work", type=Path, default=Path("D:/dna_decode_cache/pneumo_gps"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "pneumo_serotype_cohort_validation.json")
    a = ap.parse_args(argv)

    results = a.work / "quellung_results.jsonl"
    results.parent.mkdir(parents=True, exist_ok=True)
    done = ({json.loads(x)["err"] for x in results.read_text().splitlines() if x.strip()}
            if results.exists() else set())
    cohort = load_cohort(a.cohort_tsv)
    todo = [c for c in cohort if c["err"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"[pneumo-quellung] {len(cohort)} isolates | {len(done)} done | running {len(todo)}")

    asm_dir = a.work / "asm"
    asm_dir.mkdir(parents=True, exist_ok=True)
    for n, iso in enumerate(todo, 1):
        fa = asm_dir / f"{iso['ers']}.fna"
        if not _fetch_assembly(iso["ers"], fa):
            rec = {"err": iso["err"], "ers": iso["ers"], "measured": iso["measured"],
                   "status": "assembly_unavailable", "predicted": None, "exact": None, "serogroup": None}
        else:
            r = call_pneumo_serotype(fa, a.db_dir)
            pred = r.get("serotype")
            predn = _norm(pred) if pred else None
            measn = _norm(iso["measured"])
            rec = {"err": iso["err"], "ers": iso["ers"], "measured": iso["measured"], "status": r["status"],
                   "predicted": pred, "predicted_norm": predn,
                   "exact": (predn == measn) if pred else None,
                   "serogroup": (_serogroup(predn) == _serogroup(measn)) if pred else None,
                   "pid": r.get("percent_identity"), "cov": r.get("percent_coverage")}
            fa.unlink(missing_ok=True)
        with results.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"  [{n}/{len(todo)}] {iso['err']}: {rec['status']} pred={rec.get('predicted')} "
              f"meas={iso['measured']} exact={rec.get('exact')}")

    recs = [json.loads(x) for x in results.read_text().splitlines() if x.strip()]
    scored = [r for r in recs if r.get("exact") is not None]
    ex = sum(1 for r in scored if r["exact"])
    sg = sum(1 for r in scored if r.get("serogroup"))
    agg = {"n_total": len(recs), "n_scored": len(scored),
           "exact_concordance": round(ex / len(scored), 3) if scored else None,
           "serogroup_concordance": round(sg / len(scored), 3) if scored else None,
           "n_assembly_unavailable": sum(1 for r in recs if r.get("status") == "assembly_unavailable")}
    a.out.write_text(json.dumps({
        "schema": "pneumo-serotype-cohort-validation-v1",
        "label": "wet-lab PHENOTYPIC serotype (GPS Poland cohort; Quellung/phenotypic; independent of genome)",
        "source": "Nat Commun 2025 GPS pipeline paper, Supplementary Data 1 (Phenotypic_serotype); ENA GPS assemblies",
        "caller": "dna_decode cps-reference blastn v0 (PneumoCaT Stage-1 DB, 95 serotypes)",
        "aggregate": agg,
        "caveat": ("Independent measured-label validation: phenotypic serology label + GPS-deposited assembly, "
                   "both independent of OUR caller. exact + serogroup reported separately (v0 within-serogroup "
                   "ceiling, e.g. 6A/6B). Published GPS-pipeline in-silico-vs-Quellung ceiling ~89%."),
    }, indent=2), encoding="utf-8")
    print(f"[pneumo-quellung] AGG: {agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
