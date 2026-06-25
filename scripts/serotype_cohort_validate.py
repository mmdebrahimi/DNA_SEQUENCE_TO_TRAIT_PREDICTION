"""Full-cohort validation for the Salmonella-serovar + pneumococcus-serotype cells vs a MEASURED label.

The GREEN-VALIDATED number: per isolate, call the deterministic caller on its assembly and score concordance
vs the WET-LAB measured label (Quellung serotype / traditional Kauffmann-White serovar) -- NOT vs an
in-silico tool (circularity rail). Mirrors scripts/ktype_cohort_validate.py but uses NATIVE blastn (no
Docker) since serotyping runs from an assembly FASTA directly.

INPUT cohort TSV: `accession<tab>measured_label` (one per line; header optional). `measured_label` is the
wet-lab serotype/serovar. Assemblies are read from `--assembly-dir/<accession>.fna` if present, else fetched
via the NCBI `datasets` CLI if it is on PATH (else the isolate is reported as `assembly_missing`, not hidden).

Checkpointed to a jsonl (restartable; skips done accessions). `--limit N` runs a pilot.

    uv run python scripts/serotype_cohort_validate.py --cell pneumo \
        --cohort-tsv data/pneumo_quellung_cohort.tsv --db-dir data/pneumoserotype_db \
        --assembly-dir D:/dna_decode_cache/pneumo_asm --out wiki/pneumo_serotype_cohort_validation.json
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _serogroup(serotype: str) -> str:
    """'19F' -> '19'; '6B' -> '6'; '4' -> '4'. Leading digits = serogroup."""
    i = 0
    while i < len(serotype) and serotype[i].isdigit():
        i += 1
    return serotype[:i] or serotype


def _fetch_assembly(acc: str, dest: Path) -> bool:
    """Use the NCBI `datasets` CLI if available; return True if dest now exists."""
    if dest.exists():
        return True
    if not shutil.which("datasets"):
        return False
    import tempfile
    import zipfile
    with tempfile.TemporaryDirectory() as td:
        zp = Path(td) / "g.zip"
        try:
            subprocess.run(["datasets", "download", "genome", "accession", acc,
                            "--include", "genome", "--filename", str(zp)],
                           check=True, capture_output=True, timeout=600)
            with zipfile.ZipFile(zp) as z:
                fna = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
                if not fna:
                    return False
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(z.read(fna[0]))
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError, zipfile.BadZipFile):
            return False


def _call(cell: str, fasta: Path, db_dir: Path) -> dict:
    if cell == "pneumo":
        from dna_decode.pneumoserotype.runner import call_pneumo_serotype
        r = call_pneumo_serotype(fasta, db_dir)
        return {"predicted": r.get("serotype"), "status": r["status"]}
    from dna_decode.salmserovar.runner import call_serovar
    r = call_serovar(fasta, db_dir)
    return {"predicted": r.get("serovar"), "formula": r.get("antigenic_formula"), "status": r["status"]}


def _concordant(cell: str, predicted: str | None, measured: str) -> dict:
    if not predicted:
        return {"exact": None, "serogroup": None}
    measured = measured.strip()
    exact = predicted.strip().lower() == measured.lower()
    if cell == "pneumo":
        return {"exact": exact, "serogroup": _serogroup(predicted) == _serogroup(measured)}
    return {"exact": exact, "serogroup": None}


def load_cohort(tsv: Path) -> list[dict]:
    out = []
    for line in tsv.read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 2 or parts[0].lower() in ("accession", "id"):
            continue
        out.append({"accession": parts[0], "measured": parts[1]})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", choices=["pneumo", "salm"], required=True)
    ap.add_argument("--cohort-tsv", type=Path, required=True, help="accession<tab>measured_label")
    ap.add_argument("--db-dir", type=Path, required=True)
    ap.add_argument("--assembly-dir", type=Path, required=True, help="<accession>.fna cache (fetched if missing + datasets CLI present)")
    ap.add_argument("--results-jsonl", type=Path, default=None, help="checkpoint jsonl (default <assembly-dir>/results.jsonl)")
    ap.add_argument("--out", type=Path, required=True, help="aggregate result JSON")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)

    results = a.results_jsonl or (a.assembly_dir / f"{a.cell}_results.jsonl")
    results.parent.mkdir(parents=True, exist_ok=True)
    done = ({json.loads(l)["accession"] for l in results.read_text().splitlines() if l.strip()}
            if results.exists() else set())
    cohort = load_cohort(a.cohort_tsv)
    todo = [c for c in cohort if c["accession"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"[serotype-cohort:{a.cell}] {len(cohort)} isolates | {len(done)} done | running {len(todo)}")

    for n, iso in enumerate(todo, 1):
        acc = iso["accession"]
        fasta = a.assembly_dir / f"{acc}.fna"
        if not _fetch_assembly(acc, fasta):
            rec = {"accession": acc, "measured": iso["measured"], "status": "assembly_missing",
                   "predicted": None, "exact": None, "serogroup": None}
        else:
            c = _call(a.cell, fasta, a.db_dir)
            conc = _concordant(a.cell, c.get("predicted"), iso["measured"])
            rec = {"accession": acc, "measured": iso["measured"], "status": c["status"],
                   "predicted": c.get("predicted"), "formula": c.get("formula"), **conc}
        with results.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"  [{n}/{len(todo)}] {acc}: {rec.get('status')} pred={rec.get('predicted')} "
              f"meas={iso['measured']} exact={rec.get('exact')}")

    recs = [json.loads(l) for l in results.read_text().splitlines() if l.strip()]
    scored = [r for r in recs if r.get("exact") is not None]
    ex = sum(1 for r in scored if r["exact"])
    sg = [r for r in scored if r.get("serogroup") is not None]
    sgc = sum(1 for r in sg if r["serogroup"])
    agg = {"n_total": len(recs), "n_scored": len(scored),
           "exact_concordance": round(ex / len(scored), 3) if scored else None,
           "serogroup_concordance": round(sgc / len(sg), 3) if sg else None,
           "n_assembly_missing": sum(1 for r in recs if r.get("status") == "assembly_missing")}
    a.out.write_text(json.dumps({
        "schema": "serotype-cohort-validation-v1", "cell": a.cell,
        "label": ("wet-lab MEASURED serotype/serovar (NOT an in-silico tool -- circularity rail)"),
        "aggregate": agg,
        "caveat": ("Faithful-to-tool caller scored vs the wet-lab measured label. pneumo: exact + serogroup "
                   "reported separately (v0 within-serogroup ceiling). salm: serovar only when the formula "
                   "resolves uniquely."),
    }, indent=2), encoding="utf-8")
    print(f"[serotype-cohort:{a.cell}] AGG: {agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
