"""Feasibility census for the two DNA-LLM forward levers (read-only; no downloads, no Docker).

After the functional-alphabet probe returned TIES (cipro within-lineage underpowered: 6 shared lineages /
43 pairs, p=0.0565), the two non-foreclosed levers were both DATA questions:
  1. Cohort expansion -- more R+S-sharing MLST lineages to power the within-lineage metric (cipro).
  2. A distributed-mechanism drug (tetracycline) where curated determinants are incomplete.

This census computes the *theoretical ceiling* for each, straight from the raw BV-BRC tables already on
disk, using the RELAXED binary R/S phenotype (not strict-MIC -- the within-lineage metric only needs MLSTs
carrying BOTH R and S with a DOWNLOADABLE assembly). It answers: is either lever even attemptable, and
how many NEW genomes would have to be downloaded + AMRFinder'd to power it?

Output: wiki/dna_llm_lever_feasibility_<date>.{md,json}. Exit 0 always (a census, not a gate).
"""
from __future__ import annotations

import datetime
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
AST_CSV = Path("C:/Users/Farshad/Downloads/BVBRC_genome_amr.csv")
GENOME_CSV = Path("C:/Users/Farshad/Downloads/BVBRC_genome (1).csv")
DRUGS = ["ciprofloxacin", "tetracycline"]
# current cipro probe substrate, for comparison
CURRENT = {"ciprofloxacin": {"shared_lineages": 6, "within_lineage_pairs": 43, "n": 147}}


def _norm(s):
    return str(s).strip().lower()


def census() -> dict:
    if not AST_CSV.exists() or not GENOME_CSV.exists():
        return {"status": "BLOCKED_NO_RAW_DATA", "ast_csv": str(AST_CSV), "genome_csv": str(GENOME_CSV)}

    ast = pd.read_csv(AST_CSV, dtype=str,
                      usecols=["Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype"])
    ast = ast[ast["Genome Name"].fillna("").str.contains("Escherichia coli", case=False)]
    gen = pd.read_csv(GENOME_CSV, dtype=str, usecols=["Genome ID", "MLST", "Assembly Accession"])
    gmap = gen.set_index("Genome ID")[["MLST", "Assembly Accession"]].to_dict("index")

    out = {}
    for drug in DRUGS:
        d = ast[ast["Antibiotic"].fillna("").str.lower() == drug].copy()
        # relaxed binary label straight from the phenotype string
        d["lab"] = d["Resistant Phenotype"].map(
            lambda p: 1 if _norm(p) == "resistant" else (0 if _norm(p) == "susceptible" else None))
        d = d[d["lab"].notna()]
        # one label per genome (drop genomes with conflicting calls)
        by_genome = defaultdict(set)
        for gid, lab in zip(d["Genome ID"], d["lab"]):
            by_genome[gid].add(int(lab))
        clean = {g: next(iter(s)) for g, s in by_genome.items() if len(s) == 1}

        total_R = sum(1 for v in clean.values() if v == 1)
        total_S = sum(1 for v in clean.values() if v == 0)

        # restrict to DOWNLOADABLE (has Assembly Accession) + has MLST
        downloadable = {}
        for g, lab in clean.items():
            meta = gmap.get(g)
            if not meta:
                continue
            acc = "" if pd.isna(meta.get("Assembly Accession")) else str(meta.get("Assembly Accession")).strip()
            mlst = "" if pd.isna(meta.get("MLST")) else str(meta.get("MLST")).strip()
            if acc and mlst:
                downloadable[g] = (lab, mlst)
        dl_R = sum(1 for lab, _ in downloadable.values() if lab == 1)
        dl_S = sum(1 for lab, _ in downloadable.values() if lab == 0)

        # shared-lineage ceiling among downloadable
        by_mlst = defaultdict(lambda: [0, 0])  # mlst -> [S, R]
        for lab, mlst in downloadable.values():
            by_mlst[mlst][lab] += 1
        shared = {m: (s, r) for m, (s, r) in by_mlst.items() if s > 0 and r > 0}
        pairs = sum(s * r for s, r in shared.values())
        strains_in_shared = sum(s + r for s, r in shared.values())

        out[drug] = {
            "ast_R": total_R, "ast_S": total_S,
            "downloadable_R": dl_R, "downloadable_S": dl_S,
            "shared_RS_lineages_ceiling": len(shared),
            "within_lineage_pairs_ceiling": pairs,
            "strains_in_shared_lineages": strains_in_shared,
            "current_probe": CURRENT.get(drug),
        }
    return {"status": "OK", "drugs": out}


def build_manifest(drug: str, n_lineages: int = 20, per_class: int = 3) -> dict:
    """Pick a concrete, balanced shared-lineage cohort for `drug` -> a fetch manifest (accessions to download).

    Greedy: take the shared MLSTs with the largest min(R,S), up to `per_class` R + `per_class` S each, until
    `n_lineages` lineages are covered. This is the cohort a powered within-lineage probe would need fetched."""
    ast = pd.read_csv(AST_CSV, dtype=str, usecols=["Genome ID", "Genome Name", "Antibiotic", "Resistant Phenotype"])
    ast = ast[ast["Genome Name"].fillna("").str.contains("Escherichia coli", case=False)]
    gen = pd.read_csv(GENOME_CSV, dtype=str, usecols=["Genome ID", "MLST", "Assembly Accession"])
    gmap = gen.set_index("Genome ID")[["MLST", "Assembly Accession"]].to_dict("index")
    d = ast[ast["Antibiotic"].fillna("").str.lower() == drug]
    by_genome = defaultdict(set)
    for gid, p in zip(d["Genome ID"], d["Resistant Phenotype"]):
        lab = 1 if _norm(p) == "resistant" else (0 if _norm(p) == "susceptible" else None)
        if lab is not None:
            by_genome[gid].add(lab)
    by_mlst = defaultdict(lambda: {0: [], 1: []})
    for g, s in by_genome.items():
        if len(s) != 1:
            continue
        lab = next(iter(s))
        meta = gmap.get(g) or {}
        acc = "" if pd.isna(meta.get("Assembly Accession")) else str(meta.get("Assembly Accession")).strip()
        mlst = "" if pd.isna(meta.get("MLST")) else str(meta.get("MLST")).strip()
        if acc and mlst:
            by_mlst[mlst][lab].append(acc)
    shared = [(m, v) for m, v in by_mlst.items() if v[0] and v[1]]
    shared.sort(key=lambda kv: min(len(kv[1][0]), len(kv[1][1])), reverse=True)
    picks = []
    for m, v in shared[:n_lineages]:
        for lab in (1, 0):
            for acc in v[lab][:per_class]:
                picks.append({"assembly_accession": acc, "mlst": m, "label": "R" if lab else "S"})
    return {"drug": drug, "n_lineages": min(n_lineages, len(shared)), "per_class": per_class,
            "n_strains": len(picks), "n_R": sum(1 for p in picks if p["label"] == "R"),
            "n_S": sum(1 for p in picks if p["label"] == "S"), "strains": picks}


def main() -> int:
    import sys
    if "--manifest" in sys.argv:
        drug = sys.argv[sys.argv.index("--manifest") + 1]
        man = build_manifest(drug)
        man["analysis_date"] = datetime.date.today().isoformat()
        out = ROOT / "wiki" / f"dna_llm_shared_lineage_manifest_{drug}_{man['analysis_date']}.json"
        out.write_text(json.dumps(man, indent=2), encoding="utf-8")
        print(f"manifest {drug}: {man['n_lineages']} shared lineages, {man['n_strains']} strains "
              f"({man['n_R']}R/{man['n_S']}S) to fetch -> {out.name}")
        return 0
    rep = census()
    rep["analysis_date"] = datetime.date.today().isoformat()
    base = ROOT / "wiki" / f"dna_llm_lever_feasibility_{rep['analysis_date']}"
    Path(f"{base}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    L = [f"# DNA-LLM forward-lever feasibility census ({rep['analysis_date']})", "",
         "_Read-only ceiling from the raw BV-BRC tables on disk; relaxed binary R/S phenotype. The "
         "within-lineage metric needs MLSTs carrying BOTH R and S with a downloadable assembly._", ""]
    if rep["status"] != "OK":
        L.append(f"**STATUS: {rep['status']}** -- raw BV-BRC CSVs not found; cannot census.")
    else:
        L += ["| drug | AST R/S | downloadable R/S | shared R+S lineages (ceiling) | within-lineage pairs (ceiling) | strains in shared | current probe |",
              "|---|---|---|---|---|---|---|"]
        for drug, d in rep["drugs"].items():
            cur = d["current_probe"]
            curs = f"{cur['shared_lineages']} lin / {cur['within_lineage_pairs']} pairs (N={cur['n']})" if cur else "—"
            L.append(f"| {drug} | {d['ast_R']}/{d['ast_S']} | {d['downloadable_R']}/{d['downloadable_S']} | "
                     f"{d['shared_RS_lineages_ceiling']} | {d['within_lineage_pairs_ceiling']} | "
                     f"{d['strains_in_shared_lineages']} | {curs} |")
    L += ["", "_Ceiling = the most a re-selected cohort could achieve IF every downloadable strain in a "
          "shared lineage were fetched + AMRFinder'd. The gap to the current probe is the new-genome "
          "download + AMRFinder cost (Docker), NOT a compute/GPU question._"]
    Path(f"{base}.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"[census -> {base}.{{md,json}}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
