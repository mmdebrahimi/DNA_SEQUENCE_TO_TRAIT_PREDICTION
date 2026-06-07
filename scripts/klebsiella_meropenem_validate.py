"""Klebsiella pneumoniae meropenem — 2nd-organism, NEW mechanism class (carbapenem). Phase 3, slice 2.

Meropenem is the highest-value Klebsiella drug: carbapenem resistance (KPC/NDM/OXA-48 carbapenemases) is
the defining K. pneumoniae clinical threat and a mechanism class E. coli AMR never covered. The deterministic
rule: acquired carbapenemase determinant (AMRFinder Subclass CARBAPENEM) ≥1 → R. Same shape as cef
(acquired bla + Subclass refinement); excludes ESBL/AmpC that raise meropenem MIC without hydrolyzing it.

Labels: NCBI Pathogen Detection AST (independent source). Cohort: balanced meropenem R/S, K. pneumoniae,
zero-tuning. AMRFinder `-O Klebsiella_pneumoniae`. Reuses cached cipro-cohort AMRFinder runs where strains
overlap (an AMRFinder run is drug-agnostic — it detects all determinants).

MVP bar: acc >= 0.80 AND sens >= 0.80 on the independent labels (matches the cipro-transfer bar).

Restartable. AMRFinder ~95s/uncached strain.
    uv run python scripts/klebsiella_meropenem_validate.py
    uv run python scripts/klebsiella_meropenem_validate.py --eval-only
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from datetime import date as _date
from pathlib import Path

from dna_decode.data import refseq
from dna_decode.eval.amr_rules import (
    AMRFINDER_IMAGE_PINNED, call_resistance, cipro_determinants_from_main, discordance_bucket,
)

DRUG = "meropenem"
ORG = "Klebsiella_pneumoniae"
KDIR = Path("data/raw/klebsiella_meropenem")
SELECTED = KDIR / "selected.tsv"
GENOME_CACHE = Path("data/raw/klebsiella_cipro/refseq")     # shared genome cache (drug-agnostic)
RUNS = KDIR / "amrfinder_runs"
REUSE_RUNS = Path("data/raw/klebsiella_cipro/amrfinder_runs")  # reuse cipro-cohort runs on overlap
PDG = "PDG000000012.2431"
META_URL = (f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Metadata/"
            f"{PDG}.metadata.tsv")
TARGET_PER_CLASS = 15


def _run_dir(acc: str) -> Path | None:
    """Locate an AMRFinder run for acc — prefer the meropenem dir, else reuse the cipro-cohort run."""
    for root in (RUNS, REUSE_RUNS):
        m = root / acc / "main.tsv"
        if m.exists() and m.stat().st_size > 0:
            return root / acc
    return None


def select_cohort() -> dict[str, int]:
    """Stream NCBI metadata → balanced meropenem R/S cohort, preferring already-cached accessions."""
    if SELECTED.exists():
        out = {}
        for ln in SELECTED.read_text().splitlines():
            if ln.strip():
                a, rs = ln.split("\t"); out[a] = 1 if rs == "R" else 0
        return out
    cached = {p.name for p in REUSE_RUNS.iterdir()} if REUSE_RUNS.exists() else set()
    cand_cached, cand_new = [], []   # (acc, label)
    with urllib.request.urlopen(META_URL, timeout=400) as resp:
        header = resp.readline().decode().rstrip("\n").split("\t")
        ai = header.index("asm_acc"); pi = header.index("AST_phenotypes")
        seen = set()
        for raw in resp:
            cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(cells) <= max(ai, pi):
                continue
            acc, ast = cells[ai], cells[pi]
            if not acc.startswith("GCA_") or acc in seen or not ast or ast == "NULL":
                continue
            for kv in ast.split(","):
                if kv == "meropenem=R" or kv == "meropenem=S":
                    lab = 1 if kv.endswith("=R") else 0
                    (cand_cached if acc in cached else cand_new).append((acc, lab))
                    seen.add(acc)
                    break
            # enough to fill both classes preferring cached
            if sum(1 for _, l in cand_cached if l == 1) >= TARGET_PER_CLASS and \
               sum(1 for _, l in cand_cached if l == 0) >= TARGET_PER_CLASS:
                break
    chosen = {}
    for pool in (cand_cached, cand_new):           # cached first (reuse runs)
        for lab in (1, 0):
            have = sum(1 for v in chosen.values() if v == lab)
            for acc, l in pool:
                if l == lab and acc not in chosen and have < TARGET_PER_CLASS:
                    chosen[acc] = lab; have += 1
    KDIR.mkdir(parents=True, exist_ok=True)
    SELECTED.write_text("".join(f"{a}\t{'R' if v else 'S'}\n" for a, v in chosen.items()), encoding="utf-8")
    return chosen


def ensure_run(acc: str) -> Path | None:
    rd = _run_dir(acc)
    if rd:
        return rd / "main.tsv"
    try:
        refseq.download_genome(acc, GENOME_CACHE)
        fasta = refseq.fasta_path(acc, GENOME_CACHE)
    except Exception as e:
        print(f"  [{acc}] download FAILED: {type(e).__name__}: {e}"); return None
    if not Path(fasta).exists():
        return None
    try:
        import scripts.drug_mechanism_audit as dma
        from scripts.drug_mechanism_audit import _run_amrfinder
        dma.AMRFINDER_DB = str(Path("data/amrfinder_db").resolve())
        out_dir = RUNS / acc; out_dir.mkdir(parents=True, exist_ok=True)
        _run_amrfinder(Path(fasta), out_dir, organism=ORG)
    except Exception as e:
        print(f"  [{acc}] AMRFinder FAILED: {type(e).__name__}: {e}"); return None
    m = RUNS / acc / "main.tsv"
    return m if m.exists() else None


def _conf(pairs):
    tp = sum(1 for p, y in pairs if p == 1 and y == 1); fp = sum(1 for p, y in pairs if p == 1 and y == 0)
    tn = sum(1 for p, y in pairs if p == 0 and y == 0); fn = sum(1 for p, y in pairs if p == 0 and y == 1)
    n = tp + fp + tn + fn
    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None}


def evaluate(sel):
    dna, naive, na = [], [], 0
    disc = {"FN_undetected_mechanism": 0, "FP_determinant_without_phenotype": 0}
    for acc, y in sel.items():
        rd = _run_dir(acc)
        if not rd:
            na += 1; continue
        c = call_resistance(rd / "main.tsv", DRUG)
        if c["prediction"] == "INDETERMINATE":
            na += 1; continue
        dna.append((1 if c["prediction"] == "R" else 0, y))
        b = discordance_bucket(c["prediction"], y)
        if b:
            disc[b] += 1
        naive.append((1 if len(cipro_determinants_from_main(rd / "main.tsv", DRUG)) >= 1 else 0, y))
    d, nv = _conf(dna), _conf(naive)
    transfers = bool(d["accuracy"] and d["accuracy"] >= 0.80 and d["sensitivity"] and d["sensitivity"] >= 0.80)
    return {"dna_amr": d, "naive_amrfinder": nv, "discordance": disc, "na": na,
            "verdict": "VALIDATED" if transfers else "FAILS_BAR"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-only", action="store_true")
    args = ap.parse_args(argv)
    sel = select_cohort()
    print(f"Klebsiella meropenem cohort: {len(sel)} ({sum(sel.values())}R/{len(sel)-sum(sel.values())}S)")
    if not args.eval_only:
        for i, acc in enumerate(sel, 1):
            print(f"[{i}/{len(sel)}] {acc} ...")
            ensure_run(acc)
    n_run = sum(1 for a in sel if _run_dir(a))
    print(f"with AMRFinder runs: {n_run}/{len(sel)}")
    res = evaluate(sel)
    d = _date.today().isoformat()
    a, nv = res["dna_amr"], res["naive_amrfinder"]
    md = [
        f"# Klebsiella pneumoniae meropenem — 2nd-organism, carbapenem mechanism — {d}",
        "", "> Phase 3 slice 2. Carbapenem (KPC/NDM/OXA-48) — a mechanism class E. coli AMR never covered.",
        f"- Source: NCBI Pathogen Detection `{PDG}`; cohort {len(sel)} K. pneumoniae ({sum(sel.values())}R/{len(sel)-sum(sel.values())}S), {n_run} with runs",
        f"- AMRFinder `{AMRFINDER_IMAGE_PINNED}` `-O {ORG}`; rule: carbapenemase (CARBAPENEM-subclass) >=1",
        "", f"## VERDICT: {res['verdict']}", "",
        "| caller | N | acc | sens | spec |", "|---|---:|---:|---:|---:|",
        f"| **dna-amr (CARBAPENEM-subclass)** | {a['n']} | **{a['accuracy']}** | {a['sensitivity']} | {a['specificity']} |",
        f"| naive AMRFinder (any beta-lactam/carbapenem determinant) | {nv['n']} | {nv['accuracy']} | {nv['sensitivity']} | {nv['specificity']} |",
        "", "## Discordance",
        f"- FN (R missed — porin-loss/ESBL+impermeability/low-level): {res['discordance']['FN_undetected_mechanism']}",
        f"- FP (called R, susceptible): {res['discordance']['FP_determinant_without_phenotype']}",
        "", "## Honest scope",
        "1 organism, 1 drug, NCBI labels (different source/curation, not a different-lab study). The rule is",
        "blind to porin-loss-mediated carbapenem resistance (no carbapenemase gene) — expected FN mode.",
    ]
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"klebsiella_meropenem_validate_{d}.md").write_text("\n".join(md), encoding="utf-8")
    (Path("wiki") / f"klebsiella_meropenem_validate_{d}.json").write_text(json.dumps(
        {"organism": ORG, "drug": DRUG, "pdg": PDG, "date": d, "n": len(sel), "n_run": n_run,
         "accessions": sorted(sel), **res}, indent=2), encoding="utf-8")
    print(f"VERDICT: {res['verdict']}  dna-amr acc={a['accuracy']} sens={a['sensitivity']} spec={a['specificity']} | naive acc={nv['accuracy']}")
    print(f"Wrote wiki/klebsiella_meropenem_validate_{d}.{{md,json}}")
    return 0 if res["verdict"] == "VALIDATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
