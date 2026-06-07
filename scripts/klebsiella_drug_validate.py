"""Klebsiella pneumoniae per-drug validation — drug-agnostic (completes the 2nd-organism drug matrix).

Generalizes the cipro/meropenem Klebsiella validators: `--drug <cef|gent|tet|...>` validates the deployed
dna-amr per-drug rule (DRUG_RULE, applied UNCHANGED from E. coli) on a balanced K. pneumoniae cohort from
NCBI Pathogen Detection. call_resistance dispatches per-drug, so the eval is drug-agnostic.

Reuses ALL cached Klebsiella AMRFinder runs (data/raw/klebsiella_*/amrfinder_runs) — runs are
drug-agnostic (detect all determinants), so strains shared across drug cohorts need no re-run.

MVP bar per drug: acc >= 0.80 AND sens >= 0.80 on the independent labels (VALIDATED), else FAILS_BAR with
the determinant breakdown so the organism-specific failure mode is documented (the Phase-3 falsifier).

Restartable. AMRFinder ~95s/uncached strain, `-O Klebsiella_pneumoniae`.
    uv run python scripts/klebsiella_drug_validate.py --drug ceftriaxone
    uv run python scripts/klebsiella_drug_validate.py --drug ceftriaxone --eval-only
"""
from __future__ import annotations

import argparse
import glob
import json
import urllib.request
from datetime import date as _date
from pathlib import Path

from dna_decode.data import refseq
from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.eval.amr_rules import (
    AMRFINDER_IMAGE_PINNED, call_resistance, cipro_determinants_from_main, discordance_bucket,
)

ORG = "Klebsiella_pneumoniae"
GENOME_CACHE = Path("data/raw/klebsiella_cipro/refseq")     # shared genome cache (drug-agnostic)
PDG = "PDG000000012.2431"
META_URL = (f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Metadata/"
            f"{PDG}.metadata.tsv")
TARGET_PER_CLASS = 15
REUSE_GLOB = "data/raw/klebsiella_*/amrfinder_runs"


def _reuse_roots() -> list[Path]:
    return [Path(p) for p in glob.glob(REUSE_GLOB)]


def _run_dir(acc: str, own: Path) -> Path | None:
    for root in [own, *_reuse_roots()]:
        m = root / acc / "main.tsv"
        if m.exists() and m.stat().st_size > 0:
            return root / acc
    return None


def select_cohort(drug: str, kdir: Path, selected: Path) -> dict[str, int]:
    if selected.exists():
        out = {}
        for ln in selected.read_text().splitlines():
            if ln.strip():
                a, rs = ln.split("\t"); out[a] = 1 if rs == "R" else 0
        return out
    cached = set()
    for root in _reuse_roots():
        cached |= {p.name for p in root.iterdir()} if root.exists() else set()
    cand_cached, cand_new = [], []
    tok_r, tok_s = f"{drug}=R", f"{drug}=S"
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
                if kv == tok_r or kv == tok_s:
                    (cand_cached if acc in cached else cand_new).append((acc, 1 if kv == tok_r else 0))
                    seen.add(acc); break
            if (sum(1 for _, l in cand_cached if l == 1) >= TARGET_PER_CLASS and
                    sum(1 for _, l in cand_cached if l == 0) >= TARGET_PER_CLASS):
                break
    chosen = {}
    for pool in (cand_cached, cand_new):
        for lab in (1, 0):
            have = sum(1 for v in chosen.values() if v == lab)
            for acc, l in pool:
                if l == lab and acc not in chosen and have < TARGET_PER_CLASS:
                    chosen[acc] = lab; have += 1
    kdir.mkdir(parents=True, exist_ok=True)
    selected.write_text("".join(f"{a}\t{'R' if v else 'S'}\n" for a, v in chosen.items()), encoding="utf-8")
    return chosen


def ensure_run(acc: str, own_runs: Path) -> Path | None:
    rd = _run_dir(acc, own_runs)
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
        out_dir = own_runs / acc; out_dir.mkdir(parents=True, exist_ok=True)
        _run_amrfinder(Path(fasta), out_dir, organism=ORG)
    except Exception as e:
        print(f"  [{acc}] AMRFinder FAILED: {type(e).__name__}: {e}"); return None
    m = own_runs / acc / "main.tsv"
    return m if m.exists() else None


def _conf(pairs):
    tp = sum(1 for p, y in pairs if p == 1 and y == 1); fp = sum(1 for p, y in pairs if p == 1 and y == 0)
    tn = sum(1 for p, y in pairs if p == 0 and y == 0); fn = sum(1 for p, y in pairs if p == 0 and y == 1)
    n = tp + fp + tn + fn
    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None}


def evaluate(sel, drug, own_runs):
    dna, naive, na = [], [], 0
    disc = {"FN_undetected_mechanism": 0, "FP_determinant_without_phenotype": 0}
    for acc, y in sel.items():
        rd = _run_dir(acc, own_runs)
        if not rd:
            na += 1; continue
        c = call_resistance(rd / "main.tsv", drug)
        if c["prediction"] == "INDETERMINATE":
            na += 1; continue
        dna.append((1 if c["prediction"] == "R" else 0, y))
        b = discordance_bucket(c["prediction"], y)
        if b:
            disc[b] += 1
        naive.append((1 if len(cipro_determinants_from_main(rd / "main.tsv", drug)) >= 1 else 0, y))
    d, nv = _conf(dna), _conf(naive)
    ok = bool(d["accuracy"] and d["accuracy"] >= 0.80 and d["sensitivity"] and d["sensitivity"] >= 0.80)
    return {"dna_amr": d, "naive_amrfinder": nv, "discordance": disc, "na": na,
            "verdict": "VALIDATED" if ok else "FAILS_BAR"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", required=True, choices=supported_drugs())
    ap.add_argument("--eval-only", action="store_true")
    args = ap.parse_args(argv)
    drug = args.drug
    kdir = Path(f"data/raw/klebsiella_{drug}")
    own_runs = kdir / "amrfinder_runs"
    selected = kdir / "selected.tsv"
    sel = select_cohort(drug, kdir, selected)
    print(f"Klebsiella {drug} cohort: {len(sel)} ({sum(sel.values())}R/{len(sel)-sum(sel.values())}S)")
    if not args.eval_only:
        for i, acc in enumerate(sel, 1):
            print(f"[{i}/{len(sel)}] {acc} ...")
            ensure_run(acc, own_runs)
    n_run = sum(1 for a in sel if _run_dir(a, own_runs))
    print(f"with AMRFinder runs: {n_run}/{len(sel)}")
    res = evaluate(sel, drug, own_runs)
    d = _date.today().isoformat()
    a, nv = res["dna_amr"], res["naive_amrfinder"]
    md = [
        f"# Klebsiella pneumoniae {drug} — 2nd-organism drug-matrix validation — {d}",
        "", f"> Phase 3 matrix. Deployed dna-amr {drug} rule applied UNCHANGED from E. coli.",
        f"- Source: NCBI Pathogen Detection `{PDG}`; cohort {len(sel)} K. pneumoniae "
        f"({sum(sel.values())}R/{len(sel)-sum(sel.values())}S), {n_run} with runs; `-O {ORG}`; "
        f"AMRFinder `{AMRFINDER_IMAGE_PINNED}`",
        "", f"## VERDICT: {res['verdict']}", "",
        "| caller | N | acc | sens | spec |", "|---|---:|---:|---:|---:|",
        f"| **dna-amr ({drug} rule, unchanged)** | {a['n']} | **{a['accuracy']}** | {a['sensitivity']} | {a['specificity']} |",
        f"| naive AMRFinder (any drug-class determinant) | {nv['n']} | {nv['accuracy']} | {nv['sensitivity']} | {nv['specificity']} |",
        "", "## Discordance",
        f"- FN (R missed): {res['discordance']['FN_undetected_mechanism']}",
        f"- FP (called R, susceptible): {res['discordance']['FP_determinant_without_phenotype']}",
        "", "## Honest scope",
        f"1 organism, 1 drug, N={len(sel)}, NCBI labels (different source/curation, not a different-lab study).",
    ]
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"klebsiella_{drug}_validate_{d}.md").write_text("\n".join(md), encoding="utf-8")
    (Path("wiki") / f"klebsiella_{drug}_validate_{d}.json").write_text(json.dumps(
        {"organism": ORG, "drug": drug, "pdg": PDG, "date": d, "n": len(sel), "n_run": n_run,
         "accessions": sorted(sel), **res}, indent=2), encoding="utf-8")
    print(f"VERDICT: {res['verdict']}  dna-amr acc={a['accuracy']} sens={a['sensitivity']} spec={a['specificity']} | naive acc={nv['accuracy']}")
    print(f"Wrote wiki/klebsiella_{drug}_validate_{d}.{{md,json}}")
    return 0 if res["verdict"] == "VALIDATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
