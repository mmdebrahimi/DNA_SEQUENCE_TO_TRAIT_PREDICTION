"""Klebsiella pneumoniae cipro — cross-ORGANISM transferability test (roadmap Phase 3, slice 1).

The key open platform question: does the deterministic mechanism-feature AMR caller — built + validated
on E. coli — TRANSFER to a different organism with NO re-tuning? Klebsiella is the right first test
(gram-negative, similar plasmid/QRDR biology → roadmap rates transfer MEDIUM). cipro is the cleanest drug
(QRDR point-mutations in gyrA/parC; the cipro DRUG_RULE = threshold 2, organism-agnostic QUINOLONE classes).

The ONLY adaptation vs E. coli: AMRFinder runs with `-O Klebsiella_pneumoniae` (QRDR point-mutation
detection is organism-specific in AMRFinder's curated set; the wrong -O silently misses Klebsiella QRDR).
The dna-amr rule itself is applied UNCHANGED — that's the transferability test.

Labels: NCBI Pathogen Detection AST (independent source, same as the E. coli cross-source check).
Reports dna-amr cipro rule vs the Klebsiella labels AND vs naive AMRFinder ("any quinolone determinant").

Verdict: TRANSFERS if dna-amr acc >= 0.80 AND sens >= 0.80 (matching the E. coli operating range) — the
deterministic method generalizes across organisms with only the -O swap. Otherwise documents the
organism-specific failure mode (per the Phase 3 falsifier).

Restartable (skips downloaded genomes + existing AMRFinder runs). AMRFinder ~95s/strain.
    uv run python scripts/klebsiella_cipro_transfer.py
    uv run python scripts/klebsiella_cipro_transfer.py --eval-only
"""
from __future__ import annotations

import argparse
import json
from datetime import date as _date
from pathlib import Path

from dna_decode.data import refseq
from dna_decode.eval.amr_rules import (
    AMRFINDER_IMAGE_PINNED, call_resistance, cipro_determinants_from_main, discordance_bucket,
)

DRUG = "ciprofloxacin"
ORGANISM_AMRFINDER = "Klebsiella_pneumoniae"
KDIR = Path("data/raw/klebsiella_cipro")
SELECTED = KDIR / "selected.tsv"
GENOME_CACHE = KDIR / "refseq"
RUNS_ROOT = KDIR / "amrfinder_runs"          # separate from E. coli runs (different -O organism)
PDG = "PDG000000012.2431"
SOURCE_URL = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Klebsiella/latest_snps/Metadata/"
              "PDG000000012.2431.metadata.tsv")


def load_selected() -> dict[str, int]:
    out = {}
    for ln in SELECTED.read_text().splitlines():
        if not ln.strip():
            continue
        acc, rs = ln.split("\t")
        out[acc] = 1 if rs.strip() == "R" else 0
    return out


def ensure_genome_and_run(acc: str) -> Path | None:
    main_tsv = RUNS_ROOT / acc / "main.tsv"
    if main_tsv.exists() and main_tsv.stat().st_size > 0:
        return main_tsv
    try:
        refseq.download_genome(acc, GENOME_CACHE)
        fasta = refseq.fasta_path(acc, GENOME_CACHE)
    except Exception as e:
        print(f"  [{acc}] download FAILED: {type(e).__name__}: {e}")
        return None
    if not Path(fasta).exists():
        print(f"  [{acc}] no FASTA after download")
        return None
    try:
        import scripts.drug_mechanism_audit as dma
        from scripts.drug_mechanism_audit import _run_amrfinder
        dma.AMRFINDER_DB = str(Path("data/amrfinder_db").resolve())   # Docker-readable C: DB
        out_dir = RUNS_ROOT / acc
        out_dir.mkdir(parents=True, exist_ok=True)
        _run_amrfinder(Path(fasta), out_dir, organism=ORGANISM_AMRFINDER)
    except Exception as e:
        print(f"  [{acc}] AMRFinder FAILED: {type(e).__name__}: {e}")
        return None
    return main_tsv if main_tsv.exists() else None


def _confusion(pairs):
    tp = sum(1 for p, y in pairs if p == 1 and y == 1)
    fp = sum(1 for p, y in pairs if p == 1 and y == 0)
    tn = sum(1 for p, y in pairs if p == 0 and y == 0)
    fn = sum(1 for p, y in pairs if p == 0 and y == 1)
    n = tp + fp + tn + fn
    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None}


def evaluate(selected: dict[str, int]) -> dict:
    dna_pairs, naive_pairs, na = [], [], 0
    disc = {"FN_undetected_mechanism": 0, "FP_determinant_without_phenotype": 0}
    for acc, y in selected.items():
        main_tsv = RUNS_ROOT / acc / "main.tsv"
        if not main_tsv.exists():
            na += 1
            continue
        call = call_resistance(main_tsv, DRUG)
        if call["prediction"] == "INDETERMINATE":
            na += 1
            continue
        dna_r = 1 if call["prediction"] == "R" else 0
        dna_pairs.append((dna_r, y))
        b = discordance_bucket(call["prediction"], y)
        if b:
            disc[b] += 1
        naive_r = 1 if len(cipro_determinants_from_main(main_tsv, DRUG)) >= 1 else 0
        naive_pairs.append((naive_r, y))
    dna = _confusion(dna_pairs)
    naive = _confusion(naive_pairs)
    transfers = bool(dna["accuracy"] is not None and dna["accuracy"] >= 0.80
                     and dna["sensitivity"] is not None and dna["sensitivity"] >= 0.80)
    return {"dna_amr": dna, "naive_amrfinder": naive, "discordance": disc, "na": na,
            "verdict": "TRANSFERS" if transfers else "DOES_NOT_TRANSFER"}


def render_md(res: dict, n: int, n_run: int) -> str:
    d = _date.today().isoformat()
    a, b = res["dna_amr"], res["naive_amrfinder"]
    da = (a["accuracy"] - b["accuracy"]) if (a["accuracy"] is not None and b["accuracy"] is not None) else None
    L = [
        f"# Klebsiella pneumoniae cipro — cross-organism transferability — {d}",
        "",
        "> Roadmap Phase 3, slice 1. Does the E. coli-built deterministic AMR caller transfer to Klebsiella",
        "> with NO re-tuning (only AMRFinder -O Klebsiella_pneumoniae)? cipro DRUG_RULE applied UNCHANGED.",
        "",
        f"- Source: NCBI Pathogen Detection `{PDG}` ({SOURCE_URL})",
        f"- Cohort: {n} K. pneumoniae (balanced 15R/15S), {n_run} with AMRFinder runs",
        f"- AMRFinder: `{AMRFINDER_IMAGE_PINNED}`, `-O {ORGANISM_AMRFINDER}`",
        f"- dna-amr rule: cipro DRUG_RULE (threshold 2, QUINOLONE classes) — unchanged from E. coli",
        "",
        f"## VERDICT: {res['verdict']}",
        "",
        "| caller | N | accuracy | sensitivity | specificity |",
        "|---|---:|---:|---:|---:|",
        f"| **dna-amr (transferred unchanged)** | {a['n']} | **{a['accuracy']}** | {a['sensitivity']} | {a['specificity']} |",
        f"| naive AMRFinder (any quinolone determinant) | {b['n']} | {b['accuracy']} | {b['sensitivity']} | {b['specificity']} |",
        f"| Δ accuracy (dna-amr − naive) | | {'+' if (da or 0)>=0 else ''}{round(da,3) if da is not None else '—'} | | |",
        "",
        "## Discordance (failure modes vs independent Klebsiella labels)",
        f"- FN (R missed — efflux/porin/regulatory/low-level): {res['discordance']['FN_undetected_mechanism']}",
        f"- FP (called R, susceptible — label/expression/borderline): {res['discordance']['FP_determinant_without_phenotype']}",
        "",
        "## Interpretation",
        "- TRANSFERS (acc>=0.80 AND sens>=0.80) ⇒ the deterministic mechanism-feature method generalizes",
        "  across organisms with only the AMRFinder -O swap — the platform is cross-organism, not E.coli-specific.",
        "- Honest scope: 1 organism, 1 drug, N=30, NCBI Pathogen Detection labels (different source/curation",
        "  than BV-BRC, not a controlled different-lab study). First transferability data point, not a benchmark.",
    ]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-only", action="store_true")
    args = ap.parse_args(argv)
    selected = load_selected()
    print(f"Klebsiella cipro cohort: {len(selected)} strains")
    if not args.eval_only:
        for i, acc in enumerate(selected, 1):
            print(f"[{i}/{len(selected)}] {acc} ...")
            ensure_genome_and_run(acc)
    n_run = sum(1 for acc in selected if (RUNS_ROOT / acc / "main.tsv").exists())
    print(f"with AMRFinder runs: {n_run}/{len(selected)}")
    res = evaluate(selected)
    d = _date.today().isoformat()
    Path("wiki").mkdir(exist_ok=True)
    out_md = Path("wiki") / f"klebsiella_cipro_transfer_{d}.md"
    out_json = Path("wiki") / f"klebsiella_cipro_transfer_{d}.json"
    out_md.write_text(render_md(res, len(selected), n_run), encoding="utf-8")
    out_json.write_text(json.dumps({
        "organism": "Klebsiella_pneumoniae", "drug": DRUG, "pdg": PDG, "source_url": SOURCE_URL,
        "amrfinder_image": AMRFINDER_IMAGE_PINNED, "amrfinder_organism": ORGANISM_AMRFINDER,
        "date": d, "n_strains": len(selected), "n_with_run": n_run,
        "accessions": sorted(selected), **res}, indent=2), encoding="utf-8")
    a, b = res["dna_amr"], res["naive_amrfinder"]
    print(f"VERDICT: {res['verdict']}")
    print(f"  dna-amr: acc={a['accuracy']} sens={a['sensitivity']} spec={a['specificity']}")
    print(f"  naive  : acc={b['accuracy']} sens={b['sensitivity']} spec={b['specificity']}")
    print(f"Wrote {out_md} + {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
