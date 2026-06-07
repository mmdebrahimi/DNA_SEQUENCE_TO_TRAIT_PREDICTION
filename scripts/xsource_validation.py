"""Cross-source validation of dna-amr against an INDEPENDENT label source (EP follow-on).

Recommendation #3 from the 2026-06-06 brainstorm: the in-cohort + held-out validations both used BV-BRC
labels. This validates the deterministic AMR caller against NCBI Pathogen Detection AST phenotypes — a
DIFFERENT source/curation than BV-BRC — on a cohort with ZERO accession overlap with any BV-BRC cohort.

It reports, per drug, BOTH:
  - dna-amr per-drug rule (DRUG_RULE: threshold + Subclass refinement) vs the independent label
  - naive AMRFinder ("any drug-class determinant present -> R", threshold=1, no refinement) vs the same label
so the headline answers the product question: does the per-drug policy add value over vanilla AMRFinder
on data it was NOT tuned on? Plus the discordance taxonomy (FN_undetected_mechanism / FP_determinant_
without_phenotype).

Honest scope: NCBI Pathogen Detection aggregates public AST submissions, so this is a DIFFERENT-SOURCE /
different-curation check, NOT a controlled different-lab study. Labels are still broth-microdilution-class
clinical AST. It closes the same-database gap (BV-BRC), not the same-methodology gap.

Restartable: skips genomes already downloaded + accessions already AMRFinder-run. AMRFinder is the long
pole (~95s/strain via Docker). Run end-to-end (download -> AMRFinder -> evaluate -> artifact):

    uv run python scripts/xsource_validation.py
    uv run python scripts/xsource_validation.py --eval-only   # re-evaluate from cached runs, no Docker
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

DRUGS = ["ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline"]
XSRC = Path("data/raw/ncbi_pathogen_xsource")
SELECTED = XSRC / "selected.tsv"
GENOME_CACHE = XSRC / "refseq"
RUNS_ROOT = Path("data/amrfinder_runs")          # shared cache (non-overlap guaranteed by selection)
PDG = "PDG000000004.6094"
SOURCE_URL = ("https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Escherichia_coli_Shigella/"
              "latest_snps/Metadata/PDG000000004.6094.metadata.tsv")


def load_selected() -> dict[str, dict[str, int]]:
    """accession -> {drug: 1(R)/0(S)} for drugs with an R/S call."""
    out: dict[str, dict[str, int]] = {}
    lines = SELECTED.read_text().splitlines()
    header = lines[0].split("\t")
    didx = {d: header.index(d) for d in DRUGS}
    for ln in lines[1:]:
        cells = ln.split("\t")
        acc = cells[0]
        labels = {}
        for d in DRUGS:
            v = cells[didx[d]].strip()
            if v == "R":
                labels[d] = 1
            elif v == "S":
                labels[d] = 0
        if labels:
            out[acc] = labels
    return out


def ensure_genome_and_run(acc: str) -> Path | None:
    """Download genome + run AMRFinder (both restartable). Returns main.tsv path or None on failure."""
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
        # Force the repo-local prepared DB on C: — the home DB's `latest` symlink resolves to D:,
        # which Docker Desktop cannot share (silent empty mount -> "AMRProt.fa.phr not found").
        dma.AMRFINDER_DB = str(Path("data/amrfinder_db").resolve())
        out_dir = RUNS_ROOT / acc
        out_dir.mkdir(parents=True, exist_ok=True)
        _run_amrfinder(Path(fasta), out_dir)
    except Exception as e:
        print(f"  [{acc}] AMRFinder FAILED: {type(e).__name__}: {e}")
        return None
    return main_tsv if main_tsv.exists() else None


def _confusion(pairs: list[tuple[int, int]]) -> dict:
    """pairs = [(pred 0/1, true 0/1)]."""
    tp = sum(1 for p, y in pairs if p == 1 and y == 1)
    fp = sum(1 for p, y in pairs if p == 1 and y == 0)
    tn = sum(1 for p, y in pairs if p == 0 and y == 0)
    fn = sum(1 for p, y in pairs if p == 0 and y == 1)
    n = tp + fp + tn + fn
    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "accuracy": round((tp + tn) / n, 3) if n else None,
            "sensitivity": round(tp / (tp + fn), 3) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 3) if (tn + fp) else None}


def evaluate(selected: dict[str, dict[str, int]]) -> dict:
    per_drug = {}
    for drug in DRUGS:
        dna_pairs, naive_pairs, na = [], [], 0
        disc = {"FN_undetected_mechanism": 0, "FP_determinant_without_phenotype": 0}
        for acc, labels in selected.items():
            if drug not in labels:
                continue
            main_tsv = RUNS_ROOT / acc / "main.tsv"
            if not main_tsv.exists():
                na += 1
                continue
            y = labels[drug]
            call = call_resistance(main_tsv, drug)
            if call["prediction"] == "INDETERMINATE":
                na += 1
                continue
            dna_r = 1 if call["prediction"] == "R" else 0
            dna_pairs.append((dna_r, y))
            b = discordance_bucket(call["prediction"], y)
            if b:
                disc[b] += 1
            # naive: any broad drug-class determinant present (threshold 1, no subclass refinement)
            naive_r = 1 if len(cipro_determinants_from_main(main_tsv, drug)) >= 1 else 0
            naive_pairs.append((naive_r, y))
        per_drug[drug] = {
            "dna_amr": _confusion(dna_pairs),
            "naive_amrfinder": _confusion(naive_pairs),
            "discordance": disc,
            "na": na,
        }
    return per_drug


def render_md(per_drug: dict, n_strains: int, n_with_run: int) -> str:
    d = _date.today().isoformat()
    L = [
        f"# dna-amr cross-source validation — NCBI Pathogen Detection — {d}",
        "",
        "> Recommendation #3 (2026-06-06 brainstorm): validate the deterministic AMR caller against an",
        "> INDEPENDENT label source (NCBI Pathogen Detection), zero accession overlap with BV-BRC cohorts.",
        "> Headline question: does the per-drug rule add value over vanilla AMRFinder on UN-tuned data?",
        "",
        f"- Source: NCBI Pathogen Detection `{PDG}` ({SOURCE_URL})",
        f"- Cohort: {n_strains} E. coli, {n_with_run} with AMRFinder runs; non-overlap with BV-BRC enforced by construction",
        f"- AMRFinder image (pinned): `{AMRFINDER_IMAGE_PINNED}`",
        "",
        "## dna-amr per-drug rule vs naive AMRFinder (independent labels)",
        "",
        "| Drug | N | dna-amr acc | sens | spec | naive acc | sens | spec | Δacc |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for drug in DRUGS:
        a = per_drug[drug]["dna_amr"]
        b = per_drug[drug]["naive_amrfinder"]
        da = (a["accuracy"] - b["accuracy"]) if (a["accuracy"] is not None and b["accuracy"] is not None) else None
        L.append(f"| {drug} | {a['n']} | **{a['accuracy']}** | {a['sensitivity']} | {a['specificity']} | "
                 f"{b['accuracy']} | {b['sensitivity']} | {b['specificity']} | "
                 f"{'+' if (da or 0) >= 0 else ''}{round(da,3) if da is not None else '—'} |")
    L += ["", "## Discordance taxonomy (dna-amr failure modes vs independent labels)", "",
          "| Drug | FN (R missed: efflux/porin/regulatory/low-level) | FP (called R, susceptible: label/expression/borderline) |",
          "|---|---:|---:|"]
    for drug in DRUGS:
        disc = per_drug[drug]["discordance"]
        L.append(f"| {drug} | {disc['FN_undetected_mechanism']} | {disc['FP_determinant_without_phenotype']} |")
    L += ["", "## Interpretation", "",
          "- Δacc > 0 means the per-drug policy (threshold + Subclass refinement) beats vanilla "
          "\"any drug-class determinant -> R\" on data the rule was NOT tuned on — the product-value headline.",
          "- Honest scope: NCBI Pathogen Detection aggregates public AST submissions (different source/curation "
          "than BV-BRC, NOT a controlled different-lab study). Closes the same-database gap, not the same-methodology gap.",
          "- Non-overlap with all BV-BRC cohorts is guaranteed by construction (selection excluded the 176 BV-BRC accessions)."]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-only", action="store_true", help="skip download+AMRFinder; evaluate cached runs")
    args = ap.parse_args(argv)

    selected = load_selected()
    print(f"selected cohort: {len(selected)} strains")

    if not args.eval_only:
        for i, acc in enumerate(selected, 1):
            print(f"[{i}/{len(selected)}] {acc} ...")
            ensure_genome_and_run(acc)

    n_with_run = sum(1 for acc in selected if (RUNS_ROOT / acc / "main.tsv").exists())
    print(f"strains with AMRFinder runs: {n_with_run}/{len(selected)}")
    per_drug = evaluate(selected)

    d = _date.today().isoformat()
    XSRC.mkdir(parents=True, exist_ok=True)
    out_md = Path("wiki") / f"dna_amr_xsource_validation_{d}.md"
    out_json = Path("wiki") / f"dna_amr_xsource_validation_{d}.json"
    manifest = {
        "source": "NCBI Pathogen Detection", "pdg": PDG, "source_url": SOURCE_URL,
        "amrfinder_image": AMRFINDER_IMAGE_PINNED, "date": d,
        "n_strains": len(selected), "n_with_run": n_with_run,
        "non_overlap_with": "all BV-BRC cohorts (176 accessions excluded at selection)",
        "accessions": sorted(selected),
        "per_drug": per_drug,
    }
    out_md.write_text(render_md(per_drug, len(selected), n_with_run), encoding="utf-8")
    out_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    for drug in DRUGS:
        a = per_drug[drug]["dna_amr"]; b = per_drug[drug]["naive_amrfinder"]
        print(f"  {drug}: dna-amr acc={a['accuracy']} (sens {a['sensitivity']}/spec {a['specificity']}) "
              f"vs naive acc={b['accuracy']}")
    print(f"Wrote {out_md} + {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
