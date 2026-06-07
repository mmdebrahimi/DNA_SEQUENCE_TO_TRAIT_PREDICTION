"""Generalized cross-organism AMR validation — any NCBI Pathogen Detection organism × any supported drug.

Generalizes scripts/klebsiella_drug_validate.py to an arbitrary organism. The deployed dna-amr per-drug
rule (DRUG_RULE) is applied UNCHANGED; the only organism-specific input is AMRFinder's `-O`. Tests whether
the deterministic "count the mechanism, not the broad class bag" approach transfers to a new organism.

    uv run python scripts/organism_drug_validate.py --ncbi-group Pseudomonas_aeruginosa \
        --amrfinder-organism Pseudomonas_aeruginosa --drug ciprofloxacin
    ... --eval-only

MVP bar per (organism,drug): acc >= 0.80 AND sens >= 0.80 → VALIDATED; else FAILS_BAR (documents the
organism-specific failure mode — often an intrinsic-efflux blind spot, per the Klebsiella tet finding).
Restartable. Labels: NCBI Pathogen Detection AST (independent source). AMRFinder ~95s/strain.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import urllib.request
from datetime import date as _date
from pathlib import Path

from dna_decode.data import refseq
from dna_decode.data.mic_tiers import supported_drugs
from dna_decode.eval.amr_rules import (
    AMRFINDER_IMAGE_PINNED, call_resistance, cipro_determinants_from_main, discordance_bucket,
)

TARGET_PER_CLASS = 15


def latest_metadata_url(group: str) -> str:
    base = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{group}/latest_snps/Metadata/"
    with urllib.request.urlopen(base, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    m = re.findall(r'href="(PDG[0-9.]+\.metadata\.tsv)"', html)
    if not m:
        raise RuntimeError(f"no PDG metadata file found at {base}")
    return base + sorted(m)[-1]


def select_cohort(group, drug, selected: Path, reuse_glob: str) -> dict[str, int]:
    if selected.exists():
        out = {}
        for ln in selected.read_text().splitlines():
            if ln.strip():
                a, rs = ln.split("\t"); out[a] = 1 if rs == "R" else 0
        return out
    cached = set()
    for root in glob.glob(reuse_glob):
        cached |= {p.name for p in Path(root).iterdir()} if Path(root).exists() else set()
    cand_cached, cand_new = [], []
    tok_r, tok_s = f"{drug}=R", f"{drug}=S"
    with urllib.request.urlopen(latest_metadata_url(group), timeout=400) as resp:
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
    selected.parent.mkdir(parents=True, exist_ok=True)
    selected.write_text("".join(f"{a}\t{'R' if v else 'S'}\n" for a, v in chosen.items()), encoding="utf-8")
    return chosen


def _run_dir(acc, own: Path, reuse_glob: str):
    roots = [own, *[Path(p) for p in glob.glob(reuse_glob)]]
    for root in roots:
        m = root / acc / "main.tsv"
        if m.exists() and m.stat().st_size > 0:
            return root / acc
    return None


def ensure_run(acc, own_runs: Path, gcache: Path, amr_org: str, reuse_glob: str):
    rd = _run_dir(acc, own_runs, reuse_glob)
    if rd:
        return rd / "main.tsv"
    try:
        refseq.download_genome(acc, gcache)
        fasta = refseq.fasta_path(acc, gcache)
    except Exception as e:
        print(f"  [{acc}] download FAILED: {type(e).__name__}: {e}"); return None
    if not Path(fasta).exists():
        return None
    try:
        import scripts.drug_mechanism_audit as dma
        from scripts.drug_mechanism_audit import _run_amrfinder
        dma.AMRFINDER_DB = str(Path("data/amrfinder_db").resolve())
        out_dir = own_runs / acc; out_dir.mkdir(parents=True, exist_ok=True)
        _run_amrfinder(Path(fasta), out_dir, organism=amr_org)
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ncbi-group", required=True, help="NCBI Pathogen Detection organism group dir")
    ap.add_argument("--amrfinder-organism", required=True, help="AMRFinder -O value")
    ap.add_argument("--drug", required=True, choices=supported_drugs())
    ap.add_argument("--eval-only", action="store_true")
    args = ap.parse_args(argv)
    slug = args.ncbi_group.lower().replace(" ", "_")
    base = Path(f"data/raw/{slug}_{args.drug}")
    own_runs = base / "amrfinder_runs"
    gcache = base / "refseq"
    reuse_glob = f"data/raw/{slug}_*/amrfinder_runs"
    sel = select_cohort(args.ncbi_group, args.drug, base / "selected.tsv", reuse_glob)
    print(f"{args.ncbi_group} {args.drug} cohort: {len(sel)} ({sum(sel.values())}R/{len(sel)-sum(sel.values())}S)")
    if not args.eval_only:
        for i, acc in enumerate(sel, 1):
            print(f"[{i}/{len(sel)}] {acc} ...")
            ensure_run(acc, own_runs, gcache, args.amrfinder_organism, reuse_glob)
    n_run = sum(1 for a in sel if _run_dir(a, own_runs, reuse_glob))
    print(f"with AMRFinder runs: {n_run}/{len(sel)}")
    dna, naive, na = [], [], 0
    disc = {"FN_undetected_mechanism": 0, "FP_determinant_without_phenotype": 0}
    for acc, y in sel.items():
        rd = _run_dir(acc, own_runs, reuse_glob)
        if not rd:
            na += 1; continue
        c = call_resistance(rd / "main.tsv", args.drug)
        if c["prediction"] == "INDETERMINATE":
            na += 1; continue
        dna.append((1 if c["prediction"] == "R" else 0, y))
        b = discordance_bucket(c["prediction"], y)
        if b:
            disc[b] += 1
        naive.append((1 if len(cipro_determinants_from_main(rd / "main.tsv", args.drug)) >= 1 else 0, y))
    a, nv = _conf(dna), _conf(naive)
    ok = bool(a["accuracy"] and a["accuracy"] >= 0.80 and a["sensitivity"] and a["sensitivity"] >= 0.80)
    verdict = "VALIDATED" if ok else "FAILS_BAR"
    d = _date.today().isoformat()
    md = [
        f"# {args.ncbi_group} {args.drug} — cross-organism validation — {d}",
        "", f"> Deployed dna-amr {args.drug} rule applied UNCHANGED. AMRFinder `-O {args.amrfinder_organism}`.",
        f"- NCBI group `{args.ncbi_group}`; cohort {len(sel)} ({sum(sel.values())}R/{len(sel)-sum(sel.values())}S), {n_run} runs; `{AMRFINDER_IMAGE_PINNED}`",
        "", f"## VERDICT: {verdict}", "",
        "| caller | N | acc | sens | spec |", "|---|---:|---:|---:|---:|",
        f"| **dna-amr ({args.drug})** | {a['n']} | **{a['accuracy']}** | {a['sensitivity']} | {a['specificity']} |",
        f"| naive AMRFinder | {nv['n']} | {nv['accuracy']} | {nv['sensitivity']} | {nv['specificity']} |",
        "", f"## Discordance", f"- FN (R missed): {disc['FN_undetected_mechanism']}",
        f"- FP (called R, susceptible): {disc['FP_determinant_without_phenotype']}",
        "", "## Honest scope",
        f"1 organism, 1 drug, N={len(sel)}, NCBI labels (different source/curation, not a different-lab study).",
    ]
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"{slug}_{args.drug}_validate_{d}.md").write_text("\n".join(md), encoding="utf-8")
    (Path("wiki") / f"{slug}_{args.drug}_validate_{d}.json").write_text(json.dumps(
        {"ncbi_group": args.ncbi_group, "amrfinder_organism": args.amrfinder_organism, "drug": args.drug,
         "date": d, "n": len(sel), "n_run": n_run, "accessions": sorted(sel),
         "dna_amr": a, "naive_amrfinder": nv, "discordance": disc, "na": na, "verdict": verdict}, indent=2),
        encoding="utf-8")
    print(f"VERDICT: {verdict}  dna-amr acc={a['accuracy']} sens={a['sensitivity']} spec={a['specificity']} | naive acc={nv['accuracy']}")
    print(f"Wrote wiki/{slug}_{args.drug}_validate_{d}.{{md,json}}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
