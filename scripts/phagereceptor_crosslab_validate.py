"""CROSS-LAB independent RBP number: score the phage RBP caller against the phageReceptor database.

Corrects the 2026-07-24 /research "data-blocked" verdict (which missed phageReceptor). phageReceptor
(Zhang et al., Bioinformatics 2020, 36(10):2975; Peng lab, Hunan U — INDEPENDENT of LBNL/Arkin-Mutalik
and BASEL/Maffei) has 37 E. coli phages with outer-membrane-protein receptors in the caller's vocabulary.
Excluding the classic model phages (T4/T7/T5/lambda/N4) + any LBNL-set overlap gives an INDEPENDENT test set.

Protocol: for each independent phageReceptor E. coli phage with an OMP-class measured receptor, fetch its
genome from GenBank, extract every tail-fiber CDS protein, take the best nearest-neighbour against the
committed LBNL RBP reference (call_rbp_from_protein), and compare the predicted receptor to phageReceptor's
MEASURED receptor. This is genuinely cross-lab: the RBP reference is LBNL, the test labels are phageReceptor.

Usage:  uv run --with biopython python scripts/phagereceptor_crosslab_validate.py \
          --api http://www.computationalbiology.cn:18887/viralRecepetor  [--json cached.json]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from dna_decode.phage.rbp_caller import call_rbp_from_protein

# descriptive receptor name (phageReceptor) -> our OMP class
def to_class(desc: str) -> str | None:
    s = (desc or "").lower()
    if "ferrichrome" in s: return "FhuA"
    if "porin c" in s: return "OmpC"
    if "porin f" in s: return "OmpF"
    if "porin a" in s or ("protein a" in s and "outer membrane" in s): return "OmpA"
    if "maltose" in s: return "LamB"
    if "cobalamin" in s or "cobinamide" in s: return "BtuB"
    if "long-chain fatty" in s: return "FadL"
    if "nucleoside" in s: return "Tsx"
    if "ferric enterobactin" in s: return "FepA"
    if "tolc" in s: return "TolC"
    if "n4 receptor" in s or ("n4" in s and "receptor" in s): return "NfrA"
    return None

CLASSIC = {"t4", "t7", "t5", "t2", "t6", "t3", "n4", "lambda"}
LBNL_OVERLAP = {"m1", "t1"}   # confirmed present in the LBNL Table_S1 by exact name
_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _is_independent(name: str) -> bool:
    toks = set(name.lower().replace("phage", " ").replace("enterobacteria", " ")
               .replace("escherichia", " ").replace("virus", " ").split())
    return not (toks & CLASSIC) and not (toks & LBNL_OVERLAP)


def build_test_set(api_json: dict) -> dict[str, str]:
    rows = [r for t in api_json.values() for r in t[1:]]
    ecoli = [r for r in rows if len(r) > 6 and "escherichia coli" in (r[4] or "").lower()]
    test: dict[str, str] = {}
    for r in ecoli:
        cls = to_class(r[6])
        name = (r[2] or "").strip()
        if cls and name and _is_independent(name):
            test.setdefault(name, cls)   # first receptor per phage
    return test


def _efetch_genbank(name: str) -> str | None:
    q = name.replace(" ", "+") + "+complete+genome"
    try:
        s = subprocess.run(["curl", "-s", "-m", "30",
                            f"{_EUTILS}/esearch.fcgi?db=nucleotide&term={q}&retmax=1"],
                           capture_output=True, text=True, timeout=40).stdout
        import re
        m = re.search(r"<Id>(\d+)</Id>", s)
        if not m:
            return None
        uid = m.group(1)
        gb = subprocess.run(["curl", "-s", "-m", "60",
                            f"{_EUTILS}/efetch.fcgi?db=nucleotide&id={uid}&rettype=gb&retmode=text"],
                           capture_output=True, text=True, timeout=80).stdout
        return gb if gb.startswith("LOCUS") else None
    except Exception:
        return None


def _tail_fiber_proteins(gb_text: str, cache: Path) -> list[str]:
    from Bio import SeqIO
    import io
    prots = []
    try:
        for rec in SeqIO.parse(io.StringIO(gb_text), "genbank"):
            for f in rec.features:
                if f.type != "CDS":
                    continue
                prod = " ".join(f.qualifiers.get("product", [])).lower()
                if any(k in prod for k in ("tail fiber", "tail spike", "tailspike",
                                           "receptor binding", "receptor-binding", "host specificity")):
                    tr = (f.qualifiers.get("translation") or [""])[0]
                    if tr:
                        prots.append(tr)
    except Exception:
        pass
    return prots


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--api", default="http://www.computationalbiology.cn:18887/viralRecepetor")
    ap.add_argument("--json", help="cached phageReceptor JSON (skip the API fetch)")
    ap.add_argument("--cache-dir", default="data/phage_ref/phagereceptor_genomes")
    ap.add_argument("--min-similarity", type=float, default=0.05)
    ap.add_argument("--date", default="2026-07-24")
    args = ap.parse_args()

    if args.json:
        api = json.load(open(args.json, encoding="utf-8"))
    else:
        raw = subprocess.run(["curl", "-s", "-m", "40", args.api], capture_output=True, text=True, timeout=60).stdout
        api = json.loads(raw)
    test = build_test_set(api)
    print(f"independent phageReceptor E. coli OMP test set: {len(test)} phages")

    cache = Path(args.cache_dir); cache.mkdir(parents=True, exist_ok=True)
    preds = []
    for name, true_rec in sorted(test.items()):
        gb_path = cache / (name.replace(" ", "_").replace("/", "_") + ".gb")
        if gb_path.exists():
            gb = gb_path.read_text(encoding="utf-8", errors="replace")
        else:
            gb = _efetch_genbank(name)
            time.sleep(0.4)
            if gb:
                gb_path.write_text(gb, encoding="utf-8")
        if not gb:
            preds.append({"phage": name, "true": true_rec, "predicted": None,
                          "status": "NO_GENOME", "correct": False}); continue
        tfs = _tail_fiber_proteins(gb, cache)
        if not tfs:
            preds.append({"phage": name, "true": true_rec, "predicted": None,
                          "status": "NO_TAIL_FIBER", "correct": False}); continue
        best = None
        for prot in tfs:
            call = call_rbp_from_protein(prot, min_similarity=args.min_similarity)
            if call.status == "CALLED" and (best is None or (call.similarity or 0) > (best.similarity or 0)):
                best = call
        if best is None:
            preds.append({"phage": name, "true": true_rec, "predicted": None,
                          "status": "INDETERMINATE", "correct": False, "n_tail_fiber": len(tfs)}); continue
        correct = best.predicted_receptor == true_rec
        preds.append({"phage": name, "true": true_rec, "predicted": best.predicted_receptor,
                      "nearest": best.nearest_phage, "similarity": round(best.similarity or 0, 3),
                      "status": "CALLED", "correct": bool(correct), "n_tail_fiber": len(tfs)})

    from collections import Counter
    called = [p for p in preds if p["status"] == "CALLED"]
    correct = sum(1 for p in called if p["correct"])
    per = {}
    for p in called:
        b = per.setdefault(p["true"], [0, 0]); b[1] += 1; b[0] += int(p["correct"])
    out = {
        "cell": "phage_receptor_class_RBP_crosslab", "date": args.date,
        "reference": "committed LBNL RBP reference (data/phage_ref/rbp_reference.faa)",
        "test_source": "phageReceptor DB (Zhang et al. Bioinformatics 2020; Peng lab, Hunan U) — INDEPENDENT of LBNL + BASEL",
        "independence": "different lab, literature-curated measured receptors; classic model phages + LBNL-name overlaps (M1,T1) excluded",
        "n_test": len(test), "n_called": len(called), "n_correct": correct,
        "n_no_genome": sum(1 for p in preds if p["status"] == "NO_GENOME"),
        "n_no_tail_fiber": sum(1 for p in preds if p["status"] == "NO_TAIL_FIBER"),
        "n_indeterminate": sum(1 for p in preds if p["status"] == "INDETERMINATE"),
        "crosslab_accuracy": (correct / len(called)) if called else None,
        "per_receptor_correct_called": per,
        "predictions": preds,
        "honest_scope": "CROSS-LAB (LBNL reference vs phageReceptor test labels). RBP extracted by 'tail fiber' "
                        "product annotation (imperfect — a phage may carry several tail fibers; best match taken). "
                        "phageReceptor receptors mapped from descriptive names to OMP classes. This is the "
                        "genuinely cross-lab RBP number the within-LBNL LOO (0.975) could not be.",
    }
    Path("wiki").mkdir(exist_ok=True)
    (Path("wiki") / f"phage_rbp_crosslab_result_{args.date}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_test", "n_called", "n_correct", "crosslab_accuracy",
          "n_no_genome", "n_no_tail_fiber", "n_indeterminate", "per_receptor_correct_called")}, indent=2))


if __name__ == "__main__":
    main()
