"""Does the shipped ResFinder caller agree with an INDEPENDENT curated caller? Measure it, before/after.

THE DEFECT. `dna-resfinder` was registered FAITHFUL_TO_TOOL and never measured. Pointing it at one
cached E. coli genome returned **190 beta-lactam genes**, against AMRFinder's 3, with AMRFinder's calls
a SUBSET. The mechanism is exact and it is the third instance of one failure in this repo: TEM
beta-lactamase variants differ by one to three point mutations, so a single blaTEM locus matches ~180
catalog TEM alleles above the 90% identity bar, and the caller keyed its output on the ALLELE name --
so every variant cleared independently and all were reported present. The genome carries blaTEM-1B
(narrow-spectrum penicillinase, 100.0/100.0); the caller ALSO reported blaTEM-52B, blaTEM-52C,
blaTEM-12, blaTEM-10 and blaTEM-24 at 99.0-99.8% identity. Those are ESBLs. That is not a counting
problem -- it is a wrong clinical interpretation of a genome.

THE FIX. Collapse called alleles to one entry per genomic LOCUS (position, not name -- blaOXA-1 and
blaOXA-48 share a prefix, are functionally unrelated, and can genuinely co-occur), picking the
best-matching allele there IDENTITY-PRIMARY. Greedy-representative clustering, not single-linkage,
so a tandem array cannot chain into a single call.

WHAT THIS SCRIPT MEASURES. Both callers on the same cached genomes, scored against AMRFinder's own
committed output for the same accession -- an independent, curated, widely-used tool that this project
already depends on. No new labels are needed: the question is whether the caller agrees with an
independent implementation, and it is answered before and after the fix on identical inputs.

Offline: cached assemblies + the committed ResFinder DB + native blastn + committed AMRFinder runs.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.resfinder.runner import (  # noqa: E402
    cluster_alleles_by_locus, gene_of,
)
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402
from gentamicin_rmt_candidate import amrfinder_index, gene_symbol  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
CLASSES = {"aminoglycoside": "AMINOGLYCOSIDE", "beta-lactam": "BETA-LACTAM"}

# A trailing single letter after the allele number is a synonymous nucleotide variant of the same
# protein (blaTEM-1A/1B/1C/1D all encode TEM-1). AMRFinder reports `blaTEM-1`, ResFinder's allele is
# `blaTEM-1B`, so an EXACT symbol comparison would score a correct agreement as a miss. Both the exact
# and the normalized number are reported; the normalized one is the LENIENT reading and is labelled so.
_TRAILING_VARIANT_LETTER = re.compile(r"^(?P<stem>.+-\d+)[A-Za-z]$")


def normalize_symbol(sym: str) -> str:
    m = _TRAILING_VARIANT_LETTER.match(sym.strip())
    return m.group("stem") if m else sym.strip()


def old_rule_genes(called: list[tuple[str, dict]]) -> set[str]:
    """The PRE-FIX behaviour, reproduced locally: one entry per distinct ALLELE NAME.

    Re-implemented rather than imported so the comparison cannot silently change when the module does.
    """
    return {gene_of(aid) for aid, _ in called}


def new_rule_genes(called: list[tuple[str, dict]]) -> set[str]:
    """The POST-FIX behaviour: one entry per LOCUS, identity-primary winner."""
    out = set()
    for cluster in cluster_alleles_by_locus(called):
        allele_id, _ = max(cluster, key=lambda kv: (kv[1]["percent_identity"],
                                                    kv[1]["percent_coverage"]))
        out.add(gene_of(allele_id))
    return out


def amrfinder_genes(main_tsv: Path, amr_class: str) -> set[str]:
    """Acquired genes AMRFinder called in this class. Point mutations are EXCLUDED.

    A ResFinder allele DB contains acquired GENES only, so scoring it against AMRFinder rows that are
    point mutations (`Element type` == POINT, e.g. the promoter variant blaTEMp_G162T) would charge the
    caller with missing something it cannot represent.
    """
    out = set()
    with open(main_tsv, encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if amr_class not in (row.get("Class") or "").upper():
                continue
            etype = (row.get("Element type") or row.get("Type") or "").strip().upper()
            if etype == "POINT":
                continue
            sym = gene_symbol(row)
            if sym:
                out.add(sym)
    return out


def jaccard(a: set[str], b: set[str]) -> float | None:
    return None if not (a | b) else len(a & b) / len(a | b)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    ap.add_argument("--db-dir", type=Path, default=ROOT / "data" / "resfinder_db")
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/resfinder_collapse/rows.jsonl"))
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"resfinder_locus_collapse_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    idx = amrfinder_index()
    fastas = {}
    for d in sorted(a.refseq_cache.iterdir()) if a.refseq_cache.exists() else []:
        if d.is_dir() and d.name in idx:
            fa = next((p for p in d.glob("*.fna")), None)
            if fa:
                fastas[d.name] = fa
    accs = sorted(fastas)
    if a.limit:
        accs = accs[:a.limit]
    if not accs:
        print("no accession has BOTH a cached assembly and an AMRFinder run", file=sys.stderr)
        return 2
    print(f"{len(accs)} accessions with both a cached assembly and a committed AMRFinder run\n")

    a.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for line in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["acc"]] = r
        print(f"  resuming: {len(done)} already scored")
    fh = open(a.checkpoint, "a", encoding="utf-8")

    for n, acc in enumerate(accs, 1):
        if acc in done:
            continue
        rec: dict = {"acc": acc, "classes": {}}
        for cls, amr_class in CLASSES.items():
            db = a.db_dir / f"{cls}.fsa"
            if not db.exists():
                continue
            try:
                res = call_alleles(fastas[acc], db, identity_threshold=90.0, coverage_threshold=60.0,
                                   blastn_bin=a.blastn, timeout=600, with_positions=True)
            except Exception as e:                       # noqa: BLE001
                rec["classes"][cls] = {"status": f"error:{type(e).__name__}"}
                continue
            if res.get("status") != "ok":
                rec["classes"][cls] = {"status": res.get("status")}
                continue
            called = [(aid, h) for aid, h in res["per_allele"].items() if h["called"]]
            old, new = old_rule_genes(called), new_rule_genes(called)
            amr = amrfinder_genes(idx[acc], amr_class)
            rec["classes"][cls] = {
                "status": "ok", "n_alleles_called": len(called),
                "n_old": len(old), "n_new": len(new), "n_amrfinder": len(amr),
                "old": sorted(old), "new": sorted(new), "amrfinder": sorted(amr),
                "j_old_exact": jaccard(old, amr), "j_new_exact": jaccard(new, amr),
                "j_old_norm": jaccard({normalize_symbol(g) for g in old},
                                      {normalize_symbol(g) for g in amr}),
                "j_new_norm": jaccard({normalize_symbol(g) for g in new},
                                      {normalize_symbol(g) for g in amr}),
            }
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        done[acc] = rec
        if n % 50 == 0:
            print(f"  [{n}/{len(accs)}] {acc}", flush=True)
    fh.close()

    rows = [done[x] for x in accs if x in done]
    summary = {}
    for cls in CLASSES:
        ok = [r["classes"][cls] for r in rows
              if r.get("classes", {}).get(cls, {}).get("status") == "ok"]
        if not ok:
            continue
        def mean(key):
            v = [c[key] for c in ok if c.get(key) is not None]
            return sum(v) / len(v) if v else None
        summary[cls] = {
            "n_genomes": len(ok),
            "total_alleles_called": sum(c["n_alleles_called"] for c in ok),
            "mean_genes_old": mean("n_old"), "mean_genes_new": mean("n_new"),
            "mean_genes_amrfinder": mean("n_amrfinder"),
            "jaccard_vs_amrfinder_exact_old": mean("j_old_exact"),
            "jaccard_vs_amrfinder_exact_new": mean("j_new_exact"),
            "jaccard_vs_amrfinder_normalized_old": mean("j_old_norm"),
            "jaccard_vs_amrfinder_normalized_new": mean("j_new_norm"),
            "n_genomes_where_old_multireported_a_locus":
                sum(1 for c in ok if c["n_old"] > c["n_new"]),
        }

    # NON-VACUITY. If blastn never called anything, both rules agree on the empty set and every metric
    # above is a statement about a broken pipeline wearing the costume of a clean result.
    total_alleles = sum(s["total_alleles_called"] for s in summary.values())
    n_scored = sum(s["n_genomes"] for s in summary.values())
    if not summary or total_alleles == 0 or n_scored == 0:
        print(f"\nREFUSING to report: {n_scored} genome-classes scored, {total_alleles} alleles called "
              "in total. Neither rule was ever exercised, so agreement here is a plumbing result.",
              file=sys.stderr)
        return 3

    print(f"\n=== {len(rows)} genomes ===")
    for cls, s in summary.items():
        print(f"\n{cls}  (n={s['n_genomes']} genomes, {s['total_alleles_called']} alleles called)")
        print(f"  mean genes/genome   old={s['mean_genes_old']:.2f}  "
              f"new={s['mean_genes_new']:.2f}  AMRFinder={s['mean_genes_amrfinder']:.2f}")
        print(f"  Jaccard vs AMRFinder (exact)      old={s['jaccard_vs_amrfinder_exact_old']:.4f}  "
              f"new={s['jaccard_vs_amrfinder_exact_new']:.4f}")
        print(f"  Jaccard vs AMRFinder (normalized) old={s['jaccard_vs_amrfinder_normalized_old']:.4f}  "
              f"new={s['jaccard_vs_amrfinder_normalized_new']:.4f}")
        print(f"  genomes where the OLD rule multi-reported a single locus: "
              f"{s['n_genomes_where_old_multireported_a_locus']}/{s['n_genomes']}")

    gains = [s["jaccard_vs_amrfinder_normalized_new"] - s["jaccard_vs_amrfinder_normalized_old"]
             for s in summary.values()]
    if all(g > 0 for g in gains):
        verdict = "LOCUS_COLLAPSE_IMPROVES_AGREEMENT_WITH_AN_INDEPENDENT_CALLER"
        why = ("collapsing called alleles to one entry per genomic locus raises agreement with "
               "AMRFinder in EVERY drug class measured, on identical inputs. The pre-fix caller "
               "reported a single locus as many separately-present genes, which for beta-lactamases "
               "means reporting ESBLs in a genome that carries only a narrow-spectrum enzyme.")
    elif any(g > 0 for g in gains):
        verdict = "LOCUS_COLLAPSE_IMPROVES_SOME_CLASSES_ONLY"
        why = ("agreement improves in some drug classes and not others, so the fix is not uniformly "
               "beneficial and the per-class numbers must be quoted rather than a pooled figure.")
    else:
        verdict = "LOCUS_COLLAPSE_DOES_NOT_IMPROVE_AGREEMENT"
        why = ("collapsing to one entry per locus does NOT raise agreement with AMRFinder. The "
               "over-reporting is real but the collapse is not the right remedy on this evidence.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "resfinder-locus-collapse-v1", "date": _date.today().isoformat(),
           "question": ("does collapsing ResFinder allele hits to one call per genomic LOCUS improve "
                        "agreement with an independent curated caller (AMRFinder)?"),
           "comparator": ("AMRFinder committed per-genome output for the SAME accession. Independent "
                          "curated implementation; acquired genes only (POINT rows excluded, since a "
                          "ResFinder allele DB cannot represent a point mutation)."),
           "n_genomes": len(rows), "summary": summary, "verdict": verdict, "why": why,
           "honest_limits": [
               "The comparator is a TOOL, not a wet-lab label. This measures agreement with an "
               "independent implementation, NOT correctness -- both callers could be wrong together. "
               "The cell stays FAITHFUL_TO_TOOL; what changes is WHICH tool it is faithful to and how "
               "closely.",
               "Only two ResFinder class DBs are committed (aminoglycoside, beta-lactam), so this says "
               "nothing about the other resistance classes the real ResFinder DB covers.",
               "The normalized Jaccard strips a trailing variant letter (blaTEM-1B -> blaTEM-1) to "
               "avoid scoring a correct agreement as a miss. It is the LENIENT reading; the exact "
               "number is reported beside it.",
               "Genomes are whatever this project already had cached, drawn from AMR cohorts -- they "
               "are enriched for resistance and are NOT a random sample of the species.",
               "Locus clustering uses a 50% reciprocal-overlap bar on the shorter interval. A tandem "
               "array whose copies overlap more than that would still collapse into one call.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
