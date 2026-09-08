"""Does the ResFinder locus-collapse fix change anything for DisinFinder? Measure it, don't assume it.

DisinFinder shares the defective pattern LITERALLY: `dna_decode/disinfinder/runner.py` imports
`resfinder.gene_of` and keys its output on the allele name exactly as ResFinder did before the
2026-09-05 locus-collapse fix. The obvious move is to propagate the fix. This project's standing rule
says otherwise -- a shared code pattern is a LEAD, not a diagnosis, and `plasmid` and `pneumoserotype`
were each measured separately rather than patched by analogy.

THE STRUCTURAL REASON IT MIGHT NOT MATTER. The ResFinder failure needed a DENSE variant family: ~180
catalogued blaTEM alleles all clearing a 90% identity bar against one blaTEM locus. The disinfectant DB
is tiny by comparison, so if it holds no dense families there is nothing to over-report and the rule is
inert -- the same shape as the `plasmid` finding.

This runs BOTH rules from ONE blastn pass per genome, so only the grouping differs, and reports whether
the reported GENE SET ever moves.

Offline: cached assemblies + the committed DisinFinder DB + native blastn.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.resfinder.runner import cluster_alleles_by_locus, gene_of  # noqa: E402
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"


def old_rule(called: list[tuple[str, dict]]) -> set[str]:
    """Pre-fix: one entry per distinct ALLELE NAME (what disinfinder still does)."""
    return {gene_of(aid) for aid, _ in called}


def new_rule(called: list[tuple[str, dict]]) -> set[str]:
    """Post-fix: one entry per genomic LOCUS, identity-primary winner."""
    def winner(cluster):
        return max(cluster, key=lambda kv: (kv[1]["percent_identity"], kv[1]["percent_coverage"]))[0]
    return {gene_of(winner(c)) for c in cluster_alleles_by_locus(called)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data" / "disinfinder_db" / "disinfectants.fsa")
    ap.add_argument("--refseq-cache", type=Path, default=Path("D:/dna_decode_cache/refseq"))
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"disinfinder_locus_collapse_probe_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.db.exists():
        print(f"DisinFinder DB absent at {a.db}", file=sys.stderr)
        return 2
    n_alleles = sum(1 for ln in a.db.read_text(encoding="utf-8", errors="replace").splitlines()
                    if ln.startswith(">"))

    fastas = []
    for d in sorted(a.refseq_cache.iterdir()) if a.refseq_cache.exists() else []:
        if len(fastas) >= a.limit:
            break
        if d.is_dir():
            fa = next((p for p in d.glob("*.fna")), None)
            if fa:
                fastas.append((d.name, fa))
    if not fastas:
        print("no cached assemblies found", file=sys.stderr)
        return 2
    print(f"DisinFinder DB holds {n_alleles} alleles; probing {len(fastas)} genomes\n")

    rows, n_ok, set_diffs, total_old, total_new, multi = [], 0, 0, 0, 0, 0
    for n, (acc, fa) in enumerate(fastas, 1):
        try:
            res = call_alleles(fa, a.db, identity_threshold=90.0, coverage_threshold=60.0,
                               blastn_bin=a.blastn, timeout=600, with_positions=True)
        except Exception as e:                       # noqa: BLE001
            rows.append({"acc": acc, "status": f"error:{type(e).__name__}"})
            continue
        if res.get("status") != "ok":
            rows.append({"acc": acc, "status": res.get("status")})
            continue
        n_ok += 1
        called = [(aid, h) for aid, h in res["per_allele"].items() if h["called"]]
        o, nw = old_rule(called), new_rule(called)
        total_old += len(o); total_new += len(nw)
        if o != nw:
            set_diffs += 1
        if len(o) > len(nw):
            multi += 1
        rows.append({"acc": acc, "status": "ok", "n_alleles_called": len(called),
                     "n_old": len(o), "n_new": len(nw), "identical": o == nw,
                     "only_old": sorted(o - nw), "only_new": sorted(nw - o)})
        if n % 10 == 0:
            print(f"  [{n}/{len(fastas)}] {acc} old={len(o)} new={len(nw)}", flush=True)

    # NON-VACUITY: if nothing was ever called, both rules trivially agree on the empty set and the
    # "inert" verdict is a statement about a pipeline that did nothing.
    if n_ok == 0 or total_old == 0:
        print(f"\nREFUSING: {n_ok} genomes scored and {total_old} genes called in total, so neither "
              "rule was ever exercised. A null here is a plumbing result, not a finding.",
              file=sys.stderr)
        return 3

    print(f"\n=== {n_ok} genomes, {total_old} gene calls under the old rule ===")
    print(f"  genomes whose reported GENE SET differs : {set_diffs}")
    print(f"  genomes where old reported MORE than new: {multi}")

    if set_diffs == 0:
        verdict = "LOCUS_COLLAPSE_IS_INERT_FOR_DISINFINDER"
        why = (f"the reported gene set is IDENTICAL under both rules on all {n_ok} genomes "
               f"({total_old} calls). The DisinFinder DB holds only {n_alleles} alleles with no dense "
               "variant family, so there is nothing for the allele-name keying to over-report. The "
               "ResFinder fix does NOT transfer here and the caller is left UNCHANGED -- but re-measure "
               "if this DB ever grows.")
    else:
        verdict = "LOCUS_COLLAPSE_CHANGES_DISINFINDER_OUTPUT"
        why = (f"the reported gene set differs on {set_diffs} of {n_ok} genomes, so DisinFinder DOES "
               "carry the same over-reporting defect and the fix SHOULD be propagated.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "disinfinder-locus-collapse-probe-v1", "date": _date.today().isoformat(),
           "question": ("does the ResFinder locus-collapse fix change DisinFinder's reported gene "
                        "set? DisinFinder imports resfinder.gene_of and keys on the allele name the "
                        "same way, so the pattern is present -- but a shared pattern is a LEAD, not a "
                        "diagnosis"),
           "n_db_alleles": n_alleles, "n_genomes_scored": n_ok,
           "total_gene_calls_old_rule": total_old, "total_gene_calls_new_rule": total_new,
           "n_genomes_with_different_gene_set": set_diffs,
           "n_genomes_where_old_reported_more": multi,
           "rows": rows, "verdict": verdict, "why": why,
           "honest_limits": [
               "This answers WHETHER the rule changes the reported set, NOT which rule is better. That "
               "would need wet-lab biocide-resistance labels, and it is only worth asking if the set "
               "actually moves.",
               f"The DB holds {n_alleles} alleles. The finding is a property of THIS DB build -- a "
               "larger disinfectant catalogue with dense variant families could behave like the "
               "beta-lactamases and would need re-measuring.",
               "Genomes are cached AMR-cohort assemblies, not a biocide-focused set; disinfectant-gene "
               "content is whatever they happen to carry.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
