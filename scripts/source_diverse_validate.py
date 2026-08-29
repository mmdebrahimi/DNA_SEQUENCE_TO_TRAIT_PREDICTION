"""Score the FROZEN decoder on a source-DIVERSE cohort, and refuse to report if the cohort isn't.

WHY THIS ARM EXISTS. The 10 provenance-disjoint SCORED cells are disjoint from the tuning data by
construction — that was the goal and it was met. Measuring their SOURCE diversity showed 3 of 10 rest on a
single BioProject, and one consequence is concrete: `escherichia_coli_shigella x gentamicin` is 95% one
BioProject containing ZERO `rmt`-family carriers, so it reports sens 0.893 while source-diverse
measurements of the same cell report 0.429 and 0.523. A cohort with no carriers of a determinant family
cannot detect a rule blind to that family.

That measurement was made in a throwaway snippet. This is the same measurement as REPRODUCIBLE TOOLING.

THE SELF-APPLIED STANDARD, which is the point. An arm whose whole argument is "your cohort was too
concentrated" must not ship a concentrated cohort of its own. `MIN_BIOPROJECTS` / `MAX_SOURCE_SHARE` gate
every cell: a scoring set that fails them emits `status: source_concentrated` and NO metrics. Reporting a
number from a single-source cohort here would be the exact error this arm was built to expose.

NAMESPACE-SEPARATE, like every other disclosure layer. Results land in
`wiki/source_diverse_validation_<organism>_<drug>.json` — never in `provenance_disjoint_validation_*`,
whose filename pattern the report card's `load_scored()` globs. Writing there would silently overwrite a
frozen cell with a different number: the documented shared-key trap.

LEAKAGE is the fail-closed accession manifest, not a hand-rolled check. The cheap "never appeared in a
selected.tsv" filter was measured to under-cover by two thirds (956 -> 311).

Run: uv run python scripts/source_diverse_validate.py [--min-per-class 20]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
WIKI = ROOT / "wiki"

# The self-applied bar. Derived from the measurement that motivated the arm: the cells it critiques sit at
# 1-4 BioProjects with 95-100% concentration, and the set that exposed them had 8 at 31%.
MIN_BIOPROJECTS = 5
MAX_SOURCE_SHARE = 0.60
MIN_PER_CLASS = 20

REGISTRY_ORGANISM = {"Escherichia_coli_Shigella": "Escherichia_coli_Shigella",
                     "Klebsiella": "Klebsiella_pneumoniae"}


def load_pool() -> dict:
    """The disjoint, PD-labelled pool produced by `unscored_genome_label_census.py`."""
    f = WIKI / "unscored_genome_label_census.json"
    if not f.exists():
        return {}
    d = json.loads(f.read_text(encoding="utf-8"))
    if not d.get("complete"):
        return {}
    return d.get("labels") or {}


def source_profile(accessions: list[str], provenance: dict) -> dict:
    """Distinct BioProjects + the largest one's share, over the SCORED set. PURE.

    Computed on the accessions actually scored, not on the pool they were drawn from -- a diverse pool
    can still yield a concentrated cell once the drug's labels are applied.
    """
    bps = [provenance.get(a, {}).get("bioproject_acc", "").strip() for a in accessions]
    known = [b for b in bps if b]
    c = Counter(known)
    return {"n": len(accessions), "n_known": len(known), "distinct": len(c),
            "largest_share": round(c.most_common(1)[0][1] / len(known), 3) if known else None,
            "dominant": c.most_common(1)[0][0] if c else None}


def diversity_verdict(prof: dict) -> tuple[bool, str]:
    """Does the SCORING SET clear the arm's own bar? PURE."""
    if prof["n_known"] < prof["n"] * 0.5:
        return False, "provenance unknown for more than half the cohort"
    if prof["distinct"] < MIN_BIOPROJECTS:
        return False, f"only {prof['distinct']} BioProject(s), bar is {MIN_BIOPROJECTS}"
    if prof["largest_share"] is not None and prof["largest_share"] > MAX_SOURCE_SHARE:
        return False, f"largest source holds {prof['largest_share']:.0%}, bar is {MAX_SOURCE_SHARE:.0%}"
    return True, "source-diverse"


def score_cell(group: str, drug: str, labels: dict, provenance: dict, min_per_class: int) -> dict:
    from gentamicin_rmt_candidate import amrfinder_index
    from dna_decode.eval.amr_rules import call_resistance

    idx = amrfinder_index()
    accs = [a for a, r in labels.items()
            if r.get("group") == group
            and (r.get("calls", {}).get(drug) or "").upper() in ("R", "S")
            and a in idx]
    out = {"schema": "source-diverse-validation-v1", "organism": group, "drug": drug,
           "generated": date.today().isoformat(), "n_candidates": len(accs)}
    if not accs:
        out["status"] = "no_candidates"
        return out

    prof = source_profile(accs, provenance)
    out["source_profile"] = prof
    ok, why = diversity_verdict(prof)
    out["diversity_verdict"] = why
    if not ok:
        # The arm refuses to publish a number from a cohort as concentrated as the ones it critiques.
        out["status"] = "source_concentrated"
        return out

    tp = fp = tn = fn = ind = 0
    for a in accs:
        lab = labels[a]["calls"][drug].upper()
        pred = call_resistance(idx[a], drug, organism=REGISTRY_ORGANISM.get(group))["prediction"]
        if pred == "INDETERMINATE":
            ind += 1
            continue
        r = pred == "R"
        if lab == "R":
            tp += r
            fn += not r
        else:
            fp += r
            tn += not r
    n = tp + fp + tn + fn
    out["confusion"] = {"n_scored": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn, "indeterminate": ind}
    if (tp + fn) < min_per_class or (tn + fp) < min_per_class:
        out["status"] = "underpowered"
        out["powering"] = {"R": tp + fn, "S": tn + fp, "bar": min_per_class}
        return out
    out["status"] = "scored"
    out["acc"] = round((tp + tn) / n, 3)
    out["sens"] = round(tp / (tp + fn), 3)
    out["spec"] = round(tn / (tn + fp), 3)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-per-class", type=int, default=MIN_PER_CLASS)
    args = ap.parse_args()

    labels = load_pool()
    if not labels:
        print("no complete disjoint pool on disk -- run scripts/unscored_genome_label_census.py first")
        return 1

    prov_file = WIKI / "source_diverse_provenance.json"
    if not prov_file.exists():
        print(f"provenance sidecar missing: {prov_file.name}")
        print("  the arm cannot verify its OWN cohort's diversity without it, and will not")
        print("  publish a number it cannot check. Generate it, then re-run.")
        return 2
    provenance = json.loads(prov_file.read_text(encoding="utf-8"))

    groups = sorted({r.get("group") for r in labels.values() if r.get("group")})
    drugs = sorted({d for r in labels.values() for d in (r.get("calls") or {})})
    wrote = concentrated = 0
    for g in groups:
        for drug in drugs:
            res = score_cell(g, drug, labels, provenance, args.min_per_class)
            if res.get("status") == "no_candidates":
                continue
            tag = f"{g}_{drug}".lower()
            (WIKI / f"source_diverse_validation_{tag}.json").write_text(
                json.dumps(res, indent=2), encoding="utf-8")
            wrote += 1
            p = res.get("source_profile", {})
            line = (f"{g} x {drug}: {res['status']:20} "
                    f"n={res.get('confusion', {}).get('n_scored', p.get('n', 0))} "
                    f"sources={p.get('distinct', '?')} share={p.get('largest_share')}")
            if res["status"] == "scored":
                line += f"  acc {res['acc']} sens {res['sens']} spec {res['spec']}"
            if res["status"] == "source_concentrated":
                concentrated += 1
                line += f"  ({res['diversity_verdict']})"
            print(" ", line)

    print()
    print(f"wrote {wrote} cell artifact(s); {concentrated} REFUSED as source-concentrated.")
    print("Namespace: source_diverse_validation_* -- never provenance_disjoint_validation_*, whose glob")
    print("the report card reads. These AUGMENT the frozen cells; they never replace them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
