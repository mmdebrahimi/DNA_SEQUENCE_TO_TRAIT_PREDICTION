"""Leave-one-out validation of the organism router + the mandatory out-of-set abstention control.

    uv run python scripts/identify_validate.py

LEAVE-ONE-OUT WITHOUT REBUILDING THE SKETCH. For a nearest-neighbour router, holding a genome out of
the reference is exactly "ignore its own row in the all-pairs matrix". The plan specified rebuilding
the sketch per held-out genome -- ~212 container spin-ups for identical information. One cached
`mash dist reference.msh reference.msh` serves every fold. Same strictness: the self-hit is excluded,
so no genome is ever scored against itself.

TWO REFUSALS, both deliberate:
  * an EMPTY out-of-set group -> no verdict (exit 3). A router validated only on organisms it knows
    has not been tested on the thing it exists to refuse.
  * out-of-set abstention of 0% -> the run FAILS (exit 1). A router that never abstains is not usable
    as a router, regardless of how good its accuracy looks.

PER-ORGANISM, NEVER POOLED ALONE. The pooled figure is emitted only beside the composition table that
explains it. On the earlier unbalanced corpus E. coli was 68% of the reference, so a pooled number
measured cohort composition; the reference is now balanced (largest share ~0.12) but the discipline
stands -- a pooled rate over unequal groups is not a property of the router.

ERROR != ABSTAIN. An identification that raised is counted separately from one that declined, so a
crash can never read as caution.
"""
from __future__ import annotations

import argparse
import collections
import json
import subprocess
from datetime import date
from pathlib import Path

from dna_decode.data.organism_vocab import UnknownOrganism, canonical
from dna_decode.identify.core import accession_from_reference_id, decide, parse_mash_dist
from dna_decode.identify.thresholds import frozen
from scripts.build_identify_reference import genus_species, organism_from_genbank, refseq_roots

REPO = Path(__file__).resolve().parent.parent
IDENT = Path("D:/dna_decode_cache/identify")
MATRIX_CACHE = IDENT / "allpairs.tsv"
MASH_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"


def _mash(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["docker", "run", "--rm", "-v", "D:/:/d", MASH_IMAGE] + args,
        capture_output=True, text=True, timeout=3600)


def labels_from_manifest(manifest: Path) -> dict[str, str]:
    m = json.loads(manifest.read_text(encoding="utf-8"))
    return {a: o for o, accs in m["accessions_by_organism"].items() for a in accs}


def leave_one_out(matrix_text: str, labels: dict[str, str], thresholds) -> dict:
    """Score every reference genome against the reference MINUS itself."""
    by_query: dict[str, list] = collections.defaultdict(list)
    for h in parse_mash_dist(matrix_text):
        q = accession_from_reference_id(h.query_id)
        if accession_from_reference_id(h.reference_id) == q:
            continue                                   # THE held-out step: never score against self
        by_query[q].append(h)

    per_org: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    confusions: collections.Counter = collections.Counter()
    unscorable: list[str] = []

    for acc, truth in labels.items():
        hits = sorted(by_query.get(acc, []), key=lambda h: h.distance)
        # an organism with ONE genome has no same-organism neighbour left after holding it out, so it
        # is structurally unscorable -- named, never silently counted as a miss
        if sum(1 for a, o in labels.items() if o == truth and a != acc) == 0:
            unscorable.append(acc)
            per_org[truth]["unscorable"] += 1
            continue
        call = decide(hits, thresholds,
                      lambda a: labels.get(a),
                      lambda o: o)
        if call.abstained:
            per_org[truth]["abstain"] += 1
            per_org[truth][f"abstain_{call.reason.value}"] += 1
        elif call.organism == truth:
            per_org[truth]["hit"] += 1
        else:
            per_org[truth]["miss"] += 1
            confusions[f"{truth}->{call.organism}"] += 1

    scored = {o: c for o, c in per_org.items()}
    tot_hit = sum(c["hit"] for c in scored.values())
    tot_miss = sum(c["miss"] for c in scored.values())
    tot_abst = sum(c["abstain"] for c in scored.values())
    n = tot_hit + tot_miss + tot_abst
    return {
        "per_organism": {o: dict(c) for o, c in sorted(scored.items())},
        "pooled": {"n_scored": n, "hit": tot_hit, "miss": tot_miss, "abstain": tot_abst,
                   "accuracy_on_called": round(tot_hit / (tot_hit + tot_miss), 4)
                   if (tot_hit + tot_miss) else None,
                   "accuracy_over_all": round(tot_hit / n, 4) if n else None},
        "pooled_caveat": "read ONLY beside per_organism -- a pooled rate over unequal groups is a "
                         "property of the composition, not of the router",
        "confusions": dict(confusions.most_common()),
        "structurally_unscorable": unscorable,
    }


def out_of_set_control(labels: dict[str, str], thresholds) -> dict:
    """Score genomes whose organism is NOT in the supported set. The control that decides usability."""
    oos: dict[str, tuple[str, Path]] = {}
    seen: set[str] = set()
    for root in refseq_roots():
        if not root.is_dir():
            continue
        for p in root.iterdir():
            if not p.is_dir() or not (p / "genome.fna").exists() or p.name in seen:
                continue
            seen.add(p.name)
            nm = organism_from_genbank(p / "annotations.gbk")
            if not nm:
                continue
            key = genus_species(nm)
            try:
                canonical(key)
            except UnknownOrganism:
                oos[p.name] = (key, (p / "genome.fna").resolve())

    if not oos:
        return {"status": "EMPTY", "n_genomes": 0}

    fof = IDENT / "oos_inputs.txt"
    fof.write_text("\n".join("/d/" + v[1].relative_to(Path("D:/")).as_posix()
                             for v in oos.values()) + "\n", encoding="utf-8", newline="\n")
    r = _mash(["mash", "sketch", "-l", "/d/dna_decode_cache/identify/oos_inputs.txt",
               "-o", "/d/dna_decode_cache/identify/oos", "-k", "21", "-s", "1000"])
    if not (IDENT / "oos.msh").exists():
        return {"status": "SKETCH_FAILED", "stderr": (r.stderr or "")[-400:]}
    r = _mash(["mash", "dist", "/d/dna_decode_cache/identify/reference.msh",
               "/d/dna_decode_cache/identify/oos.msh"])

    by_query: dict[str, list] = collections.defaultdict(list)
    for h in parse_mash_dist(r.stdout or ""):
        by_query[accession_from_reference_id(h.query_id)].append(h)

    per_sp: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    nearest: dict[str, dict] = {}
    for acc, (sp, _) in oos.items():
        hits = sorted(by_query.get(acc, []), key=lambda h: h.distance)
        call = decide(hits, thresholds, lambda a: labels.get(a), lambda o: o)
        per_sp[sp]["abstain" if call.abstained else "false_call"] += 1
        if call.abstained:
            per_sp[sp][f"abstain_{call.reason.value}"] += 1
        else:
            per_sp[sp][f"called_{call.organism}"] += 1
        if hits:
            lab = labels.get(accession_from_reference_id(hits[0].reference_id))
            d = nearest.setdefault(sp, {"min": 9.0, "max": 0.0, "nearest_in_set": lab})
            d["min"] = min(d["min"], hits[0].distance)
            d["max"] = max(d["max"], hits[0].distance)

    n = len(oos)
    abst = sum(c["abstain"] for c in per_sp.values())
    return {
        "status": "OK",
        "n_genomes": n,
        "n_abstained": abst,
        "abstention_rate": round(abst / n, 4),
        "per_species": {s: dict(c) for s, c in sorted(per_sp.items())},
        "nearest_in_set_distance": {s: {k: (round(v, 5) if isinstance(v, float) else v)
                                        for k, v in d.items()} for s, d in sorted(nearest.items())},
        "control_strength": "dominated by CONGENERS of in-set species (K. aerogenes vs "
                            "K. pneumoniae), so this is the HARD near-neighbour case rather than a "
                            "trivially distant one",
    }


def reference_clonality(matrix_text: str, labels: dict[str, str]) -> dict:
    """How close each genome's nearest SAME-organism neighbour is -- i.e. how optimistic LOO is.

    THE MAIN THREAT TO THE IN-SET HEADLINE, quantified instead of caveated. Leave-one-out holds out a
    GENOME, not a lineage, so where a near-identical twin remains in the reference the fold is close
    to free. This repo has that distinction on record (isolate-disjoint vs lineage-disjoint); a number
    is more useful than the warning alone.
    """
    nearest: dict[str, float] = {}
    for h in parse_mash_dist(matrix_text):
        q = accession_from_reference_id(h.query_id)
        r = accession_from_reference_id(h.reference_id)
        if q == r or not labels.get(q) or labels.get(q) != labels.get(r):
            continue
        if h.distance < nearest.get(q, 9.0):
            nearest[q] = h.distance
    if not nearest:
        return {"n": 0}
    import statistics
    ds = sorted(nearest.values())
    per: dict[str, list[float]] = collections.defaultdict(list)
    for a, d in nearest.items():
        per[labels[a]].append(d)
    near_twin = sum(1 for d in ds if d < 0.005)
    return {
        "n": len(ds),
        "median_nearest_same_organism": round(statistics.median(ds), 5),
        "max_nearest_same_organism": round(ds[-1], 5),
        "n_with_near_identical_twin_lt_0.0005": sum(1 for d in ds if d < 0.0005),
        "n_with_twin_lt_0.005": near_twin,
        "fraction_with_twin_lt_0.005": round(near_twin / len(ds), 4),
        "per_organism_median": {o: round(statistics.median(v), 5) for o, v in sorted(per.items())},
        "interpretation": "leave-one-out is OPTIMISTIC in proportion to this: for the "
                          f"{round(near_twin / len(ds) * 100, 1)}% of genomes with a ~99.5%-ANI-or-"
                          "closer twin still in the reference, holding one out leaves a near-clone to "
                          "stand in for it. A lineage-disjoint split is the stronger test and was NOT "
                          "run.",
        "failure_direction": "within-organism pairs reach 0.11991, ABOVE max_distance 0.1036, so a "
                             "DIVERGENT member of a supported organism ABSTAINS rather than being "
                             "misidentified -- the safe failure direction, and a real coverage limit.",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=None)
    a = ap.parse_args(argv)

    manifest = a.manifest or sorted((REPO / "wiki").glob("identify_reference_manifest_*.json"))[-1]
    labels = labels_from_manifest(manifest)
    th = frozen()

    if not MATRIX_CACHE.exists():
        print(f"missing all-pairs matrix {MATRIX_CACHE}; run scripts/identify_threshold_analysis.py")
        return 2

    matrix = MATRIX_CACHE.read_text(encoding="utf-8", errors="replace")
    loo = leave_one_out(matrix, labels, th)
    clonality = reference_clonality(matrix, labels)
    oos = out_of_set_control(labels, th)

    today = date.today().isoformat()
    (REPO / "wiki" / f"identify_outofset_probe_{today}.json").write_text(
        json.dumps({"schema": "identify-outofset-probe-v1", "date": today,
                    "thresholds": {"max_distance": th.max_distance,
                                   "ambiguity_margin": th.ambiguity_margin},
                    **oos}, indent=2), encoding="utf-8")

    result = {
        "schema": "identify-validation-v1", "date": today, "manifest": manifest.name,
        "thresholds": {"max_distance": th.max_distance, "ambiguity_margin": th.ambiguity_margin},
        "leave_one_out": loo, "reference_clonality": clonality, "out_of_set_control": oos,
        "method": "leave-one-out by excluding each genome's OWN row from the cached all-pairs mash "
                  "matrix; no genome is ever scored against itself",
        "honest_limits": [
            "labels are each assembly's NCBI GenBank ORGANISM field -- metadata, not wet-lab "
            "identification",
            "closed-set: an organism outside the supported set must abstain, and is never identified",
            "species-level only; no strain or subspecies claim",
            "assembled genomes only (no reads, no mixtures)",
            "held out by GENOME, not by lineage -- near-identical reference genomes make this "
            "optimistic relative to a lineage-disjoint split",
        ],
    }

    # REFUSALS
    if oos.get("status") == "EMPTY":
        result["verdict"] = "NO_VERDICT_EMPTY_OUT_OF_SET"
        rc = 3
    elif oos.get("status") != "OK":
        result["verdict"] = f"OUT_OF_SET_CONTROL_FAILED_{oos.get('status')}"
        rc = 3
    elif oos["n_abstained"] == 0:
        result["verdict"] = "FAIL_ROUTER_NEVER_ABSTAINS"
        rc = 1
    else:
        result["verdict"] = "SCORED"
        rc = 0

    out = REPO / "wiki" / f"identify_validation_{today}.json"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    p = loo["pooled"]
    print(f"LEAVE-ONE-OUT  n={p['n_scored']}  hit={p['hit']}  miss={p['miss']}  "
          f"abstain={p['abstain']}  acc_on_called={p['accuracy_on_called']}")
    for o, c in loo["per_organism"].items():
        print(f"   {o:<28} hit={c.get('hit',0):<3} miss={c.get('miss',0):<3} "
              f"abstain={c.get('abstain',0):<3} unscorable={c.get('unscorable',0)}")
    if loo["confusions"]:
        print(f"   confusions: {loo['confusions']}")
    if loo["structurally_unscorable"]:
        print(f"   structurally unscorable (n<2 for that organism): {loo['structurally_unscorable']}")
    print(f"OUT-OF-SET     n={oos.get('n_genomes')}  abstained={oos.get('n_abstained')}  "
          f"rate={oos.get('abstention_rate')}")
    for s, c in (oos.get("per_species") or {}).items():
        print(f"   {s:<28} {dict(c)}")
    print(f"REFERENCE CLONALITY  median nearest-same={clonality.get('median_nearest_same_organism')}  "
          f"with ~99.5%-ANI twin={clonality.get('n_with_twin_lt_0.005')}/{clonality.get('n')} "
          f"({clonality.get('fraction_with_twin_lt_0.005')})  => LOO is optimistic by this much")
    print(f"VERDICT: {result['verdict']}   artifact: {out}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
