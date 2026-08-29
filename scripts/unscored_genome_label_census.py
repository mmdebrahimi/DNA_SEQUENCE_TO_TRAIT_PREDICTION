"""How much FREE validation surface is sitting in already-cached genomes that were never scored?

THE OPPORTUNITY. This repo has ~1,800 genomes with cached AMRFinder output, and ~950 of them have never
appeared in any labelled cohort. AMRFinder is the expensive step (~95 s/genome, Docker) and it is already
paid for. If NCBI Pathogen Detection carries an `AST_phenotypes` call for those genomes, each one is a
scoreable isolate at zero marginal cost.

WHAT THIS SCRIPT DOES: measures the size of that opportunity, per (organism, drug). It does NOT score
anything and does NOT touch the frozen surface. Sizing first is the point -- if the answer is "a handful",
the honest move is to say the door is closed rather than build a harness for it.

THE LEAKAGE GATE IS THE MANIFEST, AND THE WEAKER CHECK WAS NOT ENOUGH -- MEASURED. The first version of
this census filtered on "never appeared in a selected.tsv". Run against the real gate
(`cohort_manifest.prior_accessions`, which also scans the parquet cohorts) that pool of 956 turned out to
contain **645 leaked accessions** -- two thirds of it. The disjoint pool is 311. The manifest is now the
primary filter, and the selected.tsv set is kept only as a cheap pre-filter.

This is the hardcoded-vs-derived trap in a new costume: a hand-rolled exclusion check beside the data will
under-cover, and the fail-closed manifest exists precisely so nobody has to re-derive it.

Run: uv run python scripts/unscored_genome_label_census.py
"""
from __future__ import annotations

import csv
import glob
import gzip
import io
import json
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

GROUPS = ("Klebsiella", "Escherichia_coli_Shigella", "Acinetobacter", "Salmonella")
DRUGS = ("ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline", "meropenem")


def already_labelled() -> set[str]:
    """Accessions that appear in ANY cohort selected.tsv with an R/S label."""
    out: set[str] = set()
    for sel in glob.glob(str(ROOT / "data" / "raw" / "*" / "selected.tsv")):
        for line in Path(sel).read_text(encoding="utf-8", errors="replace").splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1].strip().upper() in ("R", "S"):
                out.add(parts[0].strip())
    return out


def stream_group(group: str, wanted: set[str]) -> dict[str, dict[str, str]]:
    """{asm_acc: {drug: R/S}} for wanted accessions in one PD group. Network."""
    from ncbi_pd_provenance_census import latest_metadata_url
    from dna_decode.data.pd_ast import parse_ast_phenotypes

    url = latest_metadata_url(group)
    found: dict[str, dict[str, str]] = {}
    with urllib.request.urlopen(url, timeout=900) as resp:
        raw = gzip.GzipFile(fileobj=resp) if url.endswith(".gz") else resp
        rd = csv.reader(io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), delimiter="\t")
        header = next(rd, None)
        if not header:
            return found
        idx = {c: i for i, c in enumerate(header)}
        ai, pi = idx.get("asm_acc"), idx.get("AST_phenotypes")
        if ai is None or pi is None:
            return found
        for row in rd:
            if len(row) <= max(ai, pi):
                continue
            acc = row[ai].strip()
            if acc not in wanted:
                continue
            calls = {d: v for d, v in parse_ast_phenotypes(row[pi], set(DRUGS)).items()
                     if v.upper() in ("R", "S")}
            if calls:
                found[acc] = calls
    return found


def main() -> int:
    from gentamicin_rmt_candidate import amrfinder_index

    from dna_decode.eval.cohort_manifest import build_manifest, prior_accessions

    cached = set(amrfinder_index())
    labelled = already_labelled()
    manifest = build_manifest()
    if getattr(manifest, "incomplete", False):
        # An incomplete manifest cannot prove disjointness. Same fail-closed rule the
        # provenance-disjoint scorer uses: refuse rather than report false independence.
        print("INCOMPLETE_MANIFEST: cannot establish disjointness; refusing to report a pool.")
        return 2
    prior = prior_accessions(manifest, exclude_cohort="__none__")
    unscored = (cached - labelled) - prior
    print(f"cached AMRFinder genomes: {len(cached)}")
    print(f"  already in a labelled cohort: {len(cached & labelled)}")
    print(f"  additionally caught by the accession manifest: {len((cached - labelled) & prior)}")
    print(f"  DISJOINT candidate pool: {len(unscored)}\n")

    per_group: dict[str, int] = {}
    errors: dict[str, str] = {}
    hits: dict[str, dict[str, str]] = {}
    group_of: dict[str, str] = {}
    for g in GROUPS:
        try:
            got = stream_group(g, unscored)
        except Exception as exc:
            errors[g] = f"{type(exc).__name__}: {exc}"
            print(f"  {g}: FAILED ({errors[g]})")
            continue
        per_group[g] = len(got)
        for acc, calls in got.items():
            hits[acc] = calls
            group_of[acc] = g
        print(f"  {g}: {len(got)} previously-unscored genome(s) carry an AST call")

    cells: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for acc, calls in hits.items():
        for drug, v in calls.items():
            cells[(group_of[acc], drug)][v.upper()] += 1

    print()
    if errors:
        print(f"INCOMPLETE: {len(errors)} group(s) failed -- every count below is a partial floor.\n")
    print(f"previously-unscored genomes with at least one AST call: {len(hits)}")
    print()
    print(f"{'organism':28}{'drug':16}{'R':>5}{'S':>5}   both-classes >= 20?")
    powered = []
    for (g, d), c in sorted(cells.items(), key=lambda kv: -(kv[1]['R'] + kv[1]['S'])):
        ok = c["R"] >= 20 and c["S"] >= 20
        if ok:
            powered.append((g, d, c["R"], c["S"]))
        print(f"{g:28}{d:16}{c['R']:>5}{c['S']:>5}   {'YES' if ok else 'no'}")

    print()
    if powered:
        print(f"=> {len(powered)} (organism, drug) cell(s) could be scored at ZERO marginal AMRFinder cost:")
        for g, d, r, s in powered:
            print(f"     {g} x {d}  ({r}R/{s}S)")
        print("   These are CANDIDATE pools. Anything actually scored must first pass the")
        print("   accession-manifest leakage gate (cohort_manifest.prior_accessions), which is a")
        print("   stricter check than 'never appeared in a selected.tsv'.")
    elif errors:
        print("=> NO CONCLUSION: the sweep did not complete, so an empty result means the sweep failed.")
    else:
        print("=> No cell reaches 20R/20S. The cached-genome pool does NOT open new validation surface,")
        print("   and that is a measured door-closed rather than an assumption.")

    (ROOT / "wiki" / "unscored_genome_label_census.json").write_text(json.dumps(
        {"n_cached": len(cached), "n_unscored": len(unscored), "n_with_ast": len(hits),
         "per_group": per_group, "errors": errors, "complete": not errors,
         "cells": {f"{g}|{d}": dict(c) for (g, d), c in cells.items()},
         # persist per-accession labels so the pool can be SCORED offline without re-streaming
         "labels": {a: {"group": group_of[a], "calls": c} for a, c in hits.items()},
         "powered_cells": [{"organism": g, "drug": d, "R": r, "S": s} for g, d, r, s in powered]},
        indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
