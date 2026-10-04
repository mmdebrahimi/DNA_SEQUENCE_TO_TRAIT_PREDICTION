"""Measure the within- vs between-organism Mash distance distributions, THEN freeze the thresholds.

    uv run python scripts/identify_threshold_analysis.py

ORDER IS THE POINT. The thresholds are written only AFTER this analysis, never before. This repo has
TWO recorded cases where a pre-registered bar was itself the error -- a NET ceiling used to police a
DIRECTIONAL failure, and a `+4 hits` bar that was unachievable by construction -- and the lesson drawn
was *do the mechanism analysis, THEN freeze*, not *do not freeze*.

ONE CONTAINER CALL. `mash dist reference.msh reference.msh` yields the full all-pairs matrix, which
serves BOTH this analysis and Step 6's leave-one-out (where holding a genome out is just excluding its
own row). The plan's per-genome rebuild would have been ~212 container spin-ups for identical
information. The matrix is cached on D: so Step 6 does not recompute it.

What this does NOT do: choose thresholds to make a number look good. If the distributions overlap, that
is the finding, and the router ships abstention-heavy.
"""
from __future__ import annotations

import argparse
import json
import statistics
from datetime import date
from pathlib import Path

from dna_decode.identify.core import accession_from_reference_id, parse_mash_dist

REPO = Path(__file__).resolve().parent.parent
HOST_DRIVE = Path("D:/")
CONTAINER_MOUNT = "/d"
REFERENCE = Path("D:/dna_decode_cache/identify/reference.msh")
MATRIX_CACHE = Path("D:/dna_decode_cache/identify/allpairs.tsv")
MASH_IMAGE = "quay.io/biocontainers/mash:2.3--hb105d93_10"


def _container(p: Path) -> str:
    return CONTAINER_MOUNT + "/" + p.resolve().relative_to(HOST_DRIVE.resolve()).as_posix()


def load_labels(manifest: Path) -> dict[str, str]:
    """accession -> canonical organism, from the committed reference manifest."""
    m = json.loads(manifest.read_text(encoding="utf-8"))
    return {acc: org for org, accs in m["accessions_by_organism"].items() for acc in accs}


def all_pairs(reference: Path = REFERENCE, cache: Path = MATRIX_CACHE,
              *, refresh: bool = False) -> str:
    """One `mash dist <ref> <ref>` over the whole sketch, cached as raw TSV on D:."""
    if cache.exists() and not refresh and cache.stat().st_size > 0:
        return cache.read_text(encoding="utf-8", errors="replace")
    from tools.docker_runner import run as docker_run

    res = docker_run(
        MASH_IMAGE,
        ["mash", "dist", _container(reference), _container(reference)],
        mounts={str(HOST_DRIVE): CONTAINER_MOUNT},
        capture_output=True, check=True, timeout=3600,
    )
    out = res.stdout or ""
    if not out.strip():
        raise RuntimeError("mash dist produced no output; exit 0 is not evidence")
    cache.write_text(out, encoding="utf-8", newline="\n")
    return out


def distributions(stdout: str, labels: dict[str, str]) -> dict:
    """Split every off-diagonal pair into within-organism vs between-organism distances."""
    within: list[float] = []
    between: list[float] = []
    # nearest BETWEEN-organism neighbour per genome: the quantity max_distance must sit below
    nearest_between: dict[str, float] = {}
    nearest_within: dict[str, float] = {}

    for h in parse_mash_dist(stdout):
        a = accession_from_reference_id(h.reference_id)
        b = accession_from_reference_id(h.query_id)
        if a == b:
            continue                                      # self-pair, the matrix diagonal
        oa, ob = labels.get(a), labels.get(b)
        if not oa or not ob:
            continue
        if oa == ob:
            within.append(h.distance)
            if h.distance < nearest_within.get(b, 2.0):
                nearest_within[b] = h.distance
        else:
            between.append(h.distance)
            if h.distance < nearest_between.get(b, 2.0):
                nearest_between[b] = h.distance

    def q(xs: list[float]) -> dict:
        if not xs:
            return {"n": 0}
        s = sorted(xs)
        return {
            "n": len(s), "min": s[0], "p05": s[int(0.05 * (len(s) - 1))],
            "median": statistics.median(s), "p95": s[int(0.95 * (len(s) - 1))], "max": s[-1],
        }

    nb = sorted(nearest_between.values())
    nw = sorted(nearest_within.values())
    return {
        "within_organism": q(within),
        "between_organism": q(between),
        "nearest_same_organism_neighbour": q(nw),
        "nearest_different_organism_neighbour": q(nb),
        "separation": {
            "within_p95": q(within).get("p95"),
            "between_p05": q(between).get("p05"),
            "overlaps": (q(within).get("p95") is not None and q(between).get("p05") is not None
                         and q(within)["p95"] >= q(between)["p05"]),
            "worst_case_nearest_within": nw[-1] if nw else None,
            "closest_cross_organism_pair": nb[0] if nb else None,
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    manifest = a.manifest or sorted((REPO / "wiki").glob("identify_reference_manifest_*.json"))[-1]
    labels = load_labels(manifest)
    stdout = all_pairs(refresh=a.refresh)
    d = distributions(stdout, labels)
    d["schema"] = "identify-distance-distribution-v1"
    d["date"] = date.today().isoformat()
    d["manifest"] = manifest.name
    d["n_labelled_accessions"] = len(labels)
    d["matrix_cache"] = str(MATRIX_CACHE)

    out = a.out or (REPO / "wiki" / f"identify_distance_distribution_{date.today().isoformat()}.json")
    out.write_text(json.dumps(d, indent=2), encoding="utf-8")

    w, b = d["within_organism"], d["between_organism"]
    print(f"pairs: within={w['n']}  between={b['n']}  (accessions labelled: {len(labels)})")
    print(f"within-organism   median={w.get('median'):.5f}  p95={w.get('p95'):.5f}  max={w.get('max'):.5f}")
    print(f"between-organism  p05={b.get('p05'):.5f}  median={b.get('median'):.5f}")
    nw = d["nearest_same_organism_neighbour"]; nb = d["nearest_different_organism_neighbour"]
    print(f"nearest SAME-organism neighbour : median={nw.get('median'):.5f}  max={nw.get('max'):.5f}")
    print(f"nearest OTHER-organism neighbour: min={nb.get('min'):.5f}  median={nb.get('median'):.5f}")
    print(f"distributions overlap: {d['separation']['overlaps']}")
    print(f"artifact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
