"""W0 census: measure the ENCODE substrate a regulatory-activity target would train on.

WHY THIS IS A W0 PROBE, NOT A DOWNLOADER. The 2026-05-12 pivot gave up an Enformer-shaped goal
trained on ENCODE / Roadmap multi-tissue data, and the compute half of that cut is now measured
(`wiki/encode_sizing_t4_2026-11-24.json`). The DATA half was never measured at all. This pins the
substrate's real cardinalities BEFORE anything is built against them, in the pattern of
`scripts/oxford_w0_probe.py` and `scripts/essentiality_ladder_w0_probe.py`: pure summarise helpers
split from the fetch, so the arithmetic unit-tests offline with no network.

THE DISTINCTION THIS EXISTS TO ENFORCE. An ENCODE "experiment" is not a trainable track. Measured
2026-11-24 on DNase-seq: 3,582 experiments against 8,854 released bigWig files, a ratio of 2.47x.
A census that reported the experiment count as the substrate size would be wrong in the FLATTERING
direction, which is why files and experiments are counted and labelled separately here and the
measured ratio ships beside both.

A CORRECTION THIS PROBE'S OWN MEASUREMENT FORCED. The plan that commissioned it asserted the two
counts "differ by roughly an order of magnitude". Measured, the DNase-seq ratio is 2.47x. The
plan's POINT stands -- the denominators are different and conflating them inflates the substrate --
but the magnitude it asserted was wrong, so the artifact records the measured ratio and this note
rather than the claim.

THE CHEAP ROUTE, verified live. `limit=0` returns `total` with an EMPTY `@graph`, and the response
still carries `facets` with per-term doc counts. So one request per assay yields the experiment
total, the distinct `biosample_ontology.term_name` count (the TISSUE axis, which is what
"multi-tissue" means) and the assembly values, with zero records paged. Walking 3,582 records to
learn a number the response already carries would be wasted bandwidth.

BYTES ARE AN ESTIMATE AND SAY SO. Summing `file_size` over every file would need full paging, so
this samples a bounded number of file sizes, takes the MEDIAN, and multiplies by the counted total.
The artifact labels that `estimated_total_bytes` with its method, sample size and the sample's
spread, so a reader can see how much to trust it. An estimate presented as a measurement is the
failure mode this wording exists to prevent.

Read-only. No credentials. Writes only the artifact passed to `--out`.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://www.encodeproject.org/search/"
OUT = ROOT / "wiki" / f"encode_substrate_census_{date.today().isoformat()}.json"

#: The assay families a tissue-specific regulatory-activity target needs, NAMED rather than implied
#: so the artifact records what was asked for and a later reader can see what was left out.
ASSAYS = ("DNase-seq", "ATAC-seq", "Histone ChIP-seq", "total RNA-seq")

#: The facet carrying the tissue axis. "Multi-tissue" is a claim about THIS number.
TISSUE_FACET = "biosample_ontology.term_name"
ASSEMBLY_FACET = "assembly"
BIGWIG_SAMPLE = 200          # bounded sample for the byte estimate
TIMEOUT = 45


class CensusError(Exception):
    """A condition under which no count may be reported."""


# --------------------------------------------------------------------------------------------
# PURE helpers -- no network, unit-tested offline
# --------------------------------------------------------------------------------------------
def facet_terms(payload: dict, field: str) -> list[dict]:
    """Return a facet's term list, or [] when the facet is absent."""
    for f in payload.get("facets") or []:
        if f.get("field") == field:
            return list(f.get("terms") or [])
    return []


def read_total(payload: dict, *, what: str) -> int:
    """REFUSE on a missing `total`.

    A census that silently reports zero is indistinguishable from a broken query, and zero is the
    answer that would quietly kill the goal. So an absent key raises rather than defaulting.
    """
    if "total" not in payload:
        raise CensusError(f"{what}: response carries no `total` key, so the query shape is wrong "
                          f"or the API changed. Refusing rather than reporting a count.")
    total = payload["total"]
    if not isinstance(total, int):
        raise CensusError(f"{what}: `total` is {total!r}, not an int")
    return total


def summarise_sizes(sizes: list[int]) -> dict:
    """Median-based volume estimate, with the spread so its reliability is visible."""
    clean = sorted(s for s in sizes if isinstance(s, int) and s > 0)
    if not clean:
        return {"n_sampled": 0, "median_bytes": None, "note": "no file_size values sampled"}
    return {
        "n_sampled": len(clean),
        "median_bytes": int(statistics.median(clean)),
        "min_bytes": clean[0],
        "max_bytes": clean[-1],
        "p25_bytes": int(clean[len(clean) // 4]),
        "p75_bytes": int(clean[(3 * len(clean)) // 4]),
    }


def estimate_volume(n_files: int, size_stats: dict) -> dict:
    """An ESTIMATE, labelled. n_files is counted; the per-file size is sampled."""
    med = size_stats.get("median_bytes")
    if med is None or n_files <= 0:
        return {"estimated_total_bytes": None, "estimated_total_gib": None,
                "method": "unavailable -- no sampled sizes or no files"}
    total = med * n_files
    return {
        "estimated_total_bytes": total,
        "estimated_total_gib": round(total / 2**30, 1),
        "method": (f"median of {size_stats['n_sampled']} sampled file_size values "
                   f"x {n_files} counted files. The COUNT is measured; the per-file size is "
                   f"SAMPLED, so this is an estimate and not a measurement."),
    }


def census_verdict(per_assay: dict) -> dict:
    """Roll the per-assay rows up without inventing a sufficiency claim."""
    exp = sum(a["n_experiments"] for a in per_assay.values())
    files = sum(a["n_bigwig_files"] for a in per_assay.values())
    tissues = {t for a in per_assay.values() for t in a["tissue_terms_sample"]}
    gib = [a["volume"]["estimated_total_gib"] for a in per_assay.values()
           if a["volume"]["estimated_total_gib"] is not None]
    return {
        "n_assays": len(per_assay),
        "total_experiments": exp,
        "total_bigwig_files": files,
        "file_to_experiment_ratio": round(files / exp, 2) if exp else None,
        "distinct_tissue_terms_union_of_samples": len(tissues),
        "estimated_total_gib_all_assays": round(sum(gib), 1) if gib else None,
        "what_this_does_NOT_say": (
            "nothing about whether this substrate SUFFICES for the goal, nothing about label "
            "quality, and nothing about whether the volume is affordable here. It is a count."),
    }


# --------------------------------------------------------------------------------------------
# fetch
# --------------------------------------------------------------------------------------------
def get(params: dict) -> dict:
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "curl"})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            return json.loads(r.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        raise CensusError(f"fetch failed for {url}: {e}") from e


def census_one(assay: str) -> dict:
    exp = get({"type": "Experiment", "assay_title": assay, "limit": 0, "format": "json"})
    n_exp = read_total(exp, what=f"{assay} experiments")
    if n_exp == 0:
        raise CensusError(f"{assay}: zero experiments. Either the assay_title is wrong or the "
                          f"substrate is absent; refusing to record a zero as a finding.")

    tissues = facet_terms(exp, TISSUE_FACET)
    assemblies = [t.get("key") for t in facet_terms(exp, ASSEMBLY_FACET)]

    fil = get({"type": "File", "file_format": "bigWig", "assay_title": assay,
               "status": "released", "limit": 0, "format": "json"})
    n_files = read_total(fil, what=f"{assay} bigWig files")

    sizes: list[int] = []
    if n_files:
        samp = get({"type": "File", "file_format": "bigWig", "assay_title": assay,
                    "status": "released", "field": "file_size",
                    "limit": BIGWIG_SAMPLE, "format": "json"})
        sizes = [r.get("file_size") for r in (samp.get("@graph") or [])
                 if isinstance(r.get("file_size"), int)]
    stats = summarise_sizes(sizes)

    return {
        "assay_title": assay,
        "n_experiments": n_exp,
        "n_bigwig_files": n_files,
        "file_to_experiment_ratio": round(n_files / n_exp, 2) if n_exp else None,
        "n_distinct_tissue_terms": len(tissues),
        "tissue_terms_sample": [t.get("key") for t in
                                sorted(tissues, key=lambda x: -x.get("doc_count", 0))[:25]],
        "assemblies": assemblies,
        "file_size_sample": stats,
        "volume": estimate_volume(n_files, stats),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--assays", nargs="*", default=list(ASSAYS))
    a = ap.parse_args(argv)

    per: dict[str, dict] = {}
    try:
        for assay in a.assays:
            row = census_one(assay)
            per[assay] = row
            print(f"{assay:22s} exp={row['n_experiments']:>6,} "
                  f"bigwig={row['n_bigwig_files']:>7,} "
                  f"ratio={row['file_to_experiment_ratio']} "
                  f"tissues={row['n_distinct_tissue_terms']:>4} "
                  f"~{row['volume']['estimated_total_gib']} GiB", flush=True)
    except CensusError as e:
        print(f"REFUSING: {e}", file=sys.stderr)
        return 3

    rec = {
        "record": "encode-substrate-census-v1",
        "date": date.today().isoformat(),
        "question": ("what multi-tissue regulatory substrate does ENCODE actually expose, counted "
                     "rather than estimated, for the goal abandoned 2026-05-12"),
        "api": {"base": BASE, "verified_live": True,
                "cheap_route": "limit=0 returns `total` with an empty @graph while still carrying "
                               "`facets`, so one request per assay needs zero record paging"},
        "per_assay": per,
        "rollup": census_verdict(per),
        "honest_limits": [
            "Experiment counts and FILE counts are different denominators. Both are recorded with "
            "their measured ratio; neither may be quoted as the other.",
            "The volume figure is an ESTIMATE: the file count is measured, the per-file size is a "
            "bounded median sample. Its spread is recorded so its reliability is visible.",
            "Only ENCODE was probed. The Roadmap Epigenomics download surface is UNVERIFIED and is "
            "named in the commissioning plan as such -- this artifact does not speak for it.",
            "A count is not a label-quality judgement and not a feasibility verdict. The gate "
            "screen and the sizing memo are separate artifacts for separate questions.",
        ],
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    r = rec["rollup"]
    print(f"\nrollup: {r['total_experiments']:,} experiments, {r['total_bigwig_files']:,} bigWig "
          f"files (ratio {r['file_to_experiment_ratio']}), "
          f"~{r['estimated_total_gib_all_assays']} GiB estimated")
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
