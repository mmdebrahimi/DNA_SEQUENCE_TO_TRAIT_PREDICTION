"""Score sourmash FracMinHash clustering against the FROZEN Mash crosscheck bar.

WHAT IS COMPARED, AND WHY IT IS PARTITIONS
------------------------------------------
This compares PARTITIONS, not distance matrices. Comparing raw distances would
require running Mash -- the blocked thing this whole route exists to work around.
The committed `wiki/provdisjoint_lineage_metrics.json` already holds Mash's product
(per-threshold effective-lineage counts), so re-clustering the same assemblies with
sourmash and comparing the resulting counts and G8 verdicts is the comparison that
is actually available. That resolves the circularity rather than hiding it.

EXIT CODES
----------
  0  every clause met
  1  a clause failed  (a real finding: sourmash does NOT substitute for Mash)
  2  INSUFFICIENT_DATA -- a cohort's assemblies could not all be resolved
  3  sourmash unavailable
  4  the bar refused to load (not frozen, or its derived-from table has changed)

INSUFFICIENT_DATA IS A VALID RESULT, NOT A FAILURE. A cohort whose genomes do not
all resolve -- or whose labels parse to a single class, which every reference cohort
rules out -- yields NO agreement number -- the same incompleteness rule the G8 bar
commits to. A partial cohort's effective-lineage counts are lower bounds, and a
lower bound cannot be compared for EXACT EQUALITY against clause 2.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval import rejection_gates  # noqa: E402

BAR_PATH = Path("wiki/sourmash_vs_mash_crosscheck_bar.json")
OUT_STEM = Path("wiki/sourmash_vs_mash_crosscheck_2026-10-19")
COHORT_ROOT = Path("data/raw")
REFSEQ_CACHE_CANDIDATES = (
    Path("data/refseq_cache"),
    Path("D:/dna_decode_cache/refseq"),
)
THRESHOLDS = ("0.001", "0.005")

OK = "OK"
INSUFFICIENT_DATA = "INSUFFICIENT_DATA"
SOURMASH_UNAVAILABLE = "SOURMASH_UNAVAILABLE"


class BarRefused(RuntimeError):
    """The frozen bar is not usable as a commitment."""


@dataclass
class ClauseOutcome:
    clause_id: int
    name: str
    met: bool | None
    detail: str
    failures: list = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"clause_id": self.clause_id, "name": self.name, "met": self.met,
                "detail": self.detail, "failures": self.failures}


# --------------------------------------------------------------------------- #
# the bar
# --------------------------------------------------------------------------- #

def eol_renderings(raw: bytes) -> list[tuple[str, bytes]]:
    """Every line-ending rendering of `raw`, as `[(label, bytes)]`.

    `raw` first, so an exact match never pays for a conversion. CRLF is produced by
    normalising to LF and then expanding -- never by expanding `raw` directly, which
    would turn an existing CRLF into CRCRLF. A rendering identical to `raw` is
    dropped, so the list carries no duplicate digests.
    """
    lf = raw.replace(b"\r\n", b"\n")
    crlf = lf.replace(b"\n", b"\r\n")
    out: list[tuple[str, bytes]] = [("raw", raw)]
    if lf != raw:
        out.append(("lf", lf))
    if crlf != raw:
        out.append(("crlf", crlf))
    return out


def verify_derived_from(raw: bytes, pinned: str | None) -> str | None:
    """The rendering label whose sha256 equals `pinned`, else None.

    PURE -- no filesystem, so the CRLF/LF case is testable on a host that has only
    one of them. This verifies CONTENT and deliberately does not verify LINE
    ENDINGS: a pinned digest taken on a CRLF checkout is a fact about that
    checkout's bytes, not about the table's content, and the table carries no
    `.gitattributes` `eol` rule that would make the two the same thing.

    It is NOT a weakening. A table whose content changed has a different digest
    under EVERY rendering, so the only thing newly accepted is the SAME content
    written with different line endings. A pin may itself be in either form, which
    is why both directions are tried rather than only normalising the live side.
    """
    if not pinned:
        return None
    for label, data in eol_renderings(raw):
        if hashlib.sha256(data).hexdigest() == pinned:
            return label
    return None


def load_bar(bar_path: Path, *, repo: Path = REPO) -> dict:
    """Load the frozen bar, REFUSING anything that is not a commitment.

    Two refusals, both load-bearing:
      * `frozen_before_any_agreement_number_existed` must be true -- a bar written
        after a number is not a pre-commitment.
      * the `derived_from` sha256 must still match the live lineage table. A bar
        derived from a since-changed table is not a commitment ABOUT that table,
        and silently scoring against it would be the stale-pin failure this repo
        has already recorded for retired prospective locks.

    The second refusal is EOL-TOLERANT, and that is a fix, not a loosening. The
    pinned digest is of the table as it stood on a CRLF checkout, and
    `wiki/provdisjoint_lineage_metrics.json` is NOT covered by any `.gitattributes`
    `eol` rule (measured: `text: unspecified`, `eol: unspecified`, 715 CRLF). So on
    an LF checkout -- CI, Linux, a clone with different `core.autocrlf` -- the raw
    digest is `42f944cc...` against a pinned `785dcc8a...` and the bar would refuse
    with "the table has CHANGED" when not one character had. See
    `verify_derived_from`: the comparison is over EOL RENDERINGS of the live bytes,
    so content is still verified exactly while line endings are not. The frozen bar
    is NOT edited to fix this -- a pre-commitment is not rewritten to make a later
    run pass.
    """
    bar = json.loads(Path(bar_path).read_text(encoding="utf-8"))
    if not bar.get("frozen_before_any_agreement_number_existed"):
        raise BarRefused(
            f"{bar_path} does not claim to predate any agreement number; it cannot "
            "be used as an acceptance bar"
        )
    derived = bar.get("derived_from") or {}
    src = repo / derived.get("path", "")
    if not src.exists():
        raise BarRefused(f"the bar's derived_from table is missing: {src}")
    raw = src.read_bytes()
    pinned = derived.get("sha256")
    match = verify_derived_from(raw, pinned)
    if match is None:
        actual = hashlib.sha256(raw).hexdigest()
        raise BarRefused(
            f"{derived.get('path')} has CHANGED since the bar was frozen "
            f"(sha256 {actual[:12]}... != pinned {str(pinned)[:12]}..., and no "
            "line-ending rendering of the live bytes matches either). "
            "The bar is no longer a commitment about this table; re-freeze deliberately."
        )
    # Surfaced, never silent: an EOL-normalised match means this checkout's line
    # endings differ from the freeze. The CONTENT is verified; saying so lets a
    # reader tell a clean verify from a renormalised one.
    bar["_derived_from_eol_match"] = match
    return bar


def mash_reference(bar: dict) -> dict[str, dict]:
    """`{pair_label: {...}}` from the bar's own reference table."""
    out = {}
    for row in bar["mash_reference_table"]:
        label = f"{row['organism']} x {row['drug']} @ {row['threshold']}"
        out[label] = row
    return out


def near_bar_labels(bar: dict) -> set[str]:
    clause = next(c for c in bar["clauses"] if c["id"] == 2)
    return {p["label"] for p in clause["near_bar_pairs"]}


def trip_labels(bar: dict) -> set[str]:
    clause = next(c for c in bar["clauses"] if c["id"] == 3)
    return {p["label"] for p in clause["trip_pairs"]}


# --------------------------------------------------------------------------- #
# clause evaluation -- PURE, offline-testable
# --------------------------------------------------------------------------- #

def evaluate_clauses(bar: dict, observed: dict[str, dict]) -> list[ClauseOutcome]:
    """Score `{pair_label: {"minority": int, "verdict": str}}` against the bar.

    NOTE ON CLAUSE 3 (see the bar's own `clause_logic` block): clause 3 is LOGICALLY
    IMPLIED by clause 1. If a TRIP pair flips to PASS, verdict identity fails too, so
    clause 3 cannot fail while clause 1 holds. It is evaluated and reported
    separately anyway because it is a SEVERITY pre-commitment -- it names the most
    dangerous direction of error in advance and forces the TRIP outcome to be stated
    rather than averaged into a verdict-identity rate.
    """
    reference = mash_reference(bar)
    outcomes: list[ClauseOutcome] = []

    # clause 1 -- verdict identity on every pair
    mismatches = []
    for label, row in reference.items():
        obs = observed.get(label)
        if obs is None:
            mismatches.append({"pair": label, "reason": "not observed"})
            continue
        if obs["verdict"] != row["mash_verdict"]:
            mismatches.append({"pair": label, "mash": row["mash_verdict"],
                               "sourmash": obs["verdict"]})
    outcomes.append(ClauseOutcome(
        1, "verdict_identity", not mismatches,
        f"{len(reference) - len(mismatches)}/{len(reference)} verdicts identical "
        "(NECESSARY BUT NEARLY VACUOUS -- 19/20 PASS by margin alone)",
        mismatches,
    ))

    # clause 2 -- exact minority count on the near-bar pairs
    near = near_bar_labels(bar)
    off = []
    for label in sorted(near):
        obs = observed.get(label)
        if obs is None:
            off.append({"pair": label, "reason": "not observed"})
            continue
        if obs["minority"] != reference[label]["minority"]:
            off.append({"pair": label, "mash": reference[label]["minority"],
                        "sourmash": obs["minority"]})
    outcomes.append(ClauseOutcome(
        2, "exact_minority_count_near_the_bar", not off,
        f"{len(near) - len(off)}/{len(near)} near-bar pairs matched exactly "
        "(DISCRIMINATING)",
        off,
    ))

    # clause 3 -- the TRIP must reproduce
    trips = trip_labels(bar)
    lost = []
    for label in sorted(trips):
        obs = observed.get(label)
        if obs is None:
            lost.append({"pair": label, "reason": "not observed"})
        elif obs["verdict"] != rejection_gates.TRIP:
            lost.append({"pair": label, "sourmash": obs["verdict"],
                         "severity": "MOST DANGEROUS -- manufactures lineage diversity"})
    outcomes.append(ClauseOutcome(
        3, "the_trip_must_reproduce", not lost,
        f"{len(trips) - len(lost)}/{len(trips)} TRIP pair(s) reproduced "
        "(implied by clause 1; retained as a severity pre-commitment)",
        lost,
    ))
    return outcomes


# --------------------------------------------------------------------------- #
# cohort resolution
# --------------------------------------------------------------------------- #

def read_selected(path: Path) -> dict[str, int]:
    """`{accession: 1 if R else 0}` from a cohort's bare `selected.tsv`.

    DELEGATES to `scripts/compute_lineage_metrics.read_selected` -- the function that
    PRODUCED the Mash effective-lineage numbers in
    `wiki/provdisjoint_lineage_metrics.json` this bar is derived from. The labels must
    come from the same definition on both sides of the comparison, for exactly the
    reason the bar computes its verdicts THROUGH the imported gate instead of
    restating `MIN_EFFECTIVE_LINEAGES`: a second copy drifts, and a drift here is not
    cosmetic -- it silently moves the quantity being compared.

    `selected.tsv` IS HEADERLESS -- bare `<accession>\t<label>`, exact `R` for
    resistant. That is asserted by the frozen `eval/cohort_manifest.py`, by the
    producer above, and by `tests/test_ktype_kaptive_concordance.py`, whose docstring
    records a real past bug from assuming otherwise ("a DictReader silently eats the
    first genome").

    An earlier version of this function parsed a HEADER and looked for a
    `phenotype`/`label`/`resistance` column. On the real format that (a) consumed the
    first accession as a header row and (b) found no phenotype column, so every
    isolate fell through to susceptible -- `effective_lineage_N_R = 0`, minority 0, and
    G8 TRIPs on all 20 pairs. The script would then have reported its own documented
    exit-1 "a real finding: sourmash does NOT substitute for Mash" when the only thing
    measured was a label parser disagreeing with the producer of the numbers it was
    being compared against.
    """
    from scripts.compute_lineage_metrics import read_selected as canonical

    return canonical(Path(path))


def find_refseq_cache(repo: Path = REPO) -> Path | None:
    for candidate in REFSEQ_CACHE_CANDIDATES:
        path = candidate if candidate.is_absolute() else repo / candidate
        if path.is_dir():
            return path
    return None


def resolve_genomes(accessions: dict[str, int], cache: Path | None) -> tuple[dict, list]:
    """`({accession: fasta_path}, missing)`."""
    found: dict[str, Path] = {}
    missing: list[str] = []
    for acc in accessions:
        hit = None
        if cache is not None:
            for name in ("genome.fna", f"{acc}.fna", f"{acc}_genomic.fna"):
                candidate = cache / acc / name
                if candidate.is_file():
                    hit = candidate
                    break
        if hit is None:
            missing.append(acc)
        else:
            found[acc] = hit
    return found, missing


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bar", default=str(BAR_PATH))
    ap.add_argument("--cohort-root", default=str(COHORT_ROOT))
    ap.add_argument("--out-stem", default=str(OUT_STEM))
    ap.add_argument("--ksize", type=int, default=21)
    ap.add_argument("--scaled", type=int, default=1000)
    args = ap.parse_args(argv)

    try:
        bar = load_bar(REPO / args.bar)
    except BarRefused as exc:
        print(f"BAR REFUSED: {exc}")
        return 4

    reference = mash_reference(bar)
    cohort_root = REPO / args.cohort_root
    cache = find_refseq_cache()

    cells: list[dict] = []
    observed: dict[str, dict] = {}
    blocked: list[dict] = []

    try:
        from dna_decode.eval.sourmash_distance import (
            SourmashUnavailable, compute_sourmash_distances,
        )
        from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
        from dna_decode.eval.cabbage_lineage import cluster_structure
    except ImportError as exc:
        print(f"SOURMASH_UNAVAILABLE: {exc}")
        return 3

    lineage = json.loads(
        (REPO / bar["derived_from"]["path"]).read_text(encoding="utf-8")
    )

    for cell in lineage["cells"]:
        cohort = cell["cohort"]
        cohort_dir = cohort_root / cohort
        selected = cohort_dir / "selected.tsv"
        if not selected.is_file():
            blocked.append({"cohort": cohort, "reason": "no selected.tsv",
                            "path": str(selected).replace("\\", "/")})
            cells.append({"cohort": cohort, "status": INSUFFICIENT_DATA,
                          "reason": "cohort manifest (selected.tsv) not present",
                          "n_genomes_missing": None})
            continue

        labels = read_selected(selected)
        # DEFENCE IN DEPTH against the label-parse failure recorded in `read_selected`.
        # Every one of the ten reference cohorts carries BOTH classes (across all 20
        # pairs the minimum effective-lineage counts are N_R 2 and N_S 7), so a raw
        # parse yielding a single label class cannot be biology -- it is a parse that
        # did not match the file. Left un-refused it is invisible AND maximally
        # damaging: zero R means minority 0, G8 TRIPs on every pair, and the run
        # reports the project's headline negative. Same shape as the TB caller
        # refusing a VCF that holds data rows yet parses to zero calls, rather than
        # printing a confident S -- both states otherwise look identical downstream.
        if len(set(labels.values())) < 2:
            cells.append({
                "cohort": cohort, "status": INSUFFICIENT_DATA,
                "reason": "label parse yielded a single class",
                "n_accessions": len(labels),
                "n_resistant": sum(labels.values()),
                "n_genomes_missing": None,
                "note": ("every reference cohort has both an R and an S class, so a "
                         "single-class parse is a parser/format mismatch, not a "
                         "collapsed cohort -- scoring it would compare a parse "
                         "failure against the bar"),
            })
            blocked.append({"cohort": cohort,
                            "reason": "label parse yielded a single class",
                            "n_accessions": len(labels),
                            "n_resistant": sum(labels.values())})
            continue
        genomes, missing = resolve_genomes(labels, cache)
        if missing or not genomes:
            cells.append({
                "cohort": cohort, "status": INSUFFICIENT_DATA,
                "reason": "assemblies did not all resolve",
                "n_accessions": len(labels),
                "n_genomes_missing": len(missing),
                "note": ("a partial cohort yields LOWER-BOUND effective-lineage "
                         "counts, which cannot be compared for exact equality"),
            })
            blocked.append({"cohort": cohort, "reason": "missing assemblies",
                            "n_missing": len(missing)})
            continue

        try:
            dm = compute_sourmash_distances(genomes, ksize=args.ksize,
                                            scaled=args.scaled)
        except SourmashUnavailable as exc:
            print(f"SOURMASH_UNAVAILABLE: {exc}")
            return 3

        per_threshold = {}
        for threshold in THRESHOLDS:
            clusters = greedy_representative_clusters_from_matrix(dm, float(threshold))
            structure = cluster_structure(clusters, labels)
            minority = min(structure["effective_lineage_N_R"],
                           structure["effective_lineage_N_S"])
            verdict = rejection_gates.g8_dedup_collapses_balance(
                {"min_effective_lineages": minority}
            ).verdict
            label = f"{cell['organism']} x {cell['drug']} @ {threshold}"
            observed[label] = {"minority": minority, "verdict": verdict}
            per_threshold[threshold] = {
                "minority": minority, "verdict": verdict,
                # reported, never gated:
                "n_clusters": structure["n_clusters"],
                "largest_cluster_fraction": structure["largest_cluster_fraction"],
                "discordant_isolate_fraction": structure["discordant_isolate_fraction"],
                "effective_lineage_N_R": structure["effective_lineage_N_R"],
                "effective_lineage_N_S": structure["effective_lineage_N_S"],
            }
        cells.append({"cohort": cohort, "status": OK, "thresholds": per_threshold})

    if blocked:
        status, exit_code = INSUFFICIENT_DATA, 2
        clauses = [ClauseOutcome(c["id"], c["name"], None,
                                 "not evaluated -- INSUFFICIENT_DATA").as_dict()
                   for c in bar["clauses"]]
    else:
        outcomes = evaluate_clauses(bar, observed)
        clauses = [o.as_dict() for o in outcomes]
        all_met = all(o.met for o in outcomes)
        status, exit_code = (OK, 0) if all_met else ("CLAUSE_FAILED", 1)

    payload = {
        "_schema": "sourmash-vs-mash-crosscheck-v1",
        "generated": "2026-10-19",
        "status": status,
        "exit_code": exit_code,
        "bar": str(Path(args.bar)).replace("\\", "/"),
        "bar_sha256": hashlib.sha256(
            (REPO / args.bar).read_bytes()).hexdigest(),
        "sourmash_params": {"ksize": args.ksize, "scaled": args.scaled},
        "refseq_cache": str(cache).replace("\\", "/") if cache else None,
        "n_pairs_expected": len(reference),
        "n_pairs_observed": len(observed),
        "clauses": clauses,
        "cells": cells,
        "blocked": blocked,
        "agreement_number_withheld": bool(blocked),
        "why_withheld": (
            "A cohort whose assemblies do not all resolve produces LOWER-BOUND "
            "effective-lineage counts. Reporting an agreement rate over a partial "
            "cohort would compare a lower bound against an exact-equality clause, so "
            "no agreement number is emitted. INSUFFICIENT_DATA is a valid result."
        ) if blocked else None,
        "honest_limits": bar["honest_limits"],
    }

    out_json = REPO / f"{args.out_stem}.json"
    out_md = REPO / f"{args.out_stem}.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        f"# sourmash vs Mash crosscheck — {status}",
        "",
        f"- bar: `{args.bar}`",
        f"- sourmash: ksize={args.ksize} scaled={args.scaled}",
        f"- refseq cache: `{cache}`" if cache else "- refseq cache: **none found**",
        f"- pairs observed: {len(observed)} / {len(reference)}",
        "",
    ]
    if blocked:
        lines += [
            "## INSUFFICIENT_DATA — no agreement number",
            "",
            payload["why_withheld"],
            "",
            "| cohort | reason |",
            "|---|---|",
            *[f"| {b['cohort']} | {b['reason']} |" for b in blocked],
            "",
        ]
    lines += ["## Clauses", "", "| # | clause | met | detail |", "|---|---|---|---|"]
    for c in clauses:
        met = {True: "yes", False: "NO", None: "—"}[c["met"]]
        lines.append(f"| {c['clause_id']} | {c['name']} | {met} | {c['detail']} |")
    lines += ["", "## Honest limits", ""]
    lines += [f"- {limit}" for limit in bar["honest_limits"]]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"STATUS {status} (exit {exit_code})")
    print(f"  wrote {out_json}")
    print(f"  wrote {out_md}")
    for c in clauses:
        print(f"  clause {c['clause_id']} {c['name']}: met={c['met']}")
    if blocked:
        print(f"  blocked cohorts: {len(blocked)} / {len(lineage['cells'])}")
        for b in blocked[:4]:
            print(f"      {b['cohort']}: {b['reason']}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
