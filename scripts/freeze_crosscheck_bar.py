"""Freeze the sourmash-vs-Mash crosscheck acceptance bar -- DERIVED, and offline.

WHAT THIS IS FOR
----------------
Before `sourmash` may be accepted as an admissible genome-distance source for the
CABBAGE G8 gate, it has to reproduce what Mash already measured. This script writes
the bar that measurement will be scored against, and it writes it BEFORE any
agreement number exists. Amending a bar after seeing a number is a failure this
project has recorded twice; the only defence is to land the bar first.

Runs with NO drive, NO sourmash, NO network: everything here is computed from the
already-committed `wiki/provdisjoint_lineage_metrics.json`.

WHY THE CLAUSES ARE DERIVED RATHER THAN ASSERTED
------------------------------------------------
An earlier draft asserted clauses (an ARI floor, a cluster-count tolerance) that
were never checked against the distribution they would police. Reading the real
20-pair table changes what a meaningful bar looks like:

  * 19 of 20 pairs PASS, median minority count 14 against a bar of 3. So
    "verdict identity on all 20 pairs" is NECESSARY BUT NEARLY VACUOUS -- a badly
    wrong clustering still reproduces 19/20 by sheer margin. It is kept, and the
    artifact says in its own text that it is weak.
  * The discriminating pairs are the ones NEAR the bar, and that set is COMPUTED
    (within +/-2 of `MIN_EFFECTIVE_LINEAGES`), not hardcoded. At freeze time it
    resolves to three pairs.
  * The project's only real negative is a single TRIP. A bar that cannot fail on
    its one real negative is theatre, so reproducing that TRIP is its own clause.

ARI and cluster-count were DEMOTED to reported-not-gated. ARI is unanchored to
anything the G8 gate consumes, and a cluster-count tolerance binds only where
verdict identity already binds harder.

Naming the three near-bar pairs in the frozen artifact is NOT leaking future
information: they are computed from Mash numbers that have been committed for
months, so pinning them pins against known data. The artifact states this.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.eval import rejection_gates  # noqa: E402

LINEAGE_ARTIFACT = Path("wiki/provdisjoint_lineage_metrics.json")
BAR_PATH = Path("wiki/sourmash_vs_mash_crosscheck_bar.json")
SCHEMA = "sourmash-vs-mash-crosscheck-bar-v1"
NEAR_BAR_WINDOW = 2

# Reported in the crosscheck result, never gated on. Listed here so the demotion is
# part of the frozen commitment rather than an undocumented omission.
REPORTED_NOT_GATED = (
    "adjusted_rand_index",
    "n_clusters",
    "largest_cluster_fraction",
    "discordant_fraction",
    "full_20_pair_table",
)


class IncompleteLineageTable(RuntimeError):
    """The lineage table is partial, so no bar can be derived from it."""


class VacuousBar(RuntimeError):
    """The table carries no TRIP, so a TRIP-reproduction clause cannot be calibrated."""


@dataclass(frozen=True)
class Pair:
    organism: str
    drug: str
    threshold: str
    effective_lineage_N_R: int
    effective_lineage_N_S: int
    minority: int
    mash_verdict: str

    @property
    def label(self) -> str:
        return f"{self.organism} x {self.drug} @ {self.threshold}"


def load_lineage_table(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def assert_complete(data: dict) -> None:
    """Refuse a partial table.

    A bar derived from an incomplete table is a bar derived from nothing: a cohort
    whose genomes did not all resolve has effective-lineage counts that are lower
    bounds, and a lower bound cannot calibrate an exact-equality clause. This is the
    same incompleteness rule the G8 bar itself commits to.
    """
    bad = []
    for cell in data.get("cells", []):
        if cell.get("partial") or (cell.get("n_genomes_missing") or 0) > 0:
            bad.append({
                "cohort": cell.get("cohort"),
                "partial": cell.get("partial"),
                "n_genomes_missing": cell.get("n_genomes_missing"),
            })
    if bad:
        raise IncompleteLineageTable(
            f"{len(bad)} cell(s) are partial or missing genomes: {bad}"
        )


def pairs_from(data: dict) -> list[Pair]:
    """One Pair per (cell, threshold). The minority class is what G8 consumes."""
    out: list[Pair] = []
    for cell in data.get("cells", []):
        for threshold, block in sorted((cell.get("thresholds") or {}).items()):
            n_r = block.get("effective_lineage_N_R")
            n_s = block.get("effective_lineage_N_S")
            if n_r is None or n_s is None:
                raise IncompleteLineageTable(
                    f"{cell.get('cohort')} @ {threshold} has no effective-lineage counts"
                )
            minority = min(int(n_r), int(n_s))
            # Verdict via the IMPORTED gate, never a local reimplementation -- so a
            # change to MIN_EFFECTIVE_LINEAGES moves this bar instead of silently
            # diverging from the gate the measurement will actually be judged by.
            result = rejection_gates.g8_dedup_collapses_balance(
                {"min_effective_lineages": minority}
            )
            out.append(Pair(
                organism=cell.get("organism"),
                drug=cell.get("drug"),
                threshold=threshold,
                effective_lineage_N_R=int(n_r),
                effective_lineage_N_S=int(n_s),
                minority=minority,
                mash_verdict=result.verdict,
            ))
    return out


def near_bar_pairs(pairs: list[Pair], window: int = NEAR_BAR_WINDOW) -> list[Pair]:
    """Pairs whose Mash minority count is within `window` of the gate's bar.

    COMPUTED, not hardcoded. These are the only pairs where a small clustering
    difference can flip a verdict, which is why they carry an exact-equality clause
    while the other seventeen do not.
    """
    bar = rejection_gates.MIN_EFFECTIVE_LINEAGES
    return [p for p in pairs if abs(p.minority - bar) <= window]


def non_vacuity(pairs: list[Pair]) -> dict:
    """The base rate, stated so nobody reads a 19/20 headline as a result."""
    total = len(pairs)
    n_pass = sum(1 for p in pairs if p.mash_verdict == rejection_gates.PASS)
    n_trip = sum(1 for p in pairs if p.mash_verdict == rejection_gates.TRIP)
    if total and n_trip == 0:
        raise VacuousBar(
            "every pair PASSes, so clause 3 (the TRIP must reproduce) cannot be "
            "calibrated -- a bar with no reachable negative is theatre"
        )
    minorities = [p.minority for p in pairs]
    return {
        "n_pairs": total,
        "n_pass": n_pass,
        "n_trip": n_trip,
        "pass_base_rate": f"{n_pass}/{total}" if total else "0/0",
        "median_minority_count": statistics.median(minorities) if minorities else None,
        "min_minority_count": min(minorities) if minorities else None,
        "max_minority_count": max(minorities) if minorities else None,
        "bar": rejection_gates.MIN_EFFECTIVE_LINEAGES,
        "why_this_matters": (
            f"{n_pass} of {total} pairs PASS with a median minority count of "
            f"{statistics.median(minorities) if minorities else 'n/a'} against a bar of "
            f"{rejection_gates.MIN_EFFECTIVE_LINEAGES}. Clause 1 (verdict identity) is "
            "therefore NECESSARY BUT NEARLY VACUOUS: a badly wrong clustering still "
            "reproduces most verdicts by margin alone. Clauses 2 and 3 are what "
            "actually discriminate. Do NOT report a high verdict-identity rate as "
            "evidence that sourmash agrees with Mash."
        ),
    }


def _to_lf(raw: bytes) -> bytes:
    """`raw` with CRLF collapsed to LF. Kept escape-free and local: this exists only
    so the frozen artifact records BOTH digests of its own bytes."""
    return raw.replace(bytes((13, 10)), bytes((10,)))


def build_bar(data: dict, artifact_path: Path, raw_bytes: bytes) -> dict:
    assert_complete(data)
    pairs = pairs_from(data)
    near = near_bar_pairs(pairs)
    vac = non_vacuity(pairs)
    trips = [p for p in pairs if p.mash_verdict == rejection_gates.TRIP]

    return {
        "_schema": SCHEMA,
        "frozen": "2026-10-19",
        "frozen_before_any_agreement_number_existed": True,
        "purpose": (
            "Decide whether sourmash FracMinHash may be accepted as an admissible "
            "genome-distance source for the CABBAGE G8 effective-lineage gate, by "
            "requiring it to reproduce what Mash already measured on the committed "
            "provenance-disjoint cohorts."
        ),
        "derived_from": {
            "path": str(artifact_path).replace("\\", "/"),
            "sha256": hashlib.sha256(raw_bytes).hexdigest(),
            "sha256_lf": hashlib.sha256(_to_lf(raw_bytes)).hexdigest(),
            "note": (
                "The bar is DERIVED from this table's own arithmetic. If the table "
                "changes, this bar is no longer a commitment about it -- the scorer "
                "re-checks this sha256 and refuses rather than scoring."
            ),
            "eol_note": (
                "`sha256` is of the bytes AS THEY WERE WHEN FROZEN. This table carries "
                "no .gitattributes `eol` rule, so on a checkout with different line "
                "endings the raw digest differs while not one character has changed -- "
                "which is why BOTH digests are recorded and why the scorer verifies "
                "CONTENT via sourmash_vs_mash_crosscheck.verify_derived_from (it "
                "compares every EOL rendering) instead of comparing raw bytes."
            ),
        },
        "clause_logic": {
            "clause_3_is_implied_by_clause_1": True,
            "discovered": (
                "2026-10-19, while writing the scorer, BEFORE any agreement number existed"
            ),
            "note": (
                "Clause 1 requires verdict identity on ALL pairs; clause 3 requires the "
                "TRIP pairs to reproduce TRIP. The second is a SUBSET of the first, so it "
                "adds no logical strength -- if the TRIP flips to PASS, clause 1 fails "
                "too. Clause 3 is RETAINED deliberately, as a SEVERITY and REPORTING "
                "pre-commitment: it names in advance which direction of failure matters "
                "most, and forces the result to state the TRIP outcome explicitly instead "
                "of burying it inside a '19/20 verdicts identical' aggregate."
            ),
            "why_this_is_recorded_rather_than_fixed": (
                "The three requirements are left EXACTLY as frozen. Rewriting a clause "
                "after noticing it is weaker than it looked would be bar-shopping even in "
                "the honest direction. This annotation only makes the bar sound weaker "
                "than its first reading, never stronger, and it is stamped before any "
                "number exists."
            ),
            "why_build_bar_emits_it": (
                "It is EMITTED rather than hand-added so that `--force` re-freezing this "
                "bar cannot silently drop it. The implication is a property of the clause "
                "set this function builds, so it holds for every bar this function writes."
            ),
        },
        "gate_source": {
            "module": "dna_decode.eval.rejection_gates",
            "function": "g8_dedup_collapses_balance",
            "min_effective_lineages": rejection_gates.MIN_EFFECTIVE_LINEAGES,
            "note": (
                "Verdicts are computed THROUGH the imported gate, not restated. "
                "MIN_EFFECTIVE_LINEAGES is referenced, never copied."
            ),
        },
        "clauses": [
            {
                "id": 1,
                "name": "verdict_identity",
                "requirement": (
                    "Every one of the "
                    f"{len(pairs)} (cell, threshold) pairs must yield the same G8 "
                    "verdict under sourmash clustering as under Mash."
                ),
                "strength": "NECESSARY BUT NEARLY VACUOUS",
                "why_weak": vac["why_this_matters"],
            },
            {
                "id": 2,
                "name": "exact_minority_count_near_the_bar",
                "requirement": (
                    "On every pair whose Mash minority effective-lineage count is "
                    f"within +/-{NEAR_BAR_WINDOW} of the bar "
                    f"({rejection_gates.MIN_EFFECTIVE_LINEAGES}), the sourmash "
                    "minority count must be EXACTLY equal to the Mash count."
                ),
                "strength": "DISCRIMINATING",
                "near_bar_set_is_computed": True,
                "near_bar_window": NEAR_BAR_WINDOW,
                "near_bar_pairs": [
                    {"label": p.label, "minority": p.minority,
                     "mash_verdict": p.mash_verdict} for p in near
                ],
                "note": (
                    "Naming these pairs is not leaking future information: they are "
                    "computed from Mash numbers committed long before this freeze, so "
                    "pinning them pins against already-known data."
                ),
            },
            {
                "id": 3,
                "name": "the_trip_must_reproduce",
                "requirement": (
                    "Every pair that Mash scores TRIP must also score TRIP under "
                    "sourmash."
                ),
                "strength": "DISCRIMINATING -- the single most important clause",
                "trip_pairs": [
                    {"label": p.label, "minority": p.minority} for p in trips
                ],
                "why": (
                    "A method that converts the project's only real negative into a "
                    "PASS manufactures lineage diversity. That is the most dangerous "
                    "possible direction of error, because it would license sourmash on "
                    "CABBAGE and then report false lineage independence. A bar that "
                    "cannot fail on its one real negative is theatre."
                ),
            },
        ],
        "reported_not_gated": {
            "fields": list(REPORTED_NOT_GATED),
            "why": (
                "ARI and cluster-count were gate clauses in an earlier draft and are "
                "DEMOTED here with the reason recorded: ARI is unanchored to anything "
                "the G8 gate consumes, and a cluster-count tolerance binds only where "
                "verdict identity already binds harder. They are still emitted, so a "
                "reader can see HOW the partitions differ even when all clauses pass."
            ),
        },
        "non_vacuity": vac,
        "mash_reference_table": [asdict(p) for p in pairs],
        "honest_limits": [
            "FracMinHash is not Mash's bottom-sketch MinHash, so agreement on these "
            "cohorts does not make the two estimators interchangeable in general.",
            "The thresholds (0.001, 0.005) were calibrated against MASH distances. "
            "Passing this bar shows sourmash reproduces Mash's PARTITION VERDICTS at "
            "those thresholds on these cohorts; it does not re-derive the thresholds.",
            "Clause 3 rests on a SINGLE TRIP pair. One negative is enough to make the "
            "bar non-vacuous and far too few to characterise sourmash's behaviour "
            "near the bar in general.",
            "These are provenance-disjoint AMR cohorts, not CABBAGE. Accepting the "
            "distance source for CABBAGE is an extrapolation across cohorts.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lineage-artifact", default=str(LINEAGE_ARTIFACT))
    ap.add_argument("--out", default=str(BAR_PATH))
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing frozen bar (refused by default)")
    args = ap.parse_args(argv)

    src = (REPO / args.lineage_artifact) if not Path(args.lineage_artifact).is_absolute() \
        else Path(args.lineage_artifact)
    out = (REPO / args.out) if not Path(args.out).is_absolute() else Path(args.out)

    if out.exists() and not args.force:
        print(f"REFUSED: {out} already exists. A frozen bar is not silently rewritten; "
              f"pass --force only if you intend to re-freeze.")
        return 1

    raw = src.read_bytes()
    data = json.loads(raw.decode("utf-8"))
    bar = build_bar(data, Path(args.lineage_artifact), raw)
    out.write_text(json.dumps(bar, indent=2) + "\n", encoding="utf-8")

    vac = bar["non_vacuity"]
    print(f"WROTE {out}")
    print(f"  pairs            : {vac['n_pairs']}")
    print(f"  PASS / TRIP      : {vac['n_pass']} / {vac['n_trip']}")
    print(f"  base rate        : {vac['pass_base_rate']}  (median minority "
          f"{vac['median_minority_count']}, bar {vac['bar']})")
    print(f"  near-bar pairs   : {len(bar['clauses'][1]['near_bar_pairs'])} (computed)")
    for p in bar["clauses"][1]["near_bar_pairs"]:
        print(f"      {p['label']}  minority={p['minority']}  {p['mash_verdict']}")
    print(f"  TRIP pairs       : {len(bar['clauses'][2]['trip_pairs'])}")
    for p in bar["clauses"][2]["trip_pairs"]:
        print(f"      {p['label']}  minority={p['minority']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
