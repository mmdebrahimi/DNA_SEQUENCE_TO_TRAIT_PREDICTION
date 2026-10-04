"""Pure parsing + the abstention decision for the organism-identification router.

No I/O, no Docker, no thresholds baked in. Thresholds arrive as a `Thresholds` parameter so this
module is complete and fully testable BEFORE any value is chosen -- the seam pattern the ECOFF work
established for a value that must not be guessed.

`mash dist` OUTPUT FORMAT, verified live against the pinned image (`mash dist -h`), not inherited:

    [reference-ID, query-ID, distance, p-value, shared-hashes]      tab-separated

`reference-ID` is the genome's path as it appeared in the sketch, so its parent directory name is the
accession -- which the reference manifest maps to an organism.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Sequence


class AbstainReason(str, Enum):
    """Why no organism was called. Distinct values, never collapsed into one 'unknown'.

    A caller needs to tell "this genome is not a supported organism" (NOT_IN_SET territory, a correct
    answer) from "the reference never loaded" (REFERENCE_UNAVAILABLE, an operational fault). Collapsing
    them would let a broken container read as a confident biological finding.
    """

    NO_HIT = "no_hit"
    ABOVE_MAX_DISTANCE = "above_max_distance"
    AMBIGUOUS_TOP2 = "ambiguous_top2"
    REFERENCE_UNAVAILABLE = "reference_unavailable"
    UNLABELLED_REFERENCE = "unlabelled_reference"


@dataclass(frozen=True)
class Thresholds:
    """The two decision parameters.

    max_distance: a nearest neighbour further than this is not evidence of anything; ABSTAIN.
    ambiguity_margin: if the best hits of two DIFFERENT organisms sit within this distance of each
        other, the genome is between them; ABSTAIN rather than pick.
    """

    max_distance: float
    ambiguity_margin: float

    def __post_init__(self) -> None:
        if not (0.0 < self.max_distance <= 1.0):
            raise ValueError(f"max_distance must be in (0, 1]; got {self.max_distance}")
        if self.ambiguity_margin < 0.0:
            raise ValueError(f"ambiguity_margin must be >= 0; got {self.ambiguity_margin}")


@dataclass(frozen=True)
class Hit:
    """One `mash dist` row."""

    reference_id: str
    query_id: str
    distance: float
    p_value: float
    shared_hashes: str


@dataclass(frozen=True)
class IdentifyCall:
    """The router's answer. Either an organism, or an abstention with a reason -- never both."""

    organism: str | None
    routing_token: str | None
    abstained: bool
    reason: AbstainReason | None = None
    nearest_distance: float | None = None
    runner_up_organism: str | None = None
    runner_up_distance: float | None = None
    top_hits: tuple[Hit, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.abstained and self.organism is not None:
            raise ValueError("an abstention cannot also carry an organism")
        if not self.abstained and self.organism is None:
            raise ValueError("a non-abstention must carry an organism")
        if self.abstained and self.reason is None:
            raise ValueError("an abstention must carry a machine-readable reason")

    def as_dict(self) -> dict:
        return {
            "organism": self.organism,
            "routing_token": self.routing_token,
            "abstained": self.abstained,
            "abstain_reason": self.reason.value if self.reason else None,
            "nearest_distance": self.nearest_distance,
            "runner_up_organism": self.runner_up_organism,
            "runner_up_distance": self.runner_up_distance,
            "top_hits": [
                {"reference_id": h.reference_id, "distance": h.distance,
                 "p_value": h.p_value, "shared_hashes": h.shared_hashes}
                for h in self.top_hits
            ],
            "closed_set": True,
            "does_not_support": [
                "open-world taxonomy -- an organism outside the supported set ABSTAINS and is NOT "
                "identified",
                "strain or subspecies resolution",
                "metagenomic mixtures (one assembled genome per call)",
            ],
        }


def parse_mash_dist(stdout: str) -> list[Hit]:
    """Parse `mash dist` stdout into Hits, sorted nearest-first.

    Malformed or short rows are SKIPPED rather than raising: a partially-written stream should not
    abort a call. The caller sees an empty list and abstains with NO_HIT, which is honest -- but note
    that means a garbage stream and a genuine no-match both reach NO_HIT, so the runner is responsible
    for distinguishing a failed invocation (REFERENCE_UNAVAILABLE) before calling here.
    """
    hits: list[Hit] = []
    for line in (stdout or "").splitlines():
        parts = line.rstrip("\r").split("\t")
        if len(parts) < 5:
            continue
        try:
            dist = float(parts[2])
            pval = float(parts[3])
        except ValueError:
            continue
        hits.append(Hit(parts[0], parts[1], dist, pval, parts[4]))
    # stable sort on distance only; ties keep input order and are resolved by the AMBIGUITY rule,
    # never by sort position.
    return sorted(hits, key=lambda h: h.distance)


def accession_from_reference_id(reference_id: str) -> str:
    """`/d/dna_decode_cache/refseq/GCA_000023665.1/genome.fna` -> `GCA_000023665.1`.

    The sketch stores each genome's path, and the accession is its parent directory name. Returns ""
    when there is no parent component, so the caller treats it as unlabelled rather than guessing.
    """
    norm = reference_id.replace("\\", "/").rstrip("/")
    parts = [p for p in norm.split("/") if p]
    return parts[-2] if len(parts) >= 2 else ""


def decide(
    hits: Sequence[Hit],
    thresholds: Thresholds,
    organism_of_accession: Callable[[str], str | None],
    routing_token_of: Callable[[str], str | None],
    *,
    top_k: int = 5,
) -> IdentifyCall:
    """Turn ranked hits into an organism call or a reasoned abstention.

    Rules, in order:
      1. no hits                                     -> ABSTAIN NO_HIT
      2. no hit maps to a labelled organism           -> ABSTAIN UNLABELLED_REFERENCE
      3. nearest labelled hit > max_distance          -> ABSTAIN ABOVE_MAX_DISTANCE
      4. best hit of a DIFFERENT organism is within
         ambiguity_margin of the nearest              -> ABSTAIN AMBIGUOUS_TOP2
      5. otherwise                                    -> call the nearest organism

    Rule 4 compares the best hit of each DISTINCT organism, not the raw top two rows. Two reference
    genomes of the same organism are the common case and are not ambiguity -- comparing raw rows would
    abstain constantly on the E. coli-heavy reference, where the top two rows are almost always both
    E. coli.
    """
    ranked = list(hits)
    if not ranked:
        return IdentifyCall(None, None, True, AbstainReason.NO_HIT)

    # best (smallest) distance per organism, preserving nearest-first order of first appearance
    best: list[tuple[str, Hit]] = []
    seen: set[str] = set()
    for h in ranked:
        org = organism_of_accession(accession_from_reference_id(h.reference_id))
        if not org or org in seen:
            continue
        seen.add(org)
        best.append((org, h))

    top = tuple(ranked[:top_k])
    if not best:
        return IdentifyCall(None, None, True, AbstainReason.UNLABELLED_REFERENCE,
                            nearest_distance=ranked[0].distance, top_hits=top)

    org0, hit0 = best[0]
    if hit0.distance > thresholds.max_distance:
        return IdentifyCall(None, None, True, AbstainReason.ABOVE_MAX_DISTANCE,
                            nearest_distance=hit0.distance, top_hits=top)

    if len(best) >= 2:
        org1, hit1 = best[1]
        if (hit1.distance - hit0.distance) < thresholds.ambiguity_margin:
            return IdentifyCall(None, None, True, AbstainReason.AMBIGUOUS_TOP2,
                                nearest_distance=hit0.distance,
                                runner_up_organism=org1, runner_up_distance=hit1.distance,
                                top_hits=top)
        return IdentifyCall(org0, routing_token_of(org0), False,
                            nearest_distance=hit0.distance,
                            runner_up_organism=org1, runner_up_distance=hit1.distance,
                            top_hits=top)

    return IdentifyCall(org0, routing_token_of(org0), False,
                        nearest_distance=hit0.distance, top_hits=top)
