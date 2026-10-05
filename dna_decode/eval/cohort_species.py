"""Which SPECIES is a validation cohort actually made of? (NON-FROZEN, offline, pure.)

A cohort directory named `klebsiella_*` scored with AMRFinder `-O Klebsiella_pneumoniae` is consistent at
the GENUS level and can be silently wrong at the SPECIES level. Nothing checked that, and the cost is
measured: the `Klebsiella x meropenem` provenance-disjoint cell reports sens 0.467 with 16 false
negatives whose only recorded explanation was that the frozen rule is "blind to porin-loss-mediated R" --
an explanation that is false (the rule counts porin truncations, which AMRFinder files under
Subclass=CARBAPENEM). This module exists to answer what those isolates ARE instead.

WHAT THIS MODULE IS NOT. It does not call `call_resistance`, touch Docker, or read the network. Predictions
are SUPPLIED by the caller (from `compute_lineage_metrics.reconcile_raw_metrics`, which refuses to hand
them over unless the cohort's committed confusion matrix reproduces exactly). That split is deliberate:
attribution computed from an unreconciled prediction set is how a plausible-looking wrong story gets
published, so the gate lives upstream and this module stays a pure function of its inputs.

IT DISCLOSES, IT DOES NOT RE-SCORE. The 10 `provenance_disjoint_validation_*.json` artifacts are named
frozen units of the reproducibility freeze. Nothing here rewrites a published metric; the output is a
separate, namespace-separate composition record.

THE LABEL'S EVIDENCE CLASS. Species come from each assembly's own GenBank `ORGANISM` field -- the
SUBMITTER's assertion, not a wet-lab identification. Strong evidence about what a cohort CONTAINS; never
an independent species call. See `dna_decode/data/genbank_organism.py`.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from dna_decode.data.genbank_organism import genus_of, organism_from_genbank, species_of
from dna_decode.data.organism_vocab import UnknownOrganism, canonical

#: Above this unresolved fraction, no composition verdict is issued. A composition computed over a
#: silently-shrunken denominator is this module's own failure mode, so it refuses instead.
#: ASSERTED, not derived -- recorded so a reader can disagree with the number rather than guess it.
MAX_UNRESOLVED_FRACTION = 0.10

#: Verdict states. `INSUFFICIENT_RESOLUTION` carries NO composition numbers, by construction.
MATCHES = "matches_expected"
SAME_GENUS_OTHER_SPECIES = "same_genus_other_species"
OTHER_GENUS = "other_genus"
UNRESOLVED = "unresolved"

INSUFFICIENT_RESOLUTION = "INSUFFICIENT_RESOLUTION"
#: NAMED FOR GRANULARITY-NEUTRALITY, after a first draft called this `SINGLE_SPECIES_AS_EXPECTED` and the
#: real data proved that false: the Campylobacter cohort holds 31 *C. jejuni* + 9 *C. coli* and is scored
#: with the GENUS-level `-O Campylobacter`, so it legitimately matches what it was scored as while being
#: two species. The verdict asserts agreement with the scoring organism at whatever granularity that
#: organism has -- never that the cohort is monomorphic.
COMPOSITION_CLEAN = "MATCHES_SCORED_ORGANISM"
COMPOSITION_MIXED = "MIXED_SPECIES"


class CrosstabMismatch(Exception):
    """The per-species outcome cells do not sum to the confusion matrix they were derived from."""


@dataclass(frozen=True)
class SpeciesClassification:
    """How each accession's resolved species relates to the organism the cohort was SCORED as."""

    expected_species: str | None
    by_bucket: dict[str, list[str]]
    species_counts: dict[str, int]
    n_total: int

    def count(self, bucket: str) -> int:
        return len(self.by_bucket.get(bucket, ()))

    @property
    def unresolved_fraction(self) -> float:
        return self.count(UNRESOLVED) / self.n_total if self.n_total else 0.0


def _expected_species_from_amrfinder_organism(amrfinder_organism: str) -> str | None:
    """`Klebsiella_pneumoniae` -> `Klebsiella pneumoniae`; a GENUS-only value -> None.

    AMRFinder's `-O` vocabulary mixes granularities: `Escherichia` and `Campylobacter` are genus-level
    while `Klebsiella_pneumoniae` is species-level. A genus-level flag cannot pin a species, so this
    returns None there and the classifier falls back to a genus-only comparison rather than inventing a
    species to compare against.
    """
    if not amrfinder_organism:
        return None
    return species_of(amrfinder_organism.replace("_", " "))


def resolve_cohort_species(refseq_root: Path, accessions: list[str]) -> dict[str, str | None]:
    """{accession: organism-name-or-None} from each assembly's own GenBank header. Nones are kept."""
    out: dict[str, str | None] = {}
    for acc in accessions:
        out[acc] = organism_from_genbank(Path(refseq_root) / acc / "annotations.gbk")
    return out


def classify_species(resolved: dict[str, str | None], amrfinder_organism: str) -> SpeciesClassification:
    """Partition accessions by how their species relates to the scored organism.

    Four buckets, and they PARTITION: every accession lands in exactly one, so the counts always sum to
    the input size. A silently-shrinking denominator is the failure mode this guards.
    """
    expected_species = _expected_species_from_amrfinder_organism(amrfinder_organism)
    expected_genus = genus_of(amrfinder_organism.replace("_", " ")) if amrfinder_organism else None

    buckets: dict[str, list[str]] = {MATCHES: [], SAME_GENUS_OTHER_SPECIES: [], OTHER_GENUS: [],
                                     UNRESOLVED: []}
    species_counts: Counter[str] = Counter()

    for acc, name in sorted(resolved.items()):
        if not name:
            buckets[UNRESOLVED].append(acc)
            continue
        sp, gn = species_of(name), genus_of(name)
        species_counts[sp or (gn or "?")] += 1
        if expected_species is not None:
            if sp == expected_species:
                buckets[MATCHES].append(acc)
            elif gn is not None and gn == expected_genus:
                buckets[SAME_GENUS_OTHER_SPECIES].append(acc)
            else:
                buckets[OTHER_GENUS].append(acc)
        else:
            # genus-level `-O` value: the strongest statement available is genus agreement
            buckets[MATCHES if (gn is not None and gn == expected_genus) else OTHER_GENUS].append(acc)

    return SpeciesClassification(expected_species=expected_species, by_bucket=buckets,
                                 species_counts=dict(species_counts), n_total=len(resolved))


def vocabulary_status(species_name: str | None) -> str:
    """Whether a resolved species is one the project's organism vocabulary knows.

    Reuses `organism_vocab.canonical` and catches `UnknownOrganism` so an unsupported species is
    REPORTED rather than crashing an audit -- K. aerogenes is exactly such a species, and it is the
    finding, so it must survive to the artifact.
    """
    if not species_name:
        return "unresolved"
    try:
        canonical(species_name)
    except UnknownOrganism:
        return "not_in_organism_vocab"
    except Exception:
        return "not_in_organism_vocab"
    return "in_organism_vocab"


def compose(
    resolved: dict[str, str | None],
    preds_by_acc: dict[str, str],
    labels: dict[str, int],
    confusion: dict | None = None,
) -> dict:
    """Species x outcome cross-tab: {species: {tp,fp,tn,fn}}.

    `preds_by_acc` must come from a RECONCILED source. When `confusion` is supplied, the assembled cells
    are checked against it and `CrosstabMismatch` is raised on disagreement -- a cross-tab that does not
    reconcile with the matrix it decomposes is not evidence.
    """
    table: dict[str, dict[str, int]] = {}
    totals = Counter()
    for acc, pred in sorted(preds_by_acc.items()):
        if acc not in labels:
            continue
        name = resolved.get(acc)
        key = species_of(name) if name else "(unresolved)"
        key = key or "(unresolved)"
        truth = bool(labels[acc])
        if pred == "R":
            cell = "tp" if truth else "fp"
        elif pred == "S":
            cell = "fn" if truth else "tn"
        else:
            cell = "abstain"
        table.setdefault(key, {"tp": 0, "fp": 0, "tn": 0, "fn": 0, "abstain": 0})[cell] += 1
        totals[cell] += 1

    if confusion is not None:
        for cell in ("tp", "fp", "tn", "fn"):
            if confusion.get(cell) is not None and confusion[cell] != totals[cell]:
                raise CrosstabMismatch(
                    f"{cell}: confusion={confusion[cell]} crosstab_sum={totals[cell]} "
                    f"(refusing to report a decomposition that does not reconcile)"
                )
    return {"by_species": table, "totals": dict(totals)}


def composition_verdict(
    cls: SpeciesClassification,
    max_unresolved_fraction: float = MAX_UNRESOLVED_FRACTION,
) -> dict:
    """The composition finding, or an explicit refusal.

    On refusal the payload carries NO composition counts at all -- not zeros, not partial numbers. A
    reader must not be able to pick a figure out of a refused verdict.
    """
    if cls.n_total == 0:
        return {"verdict": INSUFFICIENT_RESOLUTION, "reason": "no accessions", "n_total": 0}
    if cls.unresolved_fraction > max_unresolved_fraction:
        return {
            "verdict": INSUFFICIENT_RESOLUTION,
            "reason": (f"{cls.count(UNRESOLVED)}/{cls.n_total} assemblies had no readable ORGANISM line "
                       f"({cls.unresolved_fraction:.1%} > {max_unresolved_fraction:.0%} bar)"),
            "n_total": cls.n_total,
            "n_unresolved": cls.count(UNRESOLVED),
            "max_unresolved_fraction": max_unresolved_fraction,
        }

    off = cls.count(SAME_GENUS_OTHER_SPECIES) + cls.count(OTHER_GENUS)
    return {
        "verdict": COMPOSITION_MIXED if off else COMPOSITION_CLEAN,
        "expected_species": cls.expected_species,
        "n_total": cls.n_total,
        "n_matches_expected": cls.count(MATCHES),
        "n_same_genus_other_species": cls.count(SAME_GENUS_OTHER_SPECIES),
        "n_other_genus": cls.count(OTHER_GENUS),
        "n_unresolved": cls.count(UNRESOLVED),
        "species_counts": cls.species_counts,
        "off_species_accessions": sorted(cls.by_bucket[SAME_GENUS_OTHER_SPECIES]
                                        + cls.by_bucket[OTHER_GENUS]),
    }
