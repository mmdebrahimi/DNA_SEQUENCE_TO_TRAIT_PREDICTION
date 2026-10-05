"""FROZEN decision thresholds for the organism-identification router.

Written ONLY AFTER the mechanism analysis, never before. Artifacts:
    wiki/identify_distance_distribution_2026-10-04.json   (within/between distributions)
    wiki/identify_outofset_probe_2026-10-04.json          (where out-of-set genomes land)

MEASURED on the balanced 212-genome / 10-organism reference (k=21, s=1000):

    within-organism pairs            median 0.01208   p95 0.09880   max 0.11991
    between-organism pairs           p05    0.19173   median 1.00000
    nearest SAME-organism neighbour  median 0.00165   MAX 0.09880   <-- the in-set ceiling
    nearest OTHER-organism neighbour MIN    0.13098   median 0.26302
    distributions overlap            False

    out-of-set nearest in-set neighbour:
        Klebsiella variicola      x3   0.05095  -> klebsiella_pneumoniae
        Klebsiella michiganensis  x1   0.06746  -> klebsiella_oxytoca
        Klebsiella aerogenes      x43  0.10849 .. 0.13386 -> klebsiella_pneumoniae
        Gemmata obscuriglobus     x1   1.00000  -> escherichia_coli

DERIVATION OF `MAX_DISTANCE`
    An in-set genome must be ACCEPTED, so the threshold must exceed the worst nearest-same-organism
    distance: 0.09880. An out-of-set genome should be REJECTED, so it should fall below the closest
    out-of-set approach that we can actually exclude: K. aerogenes at 0.10849. That leaves the window
    (0.09880, 0.10849); MAX_DISTANCE is its midpoint. Not tuned -- derived from two measured
    quantities, and the midpoint is the maximally-robust point inside a window this narrow.

A MEASURED, NAMED BLIND SPOT -- NOT A TUNING FAILURE
    Klebsiella variicola (0.051) and K. michiganensis (0.067) sit BELOW the in-set ceiling (0.0988),
    so NO threshold can reject them without also rejecting legitimate in-set genomes. Mash sketch
    distance at this granularity cannot separate a congener from a conspecific. Consequence, stated
    plainly: the router WILL call a K. variicola genome `klebsiella_pneumoniae`. Four of 48 out-of-set
    genomes are unrejectable this way; the other 44 (43 aerogenes + Gemmata) are rejected.

    Operationally that mis-call is the least-bad routing available -- AMRFinder has no `-O` value for
    K. variicola at all -- but it is a species-level falsehood and the contract records it rather than
    reporting 91.7% out-of-set abstention as if the residual were noise. Lowering MAX_DISTANCE to
    catch it is REFUSED: it would reject real in-set genomes, trading a named blind spot for a silent
    coverage loss.

DERIVATION OF `AMBIGUITY_MARGIN`
    This must not fire on a correctly-identified in-set genome. The worst-case margin for such a
    genome is (closest other-organism approach) - (worst same-organism distance) =
    0.13098 - 0.09880 = 0.03218. The margin must stay BELOW that or legitimate calls abstain.
    0.015 is under half of it -- meaningful, with headroom for a corpus whose extremes are drawn from
    two different genomes.

SCOPE. These are calibrated on THIS reference (balanced 25/organism, k=21/s=1000).
Rebuilding the reference at a different composition, k, or sketch size invalidates them; the analysis
must be re-run, which is why the artifacts above are cited by name and the sketch parameters are
stamped into the manifest.

RE-CHECKED 2026-10-04 AFTER A COMPOSITION CHANGE, and the values HELD -- verified, not assumed. The
reference grew from 212 genomes / 10 organisms to 245 / 12 (adding Candida auris 8 and Mycobacterium
tuberculosis 25). Per the warning above that invalidates the calibration, so BOTH derivation inputs
were re-measured:

    in-set ceiling (max nearest-same-organism)   0.09880  -> 0.09880   UNCHANGED
    closest cross-organism approach              0.13098  -> 0.13098   UNCHANGED
    excludable out-of-set floor (K. aerogenes)   0.10849  -> unchanged (out-of-set abstention
                                                             identical at 44/48, so the floor that
                                                             bounds max_distance did not move)

WHY it held, rather than luck: the two organisms added are the two LEAST likely to widen the in-set
ceiling. M. tuberculosis is monomorphic (pairwise median ~6.3e-4 on record in this repo), so its
genomes crowd the low end; C. auris is a FUNGUS and sits at distance ~1.0 from every bacterium, so it
cannot pull the cross-organism floor down. A composition change that added a DIVERSE bacterial species
could still close the window -- re-run the analysis, do not inherit this result.
"""
from __future__ import annotations

from dna_decode.identify.core import Thresholds

#: Measured in-set ceiling: the largest nearest-same-organism distance over the reference.
IN_SET_CEILING = 0.09880

#: Measured closest approach of the out-of-set class we CAN exclude (K. aerogenes).
EXCLUDABLE_OUT_OF_SET_FLOOR = 0.10849

#: Midpoint of the measured window (IN_SET_CEILING, EXCLUDABLE_OUT_OF_SET_FLOOR).
MAX_DISTANCE = 0.1036

#: Worst-case margin available to a correctly-identified in-set genome: 0.13098 - 0.09880.
IN_SET_WORST_CASE_MARGIN = 0.03218

#: Comfortably under half of it, so a legitimate call does not abstain.
AMBIGUITY_MARGIN = 0.015

#: Out-of-set genomes this reference CANNOT reject at any threshold, with their measured distances.
#: Named so the limit travels with the cell instead of being averaged into an abstention rate.
UNREJECTABLE_CONGENERS = {
    "Klebsiella_variicola": 0.05095,
    "Klebsiella_michiganensis": 0.06746,
}

#: The artifacts these numbers were derived from. Cited so a reader can re-derive rather than trust.
DERIVED_FROM = (
    "wiki/identify_distance_distribution_2026-10-04.json",
    "wiki/identify_outofset_probe_2026-10-04.json",
)


def frozen() -> Thresholds:
    """The single Thresholds object every caller uses. One definition, no per-caller literals."""
    return Thresholds(max_distance=MAX_DISTANCE, ambiguity_margin=AMBIGUITY_MARGIN)
