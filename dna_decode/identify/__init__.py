"""Closed-set organism identification: given a genome, say which SUPPORTED organism it is, or abstain.

This is step 1 of the project's north star -- the piece that was missing. Every other decoder REQUIRES
the caller to name the organism (`--organism`, `-O`), so the pipeline's first step was a human typing
the answer in.

CLOSED-SET, NOT TAXONOMY. The question answered is "which of the organisms this tool supports is this,
if any?" -- not "what is this?". Anything outside the supported set must ABSTAIN, because the output
selects which downstream rule runs and a confident wrong organism routes a genome to the wrong rule.
Abstention is the safe failure.

The decision layer (`core`) is pure and takes its thresholds as parameters; the Mash invocation lives
behind the runner so the whole decision surface is testable offline.
"""
from dna_decode.identify.core import (
    AbstainReason,
    Hit,
    IdentifyCall,
    Thresholds,
    accession_from_reference_id,
    decide,
    parse_mash_dist,
)

__all__ = [
    "AbstainReason",
    "Hit",
    "IdentifyCall",
    "Thresholds",
    "accession_from_reference_id",
    "decide",
    "parse_mash_dist",
]
