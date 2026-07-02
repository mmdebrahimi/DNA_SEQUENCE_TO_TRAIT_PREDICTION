"""Thin re-export — the canonical dataset-candidate scorecard is now in the installable package
`dna_decode.deconfound.scorecard` (promoted 2026-07-02). Kept so existing `scripts.` + test imports resolve.
See dna_decode/deconfound/scorecard.py for the implementation."""
from dna_decode.deconfound.scorecard import (  # noqa: F401
    GATE_KEYS,
    GATES,
    Candidate,
    decoder_type,
    rank,
    score,
)
