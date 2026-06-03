"""v0 E. coli pathotype compatibility resolver (EP-4).

FASTA -> marker detection (detect.py) -> 11-class decision table (resolve.py) ->
audit-grade provenance JSON (cli.py). Deterministic, dependency-light, abstention-first.
Compatibility resolver, NOT a clinical predictor (ledger v5 framing).
"""
from dna_decode.pathotype.detect import detect, parse_fasta, build_vf_index, assembly_qc
from dna_decode.pathotype.resolve import resolve_call
from dna_decode.pathotype.markers import CLUSTER_MARKERS, RULES_VERSION

__all__ = ["detect", "parse_fasta", "build_vf_index", "assembly_qc",
           "resolve_call", "CLUSTER_MARKERS", "RULES_VERSION"]
