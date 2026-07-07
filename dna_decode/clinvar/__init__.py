"""Human Mendelian-disease decoder (ClinVar-backed) — the germline-pathogenicity track.

Deterministic curated-catalog decoder: iterate a VCF, look up each carried variant in the committed ClinVar
panel (ACMG SF v3.2 + carrier genes), return curated germline classifications (P/LP + B/LB) with gene /
disease / gold-star review level. Fail-closed: not-in-panel -> INDETERMINATE (absence != benign). NOT a
learned predictor, NOT a clinical tool. First-class CLI: `dna-clinvar`.
"""
from dna_decode.clinvar.decode import carried_alts, decode_vcf  # noqa: F401
