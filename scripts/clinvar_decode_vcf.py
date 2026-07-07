#!/usr/bin/env python
"""Back-compat shim — the Mendelian VCF-decode logic now lives in the `dna_decode.clinvar` package
(first-class `dna-clinvar` CLI). This re-exports it so existing callers + tests keep working.

    uv run python scripts/clinvar_decode_vcf.py --vcf X.vcf.gz --sample-id S --panel <panel.tsv>
is equivalent to `dna-clinvar X.vcf.gz --sample-id S --panel <panel.tsv> --out ...`.
"""
from __future__ import annotations

from dna_decode.clinvar.cli import main
from dna_decode.clinvar.decode import carried_alts, decode_vcf  # noqa: F401 (re-export)

_carried_alts = carried_alts  # legacy private alias used by tests


if __name__ == "__main__":
    raise SystemExit(main())
