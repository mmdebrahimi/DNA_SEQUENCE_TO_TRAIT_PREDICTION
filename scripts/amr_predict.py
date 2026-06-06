"""Thin shim → the in-package AMR decoder CLI (dna_decode.amr.cli).

The implementation moved into the package (dna_decode/amr/cli.py) so it ships in the wheel as the
`dna-amr` console entry. This shim preserves `python -m scripts.amr_predict` for back-compat.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.amr.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
