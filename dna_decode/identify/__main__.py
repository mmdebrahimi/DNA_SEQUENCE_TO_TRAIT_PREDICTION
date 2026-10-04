"""`python -m dna_decode.identify` -> the same entry point as `dna-identify`."""
from dna_decode.identify.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
