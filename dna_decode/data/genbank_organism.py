"""The ONE definition of "which organism does this assembly say it is?".

Reads the `  ORGANISM` line out of a GenBank header. Stdlib only, deliberately.

WHY THIS IS A SEPARATE MODULE AND NOT PART OF `annotations.py`. `annotations.py` is the natural home on
cohesion grounds -- it already owns GenBank parsing (`parse_genbank`, `load_annotation_table`) -- and that
is where the plan put it. Measured before following the plan: importing `annotations` costs **1.399s**
against **0.241s** for a stdlib-light module, because it pulls pandas for its DataFrame schema. Both
existing consumers (`scripts/build_identify_reference.py`, `scripts/identify_validate.py`) are currently
pandas-free, so hosting a one-line header reader there would have charged them ~1.16s per invocation for
nothing. Cohesion lost to import weight; the single-definition goal is unaffected.

WHY IT LIVES IN THE PACKAGE AT ALL. It was defined in `scripts/build_identify_reference.py`, so a package
module needing it (`dna_decode/eval/cohort_species.py`) would have had to import from `scripts/` --
inverting the layering. The direction is now scripts -> package, and the script keeps a re-export so its
five existing assertion sites and `identify_validate.py`'s import keep working untouched.

THE LABEL'S EVIDENCE CLASS, because it is easy to over-read. This is the SUBMITTER's assertion recorded in
NCBI's own metadata -- not a wet-lab identification, and not an independent species call. It is strong
evidence about what a cohort CONTAINS and it is the same evidence class this project's G1 circular-label
gate covers. Never present it as ground truth about an organism's identity.
"""
from __future__ import annotations

from pathlib import Path

#: Header lines to scan before giving up. The `  ORGANISM` line sits near the top of a GenBank record;
#: a file that has not produced one by here either is not GenBank or has an unusual header.
DEFAULT_MAX_HEADER_LINES = 60


def organism_from_genbank(gbk: Path, max_header_lines: int = DEFAULT_MAX_HEADER_LINES) -> str | None:
    """The assembly's own `  ORGANISM` value, or None.

    Reads only the header. Returns None rather than raising on a missing/short/unreadable file, so one
    bad genome cannot abort a 2,000-genome enumeration -- but the caller MUST count the Nones, because
    a silently-shrinking corpus is the failure mode here.
    """
    try:
        with Path(gbk).open(encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if line.startswith("  ORGANISM"):
                    got = line.split("ORGANISM", 1)[1].strip()
                    return got or None
                if i > max_header_lines:
                    return None
    except OSError:
        return None
    return None


def genus_of(organism_name: str) -> str | None:
    """The GENUS token of a GenBank organism name, or None when there isn't one.

    Separate from `build_identify_reference.genus_species` on purpose: that one joins the first TWO
    tokens into a species key (`Klebsiella_pneumoniae`), which cannot answer "same genus, different
    species?" -- the question a genus-named cohort hides. A cohort directory called `klebsiella_*` scored
    with `-O Klebsiella_pneumoniae` is consistent at the genus level and silently wrong at the species
    level, so the two levels need separate accessors.
    """
    if not organism_name:
        return None
    first = organism_name.split()
    return first[0] if first else None


def species_of(organism_name: str) -> str | None:
    """`Genus species` (first two tokens, space-joined), or None when the name has fewer than two.

    Returns None rather than a genus-only string so a caller cannot mistake `Klebsiella` for a resolved
    species. Strain suffixes (`Klebsiella aerogenes MGH 77` -> `Klebsiella aerogenes`) are dropped.
    """
    if not organism_name:
        return None
    toks = organism_name.split()
    return f"{toks[0]} {toks[1]}" if len(toks) >= 2 else None
