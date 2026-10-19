"""EMBL-EBI AMR Portal (CABBAGE) adapter — header dialect + a fail-closed compatibility guard.

WHY THIS MODULE EXISTS, measured 2026-10-09. The portal ships PRE-COMPUTED AMRFinderPlus output per
assembly (`<ACC>_amrfinderplus.tsv.gz`, ~1 KB) at
`https://ftp.ebi.ac.uk/pub/databases/amr_portal/genomes/{GCA|ERZ}/<d1>/<d2>/<ACC>/`. That is worth a
great deal here: this project runs AMRFinder through Docker at ~95 s/strain, so a few thousand genomes
is >100 hours of compute the portal has already done.

BUT ITS HEADER DIALECT DIFFERS FROM THE ONE `eval/amr_rules` PARSES, AND THE FAILURE IS SILENT.
The portal writes `Element_symbol` / `Protein_id` / `Contig_id` (UNDERSCORES); `amr_rules` reads
`r.get("Element symbol")` (SPACES). A missing key yields None -> "" -> the row contributes no
determinant, so the rule returns a CONFIDENT SUSCEPTIBLE instead of raising.

Measured on real portal files, through the deployed `call_resistance`:

    GCA_001416835.1 (Campylobacter jejuni; carries gyrA_T86I, tet(O), aph(2'')-If)
      as delivered:  ciprofloxacin S (n=0)   tetracycline S (n=0)   gentamicin R (n=1)
      normalized:    ciprofloxacin R (n=1)   tetracycline R (n=1)   gentamicin R (n=1)

TWO OF THREE DRUGS WERE SILENTLY WRONG, including `gyrA_T86I` -- the determinant the project's own
Campylobacter ciprofloxacin cell is overwhelmingly built on. **Gentamicin was right BY LUCK**: its rule
keys on `Subclass`, a column whose name contains no underscore, so it is unaffected. That is precisely
what makes the dialect dangerous -- a spot-check on the wrong drug makes the bug invisible.

Same failure CLASS as the TB regression this repo already recorded: a masked-VCF `FILTER` convention
mismatch made a file with 1,382 genuine `GT=1/1` rows parse to zero calls and print a confident S. The
remedy there was to REFUSE rather than report; this module does the same at the schema level, which is
strictly better than a result-level check because **zero determinants is a legitimate result** for a
susceptible genome and must not be conflated with a parse failure.
"""
from __future__ import annotations

import gzip
from pathlib import Path

# Columns `eval/amr_rules` reads by SPACE-form name. If a header carries the UNDERSCORE form of any of
# these, the repo parser silently sees nothing for it.
REPO_SPACE_COLUMNS = ("Element symbol", "Element name", "Protein id", "Contig id",
                      "Target length", "Reference sequence length", "Alignment length")

# Columns that happen to match in both dialects. `Subclass` is in here, which is exactly why the
# gentamicin rule appeared to work while the symbol-keyed rules did not.
DIALECT_AGNOSTIC_COLUMNS = ("Class", "Subclass", "Type", "Subtype", "Method", "Scope", "Strand",
                            "Start", "Stop")


class PortalHeaderIncompatible(ValueError):
    """The AMRFinder header is in the portal's underscore dialect, which `amr_rules` cannot read.

    Raised rather than auto-fixed at the read site, so a caller cannot accidentally score a whole
    cohort through the wrong dialect and get confident susceptible calls.
    """


def read_text(path: str | Path) -> str:
    """Read a portal AMRFinder TSV, transparently gunzipping a `.gz`."""
    p = Path(path)
    if p.suffix == ".gz":
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return p.read_text(encoding="utf-8", errors="replace")


def header_dialect(text: str) -> str:
    """"portal" (underscore), "repo" (space), or "unknown". Pure."""
    if not text.strip():
        return "unknown"
    head = text.splitlines()[0]
    cols = set(head.split("\t"))
    if any(c.replace(" ", "_") in cols for c in REPO_SPACE_COLUMNS):
        return "portal"
    if any(c in cols for c in REPO_SPACE_COLUMNS):
        return "repo"
    return "unknown"


def normalize_header(text: str) -> str:
    """Convert the portal dialect to the repo dialect, touching ONLY the header line.

    Underscores are replaced on line 1 and nowhere else: determinant symbols in the BODY legitimately
    contain underscores (`gyrA_T86I`, `rplV_A103V`), and rewriting those would corrupt exactly the
    QRDR point mutations the ciprofloxacin rule counts.
    """
    if not text.strip():
        return text
    lines = text.splitlines()
    lines[0] = lines[0].replace("_", " ")
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def assert_repo_compatible(text: str) -> None:
    """FAIL CLOSED on the portal dialect. Call before handing a file to `amr_rules`."""
    d = header_dialect(text)
    if d == "portal":
        raise PortalHeaderIncompatible(
            "AMRFinder header is in the EMBL-EBI AMR Portal dialect (underscored column names, e.g. "
            "`Element_symbol`). `eval/amr_rules` reads the space form (`Element symbol`), so the "
            "symbol-keyed rules (ciprofloxacin qrdr_point, tetracycline gene_prefixes, the gentamicin "
            "rmt symbol_rescue) would contribute ZERO determinants and the call would be a confident "
            "SUSCEPTIBLE. Run `normalize_header` first. NOTE: Subclass-keyed rules are unaffected, so a "
            "spot-check on gentamicin alone will NOT reveal this."
        )
    if d == "unknown":
        raise PortalHeaderIncompatible(
            "AMRFinder header carries none of the expected columns in either dialect; refusing rather "
            f"than scoring an unrecognised schema. Expected one of: {REPO_SPACE_COLUMNS}"
        )


def prepare(path: str | Path, out_path: str | Path) -> Path:
    """Read a portal AMRFinder file, normalize if needed, write a repo-readable TSV. Returns out_path.

    Idempotent: a file already in the repo dialect is passed through unchanged.
    """
    text = read_text(path)
    if header_dialect(text) == "portal":
        text = normalize_header(text)
    assert_repo_compatible(text)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def genome_url(accession: str, release: str = "2026-07") -> str:
    """The portal's per-assembly AMRFinderPlus URL.

    Layout (from the FTP README): indexed by assembly accession, split by accession TYPE (GCA or ERZ),
    then the first three and second three digits. `release` is accepted for call-site documentation;
    the genomes/ tree is not release-partitioned.
    """
    acc = accession.strip()
    if acc.startswith("GCA_"):
        digits, kind = acc[4:].split(".")[0], "GCA"
    elif acc.startswith("ERZ"):
        digits, kind = acc[3:], "ERZ"
    else:
        raise ValueError(f"expected a GCA_ or ERZ accession, got {acc!r}")
    if len(digits) < 6:
        raise ValueError(f"accession {acc!r} has too few digits to index ({digits!r})")
    return (f"https://ftp.ebi.ac.uk/pub/databases/amr_portal/genomes/{kind}/"
            f"{digits[:3]}/{digits[3:6]}/{acc}/{acc}_amrfinderplus.tsv.gz")
