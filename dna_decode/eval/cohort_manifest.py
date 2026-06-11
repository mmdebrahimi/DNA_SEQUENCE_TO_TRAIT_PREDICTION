"""Accession-manifest registry — the data-driven leakage source-of-truth for provenance-disjoint validation.

Replaces the hardcoded `_FLAGSHIP_PARQUET_COHORTS` list in `scripts/provenance_disjoint_validate.py`
(which covered only 3 of the 8 parquet cohorts — a silent under-exclusion trap that grows with the suite).

`build_manifest()` scans EVERY `data/raw/*/selected.tsv` AND EVERY `data/processed/*.parquet`, recording
each accession's cohort + role + source + organism + drug. `prior_accessions(manifest, exclude_cohort)`
returns the union of accessions across ALL cohorts EXCEPT the exact-named one.

LEAKAGE SAFETY (pre-exec /brainstorm C1):
- EXACT-self identity, NOT a substring. `exclude_cohort` is the cohort NAME (dir/stem), so excluding the
  current output cohort still excludes EVERY OTHER cohort — including a prior provenance-disjoint validation
  for the same organism×drug. A re-run is therefore forced onto fresh accessions (no reuse of an earlier
  provdisjoint set, which the old substring rule allowed).
- `Manifest.incomplete` is set True if ANY accession source fails to load. The caller MUST fail closed
  (refuse to score / write a clean-sounding artifact) unless explicitly overridden — an incomplete manifest
  is a possible false-independence claim, NOT a harmless offline fallback.
"""
from __future__ import annotations

import glob
import re
from dataclasses import dataclass, field
from pathlib import Path

# Drugs whose names may appear in cohort directory/parquet names (for advisory organism/drug parsing).
_DRUG_TOKENS = {
    "ciprofloxacin": "ciprofloxacin", "cipro": "ciprofloxacin",
    "ceftriaxone": "ceftriaxone", "cef": "ceftriaxone",
    "tetracycline": "tetracycline", "tet": "tetracycline",
    "gentamicin": "gentamicin", "gent": "gentamicin", "genta": "gentamicin",
    "meropenem": "meropenem", "mero": "meropenem",
    "oxacillin": "oxacillin",
}


def _classify_role(name: str) -> str:
    """Heuristic role from cohort name. Advisory only — exclusion uses exact identity, not role."""
    n = name.lower()
    if "provdisjoint" in n:
        return "validation"
    if "stage2_n150" in n:
        return "tuning"
    if "gate_a" in n or "gate_b" in n:
        return "held-out"
    if any(tok in n for tok in _DRUG_TOKENS) or n.endswith("_cipro"):
        return "calibration"
    return "unknown"


def _parse_organism_drug(name: str) -> tuple[str, str]:
    n = name.lower()
    drug = ""
    for tok, canon in _DRUG_TOKENS.items():
        if re.search(rf"(^|_){re.escape(tok)}(_|$)", n) or tok in n:
            drug = canon
            break
    organism = re.split(
        r"_(provdisjoint|cipro|ceftriaxone|cef|tetracycline|tet|gentamicin|gent|meropenem|mero|"
        r"oxacillin|stage2|gate|n40|n150|mini)",
        n,
    )[0]
    return organism, drug


@dataclass
class Cohort:
    name: str
    path: str
    role: str
    source: str          # "selected_tsv" | "parquet"
    organism: str
    drug: str
    accessions: set[str]


@dataclass
class Manifest:
    cohorts: list[Cohort] = field(default_factory=list)
    incomplete: bool = False
    warnings: list[str] = field(default_factory=list)

    def cohort_names(self) -> list[str]:
        return [c.name for c in self.cohorts]


def build_manifest(data_raw: str = "data/raw", data_processed: str = "data/processed") -> Manifest:
    """Scan every selected.tsv + every parquet cohort into a Manifest. Sets `incomplete=True` on any
    load failure (fail-closed signal for the caller)."""
    m = Manifest()

    for sel in sorted(glob.glob(f"{data_raw}/*/selected.tsv")):
        name = Path(sel).parent.name
        try:
            accs: set[str] = set()
            for ln in Path(sel).read_text(encoding="utf-8").splitlines():
                if "\t" in ln:
                    acc = ln.split("\t")[0].strip()
                    if acc:
                        accs.add(acc)
            org, drug = _parse_organism_drug(name)
            m.cohorts.append(Cohort(name, sel, _classify_role(name), "selected_tsv", org, drug, accs))
        except Exception as e:  # noqa: BLE001 — load failure must mark INCOMPLETE, not crash
            m.incomplete = True
            m.warnings.append(f"selected.tsv load failed: {sel}: {type(e).__name__}: {e}")

    for pq in sorted(glob.glob(f"{data_processed}/*.parquet")):
        name = Path(pq).stem
        try:
            from dna_decode.data.cohort import load_cohort
            accs = set()
            for s in load_cohort(pq).strains:
                acc = getattr(s, "assembly_accession", None)
                if acc:
                    accs.add(acc)
            org, drug = _parse_organism_drug(name)
            m.cohorts.append(Cohort(name, pq, _classify_role(name), "parquet", org, drug, accs))
        except Exception as e:  # noqa: BLE001 — a parquet that won't load marks INCOMPLETE (fail-closed)
            m.incomplete = True
            m.warnings.append(f"parquet load failed: {pq}: {type(e).__name__}: {e}")

    return m


def prior_accessions(manifest: Manifest, exclude_cohort: str) -> set[str]:
    """Union of accessions across ALL cohorts EXCEPT the exact-named `exclude_cohort`.

    `exclude_cohort` is a cohort NAME (the data/raw dir name or parquet stem), matched exactly — NOT a
    substring. This is the leakage-safe identity: the current output cohort is the ONLY thing excluded;
    every other cohort (tuning, calibration, held-out, AND any prior provdisjoint validation) contributes
    its accessions to the exclusion set.
    """
    out: set[str] = set()
    for c in manifest.cohorts:
        if c.name == exclude_cohort:
            continue
        out |= c.accessions
    return out
