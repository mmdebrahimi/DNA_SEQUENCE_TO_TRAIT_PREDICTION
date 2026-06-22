"""Step 7 — independent post-2023 TB gold-set ingestion (deliverable b, NON-FROZEN).

Ingests an INDEPENDENT (non-CRyPTIC) TB cohort into the same shape the Step-3 cell + Step-5 collapse
consume: per isolate a masked-VCF text (for determinant calls) + an optional regeno-VCF text (for
callability) + a measured DST label. Independence requirement (ratified E): from the WHO-catalogue
BUILD, i.e. post-2023 isolates (temporal hold-out) — WHO v2 swept most public pre-2023 TB WGS+pDST.

There is no gold set on disk yet (it is hand-curated per the runbook), so `load_goldset` returns [] when
the manifest is absent; the scorer then BLOCKED-gates. The manifest is a JSON list of:
  {"strain_id": "...", "masked_vcf": "<path>", "regeno_vcf": "<path|null>", "label": "R|S"}
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GoldsetIsolate:
    strain_id: str
    masked_vcf: str
    regeno_vcf: str | None
    label: str  # measured DST R/S


# --- CRyPTIC-leakage exclusion (the "provably not in CRyPTIC" independence check) ---------------
# CRyPTIC reuse table accession columns. An independent isolate must NOT appear in ANY of these.
_CRYPTIC_ACC_COLS = ("ENA_RUN", "ENA_SAMPLE", "UNIQUEID")


def cryptic_accessions(reuse_csv: Path) -> set[str]:
    """All CRyPTIC isolate identifiers (ENA_RUN ∪ ENA_SAMPLE ∪ UNIQUEID), normalized for membership tests.

    This is the exclusion set: a candidate gold-set isolate that matches ANY of these is NOT independent
    (it is in CRyPTIC, which the WHO catalogue was built from). Run/sample accessions are the load-bearing
    keys — an ERR/ERS/SAMEA hit means the same sequencing run/biosample."""
    acc: set[str] = set()
    with open(reuse_csv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            for col in _CRYPTIC_ACC_COLS:
                v = (row.get(col) or "").strip()
                if v:
                    acc.add(v.upper())
    return acc


@dataclass(frozen=True)
class IndependenceReport:
    clean: tuple[str, ...]        # candidate ids with NO CRyPTIC overlap (independent)
    leaked: tuple[str, ...]       # candidate ids found in CRyPTIC (must be dropped)
    n_checked: int


def assert_independent(candidate_ids, cryptic_acc: set[str]) -> IndependenceReport:
    """Partition candidate isolate ids into clean (independent) vs leaked (present in CRyPTIC).

    `candidate_ids` is any iterable of identifiers (ENA run/sample/biosample/strain). Matching is
    case-insensitive exact-token (mirrors the accession-string leakage discipline elsewhere in the repo —
    NOT substring). The caller drops `leaked` before scoring so a 'CRyPTIC-leaked' isolate can never inflate
    an 'independent' number."""
    clean, leaked = [], []
    for cid in candidate_ids:
        (leaked if str(cid).strip().upper() in cryptic_acc else clean).append(str(cid))
    return IndependenceReport(clean=tuple(clean), leaked=tuple(leaked), n_checked=len(clean) + len(leaked))


def goldset_available(manifest_path: Path) -> bool:
    return Path(manifest_path).exists()


def load_goldset(manifest_path: Path) -> list[GoldsetIsolate]:
    """Parse the gold-set manifest. Absent -> [] (scorer BLOCKED-gates). Raises on malformed JSON."""
    p = Path(manifest_path)
    if not p.exists():
        return []
    rows = json.loads(p.read_text(encoding="utf-8"))
    out: list[GoldsetIsolate] = []
    for r in rows:
        out.append(GoldsetIsolate(
            strain_id=str(r["strain_id"]),
            masked_vcf=str(r["masked_vcf"]),
            regeno_vcf=(str(r["regeno_vcf"]) if r.get("regeno_vcf") else None),
            label=str(r["label"]).upper(),
        ))
    return out
