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

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GoldsetIsolate:
    strain_id: str
    masked_vcf: str
    regeno_vcf: str | None
    label: str  # measured DST R/S


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
