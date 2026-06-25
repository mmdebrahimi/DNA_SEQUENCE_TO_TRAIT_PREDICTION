"""S. pneumoniae β-lactam AMR engine — deterministic PBP-type → MIC → R/S (NON-FROZEN).

The β-lactam cell the go/no-go deferred. Pneumococcal β-lactam resistance is PBP-mutation-MIC (not
gene-presence), so this engine is a PBP-TYPE → MIC LOOKUP using the CDC reference table
`Ref_PBPtype_MIC.csv` (from GlobalPneumoSeq/spn-pbp-amr, a fork of Ben Metcalf/CDC's AMR predictor) +
MIC → R/S via the breakpoint-context-aware `data.pneumo_breakpoints`.

Pipeline: (pbp1a, pbp2b, pbp2x) types -> APT key "pbp1a-pbp2b-pbp2x" -> lookup MIC for PEN/MER/TAX/CFT/CFX
-> `pneumo_breakpoints.classify(drug, context, mic)` -> R/I/S.

HONESTY (load-bearing):
  - Branding KNOWLEDGE_BASELINE / organism_routed — never the frozen E. coli surface.
  - The lookup table IS the CDC deterministic method (the simplest of CDC's three: lookup vs RandomForest vs
    ElasticNet). Faithful-to-tool: this reproduces the CDC PBP-MIC LOOKUP, not an independent baseline. The
    lookup is CDC's published reference (independent of any validation cohort — NOT a circular fit).
  - This engine takes PBP TYPES as input; it does NOT type the genome itself. PBP-typing-from-genome (BLAST
    the transpeptidase domains vs spn-pbp-amr's `Ref_PBP_3.faa`) is the deferred swap (mirrors the AMRFinder
    swap for the gene-presence cell). A novel/absent PBP type (contains 'NEW' or not in the table) -> None
    (no-call), never a fabricated MIC.
  - β-lactam R/S is breakpoint-AMBIGUOUS (meningitis/non-meningitis/oral) — the caller MUST pass a context.
"""
from __future__ import annotations

import csv
from pathlib import Path

from dna_decode.data.pneumo_breakpoints import classify

RULE_STATUS = "KNOWLEDGE_BASELINE"
RULE_SCOPE = "organism_routed"
ORGANISM = "Streptococcus_pneumoniae"

# CDC lookup column -> our breakpoint drug name.
DRUG_COL = {"penicillin": "PEN", "meropenem": "MER", "cefotaxime": "TAX",
            "ceftriaxone": "CFT", "cefuroxime": "CFX"}
SUPPORTED_DRUGS = tuple(DRUG_COL)


def apt_key(pbp1a, pbp2b, pbp2x) -> str:
    """(pbp1a,pbp2b,pbp2x) -> 'pbp1a-pbp2b-pbp2x' (the APT key in Ref_PBPtype_MIC.csv)."""
    return f"{str(pbp1a).strip()}-{str(pbp2b).strip()}-{str(pbp2x).strip()}"


def load_pbp_mic_table(csv_path: str | Path) -> dict[str, dict[str, float]]:
    """Parse Ref_PBPtype_MIC.csv -> {APT: {drug_col: mic_float}}; 'NA'/blank dropped."""
    table: dict[str, dict[str, float]] = {}
    with Path(csv_path).open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            apt = (row.get("APT") or "").strip()
            if not apt:
                continue
            mics = {}
            for col in DRUG_COL.values():
                v = (row.get(col) or "").strip()
                try:
                    mics[col] = float(v)
                except ValueError:
                    pass
            table[apt] = mics
    return table


def predict_mic(table: dict, pbp1a, pbp2b, pbp2x, drug: str) -> float | None:
    """Looked-up MIC for `drug`, or None (novel PBP type / not in table / drug out of scope)."""
    col = DRUG_COL.get((drug or "").strip().lower())
    if col is None:
        return None
    apt = apt_key(pbp1a, pbp2b, pbp2x)
    if "NEW" in apt.upper():                      # novel PBP allele -> no-call (never fabricate)
        return None
    return table.get(apt, {}).get(col)


def predict_rs(table: dict, pbp1a, pbp2b, pbp2x, drug: str, context: str) -> dict:
    """R/I/S call for `drug` at breakpoint `context`. Returns {prediction, mic, apt, ...} (prediction None
    on no-call). `context` in {meningitis, non_meningitis, oral} — REQUIRED (β-lactam R/S is ambiguous)."""
    apt = apt_key(pbp1a, pbp2b, pbp2x)
    mic = predict_mic(table, pbp1a, pbp2b, pbp2x, drug)
    rs = classify(drug, context, mic) if mic is not None else None
    return {"organism": ORGANISM, "rule_status": RULE_STATUS, "drug": drug, "context": context,
            "apt": apt, "mic": mic, "prediction": rs}
