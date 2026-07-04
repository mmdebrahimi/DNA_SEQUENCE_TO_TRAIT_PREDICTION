"""Deterministic ClinVar-backed human variant decoder — the Mendelian-disease analogue of the AMR catalog.

Given a variant (chrom, pos, ref, alt on GRCh38), return ClinVar's curated clinical significance + star
review level + gene + disease, with provenance. This is a CURATED-CATALOG deterministic decoder (regime-1:
curated catalog → deterministic wins), exactly like the AMR determinant catalog — NOT a learned predictor.
Fail-closed + honest: a variant NOT in the committed panel returns INDETERMINATE ("not-in-panel"), never a
guess. The panel is a committed subset (`data/clinvar/clinvar_panel.tsv`, P/LP + B/LB for a canonical gene
set, built by `scripts/capture_clinvar.py` from the full ClinVar VCF on D:); extensible to any gene panel.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

INDETERMINATE = "INDETERMINATE"
_PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
_BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
# ClinVar review status -> gold-star rating (ClinVar's own scale)
_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
}
_DEFAULT_PANEL = Path(__file__).resolve().parent.parent.parent / "data" / "clinvar" / "clinvar_panel.tsv"


def _verdict(sig: str) -> str:
    if sig in _PATHOGENIC:
        return "PATHOGENIC"
    if sig in _BENIGN:
        return "BENIGN"
    return INDETERMINATE


@dataclass(frozen=True)
class ClinVarCall:
    significance: str
    verdict: str                 # PATHOGENIC / BENIGN / INDETERMINATE
    stars: int | None
    gene: str | None
    disease: str | None
    provenance: str


class ClinVarDecoder:
    def __init__(self, table: dict):
        self.table = table       # {(chrom,pos,ref,alt): {significance, review_status, gene, disease, clinvar_id}}

    @classmethod
    def from_tsv(cls, path: str | Path = _DEFAULT_PANEL) -> "ClinVarDecoder":
        table = {}
        with open(path, encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            i = {k: hdr.index(k) for k in hdr}
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < len(hdr):
                    continue
                key = (str(c[i["chrom"]]), str(c[i["pos"]]), c[i["ref"]].upper(), c[i["alt"]].upper())
                table[key] = {"significance": c[i["significance"]], "review_status": c[i["review_status"]],
                              "gene": c[i["gene"]], "disease": c[i["disease"]], "clinvar_id": c[i["clinvar_id"]]}
        return cls(table)

    def call(self, chrom, pos, ref: str, alt: str) -> ClinVarCall:
        chrom = str(chrom).replace("chr", "")
        row = self.table.get((chrom, str(pos), ref.upper(), alt.upper()))
        if row is None:
            return ClinVarCall(INDETERMINATE, INDETERMINATE, None, None, None,
                               "not-in-panel (ClinVar committed gene panel; absence != benign — fetch full VCF to extend)")
        sig = row["significance"]
        return ClinVarCall(
            significance=sig, verdict=_verdict(sig),
            stars=_STARS.get(row["review_status"]), gene=row["gene"], disease=row["disease"],
            provenance=f"ClinVar {row['clinvar_id']} ({row['review_status']}) — curated germline classification")
