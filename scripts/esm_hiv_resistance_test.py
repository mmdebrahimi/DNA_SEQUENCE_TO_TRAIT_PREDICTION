"""FAIR test of 'does a masked protein-LM carry AMR-resistance signal' — HIV, R-vs-S at attributable variants.

Fixes the underpowered SARS-Mpro probe (n=2 cross-position benign; within-R-only Spearman). Here the label is
CLEAN and per-variant: SINGLE-MUTANT isolates (exactly one non-consensus residue vs HXB2) whose PhenoSense
fold is therefore ATTRIBUTABLE to that one variant. Balanced R and S (esp. NNRTI/RT: ~33 R / 24 S). No LD
confound (single mutant), no catalog circularity (label = measured fold, not the catalog), matched at the
per-variant level.

Metric: AUC of (-LLR) separating R (fold>=3) from S (fold<3) single-mutant variants, per (class, drug) + a
pooled AUC. -LLR = ESM deleteriousness (log p(wt) - log p(mut), WT-marginals, Meier 2021). AUC ~0.5 => ESM is
BLIND to resistance (the mechanistic concern, now FAIRLY tested); AUC >~0.7 => ESM DOES carry resistance
signal => the learned-scoring branch of the hybrid is greenlit. Weights -> D:.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from scripts.esm_catalog_extension_test import _auc  # noqa: E402
from scripts.hiv_nnrti_validate import _parse_fold, load_rows  # noqa: E402

AAS = "ACDEFGHIKLMNPQRSTVWY"
R_CUTOFF = 3.0
# (class, dataset, gene-ref, gene-len, drug-column). Drugs chosen for single-mutant R+S balance:
# NNRTI EFV + NRTI 3TC (M184V) give clean R; PI NFV + INSTI RAL are lower-barrier than LPV/DTG.
CLASSES = [
    ("NNRTI", "NNRTI_DataSet.Full.txt", "RT", 560, "EFV"),
    ("NRTI",  "NRTI_DataSet.Full.txt",  "RT", 560, "3TC"),
    ("PI",    "PI_DataSet.Full.txt",    "PR", 99,  "NFV"),
    ("INSTI", "INI_DataSet.Full.txt",   "IN", 288, "RAL"),
]
REFS = {"RT": "HIV1_RT_HXB2_cds.fna", "PR": "HIV1_PR_HXB2_cds.fna", "IN": "HIV1_IN_HXB2_cds.fna"}


def _wt_protein(ref_dir: Path, gene: str) -> str:
    s = "".join(l.strip() for l in (ref_dir / REFS[gene]).read_text().splitlines() if not l.startswith(">"))
    return str(Seq(s).translate()).rstrip("*")


def _single_mutant_variants(rows, wt: str, plen: int, drugcol: str) -> list[dict]:
    """Isolates with EXACTLY one non-consensus residue vs HXB2 WT + an attributable fold."""
    out = []
    for r in rows:
        if r.get("Method", "").strip() != "PhenoSense":
            continue
        muts = []
        for pos in range(1, plen + 1):
            cell = (r.get(f"P{pos}") or "").strip()
            if cell in ("", "-", ".", "NA"):
                continue
            wtaa = wt[pos - 1] if pos - 1 < len(wt) else "?"
            for aa in cell:
                if aa.isalpha() and aa != wtaa and aa != "-" and aa in AAS and wtaa in AAS:
                    muts.append((pos, wtaa, aa))
        if len(muts) != 1:
            continue
        fold = _parse_fold(r.get(drugcol, ""))
        if fold is None or fold <= 0:
            continue
        pos, wtaa, mut = muts[0]
        out.append({"pos": pos, "wt": wtaa, "mut": mut, "fold": fold, "label_R": fold >= R_CUTOFF})
    return out


def _wt_marginal_logprobs(wt: str, model_name: str) -> dict:
    """One ESM forward pass over the WT protein -> {pos(1-based): {aa: log-prob}} at every position."""
    import torch
    torch.hub.set_dir("D:/dna_decode_cache/torch")
    import esm
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    _, _, toks = alphabet.get_batch_converter()([("wt", wt)])
    with torch.no_grad():
        logits = model(toks)["logits"][0]                    # [L+2, vocab]; BOS at idx0 -> residue p at idx p
        lp = torch.log_softmax(logits, dim=-1)
    idx = {aa: alphabet.get_idx(aa) for aa in AAS}
    return {p: {aa: float(lp[p, idx[aa]]) for aa in AAS} for p in range(1, len(wt) + 1)}


def run(ref_dir: Path, data_dir: Path, model_name: str) -> dict:
    cache: dict[str, dict] = {}                              # gene -> wt-marginal logprobs (score each protein once)
    per_class = {}
    pooled: list[tuple[float, bool]] = []                    # (-LLR, label_R) across ALL single-mutant variants
    for cls, fname, gene, plen, drug in CLASSES:
        path = data_dir / fname
        if not path.exists():
            per_class[cls] = {"status": "NO_DATA"}
            continue
        wt = _wt_protein(ref_dir, gene)
        variants = _single_mutant_variants(load_rows(path), wt, plen, drug)
        if gene not in cache:
            cache[gene] = _wt_marginal_logprobs(wt, model_name)
        lp = cache[gene]
        for v in variants:
            v["esm_llr"] = round(lp[v["pos"]][v["mut"]] - lp[v["pos"]][v["wt"]], 3)
        R = [v for v in variants if v["label_R"]]
        S = [v for v in variants if not v["label_R"]]
        auc = _auc([-v["esm_llr"] for v in R], [-v["esm_llr"] for v in S])
        pooled.extend([(-v["esm_llr"], v["label_R"]) for v in variants])
        per_class[cls] = {
            "gene": gene, "drug": drug, "n_single_mutant": len(variants), "n_R": len(R), "n_S": len(S),
            "resistance_auc": round(auc, 3) if auc is not None else None,
            "note": "under-powered" if (len(R) < 8 or len(S) < 8) else "powered",
            "top_R": [{"var": f"{v['wt']}{v['pos']}{v['mut']}", "fold": v["fold"], "esm_llr": v["esm_llr"]}
                      for v in sorted(R, key=lambda x: -x["fold"])[:5]],
        }
    pooled_R = [s for s, lab in pooled if lab]
    pooled_S = [s for s, lab in pooled if not lab]
    pooled_auc = _auc(pooled_R, pooled_S)
    powered = {c: v for c, v in per_class.items() if v.get("note") == "powered" and v.get("resistance_auc") is not None}
    best = max((v["resistance_auc"] for v in powered.values()), default=None)
    win = best is not None and best >= 0.70
    return {
        "artifact": "esm_hiv_resistance_matched_test", "schema": "esm-hiv-resistance-v1",
        "date": _date.today().isoformat(), "model": model_name,
        "design": "single-mutant isolates (attributable per-variant fold); AUC(-LLR separates R from S); "
                  "WT-marginals; label = PhenoSense fold>=3 (independent wet-lab)",
        "verdict": "ESM_CARRIES_RESISTANCE_SIGNAL" if win else "ESM_BLIND_TO_RESISTANCE",
        "best_powered_class_auc": best, "n_pooled": len(pooled),
        "pooled_R": len(pooled_R), "pooled_S": len(pooled_S),
        "pooled_resistance_auc": round(pooled_auc, 3) if pooled_auc is not None else None,
        "per_class": per_class,
        "interpretation": (
            "This is the FAIR test the Mpro probe was not: balanced R+S, attributable per-variant fold, no LD "
            "confound, matched at the per-variant level. AUC ~0.5 = ESM (evolutionary fitness) is blind to "
            "resistance -> the ZERO-SHOT extension idea fails, and a SUPERVISED head is the only surviving "
            "learned-scoring form. AUC >~0.7 = zero-shot representation DOES carry resistance signal -> the "
            "learned-scoring hybrid branch is empirically greenlit."),
        "honest_caveats": [
            "single-mutant restriction shrinks N but is what makes the fold ATTRIBUTABLE (the fair-label price)",
            "fold>=3 is an illustrative uniform cutoff (per-drug DRMcv cutoffs would refine but not flip AUC)",
            "in-distribution (HIVDB fold); an independent-cohort fold is the eventual scale-up, as elsewhere",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-dir", type=Path, default=REPO / "data/hiv_ref")
    ap.add_argument("--data-dir", type=Path, default=REPO / "data/raw/hiv")
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / f"esm_hiv_resistance_matched_test_{_date.today().isoformat()}.json")
    a = ap.parse_args(argv)
    res = run(a.ref_dir, a.data_dir, a.model)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    print(f"\n[wrote {a.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
