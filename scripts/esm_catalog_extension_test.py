"""Phase-1 test of the masked-protein-LM 'extend the deterministic catalog' hypothesis (SARS-CoV-2 Mpro).

The user's masked-prediction idea, at the granularity where it is CAUSAL (protein, not genotype — see the
purifying-selection-vs-LD distinction). Does a masked protein-LM (ESM) rank the deterministic Mpro catalog's
OWN MISSED true-resistant variants (the documented FN set: novel high-fold mutants the mutant-level catalog
does not carry) as MORE deleterious than known-benign Mpro polymorphisms? If yes -> a learned layer that fills
the catalog's holes = the hybrid.

Method (Meier 2021 masked-marginals): mask each variant position in the WT Mpro protein, one ESM forward pass,
score LLR = log p(mut) - log p(wt) at that position (more negative = more deleterious). No training; WT
reference is the committed NC_045512 Mpro CDS (H41/C145 catalytic + E166 nirmatrelvir verified). Labels =
CoV-RDB measured fold (the SAME independent label the deterministic cell was validated against).

Verdict WIN iff BOTH: (1) ESM orders variants by resistance magnitude (Spearman(-LLR, log fold) materially
negative); (2) ESM RESCUES the catalog blindspot — the catalog-missed true-R variants score more deleterious
than the benign negatives (AUC materially > 0.5). WEIGHTS -> D: (torch.hub set to D:; C: is disk-tight).
"""
from __future__ import annotations

import json
import re
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

_VAR = re.compile(r"^(\d+)([A-Z])$")
# Known-benign Mpro polymorphisms (the catalog correctly calls S; truly low/no fold). Sourced:
#   P132H = the near-universal Omicron nsp5 polymorphism (benign, no nirmatrelvir resistance; Ullrich 2022).
#   K90R  = a documented benign Mpro polymorphism. Used as the negative set (CoV-RDB is R-enriched -> no TN).
BENIGN_NEGATIVES = [("P", 132, "H"), ("K", 90, "R")]


def _mpro_wt(ref_fna: Path) -> str:
    from Bio.Seq import Seq
    seq = "".join(l.strip() for l in ref_fna.read_text().splitlines() if not l.startswith(">"))
    return str(Seq(seq).translate()).rstrip("*")


def _load_variants(validation_json: Path, wt: str) -> list[dict]:
    import statistics
    d = json.loads(validation_json.read_text())
    out = []
    for key, folds in d["per_mutation_fold"].items():
        m = _VAR.match(key.strip())
        if not m:
            continue
        pos, mut = int(m.group(1)), m.group(2)
        if not (1 <= pos <= len(wt)) or mut not in "ACDEFGHIKLMNPQRSTVWY":
            continue
        vals = [float(x) for x in folds if _is_num(x)]
        if not vals:
            continue
        out.append({"pos": pos, "wt": wt[pos - 1], "mut": mut, "median_fold": statistics.median(vals)})
    return out


def _is_num(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def _auc(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """AUC = P(a random positive score > a random negative score), ties=0.5."""
    if not pos_scores or not neg_scores:
        return None
    wins = ties = 0
    for p in pos_scores:
        for n in neg_scores:
            if p > n:
                wins += 1
            elif p == n:
                ties += 1
    return (wins + 0.5 * ties) / (len(pos_scores) * len(neg_scores))


def score_esm_masked_marginals(wt: str, positions: list[int], model_name: str) -> dict:
    """Return {pos: {aa: log-prob}} at each requested position via masked-marginals (one fwd pass per pos)."""
    import torch
    torch.hub.set_dir("D:/dna_decode_cache/torch")            # weights -> D: (C: disk-tight)
    import esm
    model, alphabet = getattr(esm.pretrained, model_name)()
    model.eval()
    bc = alphabet.get_batch_converter()
    mask_i = alphabet.mask_idx
    _, _, toks = bc([("wt", wt)])                              # [1, L+2] (BOS + seq + EOS)
    out: dict[int, dict[str, float]] = {}
    with torch.no_grad():
        for pos in sorted(set(positions)):
            masked = toks.clone()
            masked[0, pos] = mask_i                            # token index = pos (BOS at 0 -> residue p at idx p)
            logits = model(masked)["logits"]
            lp = torch.log_softmax(logits[0, pos], dim=-1)
            out[pos] = {aa: float(lp[alphabet.get_idx(aa)]) for aa in "ACDEFGHIKLMNPQRSTVWY"}
    return out


def run(ref_fna: Path, validation_json: Path, model_name: str, r_threshold: float = 2.5) -> dict:
    import math
    from scipy.stats import spearmanr
    wt = _mpro_wt(ref_fna)
    variants = _load_variants(validation_json, wt)
    benign = [{"pos": p, "wt": w, "mut": m, "median_fold": None} for (w, p, m) in BENIGN_NEGATIVES]
    positions = [v["pos"] for v in variants] + [b["pos"] for b in benign]
    lp = score_esm_masked_marginals(wt, positions, model_name)

    def llr(v):  # log p(mut) - log p(wt) at the position; more negative = more deleterious
        return lp[v["pos"]][v["mut"]] - lp[v["pos"]][v["wt"]]

    for v in variants + benign:
        v["esm_llr"] = round(llr(v), 3)

    scored = [v for v in variants if v["median_fold"] is not None]
    rho = None
    if len(scored) >= 5:
        rho = float(spearmanr([-v["esm_llr"] for v in scored],
                              [math.log10(v["median_fold"]) for v in scored])[0])
    # Catalog blindspot rescue: the catalog-missed true-R (FN) vs benign negatives. -LLR = deleteriousness.
    fn_R = [v for v in scored if v["median_fold"] >= r_threshold]     # the catalog's FN blindspot (all truly R)
    blindspot_auc = _auc([-v["esm_llr"] for v in fn_R], [-b["esm_llr"] for b in benign])
    famous = sorted(fn_R, key=lambda v: v["esm_llr"])[:6]
    win = (rho is not None and rho >= 0.3) and (blindspot_auc is not None and blindspot_auc >= 0.7)
    return {
        "artifact": "esm_catalog_extension_sarscov2_mpro", "schema": "esm-catalog-extension-v1",
        "date": _date.today().isoformat(), "model": model_name, "protein": "SARS-CoV-2 Mpro (nsp5, 306aa)",
        "label_source": "CoV-RDB measured fold (same independent label as the deterministic Mpro cell)",
        "method": "ESM masked-marginals (Meier 2021); LLR = log p(mut) - log p(wt); -LLR = deleteriousness",
        "n_variants_scored": len(scored), "n_fn_blindspot_R": len(fn_R), "n_benign_neg": len(benign),
        "spearman_negLLR_vs_log_fold": round(rho, 3) if rho is not None else None,
        "blindspot_rescue_auc": round(blindspot_auc, 3) if blindspot_auc is not None else None,
        "verdict": "WIN_ESM_EXTENDS_CATALOG" if win else "NO_WIN",
        "famous_misses_esm_llr": [{"var": f"{v['wt']}{v['pos']}{v['mut']}", "fold": v["median_fold"],
                                   "esm_llr": v["esm_llr"]} for v in famous],
        "benign_esm_llr": [{"var": f"{b['wt']}{b['pos']}{b['mut']}", "esm_llr": b["esm_llr"]} for b in benign],
        "honest_caveats": [
            "CoV-RDB Mpro is R-ENRICHED -> the benign negative set is small (2 sourced polymorphisms); "
            "the blindspot AUC is a FIRST-PASS indicator, not a powered number (HIV RT/PR/IN, rich in R+S, "
            "is the scale-up if this wins)",
            "ESM is scored on the SAME CoV-RDB fold labels the deterministic cell used -> in-distribution "
            "(the honest scale-up is a held-out / independent fold source, mirroring the cell's own arc)",
            "masked-marginals (not fine-tuned) -> a pure zero-shot masked-LM signal, exactly the 'predict the "
            "missing part' paradigm the user proposed",
        ],
        "next": ("PHASE-1 WIN -> greenlight HIV scale-up + Phase-2 genotype world-model consideration"
                 if win else "PHASE-1 NO-WIN -> masked protein-LM does not clearly rescue the Mpro blindspot here"),
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", type=Path, default=REPO / "data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna")
    ap.add_argument("--validation-json", type=Path, default=REPO / "wiki/sarscov2_mpro_cov_rdb_validation.json")
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / f"esm_catalog_extension_mpro_{_date.today().isoformat()}.json")
    a = ap.parse_args(argv)
    res = run(a.ref, a.validation_json, a.model)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    print(f"\n[wrote {a.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
