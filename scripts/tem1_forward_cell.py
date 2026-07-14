"""TEM-1 beta-lactamase forward variant-effect cell — validate the Regime-B forward predictor
(dna_decode/forward) per-variant against the FREE wet-lab DMS ampicillin-fitness landscape (ProteinGym).

This is the FORWARD direction of the decoder on E. coli: make a minor edit (any point mutation in TEM-1
beta-lactamase) -> predict the change in phenotype (ampicillin growth fitness) -> validate against the
MEASURED effect for that exact variant. The label is free + per-variant + independent (deep mutational
scan), so the project's label wall does NOT bind here.

Deterministic BLOSUM62 baseline (no GPU / no network). Honest headline = signed Spearman rank correlation
between the predictor's continuous score and the measured DMS fitness, over all single-point variants whose
WT residue matches the reference sequence (a WT mismatch is a coordinate/frame error and is COUNTED, never
silently scored). ESM2 zero-shot is the drop-in upgrade (scripts/esm_zeroshot_dms.py, median ~0.49).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward import parse_mutation, predict_effect  # noqa: E402

PG = Path("D:/dna_decode_cache/proteingym")
DMS_DIR = PG / "pg_dms" / "DMS_ProteinGym_substitutions"


def load_target_seq(dms_id: str) -> str:
    with open(PG / "pg_reference.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("DMS_id") == dms_id:
                return r.get("target_seq", "")
    raise KeyError(f"{dms_id} not found in pg_reference.csv")


def spearman(xs: list[float], ys: list[float]) -> float:
    """Signed Spearman rank correlation (pure — Pearson on average ranks; ties averaged)."""
    n = len(xs)
    if n < 3:
        return float("nan")

    def ranks(vals):
        order = sorted(range(n), key=lambda i: vals[i])
        rk = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def run(dms_id: str, protein_label: str, phenotype_axis: str,
        method: str = "blosum62", esm_table: dict | None = None, am_table: dict | None = None) -> dict:
    seq = load_target_seq(dms_id)
    path = DMS_DIR / f"{dms_id}.csv"
    pred_scores: list[float] = []
    dms_scores: list[float] = []
    wt_mismatch = 0
    n_multi = 0
    nonsense_dms: list[float] = []
    # coarse confusion using DMS_score_bin (1 = functional-ish per ProteinGym binarization) if present
    tp = fp = tn = fn = 0
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mut = (r.get("mutant") or "").strip()
            if ":" in mut:                       # multi-mutant — this cell scores single edits only
                n_multi += 1
                continue
            try:
                dms = float(r.get("DMS_score"))
            except (TypeError, ValueError):
                continue
            try:
                p = predict_effect(seq, mut, protein=protein_label, phenotype_axis=phenotype_axis,
                                   method=method, esm_table=esm_table, am_table=am_table)
            except ValueError as e:
                if "WT mismatch" in str(e):
                    wt_mismatch += 1
                    continue
                if "not AlphaMissense-covered" in str(e):
                    continue                 # variant not AM-covered -> skip (reported via match count)
                raise
            pred_scores.append(p.raw_score)
            dms_scores.append(dms)
            _, _, alt = parse_mutation(mut)
            if alt in ("*", "X"):
                nonsense_dms.append(dms)
            b = (r.get("DMS_score_bin") or "").strip()
            if b in ("0", "1"):
                measured_fn = (b == "1")
                pred_fn = (p.predicted_effect == "preserved")
                tp += measured_fn and pred_fn
                fp += (not measured_fn) and pred_fn
                tn += (not measured_fn) and (not pred_fn)
                fn += measured_fn and (not pred_fn)

    rho = spearman(pred_scores, dms_scores)
    n = len(pred_scores)
    # polarity anchor: nonsense variants must sit at the LOW-fitness extreme (confirms sign expectation)
    nonsense_mean = sum(nonsense_dms) / len(nonsense_dms) if nonsense_dms else None
    overall_mean = sum(dms_scores) / n if n else None
    res = {
        "cell": "tem1_forward_variant_effect",
        "protein": protein_label,
        "dms_id": dms_id,
        "phenotype_axis": phenotype_axis,
        "method": {"blosum62": "blosum62_deterministic", "esm2": "esm2_zeroshot",
                   "alphamissense": "alphamissense_learned"}.get(method, method),
        "n_single_variants_scored": n,
        "n_multi_mutant_skipped": n_multi,
        "n_wt_mismatch": wt_mismatch,          # MUST be 0 — else a coordinate/frame error
        "spearman_pred_vs_dms": round(rho, 4),
        "sign_expectation": "positive (conservative substitution -> preserved ampicillin fitness)",
        "nonsense_mean_dms": (round(nonsense_mean, 4) if nonsense_mean is not None else None),
        "overall_mean_dms": (round(overall_mean, 4) if overall_mean is not None else None),
        # direction confirmed if the correlation is positive (matching sign_expectation); if the assay HAS
        # nonsense variants, they must additionally sit below the overall-fitness mean.
        "polarity_ok": (rho > 0 and (nonsense_mean is None or (overall_mean is not None
                                                               and nonsense_mean < overall_mean))),
        "nonsense_anchor_available": nonsense_mean is not None,
        "coarse_confusion_on_dms_bin": {"tp": tp, "fp": fp, "tn": tn, "fn": fn},
        "honesty": ("Regime-B molecular-fitness forward predictor; deterministic BLOSUM62 baseline. Validated "
                    "per-variant on FREE independent DMS. NOT an organism-level predictor (Regime C abstains); "
                    "for clinical-resistance calls use the Regime-A determinant catalogue."),
        "status": "FORWARD_CELL_DMS_VALIDATED" if (n >= 100 and wt_mismatch == 0) else "DEGRADED",
    }
    return res


def _build_or_load_esm_table(seq: str, dms_id: str, model: str) -> dict:
    """ESM2 masked-marginal {pos:{aa:logp}} table for `seq`, cached to D: (build is the slow CPU step)."""
    import json as _json
    cache_dir = Path("D:/dna_decode_cache/esm")
    cache_dir.mkdir(parents=True, exist_ok=True)
    tag = model.split("/")[-1]
    cache = cache_dir / f"{tag}__{dms_id}.json"
    if cache.exists():
        raw = _json.loads(cache.read_text(encoding="utf-8"))
        print(f"[tem1-forward] loaded cached ESM table {cache.name} ({len(raw)} positions)")
        return {int(k): v for k, v in raw.items()}
    from dna_decode.forward.esm_scorer import esm2_logp_table
    print(f"[tem1-forward] building ESM masked-marginal table ({len(seq)} positions, model={tag}, CPU) ...")
    table = esm2_logp_table(seq, model_name=model)
    cache.write_text(_json.dumps({str(k): v for k, v in table.items()}), encoding="utf-8")
    print(f"[tem1-forward] wrote ESM table cache -> {cache}")
    return table


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dms-id", default="BLAT_ECOLX_Stiffler_2015")
    ap.add_argument("--protein", default="TEM-1 beta-lactamase (BLAT_ECOLX)")
    ap.add_argument("--phenotype", default="ampicillin growth fitness (10-2500 ug/mL; DMS-measured)")
    ap.add_argument("--method", default="blosum62", choices=["blosum62", "esm2", "alphamissense"])
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D", help="ESM2 HF model id (--method esm2)")
    ap.add_argument("--uniprot", default=None, help="UniProt accession (--method alphamissense)")
    ap.add_argument("--offset", type=int, default=0, help="DMS->UniProt position offset (--method alphamissense)")
    a = ap.parse_args(argv)
    if not (DMS_DIR / f"{a.dms_id}.csv").exists():
        print(f"ERROR: DMS assay {a.dms_id}.csv not found under {DMS_DIR}", file=sys.stderr)
        return 2
    esm_table = am_table = None
    if a.method == "esm2":
        esm_table = _build_or_load_esm_table(load_target_seq(a.dms_id), a.dms_id, a.model)
    elif a.method == "alphamissense":
        if not a.uniprot:
            print("ERROR: --method alphamissense requires --uniprot <accession>", file=sys.stderr)
            return 2
        from dna_decode.forward.am_scorer import am_table_for_mutants, load_am_for_uniprot
        am_tsv = PG / "am_filtered.tsv"
        am_by_variant = load_am_for_uniprot(am_tsv, a.uniprot)
        mutants = [r["mutant"] for r in csv.DictReader(open(DMS_DIR / f"{a.dms_id}.csv", encoding="utf-8"))]
        am_table = am_table_for_mutants(am_by_variant, a.offset, mutants)
        print(f"[tem1-forward] AlphaMissense: {len(am_by_variant)} AM variants for {a.uniprot}; "
              f"{len(am_table)} of the assay's mutants covered")
    res = run(a.dms_id, a.protein, a.phenotype, method=a.method, esm_table=esm_table, am_table=am_table)
    suffix = "" if a.method == "blosum62" else f"_{a.method}"
    out = REPO / "wiki" / f"tem1_forward_cell_{a.dms_id.lower()}{suffix}_{_date.today().isoformat()}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[tem1-forward] {a.dms_id}: n={res['n_single_variants_scored']} "
          f"wt_mismatch={res['n_wt_mismatch']} | Spearman({res['method']},DMS)={res['spearman_pred_vs_dms']} "
          f"| polarity_ok={res['polarity_ok']} | status={res['status']}")
    print(f"  nonsense_mean_dms={res['nonsense_mean_dms']} vs overall={res['overall_mean_dms']} "
          f"(nonsense should be lower)")
    print(f"  coarse confusion on DMS_bin: {res['coarse_confusion_on_dms_bin']}")
    print(f"artifact -> {out}")
    return 0 if res["status"] == "FORWARD_CELL_DMS_VALIDATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
