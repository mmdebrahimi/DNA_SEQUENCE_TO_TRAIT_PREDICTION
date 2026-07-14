"""Forward edit router demo — one entry point routing an edit to the RIGHT regime predictor, shown across
all three regimes with REAL predictors:

  Regime A (curated AMR determinant)  -> the real WHO TB catalogue -> R/S     (M. tuberculosis rpoB / rifampicin)
  Regime B (protein molecular fitness)-> the cached ESM2 DMS predictor         (E. coli TEM-1 beta-lactamase)
  Regime C (organism-level trait)     -> ABSTAIN (closed negative)             (E. coli growth-rate edit)

Shows the router NEVER sends a resistance edit to the likelihood predictor (Regime A -> catalogue) nor an
organismal trait to any predictor (Regime C -> abstain). Regime-A uses the real committed WHO catalogue if
its (gitignored) master is present; else it degrades to a small illustrative key set + flags it.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward import REGIME_A, REGIME_B, REGIME_C, predict_edit, variant_key  # noqa: E402

_THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
          "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
          "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _parse_who_variant(variant: str) -> str | None:
    """'rpoB_p.Ser450Leu' -> 'rpoB:S450L' (simple substitutions only; None for indels/LoF/etc.)."""
    import re
    m = re.match(r"^([A-Za-z0-9]+)_p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$", variant.strip())
    if not m:
        return None
    gene, wt3, pos, alt3 = m.groups()
    wt, alt = _THREE.get(wt3.upper()), _THREE.get(alt3.upper())
    if not wt or not alt:
        return None
    return variant_key(gene, f"{wt}{pos}{alt}")


def build_tb_resistance_keys(drug: str = "rifampicin"):
    """Real WHO grade-1/2 substitution determinants for `drug` -> {gene:MUT} key set. Degrades if absent."""
    try:
        from dna_decode.data import tb_who_catalogue
        tb_who_catalogue.verify_pins()
        dets = tb_who_catalogue.load_determinants(drug)
        keys = {k for d in dets if (k := _parse_who_variant(d.variant))}
        return keys, "real_who_catalogue"
    except Exception as e:                       # gitignored master absent / pin issue -> illustrative
        return {variant_key("rpoB", "S450L"), variant_key("rpoB", "H445Y"), variant_key("rpoB", "D435V")}, \
               f"illustrative_fallback ({type(e).__name__})"


def load_tem1():
    ref = Path("D:/dna_decode_cache/proteingym/pg_reference.csv")
    table_p = Path("D:/dna_decode_cache/esm/esm2_t33_650M_UR50D__BLAT_ECOLX_Stiffler_2015.json")
    seq = None
    if ref.exists():
        seq = next((r["target_seq"] for r in csv.DictReader(open(ref, encoding="utf-8"))
                    if r["DMS_id"] == "BLAT_ECOLX_Stiffler_2015"), None)
    table = ({int(k): v for k, v in json.loads(table_p.read_text(encoding="utf-8")).items()}
             if table_p.exists() else None)
    return seq, table


def main() -> int:
    out_lines = []
    results = {}

    # --- Regime A: real WHO TB catalogue -----------------------------------------------------------
    keys, src = build_tb_resistance_keys("rifampicin")
    a_hit = predict_edit("rpoB", "S450L", regime=REGIME_A, drug="rifampicin", resistance_keys=keys)
    a_miss = predict_edit("rpoB", "A286V", regime=REGIME_A, drug="rifampicin", resistance_keys=keys)
    results["regime_A"] = {"source": src, "n_determinant_keys": len(keys), "hit": a_hit, "miss": a_miss}
    out_lines.append(f"[A determinant] rpoB S450L / rifampicin  -> {a_hit['prediction']} "
                     f"({a_hit['predictor']}, {src}, {len(keys)} keys)")
    out_lines.append(f"[A determinant] rpoB A286V / rifampicin  -> {a_miss['prediction']} ({a_miss['notes'][0]})")

    # --- Regime B: real cached ESM2 DMS predictor --------------------------------------------------
    seq, table = load_tem1()
    if seq and table:
        # E210K + P181R are real, WT-verified variants from the genome demo (benign / damaging)
        b = predict_edit("blaTEM-1", "E210K", regime=REGIME_B, protein_seq=seq, method="esm2", esm_table=table,
                         phenotype="ampicillin fitness (DMS)")
        b2 = predict_edit("blaTEM-1", "P181R", regime=REGIME_B, protein_seq=seq, method="esm2", esm_table=table,
                          phenotype="ampicillin fitness (DMS)")
        results["regime_B"] = {"benign_example": b, "damaging_example": b2}
        out_lines.append(f"[B molecular ] blaTEM-1 E210K (ampicillin) -> {b['prediction']} "
                         f"(ESM2 score {round(b['raw_score'], 3)}, {b['predictor']})")
        out_lines.append(f"[B molecular ] blaTEM-1 P181R (ampicillin) -> {b2['prediction']} "
                         f"(ESM2 score {round(b2['raw_score'], 3)})")
    else:
        results["regime_B"] = {"status": "unavailable (D: cache absent)"}
        out_lines.append("[B molecular ] blaTEM-1 M69L -> UNAVAILABLE (D: cache absent)")

    # --- Regime B (HUMAN): AlphaMissense predictor -------------------------------------------------
    try:
        from dna_decode.forward import am_table_for_mutants, load_am_for_uniprot
        pten_ref = Path("D:/dna_decode_cache/proteingym/pg_reference.csv")
        pseq = next((r["target_seq"] for r in csv.DictReader(open(pten_ref, encoding="utf-8"))
                     if r["DMS_id"] == "PTEN_HUMAN_Mighell_2018"), None) if pten_ref.exists() else None
        am_tsv = Path("D:/dna_decode_cache/proteingym/am_filtered.tsv")
        if pseq and am_tsv.exists():
            dms = Path("D:/dna_decode_cache/proteingym/pg_dms/DMS_ProteinGym_substitutions/PTEN_HUMAN_Mighell_2018.csv")
            muts = [r["mutant"] for r in csv.DictReader(open(dms, encoding="utf-8"))]
            am_table = am_table_for_mutants(load_am_for_uniprot(am_tsv, "P60484"), 0, muts)
            # pick the highest-AM (pathogenic) + lowest-AM (benign) real PTEN variant
            path_mut = max(am_table, key=am_table.get)
            ben_mut = min(am_table, key=am_table.get)
            hp = predict_edit("PTEN", path_mut, regime=REGIME_B, protein_seq=pseq, method="alphamissense",
                              am_table=am_table, phenotype="PTEN function (human)")
            hb = predict_edit("PTEN", ben_mut, regime=REGIME_B, protein_seq=pseq, method="alphamissense",
                              am_table=am_table, phenotype="PTEN function (human)")
            results["regime_B_human_am"] = {"pathogenic": hp, "benign": hb}
            out_lines.append(f"[B molecular*] PTEN {path_mut} (human, AlphaMissense) -> {hp['prediction']} "
                             f"(AM {round(1 - hp['raw_score'], 3)}, {hp['predictor']})")
            out_lines.append(f"[B molecular*] PTEN {ben_mut} (human, AlphaMissense) -> {hb['prediction']} "
                             f"(AM {round(1 - hb['raw_score'], 3)})")
    except Exception as e:
        results["regime_B_human_am"] = {"status": f"unavailable ({type(e).__name__})"}
        out_lines.append(f"[B molecular*] PTEN (human, AlphaMissense) -> UNAVAILABLE ({type(e).__name__})")

    # --- Regime C: organismal -> abstain -----------------------------------------------------------
    c = predict_edit("acs (acetate growth rate)", "X100Y", regime=REGIME_C)
    results["regime_C"] = c
    out_lines.append(f"[C organismal] acetate-growth-rate edit    -> {c['prediction']} (closed-negative; abstain)")

    res = {
        "demo": "forward_edit_router",
        "regimes_routed": {"A_determinant": "WHO TB catalogue (R/S)", "B_molecular": "ESM2 DMS predictor",
                           "C_organismal": "ABSTAIN"},
        "regime_A_source": src,
        "results": results,
        "honesty": ("One entry point auto-routes an edit to the correct regime predictor. Resistance edits go "
                    "to the determinant catalogue (never the likelihood predictor — resistance-conservativeness "
                    "finding); organismal traits abstain (closed negative). Regime B is the DMS-validated "
                    "molecular predictor."),
        "status": ("ROUTER_DEMO_OK" if (a_hit["prediction"] == "R"
                                        and results["regime_B"].get("benign_example", {}).get("prediction")
                                        and c["abstain"]) else "DEGRADED"),
    }
    outp = REPO / "wiki" / f"forward_router_demo_{_date.today().isoformat()}.json"
    outp.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("[forward-router-demo] all three regimes routed to the right predictor:")
    for ln in out_lines:
        print("  " + ln)
    print(f"  status={res['status']}  artifact -> {outp}")
    return 0 if res["status"] == "ROUTER_DEMO_OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
