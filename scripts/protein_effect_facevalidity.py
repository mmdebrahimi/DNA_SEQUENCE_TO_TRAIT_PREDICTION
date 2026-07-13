"""Face-validity anchor for the rung-2 mutation-effect predictor (2026-07-13).

Turns "how close are we to predicting a mutation's molecular effect" into a NUMBER, on real E. coli data.
Runs the predictor's ESM2-650M damage_llr over every single substitution in an E. coli ProteinGym DMS assay
and reports the Spearman correlation vs. the experimental effect (DMS_score). The honest framing (per the
design review):

  * PRIMARY anchor = ARGR_ECOLI_Tsuboyama_2023_1AOY (STABILITY assay) — clean rung-2 molecular function.
  * BLAT_ECOLX (organismal-fitness beta-lactamase, if present) is scored as a SEPARATE, explicitly-FLAGGED
    caveat example — it is under antibiotic-selection pressure and blurs toward AMR territory; it does NOT
    support the general molecular-function claim.
  * This is ZERO-SHOT BENCHMARK FACE-VALIDITY, not prospective validation: ESM has likely seen the WT
    sequence + homologs; the DMS effect LABELS are not in training, so the Spearman is a fair test of
    effect-RANKING on a known family, but not "prospective performance on a novel protein".

DMS_score convention: higher = MORE stable/fit. damage_llr: higher = MORE damaging. So the expected
correlation is NEGATIVE; we report spearman(-damage_llr, DMS_score) so a POSITIVE number = the scorer
ranks effects in the right direction.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.protein_effect import predictor as P  # noqa: E402

DMS_DIR = Path("D:/dna_decode_cache/proteingym/pg_dms/DMS_ProteinGym_substitutions")
CACHE_DIR = REPO / "data" / "processed" / "protein_effect_cache"
# (assay file stem, category, is-flagged-caveat)
ASSAYS = [
    ("ARGR_ECOLI_Tsuboyama_2023_1AOY", "Stability", False),
    ("BLAT_ECOLX_Deng_2012", "OrganismalFitness", True),
]


def _spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return round(num / (dx * dy), 4) if dx and dy else None


def _wt_sequence(rows):
    """Reconstruct WT from any single-substitution row: revert the mutation in its mutated_sequence."""
    r = rows[0]
    wt, pos, _ = P.parse_mutation(r["mutant"])
    s = r["mutated_sequence"]
    return s[: pos - 1] + wt + s[pos:]


def run_assay(stem: str, category: str, flagged: bool):
    path = DMS_DIR / f"{stem}.csv"
    if not path.exists():
        return {"assay": stem, "available": False, "note": f"DMS file absent at {path}"}
    rows = [r for r in csv.DictReader(open(path, encoding="utf-8")) if ":" not in r["mutant"]]
    seq = _wt_sequence(rows)
    # consistency: every mutant's WT residue must match the reconstructed WT
    for r in rows:
        wt, pos, _ = P.parse_mutation(r["mutant"])
        if pos > len(seq) or seq[pos - 1] != wt:
            return {"assay": stem, "available": True, "note": f"WT inconsistency at {r['mutant']}"}
    logp = P.masked_marginals(seq, cache_path=CACHE_DIR / f"{stem}.json", progress=True)
    dl, dms = [], []
    for r in rows:
        wt, pos, mut = P.parse_mutation(r["mutant"])
        dl.append(P.damage_llr(logp, pos, wt, mut))
        dms.append(float(r["DMS_score"]))
    # positive spearman(-damage_llr, DMS_score) => scorer ranks effects the right way
    rho = _spearman([-x for x in dl], dms)
    return {"assay": stem, "available": True, "category": category, "flagged_caveat": flagged,
            "n_mutations": len(rows), "seq_len": len(seq),
            "spearman_neg_damage_vs_dms": rho,
            "interpretation": ("positive = ESM damage rank agrees with the experimental effect direction; "
                               "compare to ProteinGym ESM2-650M ~0.49 median / ~0.52 stability")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    results = [run_assay(*x) for x in ASSAYS]
    primary = next((r for r in results if r.get("available") and not r.get("flagged_caveat")), None)
    res = {
        "artifact": "protein_effect_facevalidity", "schema": "protein-effect-facevalidity-v1",
        "date": str(_date.today()),
        "question": "Does the rung-2 predictor's ESM2-650M damage rank agree with real E. coli DMS effects "
                    "(zero-shot FACE-VALIDITY, not prospective validation)?",
        "primary_anchor": primary["assay"] if primary else None,
        "primary_spearman": primary.get("spearman_neg_damage_vs_dms") if primary else None,
        "honest_framing": ("ZERO-SHOT benchmark face-validity: WT + homologs likely in ESM pretraining; DMS "
                           "effect labels are NOT, so the Spearman is a fair effect-RANKING test on a known "
                           "family, but not prospective performance on a novel protein. Stability (ARGR) is "
                           "the clean rung-2 anchor; BLAT is a flagged fitness/selection caveat, NOT support "
                           "for general molecular-function scoring."),
        "results": results,
    }
    out = a.out or (REPO / "wiki" / f"protein_effect_facevalidity_{_date.today()}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\n=== rung-2 predictor face-validity on E. coli DMS (spearman(-damage_llr, DMS_score)) ===")
    for r in results:
        if not r.get("available"):
            print(f"  {r['assay']}: UNAVAILABLE ({r.get('note')})"); continue
        flag = "  [FLAGGED caveat: fitness/selection, not molecular-function evidence]" if r["flagged_caveat"] else ""
        print(f"  {r['assay']} ({r['category']}, n={r['n_mutations']}): spearman={r['spearman_neg_damage_vs_dms']}"
              f"{flag}")
    print(f"\nPRIMARY (stability) anchor: {res['primary_anchor']} spearman={res['primary_spearman']} "
          f"(vs ProteinGym ESM2-650M ~0.52 stability)")
    print(f"[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
