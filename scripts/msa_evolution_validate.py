"""Real-surface validation of the MSA->evolution-score pipeline (dna_decode/forward/msa_evolution.py).

On a REAL ProteinGym protein, using the on-D: MSA + ProteinGym's precomputed sequence weights, this proves:
  (a) CORRECTNESS -- my site-independent score reproduces ProteinGym's own `Site_Independent` column
      (Spearman >= --min-repro over the shared match-column variants). If my MSA parse / reweighting /
      log-odds are right, the two agree strongly.
  (b) END-TO-END -- the ESM2 (+) my-evolution hybrid runs through `rank_average_hybrid` and yields a Spearman
      vs the measured DMS. (The site-independent floor is NOT expected to beat ESM2 -- that's the measured
      R2 finding; this checks the PIPE is wired, not that the floor model lifts.)

Uses the ESM2-650M table on D: (from the Kaggle run) when present; else skips arm (b).

Run:  uv run python scripts/msa_evolution_validate.py --dms CASP3_HUMAN_Roychowdhury_2020
Exit: 0 = validated (repro >= bar); 1 = repro below bar; 2 = substrate unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward.msa_evolution import site_independent_table  # noqa: E402
from dna_decode.forward.variant_effect import rank_average_hybrid    # noqa: E402
from scripts.forward_blosum_proteingym_sweep import _spearman        # noqa: E402

PG = Path("D:/dna_decode_cache/proteingym")
MSA_DIR = PG / "pg_msa" / "DMS_msa_files"
ZS_DIR = PG / "pg_zeroshot"
REF = PG / "pg_reference.csv"
WEIGHTS_ZIP = PG / "pg_msa_weights.zip"
ESM_DIR = Path("D:/dna_decode_cache/esm")


def _ref_row(dms: str) -> dict:
    with REF.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["DMS_id"] == dms:
                return row
    raise SystemExit(f"{dms} not in reference")


def _extract_weight_npy(dms: str, theta: float, scratch: Path) -> Path | None:
    """Pull the protein's precomputed weight .npy out of pg_msa_weights.zip (canonical reweighting)."""
    if not WEIGHTS_ZIP.exists():
        return None
    z = zipfile.ZipFile(WEIGHTS_ZIP)
    # weight files are named like PROT_HUMAN_theta0.2_<build>.npy under DMS_msa_weights/ (NOT the MSA_Transformer subdir)
    prot = dms.rsplit("_", 2)[0] if dms.count("_") >= 2 else dms
    cands = [n for n in z.namelist()
             if n.endswith(".npy") and "MSA_Transformer" not in n and prot.split("_")[0] in n]
    if not cands:
        return None
    scratch.mkdir(parents=True, exist_ok=True)
    target = scratch / Path(cands[0]).name
    target.write_bytes(z.read(cands[0]))
    return target


def _load_esm(dms: str):
    q = ESM_DIR / f"esm2_t33_650M_UR50D__{dms}.json"
    if not q.exists():
        return None
    raw = json.loads(q.read_text(encoding="utf-8"))
    return {int(k): v for k, v in raw.items()}


def _esm_variant_table(esm_tab: dict, muts: list[str]) -> dict[str, float]:
    """ESM2 delta logP(alt)-logP(wt) per mutation, higher=preserved -- one table for rank_average_hybrid."""
    out = {}
    for m in muts:
        wt, pos, alt = m[0], int(m[1:-1]), m[-1]
        col = esm_tab.get(pos)
        if col and wt in col and alt in col:
            out[m] = col[alt] - col[wt]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dms", default="CASP3_HUMAN_Roychowdhury_2020")
    ap.add_argument("--min-repro", type=float, default=0.7,
                    help="Spearman(my site-independent, ProteinGym Site_Independent) bar for CORRECTNESS")
    ap.add_argument("--scratch", default=str(REPO / "data" / "processed" / "msa_eval_scratch"))
    args = ap.parse_args(argv)

    row = _ref_row(args.dms)
    msa = MSA_DIR / row["MSA_filename"]
    zs = ZS_DIR / f"{args.dms}.csv"
    if not msa.exists() or not zs.exists():
        print(f"SUBSTRATE MISSING: msa={msa.exists()} zeroshot={zs.exists()}", file=sys.stderr)
        return 2

    from dna_decode.forward.msa_evolution import parse_a2m
    _, _, match_cols = parse_a2m(msa)
    depth = len(match_cols)
    theta = float(row.get("MSA_theta") or 0.2)
    wnpy = _extract_weight_npy(args.dms, theta, Path(args.scratch))
    # canonical .npy is for ProteinGym's coverage-FILTERED MSA; use it only if the depth matches,
    # else fall back to uniform (the deployable-anywhere floor -- a novel protein has no precomputed weights).
    if wnpy is not None:
        import numpy as np
        if len(np.load(wnpy)) != depth:
            print(f"[{args.dms}] canonical weights depth {len(np.load(wnpy))} != MSA {depth} "
                  f"(PG coverage-filter) -> UNIFORM weighting")
            wnpy = None
    weighting = "canonical .npy" if wnpy else "uniform"
    print(f"[{args.dms}] MSA={msa.name} depth={depth} weighting={weighting} theta={theta}")

    my = (site_independent_table(msa, weights_npy=wnpy) if wnpy
          else site_independent_table(msa, weights=[1.0] * depth))

    # join to the zeroshot CSV: DMS_score + PG Site_Independent + my score, over shared variants
    my_s, pg_si, dms_v, muts = [], [], [], []
    with zs.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            m = r["mutant"]
            if m in my and r.get("Site_Independent") not in (None, "", "NA") and r.get("DMS_score"):
                try:
                    pg = float(r["Site_Independent"]); d = float(r["DMS_score"])
                except ValueError:
                    continue
                my_s.append(my[m]); pg_si.append(pg); dms_v.append(d); muts.append(m)

    if len(muts) < 20:
        print(f"TOO FEW shared variants ({len(muts)})", file=sys.stderr)
        return 2

    repro = abs(_spearman(my_s, pg_si))          # (a) correctness vs ProteinGym's own column
    my_vs_dms = abs(_spearman(my_s, dms_v))
    pg_vs_dms = abs(_spearman(pg_si, dms_v))
    print(f"(a) CORRECTNESS  Spearman(my SiteIndep, PG Site_Independent) = {repro:.4f}  (n={len(muts)})")
    print(f"    sanity: my-vs-DMS {my_vs_dms:.4f}  |  PG-SiteIndep-vs-DMS {pg_vs_dms:.4f}")

    # (b) end-to-end hybrid pipe: ESM2 (+) my-evolution
    esm = _load_esm(args.dms)
    hyb_line = "(b) END-TO-END  ESM2 table unavailable on D: -- pipe arm skipped"
    if esm:
        esm_tab = _esm_variant_table(esm, muts)
        shared = [m for m in muts if m in esm_tab and m in my]
        if len(shared) >= 20:
            combined = rank_average_hybrid([esm_tab, {m: my[m] for m in shared}])
            idx = {m: i for i, m in enumerate(muts)}
            hy = [combined[m] for m in shared]
            dd = [dms_v[idx[m]] for m in shared]
            ee = [esm_tab[m] for m in shared]
            h_rho, e_rho = abs(_spearman(hy, dd)), abs(_spearman(ee, dd))
            hyb_line = (f"(b) END-TO-END  hybrid ran on n={len(shared)}: "
                        f"|Spearman| hybrid {h_rho:.4f} vs ESM2-alone {e_rho:.4f} "
                        f"(floor model -> lift not expected; pipe is wired)")
    print(hyb_line)

    ok = repro >= args.min_repro
    print(f"\nVERDICT: {'VALIDATED' if ok else 'REPRO BELOW BAR'} "
          f"(repro {repro:.4f} {'>=' if ok else '<'} {args.min_repro})")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
