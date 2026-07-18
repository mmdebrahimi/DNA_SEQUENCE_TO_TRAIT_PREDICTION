"""Does ESM2 (+) OUR-OWN ProSST lift on ProteinGym? — closing the STRUCTURE path end-to-end.

The modality-hybrid sweep found `ESM2 + ProSST` +0.05 paired vs ESM2-650M (win 87%) using ProteinGym's
PRECOMPUTED `ProSST-2048` column — the biggest single modality lever. This runs ProSST OURSELVES
(`AI4Protein/ProSST-2048`) on AlphaFold structures (via the shared `fetch_alphafold_pdb`) or ProteinGym's
pre-quantized structure tokens, feeds the scores through OUR `rank_average_hybrid`, and checks (1) our ProSST
reproduces the column and (2) the hybrid lift reproduces — the structure analog of the MSA-Transformer run.

PRIMARY bar = reproduction of the ProSST-2048 column (correctness). SECONDARY = the per-category hybrid lift
(structure helps Stability/Expression). Honest-limit framing per the MSA-T lesson: report the FINAL n; the
powered per-category evidence is the N=95 precomputed-column sweep.

Needs the ProSST stack (transformers + prosst; torch_geometric only for the PdbQuantizer path) — runs on
Kaggle T4 / a Linux GPU. Restartable JSONL checkpoint (one line per protein).

Run:  uv run python scripts/prosst_lift.py --limit 12 [--structure-dir <pre-quantized-tokens-dir>]
Exit: 0 = ran; 2 = substrate/model unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.forward_blosum_proteingym_sweep import _spearman   # noqa: E402

PG = Path("D:/dna_decode_cache/proteingym")
ZS_DIR = PG / "pg_zeroshot"
REF = PG / "pg_reference.csv"
ESM_DIR = Path("D:/dna_decode_cache/esm")
STRUCT_CACHE = Path("D:/dna_decode_cache/alphafold")


# ---- pure analysis helpers (unit-tested; no model) -------------------------------------------------------

def paired_in_subset(records: list[dict], key: str = "hybrid_minus_esm") -> dict:
    """Paired median-delta + win-rate over OK records carrying `key`."""
    ds = [r[key] for r in records if r.get("status") == "OK" and key in r]
    wins = sum(1 for x in ds if x > 1e-9)
    return {"n": len(ds), "median": round(statistics.median(ds), 4) if ds else None,
            "win": f"{wins}/{len(ds)}" if ds else "0/0"}


def by_category(records: list[dict], categories: dict[str, str]) -> dict:
    """Per-phenotype-category paired hybrid lift, keyed by coarse_selection_type."""
    buckets: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        if r.get("status") == "OK":
            buckets[categories.get(r["dms"], "?")].append(r)
    return {c: paired_in_subset(rs) for c, rs in sorted(buckets.items(), key=lambda kv: -len(kv[1]))}


def reproduction_median(records: list[dict]) -> float | None:
    reps = [r["prosst_vs_pg_repro"] for r in records
            if r.get("status") == "OK" and r.get("prosst_vs_pg_repro") is not None]
    return round(statistics.median(reps), 4) if reps else None


# ---- per-protein scoring (needs the ProSST stack) --------------------------------------------------------

def _load_esm(dms: str):
    q = ESM_DIR / f"esm2_t33_650M_UR50D__{dms}.json"
    return {int(k): v for k, v in json.loads(q.read_text(encoding="utf-8")).items()} if q.exists() else None


def _esm_variant_table(esm: dict, muts: list[str]) -> dict[str, float]:
    out = {}
    for m in muts:
        wt, pos, alt = m[0], int(m[1:-1]), m[-1]
        col = esm.get(pos)
        if col and wt in col and alt in col:
            out[m] = col[alt] - col[wt]
    return out


def score_protein(dms: str, row: dict, structure_dir: Path | None, vocab: int) -> dict:
    from dna_decode.forward.prosst_scorer import prosst_variant_table, fetch_alphafold_pdb, quantize_structure
    from dna_decode.forward.variant_effect import rank_average_hybrid

    zs = ZS_DIR / f"{dms}.csv"
    esm = _load_esm(dms)
    uniprot = row.get("UniProt_ID", "")
    if not (zs.exists() and esm and uniprot):
        return {"dms": dms, "status": "SUBSTRATE_MISSING"}
    wt_seq = row["target_seq"].strip().upper()

    # pre-quantized tokens (transformer-only) preferred; else fetch + quantize (needs torch_geometric).
    # ProteinGym pre-quantized structures ship as {DMS_id}.fasta (header + comma-separated int tokens).
    tokens = None
    if structure_dir:
        tf = structure_dir / f"{dms}.fasta"
        if tf.exists():
            line = next(l for l in tf.read_text(encoding="utf-8").splitlines() if not l.startswith(">"))
            tokens = [int(x) for x in line.strip().split(",")]
    if tokens is None:
        pdb = fetch_alphafold_pdb(uniprot, STRUCT_CACHE)
        tokens = quantize_structure(pdb, vocab)

    # collect the DMS variants present in the zeroshot CSV + the PG ProSST column
    muts, dms_v, pg_prosst = [], [], {}
    with zs.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            m = r["mutant"]
            if len(m) < 3 or ":" in m or not m[1:-1].isdigit() or not r.get("DMS_score"):
                continue
            muts.append(m); dms_v.append(float(r["DMS_score"]))
            col = f"ProSST-{vocab}"
            if r.get(col) not in (None, "", "NA"):
                pg_prosst[m] = float(r[col])

    prosst_tab = prosst_variant_table(wt_seq, muts, structure_tokens=tokens, vocab=vocab)
    esm_tab = _esm_variant_table(esm, muts)
    shared = [m for m in muts if m in prosst_tab and m in esm_tab]
    if len(shared) < 20:
        return {"dms": dms, "status": "TOO_FEW", "n": len(shared)}

    idx = {m: i for i, m in enumerate(muts)}
    hyb = rank_average_hybrid([esm_tab, {m: prosst_tab[m] for m in shared}])
    dd = [dms_v[idx[m]] for m in shared]
    e_rho = abs(_spearman([esm_tab[m] for m in shared], dd))
    p_rho = abs(_spearman([prosst_tab[m] for m in shared], dd))
    h_rho = abs(_spearman([hyb[m] for m in shared], dd))
    repro = (abs(_spearman([prosst_tab[m] for m in shared if m in pg_prosst],
                           [pg_prosst[m] for m in shared if m in pg_prosst]))
             if len([m for m in shared if m in pg_prosst]) >= 20 else None)
    return {"dms": dms, "status": "OK", "n": len(shared),
            "esm_spearman": round(e_rho, 4), "prosst_spearman": round(p_rho, 4),
            "hybrid_spearman": round(h_rho, 4), "hybrid_minus_esm": round(h_rho - e_rho, 4),
            "prosst_vs_pg_repro": round(repro, 4) if repro is not None else None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=12)
    ap.add_argument("--vocab", type=int, default=2048)
    ap.add_argument("--structure-dir", default=None, help="dir of pre-quantized {UniProt}.json token lists")
    ap.add_argument("--only", default=None)
    ap.add_argument("--checkpoint", default=str(REPO / "data" / "processed" / "prosst_lift_checkpoint.jsonl"))
    args = ap.parse_args(argv)

    if not REF.exists():
        print("no ProteinGym reference on D:", file=sys.stderr); return 2
    ref = {r["DMS_id"]: r for r in csv.DictReader(REF.open(encoding="utf-8"))}
    categories = {d: r.get("coarse_selection_type", "?") for d, r in ref.items()}
    struct_dir = Path(args.structure_dir) if args.structure_dir else None

    picks = ([d.strip() for d in args.only.split(",")] if args.only else
             [d for d, r in ref.items()
              if (ZS_DIR / f"{d}.csv").exists() and (ESM_DIR / f"esm2_t33_650M_UR50D__{d}.json").exists()
              and r.get("UniProt_ID") and int(r.get("seq_len") or 0) <= 400][: args.limit])
    if not picks:
        print("no usable proteins (need UniProt_ID + ESM2 table + zeroshot)", file=sys.stderr); return 2

    cp = Path(args.checkpoint); cp.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if cp.exists():
        for line in cp.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(line); done[r["dms"]] = r
            except Exception:
                pass

    print(f"ProSST-{args.vocab} lift over {len(picks)} proteins "
          f"({'pre-quantized' if struct_dir else 'fetch+quantize'})")
    results = []
    for dms in picks:
        if dms in done and done[dms].get("status") == "OK":
            results.append(done[dms]); continue
        t0 = time.time()
        try:
            rec = score_protein(dms, ref[dms], struct_dir, args.vocab)
        except Exception as e:
            rec = {"dms": dms, "status": "ERROR", "error": str(e)[:200]}
        rec["seconds"] = round(time.time() - t0, 1)
        with cp.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        results.append(rec)
        if rec.get("status") == "OK":
            print(f"  {dms:40s} esm {rec['esm_spearman']} prosst {rec['prosst_spearman']} "
                  f"hyb {rec['hybrid_spearman']} d {rec['hybrid_minus_esm']:+} repro {rec['prosst_vs_pg_repro']}")
        else:
            print(f"  {dms:40s} {rec.get('status')} {rec.get('error','')}")

    ok = [r for r in results if r.get("status") == "OK"]
    if ok:
        print(f"\nN={len(ok)} | reproduction median (vs PG ProSST-{args.vocab}) = {reproduction_median(ok)}")
        print(f"overall {paired_in_subset(ok)}")
        for c, agg in by_category(ok, categories).items():
            print(f"  {c:20s} {agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
