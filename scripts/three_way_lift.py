"""The 3-way modality hybrid: does ESM2 (+) ProSST (+) GEMME beat the 2-way ESM2+ProSST? (sweep top 0.547)

The modality-hybrid sweep's TOP combination was `ESM2+GEMME+ProSST` (+0.056 vs ESM2-650M, win 90.5%). The
2-way `ESM2+ProSST` (structure) is now VALIDATED in our pipeline (`wiki/prosst_lift_2026-07-18.md`: +0.067,
93%). This checks whether ADDING the evolution modality (GEMME) on top lifts FURTHER — the decision "is the
3-way worth it over the 2-way".

Components, all real:
  - ESM2-650M: OUR precomputed masked-marginal tables (D:/dna_decode_cache/esm).
  - ProSST-2048: OUR forward pass (validated to reproduce ProteinGym's column at Spearman 1.0), scored from
    ProteinGym's pre-quantized structure tokens.
  - GEMME: ProteinGym's PRECOMPUTED column. GEMME has ZERO learned parameters (a deterministic
    evolutionary-conservation model), so its precomputed column IS canonical GEMME output — using it here is
    the SAME move as using the pre-quantized ProSST structures (a canonical deterministic input, not a
    learned score we could compute differently). The Windows-hostile GEMME install (JET2/R/Java) is needed
    ONLY to run GEMME on a NOVEL protein (deployment), NOT for this validation.

Reuses the FROZEN `rank_average_hybrid` (already N-ary). Per-category paired deltas; restartable checkpoint.

Run:  TORCH_HOME=D:/torch_hub HF_HOME=D:/hf_cache uv run python scripts/three_way_lift.py --limit 60 \
        --structure-dir D:/dna_decode_cache/prosst_structures/structure_sequence/2048
Exit: 0 = ran; 2 = substrate unavailable.
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
from scripts.prosst_lift import PG, ZS_DIR, REF, ESM_DIR, _load_esm, _esm_variant_table   # noqa: E402


def paired(records: list[dict], key: str) -> dict:
    ds = [r[key] for r in records if r.get("status") == "OK" and key in r]
    wins = sum(1 for x in ds if x > 1e-9)
    return {"n": len(ds), "median": round(statistics.median(ds), 4) if ds else None,
            "win": f"{wins}/{len(ds)}" if ds else "0/0"}


def by_category(records: list[dict], categories: dict[str, str], key: str) -> dict:
    buckets: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        if r.get("status") == "OK":
            buckets[categories.get(r["dms"], "?")].append(r)
    return {c: paired(rs, key) for c, rs in sorted(buckets.items(), key=lambda kv: -len(kv[1]))}


def score_protein(dms: str, row: dict, structure_dir: Path, vocab: int) -> dict:
    from dna_decode.forward.prosst_scorer import prosst_variant_table, _load_prosst
    from dna_decode.forward.variant_effect import rank_average_hybrid

    zs = ZS_DIR / f"{dms}.csv"
    esm = _load_esm(dms)
    if not (zs.exists() and esm):
        return {"dms": dms, "status": "SUBSTRATE_MISSING"}
    wt_seq = row["target_seq"].strip().upper()
    sf = structure_dir / f"{dms}.fasta"
    if not sf.exists():
        return {"dms": dms, "status": "NO_STRUCTURE"}
    line = next(l for l in sf.read_text(encoding="utf-8").splitlines() if not l.startswith(">"))
    tokens = [int(x) for x in line.strip().split(",")]

    muts, dms_v, gemme = [], [], {}
    with zs.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            m = r["mutant"]
            if len(m) < 3 or ":" in m or not m[1:-1].isdigit() or not r.get("DMS_score"):
                continue
            muts.append(m); dms_v.append(float(r["DMS_score"]))
            if r.get("GEMME") not in (None, "", "NA"):
                gemme[m] = float(r["GEMME"])

    prosst = prosst_variant_table(wt_seq, muts, structure_tokens=tokens, model_bundle=_load_prosst(vocab))
    esm_tab = _esm_variant_table(esm, muts)
    shared = [m for m in muts if m in prosst and m in esm_tab and m in gemme]
    if len(shared) < 20:
        return {"dms": dms, "status": "TOO_FEW", "n": len(shared)}

    idx = {m: i for i, m in enumerate(muts)}
    dd = [dms_v[idx[m]] for m in shared]
    two = rank_average_hybrid([{m: esm_tab[m] for m in shared}, {m: prosst[m] for m in shared}])
    three = rank_average_hybrid([{m: esm_tab[m] for m in shared},
                                 {m: prosst[m] for m in shared}, {m: gemme[m] for m in shared}])
    e_rho = abs(_spearman([esm_tab[m] for m in shared], dd))
    two_rho = abs(_spearman([two[m] for m in shared], dd))
    three_rho = abs(_spearman([three[m] for m in shared], dd))
    g_rho = abs(_spearman([gemme[m] for m in shared], dd))
    return {"dms": dms, "status": "OK", "n": len(shared),
            "esm": round(e_rho, 4), "gemme": round(g_rho, 4),
            "two_way": round(two_rho, 4), "three_way": round(three_rho, 4),
            "three_minus_two": round(three_rho - two_rho, 4),
            "three_minus_esm": round(three_rho - e_rho, 4)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=60)
    ap.add_argument("--vocab", type=int, default=2048)
    ap.add_argument("--structure-dir", required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--checkpoint", default=str(REPO / "data" / "processed" / "three_way_lift_checkpoint.jsonl"))
    args = ap.parse_args(argv)

    if not REF.exists():
        print("no ProteinGym reference on D:", file=sys.stderr); return 2
    ref = {r["DMS_id"]: r for r in csv.DictReader(REF.open(encoding="utf-8"))}
    categories = {d: r.get("coarse_selection_type", "?") for d, r in ref.items()}
    struct_dir = Path(args.structure_dir)

    picks = ([d.strip() for d in args.only.split(",")] if args.only else
             [d for d, r in ref.items()
              if (ZS_DIR / f"{d}.csv").exists() and (ESM_DIR / f"esm2_t33_650M_UR50D__{d}.json").exists()
              and (struct_dir / f"{d}.fasta").exists() and int(r.get("seq_len") or 0) <= 400][: args.limit])
    if not picks:
        print("no usable proteins", file=sys.stderr); return 2

    cp = Path(args.checkpoint); cp.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if cp.exists():
        for line in cp.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(line); done[r["dms"]] = r
            except Exception:
                pass

    print(f"3-way ESM2+ProSST+GEMME over {len(picks)} proteins")
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
            print(f"  {dms:40s} esm {rec['esm']} 2way {rec['two_way']} 3way {rec['three_way']} "
                  f"3-2 {rec['three_minus_two']:+} 3-esm {rec['three_minus_esm']:+}")
        else:
            print(f"  {dms:40s} {rec.get('status')} {rec.get('error','')}")

    ok = [r for r in results if r.get("status") == "OK"]
    if ok:
        print(f"\nN={len(ok)}")
        print(f"3-way vs 2-way (ESM2+ProSST): {paired(ok,'three_minus_two')}")
        print(f"3-way vs ESM2 baseline:       {paired(ok,'three_minus_esm')}")
        print(f"median |Spearman| 3-way {statistics.median(r['three_way'] for r in ok):.4f} "
              f"2-way {statistics.median(r['two_way'] for r in ok):.4f} "
              f"ESM2 {statistics.median(r['esm'] for r in ok):.4f}")
        print("per-category (3-way minus 2-way):")
        for c, agg in by_category(ok, categories, "three_minus_two").items():
            print(f"  {c:20s} {agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
