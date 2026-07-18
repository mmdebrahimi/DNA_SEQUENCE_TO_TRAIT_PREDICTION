"""Does ESM2 (+) MY-OWN MSA-Transformer LIFT on ProteinGym? -- closing the lifting hybrid end-to-end.

The modality-hybrid finding used ProteinGym's PRECOMPUTED `MSA_Transformer` column (+0.013 paired vs
ESM2-650M). This runs MSA-Transformer OURSELVES (esm_msa1b_t12_100M_UR50S, CPU) on the on-D: MSAs, feeds its
scores through OUR `rank_average_hybrid`, and checks the lift reproduces -- converting "ProteinGym's column
lifts" into "OUR deployable pipeline lifts", the last validation the run-time evolution component needed.

Per protein (with an on-D: MSA + an ESM2-650M table + a zeroshot CSV): subsample the MSA to N rows, run one
MSA-T forward, read the query row's per-position log-probs (wt-marginal -- the fast single-forward zero-shot;
ProteinGym uses the slower masked-marginal, so exact reproduction is not expected, the LIFT is the target),
score every DMS variant at a match column, and compare abs-Spearman of: ESM2 alone / MSA-T alone / PG's
MSA_Transformer column / the ESM2(+)MSA-T rank-average hybrid. Paired median delta (hybrid - ESM2) is the
verdict. Restartable JSONL checkpoint (one line per protein) -- unattended/rate-limit safe.

Run:  TORCH_HOME=D:/torch_hub HF_HOME=D:/hf_cache uv run python scripts/msa_transformer_lift.py --limit 8
Exit: 0 = ran; 2 = substrate/model unavailable.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward.msa_evolution import parse_a2m, query_pos_to_col  # noqa: E402
from dna_decode.forward.variant_effect import rank_average_hybrid          # noqa: E402
from scripts.forward_blosum_proteingym_sweep import _spearman              # noqa: E402

PG = Path("D:/dna_decode_cache/proteingym")
MSA_DIR = PG / "pg_msa" / "DMS_msa_files"
MSA_ZIP = PG / "pg_msa_files.zip"
ZS_DIR = PG / "pg_zeroshot"
REF = PG / "pg_reference.csv"
ESM_DIR = Path("D:/dna_decode_cache/esm")
MAX_SEQ_LEN = 400        # skip very long proteins (MSA-T CPU memory/time)


def _ensure_msa(filename: str) -> bool:
    """Extract the MSA from pg_msa_files.zip into MSA_DIR if not already extracted. Returns True if present."""
    target = MSA_DIR / filename
    if target.exists():
        return True
    if not MSA_ZIP.exists():
        return False
    import zipfile
    try:
        z = zipfile.ZipFile(MSA_ZIP)
        member = next((n for n in z.namelist() if n.endswith(filename)), None)
        if member is None:
            return False
        MSA_DIR.mkdir(parents=True, exist_ok=True)
        target.write_bytes(z.read(member))
        return True
    except Exception:
        return False


def _load_esm(dms: str):
    q = ESM_DIR / f"esm2_t33_650M_UR50D__{dms}.json"
    if not q.exists():
        return None
    return {int(k): v for k, v in json.loads(q.read_text(encoding="utf-8")).items()}


def _subsample(match_cols: list[str], n: int) -> list[str]:
    """Keep the focus (row 0) + an evenly-strided, deterministic subsample of the rest to depth n."""
    if len(match_cols) <= n:
        return match_cols
    rest = len(match_cols) - 1
    stride = max(1, rest // (n - 1))
    idx = [0] + list(range(1, len(match_cols), stride))[: n - 1]
    return [match_cols[i] for i in idx]


def _focus_residues(focus_raw: str) -> dict[int, str]:
    qpos, out = 0, {}
    for c in focus_raw:
        if c.isupper() or c.islower():
            qpos += 1
            out[qpos] = c.upper()
    return out


def msat_query_logprobs(msa_path: Path, model, alphabet, bc, n_seqs: int):
    """Return (focus_raw, lp) where lp[col+1] = query-row log-probs at match column col (BOS prepended)."""
    import torch
    name, focus_raw, match_cols = parse_a2m(msa_path, max_rows=max(4000, n_seqs * 10))  # bound memory (OOM fix)
    sub = _subsample(match_cols, n_seqs)
    data = [(f"s{i}", s) for i, s in enumerate(sub)]
    _, _, toks = bc(data)
    with torch.no_grad():
        logits = model(toks, repr_layers=[], return_contacts=False)["logits"]
    lp = torch.log_softmax(logits[0, 0].float(), dim=-1)   # query row 0 -> (L+1, vocab)
    return focus_raw, lp, len(sub)


def score_protein(dms: str, row: dict, model, alphabet, bc, n_seqs: int) -> dict | None:
    msa = MSA_DIR / row["MSA_filename"]
    zs = ZS_DIR / f"{dms}.csv"
    esm = _load_esm(dms)
    if not (_ensure_msa(row["MSA_filename"]) and zs.exists() and esm):
        return {"dms": dms, "status": "SUBSTRATE_MISSING"}

    focus_raw, lp, depth = msat_query_logprobs(msa, model, alphabet, bc, n_seqs)
    q2c = query_pos_to_col(focus_raw)
    qres = _focus_residues(focus_raw)
    idx = alphabet.get_idx

    msat_tab, esm_tab, muts, dms_v, pg_msat = {}, {}, [], [], {}
    with zs.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            m = r["mutant"]
            if len(m) < 3 or not m[1:-1].isdigit():
                continue
            wt, p, alt = m[0], int(m[1:-1]), m[-1]
            if p not in q2c or qres.get(p) != wt or not r.get("DMS_score"):
                continue
            col = q2c[p]
            try:
                dscore = float(r["DMS_score"])
                emap = esm.get(p)
                if not emap or wt not in emap or alt not in emap:
                    continue
                msat_tab[m] = lp[col + 1, idx(alt)].item() - lp[col + 1, idx(wt)].item()
                esm_tab[m] = emap[alt] - emap[wt]
                muts.append(m); dms_v.append(dscore)
                if r.get("MSA_Transformer_ensemble") not in (None, "", "NA"):
                    pg_msat[m] = float(r["MSA_Transformer_ensemble"])
            except (ValueError, KeyError):
                continue

    if len(muts) < 20:
        return {"dms": dms, "status": "TOO_FEW", "n": len(muts)}

    hyb = rank_average_hybrid([esm_tab, msat_tab])
    didx = {m: i for i, m in enumerate(muts)}
    e_rho = abs(_spearman([esm_tab[m] for m in muts], dms_v))
    m_rho = abs(_spearman([msat_tab[m] for m in muts], dms_v))
    h_rho = abs(_spearman([hyb[m] for m in muts], dms_v))
    repro = (abs(_spearman([msat_tab[m] for m in muts if m in pg_msat],
                           [pg_msat[m] for m in muts if m in pg_msat]))
             if len(pg_msat) >= 20 else None)
    return {"dms": dms, "status": "OK", "n": len(muts), "msa_depth": depth,
            "esm_spearman": round(e_rho, 4), "my_msat_spearman": round(m_rho, 4),
            "hybrid_spearman": round(h_rho, 4), "hybrid_minus_esm": round(h_rho - e_rho, 4),
            "my_msat_vs_pg_msat_repro": round(repro, 4) if repro is not None else None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=8)
    ap.add_argument("--n-seqs", type=int, default=128)
    ap.add_argument("--only", default=None, help="comma-separated DMS ids (else auto-pick from usable set)")
    ap.add_argument("--checkpoint", default=str(REPO / "data" / "processed" / "msat_lift_checkpoint.jsonl"))
    args = ap.parse_args(argv)

    if not REF.exists():
        print("no ProteinGym reference on D:", file=sys.stderr)
        return 2

    import os
    ref = {}
    with REF.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            ref[r["DMS_id"]] = r
    extracted = set(os.listdir(MSA_DIR)) if MSA_DIR.exists() else set()

    if args.only:
        picks = [d.strip() for d in args.only.split(",")]
    else:
        # any protein with an ESM2 table + zeroshot CSV + an MSA in the zip (extracted on demand), <= MAX_SEQ_LEN
        picks = [d for d, r in ref.items()
                 if (ZS_DIR / f"{d}.csv").exists()
                 and (ESM_DIR / f"esm2_t33_650M_UR50D__{d}.json").exists()
                 and int(r.get("seq_len") or 0) <= MAX_SEQ_LEN][: args.limit]
    if not picks:
        print("no usable proteins (need on-D: MSA + zeroshot + ESM2 table)", file=sys.stderr)
        return 2

    cp = Path(args.checkpoint)
    cp.parent.mkdir(parents=True, exist_ok=True)
    done = {}
    if cp.exists():
        for line in cp.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line); done[rec["dms"]] = rec
            except Exception:
                pass

    print(f"loading MSA-Transformer (esm_msa1b_t12_100M_UR50S)... [{len(picks)} proteins, N={args.n_seqs}]")
    import esm
    model, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    model.eval()
    bc = alphabet.get_batch_converter()

    results = []
    for dms in picks:
        if dms in done and done[dms].get("status") == "OK":
            results.append(done[dms]); print(f"  [cached] {dms}"); continue
        t0 = time.time()
        try:
            rec = score_protein(dms, ref[dms], model, alphabet, bc, args.n_seqs)
        except Exception as e:
            rec = {"dms": dms, "status": "ERROR", "error": str(e)[:200]}
        rec["seconds"] = round(time.time() - t0, 1)
        with cp.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
        results.append(rec)
        s = rec.get("status")
        if s == "OK":
            print(f"  {dms:42s} esm {rec['esm_spearman']}  my-msat {rec['my_msat_spearman']}  "
                  f"hybrid {rec['hybrid_spearman']}  Δ {rec['hybrid_minus_esm']:+}  "
                  f"repro {rec['my_msat_vs_pg_msat_repro']}  ({rec['seconds']}s)")
        else:
            print(f"  {dms:42s} {s} {rec.get('error','')}")

    ok = [r for r in results if r.get("status") == "OK"]
    if ok:
        deltas = [r["hybrid_minus_esm"] for r in ok]
        wins = sum(1 for d in deltas if d > 1e-9)
        print(f"\n== PAIRED (ESM2+my-MSA-T vs ESM2), n={len(ok)} ==")
        print(f"median Δ {statistics.median(deltas):+.4f}  win-rate {wins}/{len(ok)}  "
              f"median hybrid {statistics.median(r['hybrid_spearman'] for r in ok):.4f} "
              f"vs ESM2 {statistics.median(r['esm_spearman'] for r in ok):.4f}")
        rep = [r['my_msat_vs_pg_msat_repro'] for r in ok if r.get('my_msat_vs_pg_msat_repro') is not None]
        if rep:
            print(f"my-MSA-T reproduces PG MSA_Transformer: median Spearman {statistics.median(rep):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
