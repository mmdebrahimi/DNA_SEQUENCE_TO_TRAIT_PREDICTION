#!/usr/bin/env python
"""MaveDB PHARMACOGENE ESM2-650M (R2 molecular) scorer -- self-contained for a free Kaggle T4 (internet ON).

Scores the 86 pinned held-out MaveDB DMS assays (target gene NOT in the ProteinGym benchmark the frozen
`forward` hybrid was tuned+validated on; see wiki/mavedb_prospective_holdout_2026-07-21.md) with ESM2-650M
zero-shot masked-marginals, |Spearman| vs the measured functional score. Leakage-free by benchmark exclusion;
R2 has NO population-structure confound (designed mutant library). The masked-marginal scoring core is
byte-faithful to scripts/esm_zeroshot_dms.py / notebooks/j2_phase1_esm2_proteingym.py.

Kaggle: enable_gpu (T4) + enable_internet=True (fetches MaveDB directly; no dataset attach). Writes
/kaggle/working/mavedb_pgx_esm2_result.json. |Spearman| is direction-robust (MaveDB does not standardize
per-assay score direction -- the curation ProteinGym adds). R2 analog of the AMR prospective-lock.
"""
import csv, io, json, sys, urllib.request
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

API = "https://api.mavedb.org/api/v1"
AA = "ACDEFGHIKLMNPQRSTVWY"
MAXLEN = 1022
BATCH = 16
MODEL = "facebook/esm2_t33_650M_UR50D"

URNS = [
    'urn:mavedb:00000055-b-1',
    'urn:mavedb:00000055-a-1',
    'urn:mavedb:00000062-b-1',
    'urn:mavedb:00000062-a-1',
    'urn:mavedb:00000055-0-1',
    'urn:mavedb:00000078-a-1',
    'urn:mavedb:00000078-b-1',
    'urn:mavedb:00000095-b-1',
    'urn:mavedb:00000095-a-1',
    'urn:mavedb:00001199-a-2',
    'urn:mavedb:00001199-a-1',
    'urn:mavedb:00001237-a-2',
    'urn:mavedb:00001266-b-1',
    'urn:mavedb:00001266-a-1',
    'urn:mavedb:00001266-0-1',
    'urn:mavedb:00001266-c-1',
]

_AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
           "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
           "Tyr":"Y","Val":"V"}
_BASES = "TCAG"
_CODONS = [a+b+c for a in _BASES for b in _BASES for c in _BASES]
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = dict(zip(_CODONS, _AAS))


def translate(dna):
    dna = dna.upper().replace("U", "T")
    prot = []
    for i in range(0, len(dna) - 2, 3):
        aa = CODON.get(dna[i:i+3], "X")
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)


def parse_hgvs_pro(h):
    """p.Val82Ala -> 'V82A'; None for synonymous/nonsense/multi/indel/non-standard."""
    if not h or not h.startswith("p."):
        return None
    b = h[2:].strip()
    if b in ("=", "?") or "[" in b or "fs" in b or "del" in b or "ins" in b or "*" in b:
        return None
    wt3 = b[:3]
    if wt3 not in _AA3TO1:
        return None
    i = 3
    while i < len(b) and b[i].isdigit():
        i += 1
    if i == 3 or i >= len(b):
        return None
    pos = int(b[3:i]); alt3 = b[i:]
    if alt3 in ("Ter", "=") or alt3 not in _AA3TO1:
        return None
    return _AA3TO1[wt3] + str(pos) + _AA3TO1[alt3]


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx*rx).sum() * (ry*ry).sum())
    return float((rx*ry).sum()/d) if d > 0 else float("nan")


def _get(url, timeout=90):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def fetch_assay(urn):
    """-> ((protein_seq, {mutant: score}), None) or (None, reason). Translates DNA; single-missense only."""
    try:
        rec = json.loads(_get(API + "/score-sets/" + urn))
    except Exception as e:
        return None, "rec:" + str(e)
    tgs = rec.get("targetGenes") or []
    if not tgs:
        return None, "no-target"
    ts = (tgs[0].get("targetSequence") or {})
    seq = (ts.get("sequence") or "").strip().upper()
    if not seq:
        return None, "no-seq"
    if (ts.get("sequenceType") or "").lower() == "dna":
        seq = translate(seq)
    if not seq:
        return None, "empty-protein"
    try:
        text = _get(API + "/score-sets/" + urn + "/scores")
    except Exception as e:
        return None, "scores:" + str(e)
    dms = {}
    for row in csv.DictReader(io.StringIO(text)):
        sh = parse_hgvs_pro(row.get("hgvs_pro", ""))
        raw = row.get("score", "")
        if sh is None or raw in ("", "NA", None):
            continue
        try:
            dms[sh] = float(raw)
        except ValueError:
            continue
    return (seq, dms), None


def score_assay(seq, dms, tok, model, device, mask_id, aa_ids, torch):
    L = len(seq)
    if L > MAXLEN:
        return None, "too_long"
    variants, positions, mism = [], set(), 0
    for m, y in dms.items():
        wt, pos, mut = m[0], int(m[1:-1]), m[-1]
        if pos < 1 or pos > L or wt not in AA or mut not in AA:
            continue
        if seq[pos-1] != wt:
            mism += 1; continue
        variants.append((wt, pos, mut, y)); positions.add(pos)
    if len(variants) < 20:
        return None, "too_few"
    positions = sorted(positions); logp_at = {}
    with torch.no_grad():
        enc = tok(seq, return_tensors="pt"); ids = enc["input_ids"].to(device)
        for i in range(0, len(positions), BATCH):
            chunk = positions[i:i+BATCH]; stack = ids.repeat(len(chunk), 1)
            for r, p in enumerate(chunk):
                stack[r, p] = mask_id
            logits = model(stack).logits
            for r, p in enumerate(chunk):
                logp_at[p] = torch.log_softmax(logits[r, p], dim=-1).float().cpu()
    xs, ys = [], []
    for wt, pos, mut, y in variants:
        lp = logp_at[pos]; xs.append(float(lp[aa_ids[mut]] - lp[aa_ids[wt]])); ys.append(y)
    rho = spearman(xs, ys)
    rng = np.random.default_rng(0); yss = np.array(ys); rng.shuffle(yss)
    return {"n": len(xs), "rho": rho, "rho_shuf": spearman(xs, yss), "mism": mism}, None


def main():
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("=== MaveDB PHARMACOGENE ESM2-650M (R2 molecular) | device=" + device + " | " + str(len(URNS)) + " pinned assays ===")
    if device == "cpu":
        print("WARNING: no CUDA -- 650M on CPU is slow.")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForMaskedLM.from_pretrained(MODEL).to(device).eval()
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}

    results, skipped = [], {"fetch": 0, "too_long": 0, "too_few": 0}
    for urn in URNS:
        got, why = fetch_assay(urn)
        if got is None:
            skipped["fetch"] += 1; print("  " + urn + " SKIP fetch (" + str(why) + ")"); continue
        seq, dms = got
        r, why2 = score_assay(seq, dms, tok, model, device, mask_id, aa_ids, torch)
        if why2:
            skipped[why2] = skipped.get(why2, 0) + 1
            print("  " + urn + " SKIP " + why2 + " (L=" + str(len(seq)) + ", n=" + str(len(dms)) + ")"); continue
        r["urn"] = urn; results.append(r)
        print("  " + urn + " n=" + str(r["n"]) + "  |rho|=" + ("%.3f" % abs(r["rho"])) + "  shuf=" + ("%.3f" % abs(r["rho_shuf"])))

    absr = [abs(r["rho"]) for r in results if r["rho"] == r["rho"]]
    med = float(np.median(absr)) if absr else float("nan")
    shuf = float(np.median([abs(r["rho_shuf"]) for r in results])) if results else float("nan")
    print("\nESM2-650M on the MaveDB holdout: " + str(len(results)) + " scored | median |Spearman| = "
          + ("%.3f" % med) + " (shuffled " + ("%.3f" % shuf) + ") | skipped " + str(skipped))
    out = {"_schema": "mavedb-pgx-esm2-v1", "model": MODEL, "n_pinned": len(URNS), "n_scored": len(results),
           "median_abs_spearman": med, "median_abs_shuffled": shuf, "skipped": skipped, "results": results}
    json.dump(out, open("/kaggle/working/mavedb_pgx_esm2_result.json", "w"), indent=2)
    print("wrote /kaggle/working/mavedb_pgx_esm2_result.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
