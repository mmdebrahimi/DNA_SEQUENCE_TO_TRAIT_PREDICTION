#!/usr/bin/env python
"""Leakage-free ESM2-650M + ProSST-2048 HYBRID on the held-out MaveDB set -- self-contained Kaggle T4.

The full at-scale version of the forward cell's modality-hybrid finding, on genes NOT in the ProteinGym
benchmark the hybrid was tuned on (leakage-free by construction; R2 has no population-structure confound).
For each held-out human assay: fetch seq + DMS from MaveDB, score ESM2-650M masked-marginals (GPU) AND
ProSST-2048 conditioned on the PRE-QUANTIZED structure tokens (embedded MANIFEST -- no torch_geometric on
Kaggle; the risky quantizer ran locally), rank-average the two orthogonal modalities, |Spearman| vs DMS.

Reports median |Spearman| for ESM2, ProSST, and the HYBRID separately -> does the hybrid beat single-modality
at scale? Comparators: ESM2-alone full holdout 0.478 (wiki/mavedb_full_esm2_2026-07-22); AM 0.502
(wiki/mavedb_am_holdout_2026-07-23). Writes /kaggle/working/mavedb_holdout_hybrid_result.json.

Kaggle: enable_gpu (PIN T4) + enable_internet=True. ProSST force-ties the omitted decoder weight (critical --
else garbage logits). Structure tokens are assay-region-sliced (tokens[offset:offset+L]) upstream.
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
ESM_MODEL = "facebook/esm2_t33_650M_UR50D"
PROSST_MODEL = "AI4Protein/ProSST-2048"
# PIN the ProSST remote-code revision. Run 2 (2026-07-23) had a VERIFIED decoder tie yet still scored ~0.05:
# `trust_remote_code` with revision=None tracks main, and the Kaggle log showed "A new version of ...
# modeling_prosst.py was downloaded" -- i.e. Kaggle pulled NEWER remote code than the local snapshot that
# scores MTHFR at Spearman 0.40. Pinning to the known-good local revision is the fix.
PROSST_REVISION = "e94ffee7846d7f55c1bf5efa8ec7372a336ac4b8"

# {urn: {gene, uniprot, offset, seq_len, n_tokens, structure_tokens:[...]}} -- injected at push time.
MANIFEST = json.loads(r'''__MANIFEST_JSON__''')

_AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H",
           "Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
           "Tyr":"Y","Val":"V"}
_BASES = "TCAG"
_CODON = dict(zip([a+b+c for a in _BASES for b in _BASES for c in _BASES],
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(dna):
    dna = dna.upper().replace("U", "T"); out = []
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON.get(dna[i:i+3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def parse_hgvs_pro(h):
    if not h or not h.startswith("p."):
        return None
    b = h[2:].strip()
    if b in ("=", "?") or "[" in b or "fs" in b or "del" in b or "ins" in b or "*" in b:
        return None
    if b[:3] not in _AA3TO1:
        return None
    i = 3
    while i < len(b) and b[i].isdigit():
        i += 1
    if i == 3 or i >= len(b):
        return None
    pos = int(b[3:i]); alt3 = b[i:]
    if alt3 in ("Ter", "=") or alt3 not in _AA3TO1:
        return None
    return _AA3TO1[b[:3]] + str(pos) + _AA3TO1[alt3]


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 5:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float); ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx*rx).sum() * (ry*ry).sum())
    return float((rx*ry).sum()/d) if d > 0 else float("nan")


def midrank(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i]); r = [0.0]*len(vals); i = 0
    while i < len(order):
        j = i
        while j+1 < len(order) and vals[order[j+1]] == vals[order[i]]:
            j += 1
        mid = (i+j)/2.0
        for k in range(i, j+1):
            r[order[k]] = mid
        i = j+1
    return r


def _get(url, timeout=90):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def fetch_assay(urn):
    rec = json.loads(_get(API + "/score-sets/" + urn))
    tgs = rec.get("targetGenes") or []
    if not tgs:
        return None
    ts = (tgs[0].get("targetSequence") or {})
    seq = (ts.get("sequence") or "").strip().upper()
    if not seq:
        return None
    if (ts.get("sequenceType") or "").lower() == "dna":
        seq = translate(seq)
    if not seq:
        return None
    text = _get(API + "/score-sets/" + urn + "/scores")
    dms = {}
    for row in csv.DictReader(io.StringIO(text)):
        sh = parse_hgvs_pro(row.get("hgvs_pro", "")); raw = row.get("score", "")
        if sh is None or raw in ("", "NA", None):
            continue
        try:
            dms[sh] = float(raw)
        except ValueError:
            continue
    return seq, dms


def esm2_deltas(seq, variants, tok, model, device, mask_id, aa_ids, torch):
    positions = sorted({p for _, p, _ in variants})
    logp_at = {}
    with torch.no_grad():
        enc = tok(seq, return_tensors="pt"); ids = enc["input_ids"].to(device)
        for i in range(0, len(positions), BATCH):
            chunk = positions[i:i+BATCH]; stack = ids.repeat(len(chunk), 1)
            for r, p in enumerate(chunk):
                stack[r, p] = mask_id
            logits = model(stack).logits
            for r, p in enumerate(chunk):
                logp_at[p] = torch.log_softmax(logits[r, p], dim=-1).float().cpu()
    return {(wt, pos, mut): float(logp_at[pos][aa_ids[mut]] - logp_at[pos][aa_ids[wt]])
            for wt, pos, mut in variants}


def prosst_deltas(seq, variants, structure_tokens, model, tokenizer, torch, device="cpu"):
    ss = torch.tensor([[1, *[t + 3 for t in structure_tokens], 2]], dtype=torch.long).to(device)
    enc = tokenizer([seq], return_tensors="pt")
    with torch.no_grad():
        logits = model(input_ids=enc["input_ids"].to(device),
                       attention_mask=enc["attention_mask"].to(device), ss_input_ids=ss).logits
    lp = torch.log_softmax(logits[0, 1:-1].float(), dim=-1)
    vocab = tokenizer.get_vocab()
    return {(wt, pos, mut): float(lp[pos-1, vocab[mut]] - lp[pos-1, vocab[wt]]) for wt, pos, mut in variants}


def main():
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print("=== ESM2+ProSST hybrid | device=" + device + " | " + str(len(MANIFEST)) + " held-out assays ===")

    etok = AutoTokenizer.from_pretrained(ESM_MODEL)
    emodel = AutoModelForMaskedLM.from_pretrained(ESM_MODEL).to(device).eval()
    mask_id = etok.mask_token_id
    aa_ids = {aa: etok.convert_tokens_to_ids(aa) for aa in AA}

    pmodel = AutoModelForMaskedLM.from_pretrained(PROSST_MODEL, trust_remote_code=True,
                                                  revision=PROSST_REVISION)
    ptok = AutoTokenizer.from_pretrained(PROSST_MODEL, trust_remote_code=True, revision=PROSST_REVISION)
    print("PROSST_REVISION pinned:", PROSST_REVISION)
    # CRITICAL: the ProSST checkpoint OMITS cls.predictions.decoder.weight -> it loads RANDOM -> garbage
    # logits (median |rho| ~0.04, the 2026-07-23 run-1 failure). Plain attribute assignment did NOT stick
    # under Kaggle's newer "Materializing param" lazy-load path, so ALSO copy in-place, then VERIFY loudly.
    emb = pmodel.get_input_embeddings().weight
    dec = None
    for path in ("cls.predictions.decoder", "lm_head.decoder", "cls.predictions"):
        obj = pmodel
        try:
            for part in path.split("."):
                obj = getattr(obj, part)
            if hasattr(obj, "weight"):
                dec = obj
                break
        except AttributeError:
            continue
    if dec is None:
        print("PROSST_TIE: FAILED - no decoder module found; ProSST scores will be GARBAGE")
    else:
        try:
            dec.weight = emb                                  # reference tie
        except Exception as e:
            print("PROSST_TIE: assignment failed", type(e).__name__)
        try:
            with torch.no_grad():                             # in-place copy (works even if assignment didn't stick)
                dec.weight.data.copy_(emb.data)
        except Exception as e:
            print("PROSST_TIE: copy failed", type(e).__name__)
        same = bool((dec.weight.data - emb.data).abs().max().item() == 0.0)
        print("PROSST_TIE: decoder==embeddings ->", same)
    pmodel = pmodel.to(device).eval()

    results, skipped = [], {"fetch": 0, "align": 0, "too_few": 0}
    for urn, meta in MANIFEST.items():
        try:
            got = fetch_assay(urn)
        except Exception as e:
            skipped["fetch"] += 1; print("  " + urn + " SKIP fetch " + str(e)[:50]); continue
        if got is None:
            skipped["fetch"] += 1; print("  " + urn + " SKIP fetch"); continue
        seq, dms = got
        toks = meta["structure_tokens"]
        if len(seq) != len(toks) or len(seq) > MAXLEN:
            skipped["align"] += 1
            print("  " + urn + " SKIP align seq=" + str(len(seq)) + " tok=" + str(len(toks))); continue
        variants, ys = [], []
        for m, y in dms.items():
            wt, pos, mut = m[0], int(m[1:-1]), m[-1]
            if 1 <= pos <= len(seq) and wt in AA and mut in AA and seq[pos-1] == wt:
                variants.append((wt, pos, mut)); ys.append(y)
        if len(variants) < 20:
            skipped["too_few"] += 1; print("  " + urn + " SKIP too_few " + str(len(variants))); continue
        ed = esm2_deltas(seq, variants, etok, emodel, device, mask_id, aa_ids, torch)
        pd = prosst_deltas(seq, variants, toks, pmodel, ptok, torch, device)
        e_scores = [ed[v] for v in variants]; p_scores = [pd[v] for v in variants]
        er = midrank(e_scores); pr = midrank(p_scores)
        hyb = [er[i] + pr[i] for i in range(len(variants))]   # rank-average (both higher=preserved)
        rho_e = abs(spearman(e_scores, ys)); rho_p = abs(spearman(p_scores, ys)); rho_h = abs(spearman(hyb, ys))
        r = {"urn": urn, "gene": meta.get("gene"), "n": len(variants),
             "esm2": round(rho_e, 4), "prosst": round(rho_p, 4), "hybrid": round(rho_h, 4)}
        results.append(r)
        print("  " + urn + " " + str(meta.get("gene"))[:12] + " n=" + str(len(variants))
              + "  ESM2=" + ("%.3f" % rho_e) + " ProSST=" + ("%.3f" % rho_p) + " HYBRID=" + ("%.3f" % rho_h))

    def med(key):
        v = [r[key] for r in results if r[key] == r[key]]
        return float(np.median(v)) if v else float("nan")

    # DEGENERATE-ProSST GUARD: run-1 (2026-07-23) silently reported ProSST median 0.04 because the omitted
    # decoder weight loaded RANDOM. ProSST is the ProteinGym structure-tier leader (~0.50 zero-shot), so a
    # near-zero median means the model is BROKEN, not the biology. Flag it in the artifact, loudly.
    prosst_degenerate = bool(results) and med("prosst") < 0.10
    # paired: hybrid vs best-single, per assay
    wins = sum(1 for r in results if r["hybrid"] >= max(r["esm2"], r["prosst"]))
    me, mp, mh = med("esm2"), med("prosst"), med("hybrid")
    print("\n=== " + str(len(results)) + " held-out assays scored ===")
    print("median |Spearman|: ESM2=" + ("%.4f" % me) + "  ProSST=" + ("%.4f" % mp) + "  HYBRID=" + ("%.4f" % mh))
    print("hybrid >= best-single on " + str(wins) + "/" + str(len(results)) + " assays  | skipped " + str(skipped))
    print("comparators: ESM2-full-holdout 0.478 | AM 0.502")
    if prosst_degenerate:
        print("*** PROSST_DEGENERATE: median " + ("%.4f" % mp) + " < 0.10 -- the structure model is BROKEN "
              "(check PROSST_TIE above); the HYBRID number is INVALID, do NOT publish it. ***")
    out = {"_schema": "mavedb-holdout-hybrid-v1", "esm_model": ESM_MODEL, "prosst_model": PROSST_MODEL,
           "n_manifest": len(MANIFEST), "n_scored": len(results),
           "median_esm2": me, "median_prosst": mp, "median_hybrid": mh,
           "prosst_degenerate": prosst_degenerate,
           "hybrid_wins_vs_best_single": wins, "skipped": skipped,
           "comparators": {"esm2_full_holdout": 0.478, "am_holdout": 0.502}, "results": results}
    json.dump(out, open("/kaggle/working/mavedb_holdout_hybrid_result.json", "w"), indent=2)
    print("wrote /kaggle/working/mavedb_holdout_hybrid_result.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
