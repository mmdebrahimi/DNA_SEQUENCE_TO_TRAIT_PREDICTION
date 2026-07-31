# Kaggle F2 stratified reconstruction runner. Reproducible form: STRATA_JSON is LOADED (not embedded).
# Build the strata first with scripts/build_f2_strata.py -> f2_strata.json, then either attach it as a
# Kaggle dataset at /kaggle/input/<ds>/f2_strata.json or set F2_STRATA_PATH. The committed run embedded
# the JSON inline for delivery; this loadable form is the canonical reproducible artifact.
import os
_p = os.environ.get("F2_STRATA_PATH", "/kaggle/input/dog-f2-strata/f2_strata.json")
with open(_p) as _f:
    STRATA_JSON = _f.read()

# F2 STRUCTURED-vs-CALIBRATION test: is NT-500M's masked-reconstruction win over Markov CONCENTRATED
# in CODING regions (real structural signal) or FLAT across coding/intergenic (mere calibration)?
# Per-base marginal NLL, disjoint-trained Markov, tokens bucketed by canFam4 CDS annotation.
# Self-contained; Kaggle T4 (free). Writes /kaggle/working/dog_nt_f2_strata_result.json.
import os, json, math, subprocess, sys
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "transformers==4.44.2"], check=False)
import torch, transformers
from transformers import AutoModelForMaskedLM, AutoTokenizer
print("transformers", transformers.__version__, "torch", torch.__version__, flush=True)
S = json.loads(STRATA_JSON)
KMER = S["kmer"]; WINDOWS = S["windows"]; TRAIN = S["train_seq"]
DEV = "cuda" if torch.cuda.is_available() else "cpu"
MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
KS = list(range(1, 9))
_B = ("A", "C", "G", "T"); _BI = {b: i for i, b in enumerate(_B)}

# ---- Markov (disjoint-trained; base_distribution + per-base NLL) ----
def markov_fit(seq, k):
    tables = [dict() for _ in range(k + 1)]
    for i, b in enumerate(seq):
        if b not in _BI: continue
        for o in range(0, k + 1):
            if i < o: continue
            ctx = seq[i - o:i]
            if all(c in _BI for c in ctx): tables[o].setdefault(ctx, {x: 0 for x in _B})[b] += 1
    return tables
def base_dist(tables, k, left, alpha=1.0):
    counts = {b: 0 for b in _B}
    for o in range(min(k, len(left)), -1, -1):
        ctx = left[-o:] if o > 0 else ""
        if ctx in tables[o]: counts = dict(tables[o][ctx]); break
    sm = {b: counts[b] + alpha for b in _B}; tot = sum(sm.values())
    return {b: sm[b] / tot for b in _B}
# fit Markov at every k once on the disjoint TRAIN region
MK = {k: markov_fit(TRAIN, k) for k in KS}

# ---- NT-500M load (v2 remote code needs fp32 forward) ----
tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained(MODEL_ID, trust_remote_code=True,
                                             torch_dtype=torch.float32).to(DEV).eval()
print("loaded", MODEL_ID, flush=True)
vocab = tok.get_vocab(); special = set(tok.all_special_ids); mask_id = tok.mask_token_id
bmap = None
# accumulators per stratum: nt_nll_sum, markov_nll_sum[k], n_bases
acc = {s: {"nt": 0.0, "mk": {k: 0.0 for k in KS}, "n": 0, "nt_correct": 0} for s in ("coding", "intergenic")}

for wi, w in enumerate(WINDOWS):
    seq = w["seq"]; strata = w["token_strata"]
    enc = tok(seq, return_tensors="pt", truncation=True); ids = enc["input_ids"][0]
    strs = tok.convert_ids_to_tokens(ids.tolist())
    kmers = []; cur = 0
    for ai, tid in enumerate(ids.tolist()):
        if tid in special: continue
        kmers.append((ai, strs[ai], cur)); cur += len(strs[ai])
    # score only full-ACGT 6-mer tokens whose stratum is coding/intergenic
    targets = [t for t in range(len(kmers))
               if t < len(strata) and strata[t] in ("coding", "intergenic")
               and len(kmers[t][1]) == KMER and all(c in _BI for c in kmers[t][1])]
    for s0 in range(0, len(targets), 8):
        chunk = targets[s0:s0 + 8]
        batch = ids.unsqueeze(0).repeat(len(chunk), 1).clone()
        for r, tp in enumerate(chunk): batch[r, kmers[tp][0]] = mask_id
        with torch.no_grad():
            logits = model(input_ids=batch.to(DEV)).logits
        probs = torch.softmax(logits.float(), dim=-1).cpu()
        if bmap is None:
            V = probs.shape[-1]; bmap = []
            for j in range(KMER):
                m = torch.full((V,), -1, dtype=torch.long)
                for t, i in vocab.items():
                    if i < V and len(t) == KMER and t[j] in _BI: m[i] = _BI[t[j]]
                bmap.append(m)
        for r, tp in enumerate(chunk):
            ai, kmer, bstart = kmers[tp]; strat = strata[tp]; dist = probs[r, ai]
            for j in range(KMER):
                tb = kmer[j]
                if tb not in _BI: continue
                a4 = torch.zeros(4); mm = bmap[j]; keep = mm >= 0
                a4.index_add_(0, mm[keep], dist[keep])
                acc[strat]["nt"] += -math.log(max(float(a4[_BI[tb]]), 1e-12))
                acc[strat]["nt_correct"] += (int(torch.argmax(a4)) == _BI[tb])
                jpos = bstart + j; left = seq[max(0, jpos - 8):jpos]
                for k in KS:
                    d = base_dist(MK[k], k, left[-k:] if k else "")
                    acc[strat]["mk"][k] += -math.log(max(d[tb], 1e-12))
                acc[strat]["n"] += 1
    if wi % 5 == 0: print(f"  window {wi+1}/{len(WINDOWS)}", flush=True)

out = {"model": MODEL_ID, "gpu": torch.cuda.get_device_name(0) if DEV == "cuda" else "cpu",
       "transformers": transformers.__version__, "strata": {}}
for s in ("coding", "intergenic"):
    n = acc[s]["n"]
    if n == 0: out["strata"][s] = {"n_bases": 0}; continue
    nt_nll = acc[s]["nt"] / n
    mk_by_k = {k: acc[s]["mk"][k] / n for k in KS}
    bk = min(KS, key=lambda k: mk_by_k[k])           # hardest (lowest-NLL) Markov
    out["strata"][s] = {"n_bases": n, "nt_per_base_nll": round(nt_nll, 4),
                        "markov_per_base_nll_bestk": round(mk_by_k[bk], 4), "markov_best_k": bk,
                        "nll_delta_markov_minus_nt": round(mk_by_k[bk] - nt_nll, 4),
                        "nt_marginal_accuracy": round(acc[s]["nt_correct"] / n, 4)}
if out["strata"].get("coding", {}).get("n_bases") and out["strata"].get("intergenic", {}).get("n_bases"):
    out["structured_signal_delta_gap"] = round(
        out["strata"]["coding"]["nll_delta_markov_minus_nt"]
        - out["strata"]["intergenic"]["nll_delta_markov_minus_nt"], 4)
with open("/kaggle/working/dog_nt_f2_strata_result.json", "w") as f: json.dump(out, f, indent=2)
print("DONE", json.dumps(out, indent=2))
