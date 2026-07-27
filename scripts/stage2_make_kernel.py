"""Generate the self-contained Kaggle Enformer kernel (bundle embedded as base64)."""
import base64, json
from pathlib import Path

B = Path("D:/dna_decode_cache/stage2_kaggle")
bundle_b64 = base64.b64encode((B / "stage2_bundle.npz").read_bytes()).decode()
meta = json.loads((B / "stage2_meta.json").read_text())
out = Path("C:/Users/Farshad/stage2_kernel"); out.mkdir(exist_ok=True)

KERNEL = '''# GEUVADIS Stage-2 DNA-encoder arm: Enformer zero-shot cis-variant effects vs the linear ceiling.
import base64, io, json, subprocess, sys
subprocess.run([sys.executable,"-m","pip","install","-q","enformer-pytorch"], check=True)
import numpy as np, torch, urllib.request
from enformer_pytorch import Enformer
from scipy.stats import spearmanr

BUNDLE_B64 = """%s"""
META = %s
open("/kaggle/working/b.npz","wb").write(base64.b64decode(BUNDLE_B64))
d = np.load("/kaggle/working/b.npz")
genes = [g["gene"] for g in META["genes"]]
pops = np.array(META["pops"]); samples = META["samples"]
print("genes", len(genes), "individuals", len(samples), "pops", {p:int((pops==p).sum()) for p in set(pops)})

# --- select expression tracks: CAGE lymphoblastoid/GM12878/blood ---
tsel = []
try:
    url="https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt"
    txt=urllib.request.urlopen(url, timeout=60).read().decode().splitlines()
    hdr=txt[0].split("\t"); di=hdr.index("description")
    for i,row in enumerate(txt[1:]):
        desc=row.split("\t")[di].lower()
        if "cage" in desc and any(k in desc for k in ["gm12878","lymphoblast","b cell","b-cell","blood"]):
            tsel.append((i,desc))
    print("selected CAGE tracks:", tsel[:8], "...n=",len(tsel))
except Exception as e:
    print("targets fetch failed:", e)
if not tsel:
    track_idx = list(range(4675,5313)); print("FALLBACK: all CAGE tracks", len(track_idx))
else:
    track_idx = [i for i,_ in tsel]

model = Enformer.from_pretrained("EleutherAI/enformer-official-rough").eval().cuda()
BASE={65:0,67:1,71:2,84:3,78:4,97:0,99:1,103:2,116:3,110:4}  # ACGTN upper+lower
TSS_BIN=448
def idx_of(seqbytes):
    return torch.tensor([BASE.get(int(b),4) for b in seqbytes], dtype=torch.long)

def enf_tss(t_idx_tensor):
    with torch.no_grad():
        out = model(t_idx_tensor[None].cuda())
        h = out["human"][0] if isinstance(out, dict) else out[0]
    return h[TSS_BIN, track_idx].mean().item()

def spearman_pooled_within(y, pred, grp):
    pooled = spearmanr(y, pred).statistic
    ws=[]; ns=[]
    for g in set(grp.tolist()):
        m = grp==g
        if m.sum()>=5 and np.std(pred[m])>0 and np.std(y[m])>0:
            ws.append(spearmanr(y[m], pred[m]).statistic); ns.append(int(m.sum()))
    within = float(np.average(np.abs(ws), weights=ns)) if ws else float("nan")
    return float(abs(pooled)), within

rows=[]
for gi, gene in enumerate(genes):
    seq = d[f"{gene}__seq"]; off = d[f"{gene}__off"]; alt = d[f"{gene}__alt"]
    dos = d[f"{gene}__dos"].astype(np.float64); expr = d[f"{gene}__expr"].astype(np.float64)
    base_idx = idx_of(seq)
    ref = enf_tss(base_idx)
    eff=[]
    for o,a in zip(off, alt):
        ai = base_idx.clone(); ai[int(o)] = BASE.get(int(a),4)
        eff.append(enf_tss(ai) - ref)
    eff = np.array(eff)
    pred = eff @ dos
    if np.std(pred)==0:
        rows.append({"gene":gene,"note":"flat_pred"}); continue
    p,w = spearman_pooled_within(expr, pred, pops)
    rows.append({"gene":gene,"n_var":int(len(off)),"pooled":p,"within":w})
    print(f"[{gi+1}/{len(genes)}] {gene}: pooled={p:.3f} within={w:.3f}")

sc=[r for r in rows if "pooled" in r]
summary={"n_genes_scored":len(sc),
         "mean_abs_pooled_enformer":float(np.mean([r["pooled"] for r in sc])) if sc else None,
         "mean_abs_within_enformer":float(np.nanmean([r["within"] for r in sc])) if sc else None,
         "linear_ceiling_within":0.19, "n_tracks":len(track_idx)}
json.dump({"summary":summary,"per_gene":rows}, open("/kaggle/working/stage2_results.json","w"), indent=2)
print("STAGE-2 ENFORMER SUMMARY:", json.dumps(summary, indent=2))
'''

(out / "stage2_enformer.py").write_text(KERNEL % (bundle_b64, json.dumps(meta)), encoding="utf-8")
km = {"id":"emanueleebrahimi/geuvadis-stage2-enformer","title":"geuvadis-stage2-enformer",
      "code_file":"stage2_enformer.py","language":"python","kernel_type":"script",
      "is_private":True,"enable_gpu":True,"enable_internet":True,
      "dataset_sources":[],"competition_sources":[],"kernel_sources":[]}
(out / "kernel-metadata.json").write_text(json.dumps(km, indent=2), encoding="utf-8")
sz = (out / "stage2_enformer.py").stat().st_size
print(f"kernel written: {out}/stage2_enformer.py ({sz//1024} KB)")
