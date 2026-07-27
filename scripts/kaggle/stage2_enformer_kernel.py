# GEUVADIS Stage-2 DNA-encoder arm: Enformer zero-shot cis-variant effects vs the linear ceiling.
import json, subprocess, sys, urllib.request
subprocess.run([sys.executable,"-m","pip","install","-q","enformer-pytorch"], check=True)
import numpy as np, torch
from enformer_pytorch import Enformer
from scipy.stats import spearmanr

INP="/kaggle/input/geuvadis-stage2-enformer-bundle"
d=np.load(INP+"/stage2_bundle.npz")
META=json.load(open(INP+"/stage2_meta.json"))
genes=[g["gene"] for g in META["genes"]]; pops=np.array(META["pops"])
print("genes",len(genes),"individuals",len(META["samples"]))

TAB=chr(9); tsel=[]
try:
    url="https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt"
    txt=urllib.request.urlopen(url,timeout=60).read().decode().splitlines()
    hdr=txt[0].split(TAB); di=hdr.index("description")
    for i,row in enumerate(txt[1:]):
        desc=row.split(TAB)[di].lower()
        if "cage" in desc and any(k in desc for k in ["gm12878","lymphoblast","b cell","b-cell","blood"]):
            tsel.append((i,desc))
    print("selected CAGE tracks n=",len(tsel),tsel[:6])
except Exception as e:
    print("targets fetch failed:",e)
track_idx=[i for i,_ in tsel] if tsel else list(range(4675,5313))
print("n expression tracks used:",len(track_idx))

model=Enformer.from_pretrained("EleutherAI/enformer-official-rough").eval().cuda()
BASE={65:0,67:1,71:2,84:3,78:4,97:0,99:1,103:2,116:3,110:4}; TSS_BIN=448
def idx_of(sb): return torch.tensor([BASE.get(int(b),4) for b in sb],dtype=torch.long)
def enf_tss(t):
    with torch.no_grad():
        o=model(t[None].cuda()); h=o["human"][0] if isinstance(o,dict) else o[0]
    return h[TSS_BIN,track_idx].mean().item()
def pw(y,pred,grp):
    pooled=spearmanr(y,pred).statistic; ws=[]; ns=[]
    for g in set(grp.tolist()):
        m=grp==g
        if m.sum()>=5 and np.std(pred[m])>0 and np.std(y[m])>0:
            ws.append(abs(spearmanr(y[m],pred[m]).statistic)); ns.append(int(m.sum()))
    return float(abs(pooled)),(float(np.average(ws,weights=ns)) if ws else float("nan"))

rows=[]
for gi,gene in enumerate(genes):
    seq=d[gene+"__seq"]; off=d[gene+"__off"]; alt=d[gene+"__alt"]
    dos=d[gene+"__dos"].astype(np.float64); expr=d[gene+"__expr"].astype(np.float64)
    bi=idx_of(seq); ref=enf_tss(bi); eff=[]
    for o,a in zip(off,alt):
        ai=bi.clone(); ai[int(o)]=BASE.get(int(a),4); eff.append(enf_tss(ai)-ref)
    eff=np.array(eff); pred=eff@dos
    if np.std(pred)==0: rows.append({"gene":gene,"note":"flat"}); continue
    p,w=pw(expr,pred,pops); rows.append({"gene":gene,"n_var":int(len(off)),"pooled":p,"within":w})
    print("["+str(gi+1)+"/"+str(len(genes))+"] "+gene+" pooled="+format(p,".3f")+" within="+format(w,".3f"))
sc=[r for r in rows if "pooled" in r]
summ={"n_genes_scored":len(sc),
      "mean_abs_pooled_enformer":float(np.mean([r["pooled"] for r in sc])) if sc else None,
      "mean_abs_within_enformer":float(np.nanmean([r["within"] for r in sc])) if sc else None,
      "linear_ceiling_within":0.19,"n_tracks":len(track_idx)}
json.dump({"summary":summ,"per_gene":rows},open("/kaggle/working/stage2_results.json","w"),indent=2)
print("STAGE-2 ENFORMER SUMMARY: "+json.dumps(summ))
