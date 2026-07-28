"""Build the cross-protein / higher-order epistasis bundle for the ESM2 sweep."""
import csv, gzip, io, json, random, urllib.request
from pathlib import Path

CACHE = Path("D:/dna_decode_cache/epistasis"); CACHE.mkdir(exist_ok=True)
HF = "https://huggingface.co/datasets/OATML-Markslab/ProteinGym/resolve/main/ProteinGym_substitutions"
PROTEINS = ["GFP_AEQVI_Sarkisyan_2016", "HIS7_YEAST_Pokusaeva_2019", "GRB2_HUMAN_Faure_2021",
            "SPG1_STRSG_Olson_2014", "F7YBW8_MESOW_Aakre_2015", "BLAT_ECOLX_Firnberg_2014"]
PER_ORDER_CAP = 300; MAX_ORDER = 6

ref = {r["DMS_id"]: r["target_seq"] for r in csv.DictReader(open(CACHE/"pg_ref.csv", encoding="utf-8"))}

def fetch(dms):
    p = CACHE / f"{dms}.csv"
    if p.exists() and p.stat().st_size > 1000: return p
    url = f"{HF}/{dms}.csv"
    print(f"  fetching {dms} ...")
    urllib.request.urlretrieve(url, p)
    return p

bundle = {"proteins": []}
for dms in PROTEINS:
    seq = ref.get(dms, "")
    if not seq: print(f"  {dms}: no ref seq, skip"); continue
    p = fetch(dms)
    by_order = {}
    with open(p, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            m = r["mutant"]; 
            try: s = float(r["DMS_score"])
            except: continue
            muts = m.split(":")
            k = len(muts)
            if k < 2 or k > MAX_ORDER: continue
            rec = []
            ok = True
            for x in muts:
                wt, pos, alt = x[0], int(x[1:-1]), x[-1]
                if pos > len(seq) or seq[pos-1] != wt: ok = False; break
                rec.append([pos, alt])
            if not ok: continue
            # distinct positions
            if len({p for p,_ in rec}) != k: continue
            by_order.setdefault(k, []).append({"muts": rec, "order": k, "dms": s})
    # stratified sample per order
    sampled = []
    for k, lst in sorted(by_order.items()):
        random.Random(0).shuffle(lst)
        sampled += lst[:PER_ORDER_CAP]
    if len(sampled) < 20: print(f"  {dms}: only {len(sampled)} usable multi, skip"); continue
    bundle["proteins"].append({"dms": dms, "wt_seq": seq, "variants": sampled,
                               "order_counts": {k: min(len(v), PER_ORDER_CAP) for k, v in by_order.items()}})
    print(f"  {dms}: seqlen={len(seq)} sampled={len(sampled)} orders={ {k:len(v) for k,v in by_order.items()} }")

json.dump(bundle, open(CACHE/"epi_sweep_bundle.json", "w"))
tot = sum(len(p["variants"]) for p in bundle["proteins"])
print(f"\nBUNDLE: {len(bundle['proteins'])} proteins, {tot} total multi-mutants")
import os; print("MB:", round(os.path.getsize(CACHE/'epi_sweep_bundle.json')/1e6, 2))
