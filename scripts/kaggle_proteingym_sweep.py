#!/usr/bin/env python
"""SELF-CONTAINED Kaggle-GPU notebook: ESM2 (650M / 3B) on the FULL ProteinGym v1.1 benchmark.

Turns the n=7 humsavar anecdote into a real, field-comparable world-model quality number: scores all
~217 ProteinGym v1.1 substitution DMS assays with a large ESM2 on GPU and reports the median Spearman
vs the published leaderboard. This is the true "quality" metric for the sequence model.

HOW TO RUN ON KAGGLE
--------------------
1. kaggle.com -> New Notebook. Settings: Accelerator = GPU T4 (NOT P100); Internet = ON.
   (P100 = sm_60 < CC 7.0 floor of Kaggle's current torch; fp16 fails. Code falls back to fp32.)
2. Paste this whole file into ONE cell. Run.
3. It prints a running median every 20 assays + a final MEDIAN Spearman vs the ESM2 leaderboard.
   650M over all 217 with windowing is long (~1-3 h). For a fast first read set MAX_SEQLEN=400
   (short-protein subset, ~15-20 min) or MAX_ASSAYS=40.

Data: range-reads only the 43 MB nested substitutions zip + ref.csv from the 11 GB Zenodo record
(no 11 GB download), via `remotezip`. Reads only; writes nothing outside the session.

Reference ESM2 zero-shot numbers, recomputed 2026-07-09 from ProteinGym's own per-assay CSV
(benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv, 217 assays):

  model | per-assay MEDIAN | per-assay mean | published Average_Spearman
   35M  |      0.349       |     0.319      |     0.321
  150M  |      0.451       |     0.401      |     0.387
  650M  |    **0.484**     |     0.438      |     0.414     <- PEAK of the ESM2 family
    3B  |      0.467       |     0.432      |     0.406
   15B  |      0.438       |     0.425      |     0.400

Scale does NOT flatten after 650M -- it REGRESSES. 650M is the best ESM2 checkpoint on this
benchmark, and 3B/15B are monotonically worse on all three aggregations.

Compare THIS script's output against the per-assay MEDIAN column (it prints a plain per-assay
median + mean). The published Average_Spearman is NOT a per-assay mean: it is the mean of the 5
function-category averages (Activity/Binding/Expression/OrganismalFitness/Stability), which is why
650M reads 0.414 there but 0.438 as a plain mean. Do not compare the two directly.
"""
import io
import re
import statistics
import subprocess
import sys
import zipfile

MODEL = "facebook/esm2_t33_650M_UR50D"   # <- facebook/esm2_t36_3B_UR50D for the 3B run
BATCH = 16                                # masked positions per forward
MAXLEN = 1022                             # ESM context cap (residues); long proteins => sliding window
MAX_SEQLEN = 0                            # >0 => skip assays whose target is longer (fast subset). 0 = all
MAX_ASSAYS = 0                            # >0 => cap assay count (smoke). 0 = all

ZIP_URL = "https://zenodo.org/api/records/13936340/files/ProteinGym_v1.1.zip/content"
INNER = "ProteinGym_v1.1/DMS_ProteinGym_substitutions.zip"
REF = "ProteinGym_v1.1/DMS_substitutions.csv"
AA = list("ACDEFGHIKLMNPQRSTVWY")
_SUB = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j < len(v) and v[order[j]] == v[order[i]]:
                j += 1
            for k in range(i, j):
                r[order[k]] = (i + j - 1) / 2.0
            i = j
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def fetch_proteingym():
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "remotezip"], check=True)
    from remotezip import RemoteZip
    import csv
    assays = {}          # DMS_id -> {mutant: score}
    with RemoteZip(ZIP_URL) as z:
        ref_txt = z.read(REF).decode()
        inner = z.read(INNER)
    seqs = {}            # DMS_id -> target_seq
    rdr = csv.DictReader(io.StringIO(ref_txt))
    fn2id = {}
    for row in rdr:
        seqs[row["DMS_id"]] = row["target_seq"]
        fn2id[row["DMS_filename"]] = row["DMS_id"]
    iz = zipfile.ZipFile(io.BytesIO(inner))
    for n in iz.namelist():
        if not n.endswith(".csv"):
            continue
        base = n.rsplit("/", 1)[-1]
        did = fn2id.get(base) or base[:-4]
        d = {}
        for row in csv.DictReader(io.StringIO(iz.read(n).decode())):
            mut = row.get("mutant", "")
            if ":" in mut:            # multi-mutant — skip (single-substitution benchmark)
                continue
            try:
                d[mut] = float(row["DMS_score"])
            except (KeyError, ValueError):
                pass
        if did in seqs and d:
            assays[did] = d
    return assays, seqs


def variant_scores(seq, variants, model, tok, device, torch):
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}
    positions = sorted({int(_SUB.match(v).group(2)) for v in variants if _SUB.match(v)})
    positions = [p for p in positions if 1 <= p <= len(seq)]
    logp = {}
    with torch.no_grad():
        for s in range(0, len(positions), BATCH):
            chunk = positions[s:s + BATCH]
            wins, rels = [], []
            for pos in chunk:
                if len(seq) <= MAXLEN:
                    start = 0
                else:
                    start = min(max(0, pos - MAXLEN // 2), len(seq) - MAXLEN)
                wins.append(seq[start:start + MAXLEN])
                rels.append(pos - start)          # token index (CLS at 0) of the masked residue
            enc = tok(wins, return_tensors="pt", padding=True)
            ids = enc["input_ids"].to(device)
            am = enc["attention_mask"].to(device)
            for r, rel in enumerate(rels):
                ids[r, rel] = mask_id
            out = model(ids, attention_mask=am).logits.float().log_softmax(-1)
            for r, (pos, rel) in enumerate(zip(chunk, rels)):
                logp[pos] = {aa: out[r, rel, i].item() for aa, i in aa_ids.items()}
    sc = {}
    for v in variants:
        m = _SUB.match(v)
        if not m:
            continue
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if pos in logp and 1 <= pos <= len(seq) and seq[pos - 1] == wt:
            sc[v] = logp[pos][mut] - logp[pos][wt]
    return sc


def main():
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}  model={MODEL}")
    assays, seqs = fetch_proteingym()
    ids = sorted(assays)
    if MAX_SEQLEN:
        ids = [d for d in ids if len(seqs[d]) <= MAX_SEQLEN]
    if MAX_ASSAYS:
        ids = ids[:MAX_ASSAYS]
    print(f"scoring {len(ids)} assays (of {len(assays)} available)")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForMaskedLM.from_pretrained(MODEL).eval().to(device)
    if device == "cuda" and torch.cuda.get_device_capability()[0] >= 7:
        model = model.half()   # sm_60 (Kaggle's legacy P100) is unsupported by current torch => stay fp32
    rhos = []
    for k, did in enumerate(ids, 1):
        dms = assays[did]
        seq = seqs[did]
        try:
            esm = variant_scores(seq, list(dms), model, tok, device, torch)
        except Exception as e:                         # keep going; report the failure
            print(f"  [{k}/{len(ids)}] {did[:45]:45s} ERROR {type(e).__name__}")
            continue
        common = [v for v in dms if v in esm]
        if len(common) < 10:
            continue
        rho = spearman([esm[v] for v in common], [dms[v] for v in common])
        rhos.append(rho)
        if k % 10 == 0 or k == len(ids):
            print(f"  [{k}/{len(ids)}] {did[:40]:40s} len={len(seq):4d} n={len(common):5d} "
                  f"rho={rho:+.3f}  running median={statistics.median(rhos):.3f}")
    print(f"\nMEDIAN Spearman over {len(rhos)} assays = {statistics.median(rhos):.3f}  "
          f"(mean {statistics.mean(rhos):.3f})   model={MODEL.split('_')[2]}")
    print("ProteinGym per-assay reference (217 assays)  median | mean:")
    print("  ESM2 650M 0.484 | 0.438   <- peak      ESM2 3B 0.467 | 0.432   ESM2 15B 0.438 | 0.425")
    print("  Scale REGRESSES past 650M. Do not compare to the published Average_Spearman (0.414 for")
    print("  650M) -- that is the mean of 5 function-category averages, not a per-assay mean.")


main()
