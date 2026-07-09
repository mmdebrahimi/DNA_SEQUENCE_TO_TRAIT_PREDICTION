#!/usr/bin/env python
"""SELF-CONTAINED Kaggle-GPU notebook: ESM2 (650M / 3B) vs AlphaMissense on humsavar (7 proteins).

Answers "does a BIGGER LM close the AlphaMissense gap?" on a free Kaggle GPU. The laptop baseline
(ESM2-35M, CPU) is median AUROC 0.706 vs AlphaMissense 0.823 (+0.117). This reruns the SAME 7-protein
head-to-head with a larger ESM on GPU.

HOW TO RUN ON KAGGLE
--------------------
1. kaggle.com -> Create -> New Notebook.
2. Settings (right panel): Accelerator = GPU T4 x2; Internet = ON.
   (NOT P100 -- sm_60 is below the CC>=7.0 floor of Kaggle's current torch; fp16 dies with
    'no kernel image is available for execution on the device'. The code now falls back to fp32.)
3. Paste this whole file into ONE code cell. Run. (~5-10 min for 650M; ~15-25 min for 3B.)
4. Read the final MEDIAN line + per-protein table. Paste it back.

To try the 3B model, set MODEL = "facebook/esm2_t36_3B_UR50D" (fits a 16 GB T4 in fp16).
Nothing here writes outside the Kaggle session; it only READS UniProt + AlphaMissense over the internet.
"""
import gzip
import re
import statistics
import urllib.request

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

MODEL = "facebook/esm2_t33_650M_UR50D"   # <- switch to facebook/esm2_t36_3B_UR50D for the 3B run
BATCH = 32                                # masked positions per forward pass (GPU)
MAXLEN = 1022

AA = list("ACDEFGHIKLMNPQRSTVWY")
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PROTS = [("P40692", "MLH1"), ("P43246", "MSH2"), ("P04070", "PROC"), ("P02730", "SLC4A1"),
         ("P51608", "MECP2"), ("Q9BXM7", "PINK1"), ("Q99972", "MYOC")]
UNIS = {ac for ac, _ in PROTS}
BASELINE_35M = {"P40692": 0.666, "P43246": 0.750, "P04070": 0.531, "P02730": 0.848,
                "P51608": 0.728, "Q9BXM7": 0.706, "Q99972": 0.531}  # laptop CPU 35M AUROC


def auroc(labels, scores):
    pairs = sorted(zip(scores, labels))
    n = len(pairs)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and pairs[j][0] == pairs[i][0]:
            j += 1
        for k in range(i, j):
            ranks[k] = (i + 1 + j) / 2.0
        i = j
    npos = sum(labels)
    nneg = n - npos
    if not npos or not nneg:
        return float("nan")
    sr = sum(r for r, (_, l) in zip(ranks, pairs) if l == 1)
    return (sr - npos * (npos + 1) / 2.0) / (npos * nneg)


def ranknorm(vals):
    """Rank-normalize to [0,1] (ties get sequential ranks; fine for ensembling)."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    n = len(vals)
    for rank, i in enumerate(order):
        r[i] = rank / (n - 1) if n > 1 else 0.5
    return r


def humsavar_labels():
    urllib.request.urlretrieve(
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt",
        "humsavar.txt")
    by = {ac: {} for ac in UNIS}
    for line in open("humsavar.txt", encoding="utf-8", errors="replace"):
        if " VAR_" not in line:
            continue
        p = line.split()
        if p[1] not in UNIS:
            continue
        ft = [x for x in p if x.startswith("VAR_")]
        i = p.index(ft[0])
        m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", p[i + 1])
        cat = p[i + 2]
        if not m or cat not in ("LP/P", "LB/B"):
            continue
        wt, pos, mut = AA3.get(m.group(1)), int(m.group(2)), AA3.get(m.group(3))
        if wt and mut and wt != mut:
            by[p[1]][f"{wt}{pos}{mut}"] = 1 if cat == "LP/P" else 0
    return by


def fetch_seq(ac):
    urllib.request.urlretrieve(f"https://rest.uniprot.org/uniprotkb/{ac}.fasta", f"{ac}.fasta")
    return "".join(l.strip() for l in open(f"{ac}.fasta") if not l.startswith(">"))


def alphamissense():
    """Stream the 1.2 GB AlphaMissense gz, keep only the 7 proteins' rows: {(uniprot,variant): am}."""
    am = {}
    req = urllib.request.urlopen(
        "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz")
    with gzip.GzipFile(fileobj=req) as gz:
        for raw in gz:
            line = raw.decode()
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")   # uniprot_id, protein_variant, am_pathogenicity, am_class
            if len(p) >= 3 and p[0] in UNIS:
                am[(p[0], p[1])] = float(p[2])
    return am


def esm_scores(seq, variants, model, tok, device):
    mask_id = tok.mask_token_id
    aa_ids = {aa: tok.convert_tokens_to_ids(aa) for aa in AA}
    positions = sorted({int(re.match(r"[A-Z](\d+)[A-Z]", v).group(1)) for v in variants})
    positions = [p for p in positions if p <= min(len(seq), MAXLEN)]
    ids = tok(seq[:MAXLEN], return_tensors="pt")["input_ids"].to(device)
    logp = {}
    with torch.no_grad():
        for s in range(0, len(positions), BATCH):
            chunk = positions[s:s + BATCH]
            batch = ids.repeat(len(chunk), 1)
            for r, pos in enumerate(chunk):
                batch[r, pos] = mask_id
            out = model(batch).logits.float().log_softmax(-1)
            for r, pos in enumerate(chunk):
                logp[pos] = {aa: out[r, pos, i].item() for aa, i in aa_ids.items()}
    sc = {}
    for v in variants:
        m = re.match(r"([A-Z])(\d+)([A-Z])", v)
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if pos in logp and seq[pos - 1] == wt:
            sc[v] = logp[pos][mut] - logp[pos][wt]
    return sc


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}  model={MODEL}")
    tok = AutoTokenizer.from_pretrained(MODEL)
    model = AutoModelForMaskedLM.from_pretrained(MODEL).eval().to(device)
    if device == "cuda" and torch.cuda.get_device_capability()[0] >= 7:
        model = model.half()   # sm_60 (Kaggle's legacy P100) is unsupported by current torch => stay fp32
    labels = humsavar_labels()
    am = alphamissense()
    rows = []
    for ac, gene in PROTS:
        seq = fetch_seq(ac)
        lab = labels[ac]
        esm = esm_scores(seq, list(lab), model, tok, device)
        sh = sorted(k for k in lab if k in esm and (ac, k) in am)
        if len(sh) < 10:
            print(f"{gene}: too few shared ({len(sh)})")
            continue
        y = [lab[v] for v in sh]
        esm_p = [-esm[v] for v in sh]          # higher => more pathogenic
        am_p = [am[(ac, v)] for v in sh]
        ae = auroc(y, esm_p)
        aa = auroc(y, am_p)
        # ENSEMBLE: rank-normalize both pathogenicity-oriented scores to [0,1], average
        ens = [(a + b) / 2 for a, b in zip(ranknorm(esm_p), ranknorm(am_p))]
        aen = auroc(y, ens)
        rows.append((gene, ac, len(sh), ae, aa, aen))
        print(f"{gene:8s} n={len(sh):3d}  ESM={ae:.3f}  AM={aa:.3f}  ENSEMBLE={aen:.3f}   "
              f"(35M was {BASELINE_35M.get(ac, float('nan')):.3f})")
    esm_med = statistics.median([r[3] for r in rows])
    am_med = statistics.median([r[4] for r in rows])
    ens_med = statistics.median([r[5] for r in rows])
    # Ensemble lift MUST be paired. median(ENS) - max(median(ESM), median(AM)) takes its
    # medians over different proteins, so it reports a "lift" even when the ensemble loses
    # on most proteins (650M: unpaired +0.005, but 4/7 paired losses, mean delta -0.0006).
    lifts = [r[5] - max(r[3], r[4]) for r in rows]
    lift_med, wins = statistics.median(lifts), sum(l > 0 for l in lifts)
    print(f"\nMEDIAN over {len(rows)}: ESM({MODEL.split('_')[2]})={esm_med:.3f}  "
          f"AlphaMissense={am_med:.3f}  ENSEMBLE={ens_med:.3f}")
    print(f"  gap ESM-vs-AM (median)={am_med-esm_med:+.3f}   (35M baseline median 0.706; AM 0.823)")
    print(f"  PAIRED ensemble lift vs best-single: median={lift_med:+.4f} "
          f"mean={statistics.mean(lifts):+.4f} wins={wins}/{len(rows)}")
    print("VERDICT:", "bigger LM CLOSES the gap" if am_med - esm_med < 0.06
          else "gap PERSISTS — supervised still wins",
          "| ENSEMBLE helps" if lift_med > 0.005 and wins > len(rows) / 2
          else "| ensemble ~= best single (no paired lift)")


main()
