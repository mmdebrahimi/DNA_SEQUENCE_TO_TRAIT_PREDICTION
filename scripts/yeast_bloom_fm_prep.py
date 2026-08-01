"""Layer-2 B2 prep: build allele-context sequence windows for a FM-embedding test on the Bloom yeast
substrate. Picks LD-pruned genome-wide top-QTL markers for a trait, fetches the S288C reference window
around each (NCBI efetch; the marker name encodes chr+pos+REF+ALT), verifies REF matches the reference
base, and emits REF/ALT allele windows for Nucleotide-Transformer embedding on Kaggle.

The FM-value question: does an NT embedding of a variant's sequence CONTEXT add predictive signal beyond
the raw allele identity (which the +/-1 marker already encodes)? Prior: in a bi-parental cross the markers
are a sufficient statistic, so the FM is expected to tie -- this tests that empirically.

    uv run python scripts/yeast_bloom_fm_prep.py --trait Maltose --k 40 --window 120
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np  # noqa: E402

from scripts.yeast_bloom_gp_arm import load_genotype, load_phenotype  # noqa: E402

# S288C R64 RefSeq chromosome accessions, chr01..chr16 -> NC_...
_CHR = {f"chr{n:02d}": acc for n, acc in enumerate(
    ["NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10", "NC_001137.3", "NC_001138.5",
     "NC_001139.9", "NC_001140.6", "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4"], start=1)}
_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _efetch_window(acc, start, stop):
    import time

    q = {"db": "nuccore", "id": acc, "seq_start": str(start), "seq_stop": str(stop),
         "rettype": "fasta", "retmode": "text"}
    url = f"{_EUTILS}?{urllib.parse.urlencode(q)}"
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "dna_decode/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                text = r.read().decode("utf-8")
            time.sleep(0.4)  # stay under NCBI's 3 req/s limit
            return "".join(l.strip() for l in text.splitlines() if not l.startswith(">")).upper()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
            raise
    raise RuntimeError(f"efetch failed after retries: {acc}:{start}-{stop}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", default="Maltose")
    ap.add_argument("--k", type=int, default=40)
    ap.add_argument("--window", type=int, default=120, help="half-window each side of the variant")
    ap.add_argument("--min-dist", type=int, default=20000, help="LD-prune: min bp between selected markers")
    ap.add_argument("--out", default="D:/dna_decode_cache/bloom/fm_windows_Maltose.json")
    args = ap.parse_args(argv)

    g_ids, markers, G = load_genotype("D:/dna_decode_cache/bloom/geno_v2.txt", 1)
    p_ids, traits, P = load_phenotype("D:/dna_decode_cache/bloom/BYxRM_PhenoData.txt")
    gpos = {s: i for i, s in enumerate(g_ids)}
    common = [s for s in p_ids if s in gpos]
    X = G[[gpos[s] for s in common]]
    y = P[[p_ids.index(s) for s in common], traits.index(args.trait)]
    keep = ~np.isnan(y)
    Xk, yk = X[keep], y[keep]
    Xc = (Xk - Xk.mean(0)) / (Xk.std(0) + 1e-9)
    yc = (yk - yk.mean()) / (yk.std() + 1e-9)
    corr = np.abs(Xc.T @ yc) / len(yk)

    # LD-prune: greedily take top-corr markers >= min_dist apart (per chromosome)
    order = np.argsort(-corr)
    chosen, chosen_by_chr = [], {}
    for i in order:
        parts = markers[i].split("_")
        ch, pos = parts[1], int(parts[2])
        if any(abs(pos - p) < args.min_dist for p in chosen_by_chr.get(ch, [])):
            continue
        chosen.append(i)
        chosen_by_chr.setdefault(ch, []).append(pos)
        if len(chosen) >= args.k:
            break

    seg_ids = common
    win = args.window
    out = {"trait": args.trait, "window_half": win, "seg_ids": seg_ids, "markers": []}
    for rank, i in enumerate(chosen):
        parts = markers[i].split("_")
        ch, pos, ref, alt = parts[1], int(parts[2]), parts[3], parts[4]
        acc = _CHR[ch]
        seq = _efetch_window(acc, pos - win, pos + win)
        center = win  # 0-based index of the variant in the fetched window (pos-win .. pos+win)
        ref_ok = center < len(seq) and seq[center] == ref
        ref_win = seq[:center] + ref + seq[center + 1:]
        alt_win = seq[:center] + alt + seq[center + 1:]
        out["markers"].append({
            "marker_index": int(i), "name": markers[i], "chr": ch, "pos": pos, "ref": ref, "alt": alt,
            "abs_corr": round(float(corr[i]), 3), "ref_matches_reference": bool(ref_ok),
            "ref_window": ref_win, "alt_window": alt_win,
            # per-segregant allele: B -> REF, R -> ALT (B=BY=S288C reference strain)
            "seg_alleles": "".join("B" if v < 0 else "R" for v in G[[gpos[s] for s in seg_ids], i]),
        })
        if rank % 10 == 0:
            print(f"  fetched {rank + 1}/{len(chosen)}  {markers[i]}  ref_ok={ref_ok}", flush=True)

    Path(args.out).write_text(json.dumps(out))
    n_ok = sum(m["ref_matches_reference"] for m in out["markers"])
    chrs = sorted(set(m["chr"] for m in out["markers"]))
    print(f"\n{len(chosen)} LD-pruned QTL markers across {len(chrs)} chromosomes; "
          f"REF matches reference {n_ok}/{len(chosen)} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
