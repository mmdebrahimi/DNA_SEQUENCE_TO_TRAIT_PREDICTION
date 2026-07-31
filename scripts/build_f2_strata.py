"""Build coding-vs-intergenic stratified windows for the F2 structured-signal test (canFam4).

Fetches a canFam4 chr1 region as GenBank-with-features via NCBI efetch (UCSC's api.genome.ucsc.edu is
unreachable from this host), parses CDS exon intervals + the ORIGIN sequence, and emits stratified
1200 bp windows: coding-bearing windows (centered on the largest CDS exons) + gene-desert intergenic
windows, each token tagged coding/intergenic/mixed by whether its 6 bases fall inside a CDS. Also
fetches a DISJOINT region for leakage-free Markov training. Output feeds scripts/kaggle_dog_nt_f2_strata.py.

    uv run python scripts/build_f2_strata.py --out data/processed/f2_strata.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_BASES = set("ACGT")


def _efetch(accession: str, start: int, stop: int, rettype: str) -> str:
    q = {"db": "nuccore", "id": accession, "seq_start": str(start), "seq_stop": str(stop),
         "rettype": rettype, "retmode": "text"}
    req = urllib.request.Request(f"{_EUTILS}?{urllib.parse.urlencode(q)}",
                                 headers={"User-Agent": "dna_decode/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read().decode("utf-8")


def _cds_intervals_and_seq(gb: str):
    feats = re.findall(r"^     CDS {13}(.*?)(?=^     \S)", gb, flags=re.M | re.S)
    ivs = []
    for f in feats:
        loc = "".join(f.split())
        ivs += [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
    origin = gb.split("ORIGIN", 1)[1].split("//")[0]
    seq = "".join(re.findall(r"[acgtnACGTN]", origin)).upper()
    return ivs, seq


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build F2 coding/intergenic stratified windows.")
    ap.add_argument("--accession", default="NC_051804.1")
    ap.add_argument("--region-start", type=int, default=20000000)
    ap.add_argument("--region-stop", type=int, default=21000000)
    ap.add_argument("--train-start", type=int, default=22000000)
    ap.add_argument("--train-stop", type=int, default=22200000)
    ap.add_argument("--win", type=int, default=1200)
    ap.add_argument("--kmer", type=int, default=6)
    ap.add_argument("--max-coding-tokens", type=int, default=600)
    ap.add_argument("--n-intergenic", type=int, default=18)
    ap.add_argument("--out", default="data/processed/f2_strata.json")
    args = ap.parse_args(argv)

    gb = _efetch(args.accession, args.region_start, args.region_stop, "gbwithparts")
    ivs, seq = _cds_intervals_and_seq(gb)
    n = len(seq)
    coding = bytearray(n + 2)  # 1-based coding flag
    for a, b in ivs:
        for i in range(a, min(b, n) + 1):
            coding[i] = 1
    W, K = args.win, args.kmer

    def stratum(start0: int) -> str:
        flags = [coding[start0 + 1 + j] for j in range(K)]
        return "coding" if all(flags) else ("intergenic" if not any(flags) else "mixed")

    windows = []
    used, ctok = [], 0
    for a, b in sorted(ivs, key=lambda ab: ab[1] - ab[0], reverse=True):
        if ctok >= args.max_coding_tokens:
            break
        c = (a + b) // 2
        s = max(0, min(n - W, c - W // 2))
        if any(abs(s - u) < W // 2 for u in used):
            continue
        strata = [stratum(s + t * K) for t in range(W // K)]
        nc = strata.count("coding")
        if nc == 0:
            continue
        used.append(s); ctok += nc
        windows.append({"start": args.region_start + s, "seq": seq[s:s + W],
                        "token_strata": strata, "kind": "coding-bearing"})
    ig, pos = 0, 0
    while ig < args.n_intergenic and pos < n - W:
        strata = [stratum(pos + t * K) for t in range(W // K)]
        if strata.count("coding") == 0 and strata.count("mixed") == 0 and set(seq[pos:pos + W]) <= _BASES:
            windows.append({"start": args.region_start + pos, "seq": seq[pos:pos + W],
                            "token_strata": strata, "kind": "intergenic"})
            ig += 1; pos += 53000
        else:
            pos += W

    traw = _efetch(args.accession, args.train_start, args.train_stop, "fasta")
    train = "".join(l.strip() for l in traw.splitlines() if not l.startswith(">")).upper()

    tot_c = sum(w["token_strata"].count("coding") for w in windows)
    tot_i = sum(w["token_strata"].count("intergenic") for w in windows)
    out = {"kmer": K, "win_bp": W, "n_windows": len(windows), "coding_tokens": tot_c,
           "intergenic_tokens": tot_i, "train_seq": train, "windows": windows}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out))
    print(f"windows={len(windows)} coding_tokens={tot_c} intergenic_tokens={tot_i} "
          f"train_bp={len(train)} -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
