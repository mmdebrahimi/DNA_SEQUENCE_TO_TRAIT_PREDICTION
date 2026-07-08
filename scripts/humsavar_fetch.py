#!/usr/bin/env python
"""Fetch + parse UniProt humsavar — a large FREE human pathogenicity variant set (genotype->phenotype).

humsavar.txt (~8.7 MB) lists ~85k curated human single-missense variants across ~13k proteins, each labelled
LB/B (likely-benign/benign), LP/P (likely-pathogenic/pathogenic), or US (uncertain). It is the classic
free substrate for testing whether a protein model discriminates pathogenic from benign variants (a
world-model signal): fetch a protein's UniProt sequence, score its variants with ESM masked-marginals, and
check the pathogenic vs benign separation (see wiki/humsavar_mlh1_35M_pilot.json for a worked example).

Run:  python scripts/humsavar_fetch.py --out-dir <dir>   # stages humsavar.txt + writes summary.json
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import re
import urllib.request

URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
_MISS = re.compile(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}$")


def parse(text):
    """[(gene, ac, ftid, change, category, dbsnp)] for the VAR_ data rows + category counts."""
    rows, cat = [], collections.Counter()
    prots = set()
    for line in text.splitlines():
        if " VAR_" not in line:
            continue
        p = line.split()
        ft = [x for x in p if x.startswith("VAR_")]
        if not ft:
            continue
        i = p.index(ft[0])
        gene, ac, change, category = p[0], p[1], p[i + 1], p[i + 2]
        if not _MISS.match(change):
            continue
        rows.append((gene, ac, ft[0], change, category))
        cat[category] += 1
        prots.add(ac)
    return rows, cat, len(prots)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    dst = os.path.join(args.out_dir, "humsavar.txt")
    urllib.request.urlretrieve(URL, dst)
    text = open(dst, encoding="utf-8", errors="replace").read()
    rows, cat, n_prot = parse(text)
    summary = {"source": "UniProt humsavar", "url": URL, "n_missense": len(rows),
               "n_proteins": n_prot, "categories": dict(cat)}
    json.dump(summary, open(os.path.join(args.out_dir, "humsavar_summary.json"), "w"), indent=2)
    print(f"staged humsavar.txt ({len(text)//1_000_000}MB) + summary: {len(rows)} missense / {n_prot} "
          f"proteins, categories={dict(cat)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
