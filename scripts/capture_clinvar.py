"""Capture a curated ClinVar gene-panel catalog from the full ClinVar VCF (on D:) into a COMMITTED TSV.

ClinVar = the human Mendelian-disease analogue of the AMR determinant catalog: a free, curated
variant→clinical-significance map. This extracts a bounded, canonical gene panel (P/LP + B/LB variants only —
the deployable-claim tier; VUS excluded) into a small committed catalog that `dna_decode/data/clinvar.py`
loads. The full 192 MB VCF stays on D: (gitignored); only the committed subset ships. Same pattern as the
other curated catalogs (mic_tiers / hiv_amr / fungal_amr): committed curated data + a deterministic caller.

Demo gene panel (canonical, sourced Mendelian genes across AR/AD/X-linked; extensible via --genes):
  CFTR (cystic fibrosis) · HBB (sickle/β-thal) · LDLR (FH) · BRCA1/BRCA2 (HBOC) · TP53 (Li-Fraumeni) ·
  PAH (PKU) · G6PD (X-linked) · HFE (haemochromatosis) · F8 (haemophilia A).
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_VCF = Path("D:/dna_decode_cache/clinvar/clinvar.vcf.gz")
DEFAULT_GENES = ["CFTR", "HBB", "LDLR", "BRCA1", "BRCA2", "TP53", "PAH", "G6PD", "HFE", "F8"]
# the deployable-claim significance tiers (VUS/other excluded from the committed catalog)
KEEP_SIG = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
            "Benign", "Likely_benign", "Benign/Likely_benign"}


def _info(field: str, info: str) -> str | None:
    m = re.search(rf"(?:^|;){field}=([^;]*)", info)
    return m.group(1) if m else None


def capture(vcf: Path, genes: list[str], out: Path) -> dict:
    want = set(genes)
    rows, seen_sig = [], {}
    with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            chrom, pos, vid, ref, alt, _, _, info = c[:8]
            gi = _info("GENEINFO", info)                      # "SYMBOL:ID|SYMBOL2:ID2"
            if not gi:
                continue
            gene = gi.split(":")[0]
            if gene not in want:
                continue
            sig = (_info("CLNSIG", info) or "").replace("\\x2c", ",")
            if sig not in KEEP_SIG:
                continue
            rows.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "gene": gene,
                         "significance": sig, "review_status": _info("CLNREVSTAT", info) or "",
                         "disease": (_info("CLNDN", info) or "").replace("_", " ")[:80],
                         "clinvar_id": vid})
            seen_sig[sig] = seen_sig.get(sig, 0) + 1
    out.parent.mkdir(parents=True, exist_ok=True)
    hdr = ["chrom", "pos", "ref", "alt", "gene", "significance", "review_status", "disease", "clinvar_id"]
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]).replace("\t", " ") for k in hdr) + "\n")
    return {"n_variants": len(rows), "genes": genes, "by_significance": seen_sig, "out": str(out)}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    ap.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    ap.add_argument("--out", type=Path, default=REPO / "data" / "clinvar" / "clinvar_panel.tsv")
    a = ap.parse_args(argv)
    if not a.vcf.exists():
        print(f"ERROR: ClinVar VCF not found at {a.vcf}", file=sys.stderr)
        return 2
    res = capture(a.vcf, a.genes, a.out)
    print(f"captured {res['n_variants']} P/LP+B/LB variants across {len(a.genes)} genes -> {res['out']}")
    print(f"  by significance: {res['by_significance']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
