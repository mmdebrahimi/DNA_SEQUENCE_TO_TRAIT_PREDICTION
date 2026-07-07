"""Capture a curated ClinVar gene-panel catalog from the full ClinVar VCF (on D:) into a COMMITTED TSV.

ClinVar = the human Mendelian-disease analogue of the AMR determinant catalog: a free, curated
variant→clinical-significance map. This extracts a bounded, canonical gene panel (P/LP + B/LB variants only —
the deployable-claim tier; VUS excluded) into a small committed catalog that `dna_decode/data/clinvar.py`
loads. The full 192 MB VCF stays on D: (gitignored); only the committed subset ships. Same pattern as the
other curated catalogs (mic_tiers / hiv_amr / fungal_amr): committed curated data + a deterministic caller.

Gene panel (broadened 2026-07-06): the **ACMG SF v3.2** medically-actionable secondary-findings gene list
(81 genes; the canonical "which genes to report incidental findings" set, sourced from
ncbi.nlm.nih.gov/clinvar/docs/acmg) UNIONed with the original 5 recessive-carrier genes NOT in ACMG SF
(CFTR/HBB/PAH/G6PD/F8) -> 86 genes total. Extensible via --genes.
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_VCF = Path("D:/dna_decode_cache/clinvar/clinvar.vcf.gz")

# ACMG SF v3.2 (81 genes) — verified from ncbi.nlm.nih.gov/clinvar/docs/acmg (2026-07-06), NOT from memory.
ACMG_SF_V3_2 = [
    "APC", "MYH11", "ACTA2", "TMEM43", "DSP", "PKP2", "DSG2", "DSC2", "BTD", "BRCA1", "BRCA2", "SCN5A",
    "RYR2", "CASQ2", "CALM1", "TRDN", "FLNC", "LMNA", "TNNT2", "DES", "MYH7", "TNNC1", "RBM20", "BAG3",
    "TTN", "COL3A1", "GLA", "LDLR", "APOB", "TPM1", "MYBPC3", "PRKAG2", "TNNI3", "MYL3", "MYL2", "ACTC1",
    "RET", "PALB2", "HFE", "ENG", "ACVRL1", "SDHD", "SDHB", "TTR", "PCSK9", "BMPR1A", "SMAD4", "TP53",
    "TGFBR1", "TGFBR2", "SMAD3", "KCNQ1", "KCNH2", "CALM2", "CALM3", "MSH2", "MLH1", "PMS2", "MSH6", "RYR1",
    "CACNA1S", "FBN1", "HNF1A", "MEN1", "MUTYH", "NF2", "OTC", "SDHAF2", "SDHC", "STK11", "MAX", "TMEM127",
    "GAA", "PTEN", "RB1", "RPE65", "TSC1", "TSC2", "VHL", "WT1", "ATP7B",
]
# recessive-carrier genes from the original 10 that are NOT in ACMG SF (kept for carrier-status coverage)
_EXTRA_CARRIER = ["CFTR", "HBB", "PAH", "G6PD", "F8"]
DEFAULT_GENES = ACMG_SF_V3_2 + _EXTRA_CARRIER   # 86 genes (ACMG SF v3.2 union the 5 extra carrier genes)
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
