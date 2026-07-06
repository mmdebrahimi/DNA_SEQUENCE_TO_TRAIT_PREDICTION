#!/usr/bin/env python
"""Run the deterministic ClinVar/Mendelian decoder over a WHOLE VCF — the per-VCF wrapper the caller lacked.

`dna_decode.data.clinvar.ClinVarDecoder` is a per-variant `call(chrom,pos,ref,alt)` dict-lookup. This
iterates an individual's VCF, calls the decoder on every ALT the person actually CARRIES (GT has the ALT),
and collects the curated ClinVar annotations (Pathogenic/Likely_pathogenic + Benign/Likely_benign, with
gene / disease / gold-star review level). Fail-closed: a variant not in the committed panel is
INDETERMINATE ("not-in-panel"; absence != benign), never a guess.

BUILD-AGNOSTIC — the decoder is a pure dict lookup, so the panel and the VCF just have to be the SAME build.
For a GRCh37 PGP-UK individual, pass `--panel data/clinvar/clinvar_panel_grch37.tsv` (built by
`scripts/capture_clinvar.py --vcf <clinvar_GRCh37.vcf.gz>`) — no liftover of the 100 MB individual VCF.

HONEST FRAMING: this is a research demonstration that the deterministic Mendelian decoder runs end-to-end on
a real individual and returns curated ClinVar classifications for the variants they carry. It is NOT a
clinical diagnosis — a ClinVar "Pathogenic" annotation in a research VCF is a curated database classification
of the allele, not a clinical interpretation of this person. NOT a clinical tool.

Usage:
    uv run python scripts/clinvar_decode_vcf.py --vcf D:/.../FR07961000.vcf.gz --sample-id FR07961000 \
        --panel data/clinvar/clinvar_panel_grch37.tsv
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from dna_decode.data.clinvar import ClinVarDecoder


def _carried_alts(ref: str, alt_field: str, gt: str) -> list[str]:
    """Return the ALT alleles this genotype actually carries (>=1 copy). '.'/no-call -> none."""
    alts = alt_field.split(",")
    idx = {str(i + 1): a for i, a in enumerate(alts)}
    carried = []
    for a in gt.replace("|", "/").split("/"):
        if a in idx and idx[a] not in carried:
            carried.append(idx[a])
    return carried


def decode_vcf(vcf: Path, decoder: ClinVarDecoder, sample: str | None = None) -> dict:
    opener = gzip.open if str(vcf).endswith(".gz") else open
    hits: list[dict] = []
    n_variants = 0
    n_indeterminate = 0
    sample_idx = 0
    with opener(vcf, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                if sample and sample in samples:
                    sample_idx = samples.index(sample)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _vid, ref, alt_field = cols[0], cols[1], cols[2], cols[3], cols[4]
            gt = None
            if len(cols) >= 10 and "GT" in cols[8].split(":"):
                col = 9 + sample_idx
                if col < len(cols):
                    gt = cols[col].split(":")[cols[8].split(":").index("GT")]
            if gt is None:
                continue
            carried = _carried_alts(ref, alt_field, gt)
            for a in carried:
                n_variants += 1
                call = decoder.call(chrom, pos, ref, a)
                if call.verdict == "INDETERMINATE":
                    n_indeterminate += 1
                    continue
                hits.append({"chrom": chrom.replace("chr", ""), "pos": pos, "ref": ref, "alt": a,
                             "gt": gt, "gene": call.gene, "verdict": call.verdict,
                             "significance": call.significance, "stars": call.stars,
                             "disease": call.disease, "provenance": call.provenance})
    path = [h for h in hits if h["verdict"] == "PATHOGENIC"]
    benign = [h for h in hits if h["verdict"] == "BENIGN"]
    path.sort(key=lambda h: -(h["stars"] or 0))
    return {"sample_id": sample or vcf.stem,
            "n_carried_alts_in_panel_genes_checked": n_variants,
            "n_indeterminate_not_in_panel": n_indeterminate,
            "n_pathogenic": len(path), "n_benign": len(benign),
            "pathogenic_hits": path, "benign_hits": benign[:50]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Run the ClinVar/Mendelian decoder over a whole individual VCF.")
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--panel", type=Path, default=None, help="ClinVar panel TSV (default: committed GRCh38)")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)

    decoder = ClinVarDecoder.from_tsv(args.panel) if args.panel else ClinVarDecoder.from_tsv()
    rep = decode_vcf(args.vcf, decoder, sample=args.sample_id)
    rep["schema"] = "clinvar-vcf-decode-v0"
    rep["panel"] = str(args.panel) if args.panel else "committed GRCh38 clinvar_panel.tsv"
    rep["honest_tier"] = ("research demonstration; curated ClinVar allele classifications for carried "
                          "variants, NOT a clinical diagnosis of the individual. NOT a clinical tool.")
    print(json.dumps({k: rep[k] for k in ("sample_id", "n_pathogenic", "n_benign",
                                          "n_indeterminate_not_in_panel", "pathogenic_hits")}, indent=2))
    if args.out:
        args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
        print(f"[provenance -> {args.out}]")
    return 0 if (rep["n_pathogenic"] + rep["n_benign"]) >= 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
