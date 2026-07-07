"""Core VCF-iterating Mendelian decode logic (importable; shared by `dna-clinvar` CLI + the script wrapper).

`dna_decode.data.clinvar.ClinVarDecoder` is a per-variant `call(chrom,pos,ref,alt)` dict-lookup. This iterates
an individual's VCF, decodes every ALT the person actually CARRIES (GT has the ALT), and collects the curated
ClinVar classifications (P/LP + B/LB, with gene / disease / gold-star review level). Build-agnostic: the
decoder is a pure dict lookup, so the panel and the VCF just have to be the SAME genome build (pass a GRCh37
panel for a GRCh37 VCF — no liftover of the individual genome).
"""
from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

from dna_decode.data.clinvar import ClinVarDecoder


def carried_alts(ref: str, alt_field: str, gt: str) -> list[str]:
    """Return the ALT alleles this genotype actually carries (>=1 copy). '.'/no-call -> none."""
    alts = alt_field.split(",")
    idx = {str(i + 1): a for i, a in enumerate(alts)}
    carried: list[str] = []
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
            for a in carried_alts(ref, alt_field, gt):
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
    return {"sample_id": sample or Path(vcf).stem,
            "n_carried_alts_in_panel_genes_checked": n_variants,
            "n_indeterminate_not_in_panel": n_indeterminate,
            "n_pathogenic": len(path), "n_benign": len(benign),
            "benign_by_gene": dict(Counter(h["gene"] for h in benign).most_common()),
            "pathogenic_by_gene": dict(Counter(h["gene"] for h in path).most_common()),
            "pathogenic_hits": path, "benign_hits": benign[:50]}
