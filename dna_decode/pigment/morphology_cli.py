"""`dna-decode morphology` — dog body-size + ear-type decoder (curated + pinned + validated catalog, v0).

    dna-decode morphology --dosages IGF1=2,HMGA2=2,STC2=1,GHR=1,EAR=2      # large dog, erect ears
    dna-decode morphology --dosages IGF1=0,HMGA2=0,STC2=0,GHR=0,EAR=0      # toy/small, drop ears
    dna-decode morphology --dosages HMGA2=1,IGF1=1 --json                  # partial panel (2/4 size loci)
    dna-decode morphology --vcf dog.vcf                                    # decode a real canFam4 dog VCF

Deterministic RELATIVE-signal decoder over the loci PINNED + functionally VALIDATED on Darwin's Ark
(dna_decode.pigment.dog_body_size): body SIZE = an additive polygenic score over IGF1/HMGA2/STC2/GHR
(height r=+0.619), plus single-SNP EAR type (MSRB3, r=+0.543). Input is per-locus BIG-ALLELE DOSAGE (0/1/2) —
the natural PLINK-panel shape — OR a canFam4 dog genome VCF (--vcf), from which the pinned SNP dosages are
called by coordinate (dog_vcf_input; the dog causal SNPs have no rsIDs, so matched by chr:pos).

HONEST SCOPE: RELATIVE size RANK + ear axis, NOT calibrated absolute height/inches (Q121 is a
covariate-adjusted z-score). Coat length/curl (FGF5/KRT71) + leg length (FGF4) + the 4 covariate-adjusted
"rerun" morph traits ABSTAIN — no strong single-known-SNP mapping on this substrate. Benign companion-animal
visible-trait genetics, NOT human/forensic. Frozen AMR/forward surfaces untouched.
"""
from __future__ import annotations

import argparse
import json
import sys

# axes the cell does NOT resolve (surfaced as honest abstentions, not silent gaps)
_ABSTAIN_AXES = [
    "coat length / curl (FGF5 / KRT71) — showed no strong signal on any Darwin's Ark morph question",
    "leg length (FGF4 chondrodysplasia retrogene) — a structural insertion, absent from a SNV panel",
    "the 4 covariate-adjusted rerun morph traits (Q124/127/128/245) — no known single-SNP mapping (max |r|=0.21)",
]


def _parse_dosages(spec: str) -> dict:
    out = {}
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" not in tok:
            raise ValueError(f"bad token {tok!r}; expected LOCUS=dose (e.g. HMGA2=2)")
        loc, dose = tok.split("=", 1)
        try:
            out[loc.strip().upper()] = int(dose.strip())
        except ValueError:
            raise ValueError(f"dosage for {loc!r} must be an integer 0/1/2, got {dose!r}")
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode morphology",
        description="Dog body-size + ear-type decoder (pinned + Darwin's-Ark-validated catalog, deterministic).",
        epilog="scope: RELATIVE size rank + ear axis (not absolute inches); benign companion-animal genetics. "
               "coat length/curl/leg-length + the 4 rerun morph traits ABSTAIN. NOT a human/forensic tool.",
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--dosages",
                     help="comma-separated LOCUS=dose big-allele copies (0/1/2) for any of "
                          "IGF1/HMGA2/STC2/GHR (size) + EAR (ear type), e.g. IGF1=2,HMGA2=2,STC2=1,GHR=1,EAR=2")
    src.add_argument("--vcf",
                     help="a canFam4 (UU_Cfam_GSD_1.0) dog genome VCF; the pinned body-size + ear SNP dosages "
                          "are called by coordinate (dog causal SNPs have no rsIDs). Absent/uncallable SNPs are "
                          "skipped (partial-panel scoring)")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.dog_body_size import (
        MORPH_LOCI,
        SIZE_LOCI,
        SizeInputError,
        call_ear,
        polygenic_size_score,
    )

    try:
        if args.vcf:
            from dna_decode.pigment.dog_vcf_input import dosages_from_vcf
            dos = dosages_from_vcf(args.vcf)
            if not dos:
                print(f"error: no pinned morphology SNPs found/callable in {args.vcf} "
                      f"(expected canFam4 coords for IGF1/HMGA2/STC2/GHR/EAR)", file=sys.stderr)
                return 2
        else:
            dos = _parse_dosages(args.dosages)
        size_dos = {k: v for k, v in dos.items() if k in SIZE_LOCI}
        ear = call_ear(dos["EAR"]) if "EAR" in dos else None
        unknown = [k for k in dos if k not in SIZE_LOCI and k not in MORPH_LOCI]
        if unknown:
            raise SizeInputError(f"unknown locus/loci {unknown}; v0 loci: {list(SIZE_LOCI) + list(MORPH_LOCI)}")
        if not size_dos and ear is None:
            raise SizeInputError("supply at least one size locus (IGF1/HMGA2/STC2/GHR) or EAR")
        height = polygenic_size_score(size_dos).as_dict() if size_dos else None
    except FileNotFoundError as e:
        print(f"error: VCF not found: {e}", file=sys.stderr)
        return 2
    except (ValueError, SizeInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    result = {
        "organism": "Canis_lupus_familiaris", "trait": "morphology",
        "regime": "A_curated_pinned_validated_catalog",
        "input_source": "vcf" if args.vcf else "dosages",
        "loci_scored": sorted(dos),
        "height": height, "ear": ear, "abstains_on": _ABSTAIN_AXES,
        "measure": "RELATIVE size rank + ear axis (validated on Darwin's Ark), NOT calibrated absolute height",
        "scope_limit": "companion-animal visible-trait genetics; benign, NOT human/forensic",
    }
    if args.json:
        print(json.dumps(result, indent=2))
        return 0

    print("morphology decode — dog (body size + ear type)")
    if height:
        print(f"  body size: {height['size_rank'].upper()}   confidence: {height['confidence']}   "
              f"(polygenic score {height['polygenic_score']}/{height['max_score']}, "
              f"{height['n_loci_scored']}/4 size loci)")
    if ear:
        print(f"  ear type:  {ear['ear_type'].upper()}   confidence: {ear['confidence']}   "
              f"(MSRB3 dose {ear['high_allele_dose']}, r={ear['functional_r']})")
        print(f"  note: {ear['polarity_caveat']}")
    for ax in _ABSTAIN_AXES:
        print(f"  ABSTAIN: {ax}")
    print("  [pinned + Darwin's-Ark-validated catalog: height polygenic r=0.619, ear MSRB3 r=0.543]")
    print("  [RELATIVE rank not absolute inches; benign companion-animal visible-trait genetics]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
