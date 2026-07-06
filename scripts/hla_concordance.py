"""HLA tag-SNP validation — the wrapper-vs-truth number for the `dna-hla` cell.

Two validations, honest about what each proves:
  * POPULATION-AF concordance (default; data-in-hand): the tag SNP's per-superpopulation ALT frequency vs
    the KNOWN population frequency of the HLA allele. A PARTIAL corroboration — the tag tracks the allele's
    population structure (e.g. B*57:01 tag high in SAS/EUR, low in EAS). Necessary-not-sufficient.
  * SAMPLE-LEVEL concordance (--truth <sample<TAB>carrier(0/1)>): the real wrapper-vs-truth number — tag
    call vs a free published HLA truth (e.g. Gourraud 2014 1000G HLA types). Sens/spec/PPV of the tag proxy.
    This is the SCORED validation; until a truth TSV is supplied it is HONESTLY reported as pending (an
    external data-acquisition step, NOT a code wall).
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
from dna_decode.hla.catalog import get  # noqa: E402

PED = REPO / "data" / "pgx_1000g" / "1000G_3202_samples.ped"
# Published HLA-allele ALLELE frequency by 1000G superpopulation (order-of-magnitude, for the AF corroboration).
# B*57:01: high South-Asian/European, low East-Asian (Gonzalez-Galarza AFND; the abacavir-screen literature).
_KNOWN_ALLELE_FREQ = {
    "b5701": {"SAS": "high (~0.04-0.07)", "EUR": "moderate (~0.03-0.04)", "AMR": "low-moderate",
              "AFR": "low (~0.01-0.02)", "EAS": "very low (~0.005-0.01)"},
}


def _norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def _superpop() -> dict[str, str]:
    m = {}
    for line in PED.read_text(encoding="utf-8").splitlines():
        f = line.split()
        if len(f) < 7 or f[0] == "FamilyID":
            continue
        m[f[1]] = f[6]
    return m


def population_af(vcf: Path, allele_key: str) -> dict:
    a = get(allele_key)
    sp = _superpop()
    samples: list[str] = []
    per = {}
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if line.startswith("##"):
            continue
        cols = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            samples = cols[9:]
            continue
        if len(cols) < 8 or not cols[1].isdigit():
            continue
        if _norm_chrom(cols[0]) != a.chrom or int(cols[1]) != a.pos:
            continue
        alts = cols[4].split(",")
        ai = alts.index(a.tag_alt) + 1 if a.tag_alt in alts else -1
        if ai < 0:
            break
        for i, s in enumerate(cols[9:]):
            pop = sp.get(samples[i], "NA")
            d = per.setdefault(pop, [0, 0])
            for x in s.split(":")[0].replace("|", "/").split("/"):
                if x.isdigit():
                    d[1] += 1
                    if int(x) == ai:
                        d[0] += 1
        break
    return {pop: {"tag_alt_count": c[0], "n_hap": c[1], "tag_af": round(c[0] / c[1], 4) if c[1] else None}
            for pop, c in sorted(per.items())}


def sample_concordance(vcf: Path, allele_key: str, truth_tsv: Path) -> dict:
    from dna_decode.hla.caller import call_hla
    truth = {}
    for r in csv.DictReader(truth_tsv.open(encoding="utf-8"), delimiter="\t"):
        truth[r["sample"].strip()] = int(r["carrier"])
    # sample columns in the VCF
    samples = []
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            break
    tp = fp = tn = fn = 0
    for s in samples:
        if s not in truth:
            continue
        pred = 1 if call_hla(vcf, allele_key, sample=s)["carrier"] else 0
        t = truth[s]
        tp += pred == 1 and t == 1
        fp += pred == 1 and t == 0
        tn += pred == 0 and t == 0
        fn += pred == 0 and t == 1
    n = tp + fp + tn + fn
    return {"n_scored": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sensitivity": round(tp / (tp + fn), 4) if (tp + fn) else None,
            "specificity": round(tn / (tn + fp), 4) if (tn + fp) else None,
            "ppv": round(tp / (tp + fp), 4) if (tp + fp) else None}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="HLA tag-SNP validation (population-AF + sample-level concordance).")
    ap.add_argument("--allele", default="b5701")
    ap.add_argument("--vcf", type=Path, default=REPO / "data" / "pgx_1000g" / "hla_b5701_abacavir.vcf")
    ap.add_argument("--truth", type=Path, default=None, help="TSV sample<TAB>carrier(0/1) -> sample concordance")
    args = ap.parse_args(argv)
    a = get(args.allele)
    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf} (fetch the tag region first)")
        return 2

    rep = {
        "schema": "hla-tag-validation-v0", "allele": a.allele, "allele_key": a.key,
        "tag": f"{a.rsid} chr{a.chrom}:{a.pos} {a.ref}>{a.tag_alt}", "drug": a.drug,
        "analysis_date": datetime.date.today().isoformat(), "proxy_tier": a.proxy_tier,
        "population_af": population_af(args.vcf, args.allele),
        "known_allele_freq": _KNOWN_ALLELE_FREQ.get(a.key, "not tabulated"),
        "af_corroboration": ("PARTIAL — the tag's per-superpopulation frequency tracking the allele's known "
                             "population structure is necessary-not-sufficient. The SCORED number is the "
                             "sample-level concordance vs a real HLA truth set (--truth)."),
        "sample_concordance": None,
        "sample_concordance_status": ("PENDING — supply a free published 1000G HLA truth TSV via --truth "
                                      "(e.g. Gourraud 2014); this is an external data-acquisition step, NOT a "
                                      "code wall. The harness computes sens/spec/PPV of the tag proxy."),
    }
    if args.truth and args.truth.exists():
        rep["sample_concordance"] = sample_concordance(args.vcf, args.allele, args.truth)
        rep["sample_concordance_status"] = "SCORED"

    stem = f"hla_{a.key}_validation_{rep['analysis_date']}"
    (REPO / "wiki" / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")
    print(f"# HLA {a.allele} tag validation ({rep['analysis_date']})  drug={a.drug}  tag={rep['tag']}")
    print("per-superpopulation tag AF:")
    for pop, d in rep["population_af"].items():
        print(f"  {pop:4} tag_af={d['tag_af']}  ({d['tag_alt_count']}/{d['n_hap']})")
    print(f"known {a.allele} population freq: {rep['known_allele_freq']}")
    print(f"AF corroboration: {rep['af_corroboration']}")
    print(f"sample concordance: {rep['sample_concordance'] or rep['sample_concordance_status']}")
    print(f"[report -> wiki/{stem}.json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
