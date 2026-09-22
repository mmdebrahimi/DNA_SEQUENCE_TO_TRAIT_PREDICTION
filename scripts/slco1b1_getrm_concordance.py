"""SLCO1B1 single-SNP surface vs GeT-RM CDC consensus -- testing this cell's own unmeasurability claim.

`dna_decode/pgx/slco1b1.py` states that "rs4149056 IS the truth for a 521T>C call, so 'validation' is
genotype-readout + trio-Mendelian consistency, NEVER an independent star-concordance number." That is a
claim of UNMEASURABILITY, and it is only HALF true -- which is what this script separates.

THE SPLIT THAT MAKES THIS HONEST (see wiki/slco1b1_star_acceptance_bar.json, frozen first)
------------------------------------------------------------------------------------------
  * FUNCTION axis (521 dosage -> Normal/Decreased/Poor): GeT-RM's own star assignment derives from
    genotypes INCLUDING 521T>C, so scoring our 521-derived function call against it is close to
    self-comparison. Reported as a CONTROL that bounds plumbing. A near-perfect number is EXPECTED and is
    NOT a finding. The cell's docstring is RIGHT about this axis.
  * STAR-NAME axis: *1B/*1A/*14/*21 are defined by 388A>G (and more), which this caller CANNOT see;
    *15 = 388+521 and *17 = *15 + promoter. So the printed star name depends on information absent from
    our input, and its disagreement rate with the reference IS measurable. The docstring is WRONG that no
    star-concordance number exists -- it exists and is the point of this script.
  * SAFETY axis: does any truth reduced-function carrier (*5/*15/*17, all 521C-bearing) get called
    Normal Function? That is the dangerous direction for statin myopathy.

Reuses the confusion/dosage helpers from the UGT1A1 harness so the math is shared and already pinned.
"""
from __future__ import annotations

import csv
import datetime
import json
import re
import sys
from pathlib import Path

from dna_decode.pgx import slco1b1 as sl
from scripts.ugt1a1_getrm_concordance import _dose, confusion, load_ped

REPO = Path(__file__).resolve().parent.parent
TRUTH_TSV = REPO / "data" / "pgx_getrm" / "getrm_consolidated_truth.tsv"
VCF = REPO / "data" / "pgx_1000g" / "slco1b1_1000g.vcf"
BAR = REPO / "wiki" / "slco1b1_star_acceptance_bar.json"

REDUCED_FUNCTION = ("*5", "*15", "*17")     # every 521C-bearing allele in this vocabulary
# alleles our caller structurally CANNOT identify (need 388A>G and/or promoter variants)
INVISIBLE_TO_CALLER = ("*1A", "*1B", "*14", "*15", "*17", "*21")


def load_truth() -> dict[str, str]:
    out: dict[str, str] = {}
    with TRUTH_TSV.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] == "SLCO1B1" and r["diplotype"].strip():
                out[r["coriell_id"]] = r["diplotype"].strip()
    return out


def haplotypes(diplotype: str) -> list[str] | None:
    """The two haplotypes' allele lists, or None if the string does not parse as two haplotypes."""
    parts = diplotype.split("/")
    if len(parts) != 2:
        return None
    out = []
    for hap in parts:
        alleles = re.findall(r"\*\w+", hap)
        if not alleles:
            return None
        out.append(alleles)
    return out


def reduced_count(diplotype: str) -> int | None:
    haps = haplotypes(diplotype)
    if haps is None:
        return None
    return sum(1 for h in haps if any(a in REDUCED_FUNCTION for a in h))


def carries_invisible_allele(diplotype: str) -> bool:
    haps = haplotypes(diplotype)
    if haps is None:
        return False
    return any(a in INVISIBLE_TO_CALLER for h in haps for a in h)


def norm_star(diplotype: str) -> tuple[str, ...] | None:
    """Order-insensitive normalised haplotype names, for exact-name comparison."""
    haps = haplotypes(diplotype)
    if haps is None:
        return None
    return tuple(sorted("+".join(sorted(h)) for h in haps))


def main(argv=None) -> int:
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    truth, pop = load_truth(), load_ped()

    dosages: dict[str, int] = {}
    with VCF.open(encoding="utf-8") as fh:
        order: list[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                order = f[9:]
                continue
            if int(f[1]) == sl.POS and f[3] == sl.REF and f[4] == sl.ALT:
                for name, gt in zip(order, f[9:]):
                    d = _dose(gt)
                    if d is not None:
                        dosages[name] = d
    if not dosages:
        print(f"ERROR: rs4149056 ({sl.CHROM}:{sl.POS} {sl.REF}>{sl.ALT}) not found in {VCF.name}",
              file=sys.stderr)
        return 2

    rows, parse_failures = [], []
    for sid in sorted(set(truth) & set(dosages)):
        dip = truth[sid]
        red = reduced_count(dip)
        if red is None:
            parse_failures.append({"sample": sid, "diplotype": dip})
            continue
        dose = dosages[sid]
        _gt, star_proxy, function, risk = sl._FUNCTION[dose]
        rows.append({
            "sample": sid, "superpopulation": pop.get(sid, "UNKNOWN"),
            "truth_diplotype": dip, "truth_reduced_n": red,
            "alt_dosage": dose, "called_star_proxy": star_proxy,
            "called_function": function, "myopathy_risk": risk,
            "star_name_agrees": norm_star(dip) == norm_star(star_proxy),
            "carries_allele_invisible_to_caller": carries_invisible_allele(dip),
        })

    if not rows:
        print("ERROR: no scoreable overlap", file=sys.stderr)
        return 2

    n_red_carriers = sum(1 for r in rows if r["truth_reduced_n"] > 0)
    n_invisible = sum(1 for r in rows if r["carries_allele_invisible_to_caller"])
    if n_red_carriers == 0:
        print("REFUSED: no reduced-function carrier in the overlap -- the safety axis cannot be tested "
              "against absent positives.", file=sys.stderr)
        return 3
    if n_invisible == 0:
        print("REFUSED: no *1B/*1A/*14/*21 carrier -- the star-name axis would be uninformative.",
              file=sys.stderr)
        return 3

    # --- the three axes ---
    star_agree = sum(1 for r in rows if r["star_name_agrees"])
    star_axis = {
        "n": len(rows), "exact_name_agreement": star_agree,
        "exact_name_agreement_frac": round(star_agree / len(rows), 4),
        "n_carrying_allele_invisible_to_caller": n_invisible,
    }

    func_control = confusion([(r["truth_reduced_n"], r["alt_dosage"]) for r in rows])

    safety_misses = [r for r in rows if r["truth_reduced_n"] > 0 and r["alt_dosage"] == 0]
    safety = {
        "n_truth_reduced_function_carriers": n_red_carriers,
        "n_called_normal_function_anyway": len(safety_misses),
        "misses": [{"sample": r["sample"], "truth": r["truth_diplotype"]} for r in safety_misses],
    }

    # per-allele disagreement attribution (which reference allele drives the naming gap)
    by_allele: dict[str, dict] = {}
    for r in rows:
        haps = haplotypes(r["truth_diplotype"]) or []
        for a in {a for h in haps for a in h}:
            b = by_allele.setdefault(a, {"n": 0, "name_agrees": 0})
            b["n"] += 1
            b["name_agrees"] += 1 if r["star_name_agrees"] else 0

    frac = star_axis["exact_name_agreement_frac"]
    if safety["n_called_normal_function_anyway"] > 0:
        verdict = "DEMOTE_AND_DISCLOSE_SAFETY_MISS"
    elif frac < 0.60:
        verdict = "STAR_PROXY_NOT_REFERENCE_AGREEING"
    elif frac >= 0.90:
        verdict = "STAR_PROXY_AGREES"
    else:
        verdict = "STAR_PROXY_PARTIAL"

    out = {
        "schema": "slco1b1-getrm-concordance-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "cell": "pgx:human:slco1b1",
        "acceptance_bar": "wiki/slco1b1_star_acceptance_bar.json",
        "truth_source": bar["truth"]["source"],
        "independence": bar["truth"]["independence"],
        "n_scored": len(rows),
        "n_parse_failures": len(parse_failures),
        "parse_failures": parse_failures,
        "star_name_axis": star_axis,
        "function_axis_control": {
            **func_control,
            "IS_A_CONTROL_NOT_VALIDATION": True,
            "why": bar["the_circularity_split_that_makes_this_honest"]["function_axis_is_LARGELY_DEFINITIONAL"],
        },
        "safety_axis": safety,
        "per_reference_allele": by_allele,
        "verdict": verdict,
        "registered_predictions": bar["registered_predictions"],
        "honest_limits": bar["honest_limits_registered_in_advance"],
        "rows": rows,
    }

    stamp = datetime.date.today().isoformat()
    jf = REPO / "wiki" / f"slco1b1_getrm_concordance_{stamp}.json"
    jf.write_text(json.dumps(out, indent=2), encoding="utf-8")

    lines = [
        f"# SLCO1B1 single-SNP surface vs GeT-RM CDC consensus ({stamp})", "",
        f"Scored **{len(rows)}** 1000G samples with CDC multi-lab consensus SLCO1B1 diplotypes. Bar frozen "
        f"beforehand: `wiki/slco1b1_star_acceptance_bar.json`.", "",
        "The cell's docstring says validation is *\"never an independent star-concordance number\"*. That is "
        "half right, and the two halves are reported separately below.", "",
        "## 1. Star-NAME axis (genuinely testable -- the caller cannot see 388A>G)", "",
        f"- exact star-name agreement **{star_agree}/{len(rows)} = {frac}**",
        f"- {n_invisible}/{len(rows)} samples carry at least one allele the caller structurally cannot "
        f"identify (`*1A`/`*1B`/`*14`/`*15`/`*17`/`*21`)",
        f"- **verdict: {verdict}**", "",
        "## 2. FUNCTION axis (a CONTROL, not validation)", "",
        "GeT-RM's star assignment itself derives from genotypes including 521T>C, so this is close to "
        "self-comparison. A near-perfect number here is **expected and is not evidence**.", "",
        f"- reduced-function carriership: sens **{func_control['sensitivity']}** · spec "
        f"**{func_control['specificity']}** · PPV **{func_control['ppv']}** (n={func_control['n']}), "
        f"dosage-exact {func_control['zygosity_exact']}/{func_control['n']}", "",
        "## 3. SAFETY axis (the clinical one)", "",
        f"- truth reduced-function carriers: **{n_red_carriers}**",
        f"- of those, called Normal Function anyway: **{safety['n_called_normal_function_anyway']}**", "",
        "## Which reference allele drives the naming gap", "",
        "| allele | samples | name agrees |", "|---|---|---|",
    ]
    for a, b in sorted(by_allele.items(), key=lambda kv: -kv[1]["n"]):
        lines.append(f"| `{a}` | {b['n']} | {b['name_agrees']} |")
    lines += ["", "## Honest limits", ""] + [f"- {h}" for h in out["honest_limits"]]
    if parse_failures:
        lines.append(f"- {len(parse_failures)} truth diplotype(s) did not parse and were EXCLUDED.")
    mf = REPO / "wiki" / f"slco1b1_getrm_concordance_{stamp}.md"
    mf.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"scored={len(rows)}  reduced-fn carriers={n_red_carriers}  "
          f"carrying-invisible-allele={n_invisible}  parse_failures={len(parse_failures)}")
    print(f"STAR-NAME  exact agreement={star_agree}/{len(rows)} = {frac}")
    print(f"FUNCTION   (control) sens={func_control['sensitivity']} spec={func_control['specificity']} "
          f"exact={func_control['zygosity_exact']}/{func_control['n']}")
    print(f"SAFETY     reduced-fn carriers called Normal = {safety['n_called_normal_function_anyway']}")
    print(f"VERDICT    {verdict}")
    print(f"wrote {jf.name} + {mf.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
