"""CYP4F2 single-SNP surface vs GeT-RM CDC consensus + a real AF check of its own corroboration claim.

`dna_decode/pgx/cyp4f2.py` carries the SAME unmeasurability claim VERBATIM as slco1b1.py did -- that
validation is "genotype-readout + trio-Mendelian consistency + AF-corroboration, never an independent
star-concordance number." SLCO1B1's version was measured 2026-09-22 and found half wrong. This is the last
GeT-RM-scoreable unscored PGx cell, and it puts BOTH halves of that sentence on the record:

  1. STAR-NAME axis -- *2 is defined by a variant this caller does not genotype (it reads exactly one SNP,
     rs2108622), so the printed label uses information absent from our input and its disagreement rate with
     the reference is measurable. Bar frozen first, with a PREDICTED band of ~0.81 (52/64) that is
     deliberately NOT the 0.398 band SLCO1B1 landed in -- a result near either extreme falsifies it.
  2. AF-CORROBORATION axis -- the docstring asserts the *3 (T) allele is "~29% EUR / ~79% EAS" and the
     registry contract cites "AF-corroborated" as a validation leg. That AF has never been checked against
     the 3202-sample 1000G panel this cell ships beside. Measured here over the FULL panel (the right
     denominator for a population-frequency claim, and a DIFFERENT denominator from every other number in
     this script -- which is why it is reported in its own block).

FUNCTION axis is reported as a CONTROL, not validation: GeT-RM's own *3 assignment derives from rs2108622,
the one SNP we genotype, so a near-perfect number there is definitional and is NOT evidence.

SCOPE (frozen in the bar, repeated because it is easy to over-read): *2's FUNCTIONAL status is NOT settled
by this run. CPIC's warfarin algorithm dose-adjusts on *3; whether *2 carries an independent effect is a
separate literature question. A *2 carrier called *1 is therefore a NAMING finding here, and would only be
a SAFETY finding if *2 were reduced-function -- which this measurement does not establish and must not
assert.
"""
from __future__ import annotations

import csv
import datetime
import json
import re
import sys
from collections import Counter
from pathlib import Path

from dna_decode.pgx import cyp4f2 as c4
from scripts.slco1b1_getrm_concordance import haplotypes, norm_star
from scripts.ugt1a1_getrm_concordance import _dose, confusion, load_ped

REPO = Path(__file__).resolve().parent.parent
TRUTH_TSV = REPO / "data" / "pgx_getrm" / "getrm_consolidated_truth.tsv"
VCF = REPO / "data" / "pgx_1000g" / "cyp4f2_1000g.vcf"
BAR = REPO / "wiki" / "cyp4f2_star_acceptance_bar.json"

# the caller reads exactly ONE SNP (rs2108622 -> *3). *2 is defined by a different variant it never sees.
INVISIBLE_TO_CALLER = ("*2",)
ASSERTED_AF = {"EUR": 0.29, "EAS": 0.79}     # the docstring's claim, under test
AF_TOLERANCE = 0.10                          # frozen in the bar


def load_truth() -> dict[str, str]:
    out: dict[str, str] = {}
    with TRUTH_TSV.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] == "CYP4F2" and r["diplotype"].strip():
                out[r["coriell_id"]] = r["diplotype"].strip()
    return out


def star3_count(diplotype: str) -> int | None:
    """Haplotypes carrying *3. None if the string does not parse as two haplotypes."""
    haps = haplotypes(diplotype)
    if haps is None:
        return None
    return sum(1 for h in haps if "*3" in h)


def carries_invisible_allele(diplotype: str) -> bool:
    haps = haplotypes(diplotype)
    if haps is None:
        return False
    return any(a in INVISIBLE_TO_CALLER for h in haps for a in h)


def read_dosages() -> tuple[dict[str, int], list[str]]:
    """ALT (T) dosage at rs2108622 for every sample in the panel."""
    dosages: dict[str, int] = {}
    order: list[str] = []
    with VCF.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                order = f[9:]
                continue
            if int(f[1]) == c4.POS and f[3] == c4.REF and f[4] == c4.ALT:
                for name, gt in zip(order, f[9:]):
                    d = _dose(gt)
                    if d is not None:
                        dosages[name] = d
    return dosages, order


def allele_frequencies(dosages: dict[str, int], pop: dict[str, str]) -> dict:
    """ALT allele frequency per superpopulation over the FULL panel (not the truth overlap)."""
    alt = Counter()
    tot = Counter()
    for s, d in dosages.items():
        p = pop.get(s)
        if not p:
            continue
        alt[p] += d
        tot[p] += 2
    out = {p: {"n_samples": tot[p] // 2, "alt_alleles": alt[p],
               "alt_freq": round(alt[p] / tot[p], 4) if tot[p] else None}
           for p in sorted(tot)}
    all_alt, all_tot = sum(alt.values()), sum(tot.values())
    out["ALL"] = {"n_samples": all_tot // 2, "alt_alleles": all_alt,
                  "alt_freq": round(all_alt / all_tot, 4) if all_tot else None}
    return out


def main(argv=None) -> int:
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    truth, pop = load_truth(), load_ped()
    dosages, _order = read_dosages()
    if not dosages:
        print(f"ERROR: rs2108622 ({c4.CHROM}:{c4.POS} {c4.REF}>{c4.ALT}) not found in {VCF.name}",
              file=sys.stderr)
        return 2

    rows, parse_failures = [], []
    for sid in sorted(set(truth) & set(dosages)):
        dip = truth[sid]
        n3 = star3_count(dip)
        if n3 is None:
            parse_failures.append({"sample": sid, "diplotype": dip})
            continue
        dose = dosages[sid]
        _gt, star_proxy, function, dose_dir = c4._FUNCTION[dose]
        rows.append({
            "sample": sid, "superpopulation": pop.get(sid, "UNKNOWN"),
            "truth_diplotype": dip, "truth_star3_n": n3,
            "alt_dosage": dose, "called_star_proxy": star_proxy,
            "called_function": function, "warfarin_dose_direction": dose_dir,
            "star_name_agrees": norm_star(dip) == norm_star(star_proxy),
            "carries_allele_invisible_to_caller": carries_invisible_allele(dip),
        })

    if not rows:
        print("ERROR: no scoreable overlap", file=sys.stderr)
        return 2

    n_star3_carriers = sum(1 for r in rows if r["truth_star3_n"] > 0)
    n_invisible = sum(1 for r in rows if r["carries_allele_invisible_to_caller"])
    if n_star3_carriers == 0:
        print("REFUSED: no *3 carrier in the overlap -- the safety axis cannot be tested against absent "
              "positives.", file=sys.stderr)
        return 3
    if n_invisible == 0:
        print("REFUSED: no *2 carrier -- the star-name axis would be uninformative.", file=sys.stderr)
        return 3

    star_agree = sum(1 for r in rows if r["star_name_agrees"])
    frac = round(star_agree / len(rows), 4)
    star_axis = {"n": len(rows), "exact_name_agreement": star_agree,
                 "exact_name_agreement_frac": frac,
                 "n_carrying_allele_invisible_to_caller": n_invisible}

    func_control = confusion([(r["truth_star3_n"], r["alt_dosage"]) for r in rows])

    safety_misses = [r for r in rows if r["truth_star3_n"] > 0 and r["alt_dosage"] == 0]
    safety = {
        "n_truth_star3_carriers": n_star3_carriers,
        "n_called_normal_function_anyway": len(safety_misses),
        "misses": [{"sample": r["sample"], "truth": r["truth_diplotype"]} for r in safety_misses],
        "scope": bar["scope_limits_registered_in_advance"][0],
    }

    # --- the AF-corroboration claim, on the FULL panel ---
    af = allele_frequencies(dosages, pop)
    af_checks = {}
    for p, asserted in ASSERTED_AF.items():
        measured = af.get(p, {}).get("alt_freq")
        ok = measured is not None and abs(measured - asserted) <= AF_TOLERANCE
        af_checks[p] = {"asserted": asserted, "measured": measured,
                        "delta": None if measured is None else round(measured - asserted, 4),
                        "within_tolerance": ok}
    af_claim_holds = all(c["within_tolerance"] for c in af_checks.values())

    by_allele: dict[str, dict] = {}
    for r in rows:
        haps = haplotypes(r["truth_diplotype"]) or []
        for a in {a for h in haps for a in h}:
            b = by_allele.setdefault(a, {"n": 0, "name_agrees": 0})
            b["n"] += 1
            b["name_agrees"] += 1 if r["star_name_agrees"] else 0

    if safety["n_called_normal_function_anyway"] > 0:
        star_verdict = "DEMOTE_AND_DISCLOSE_SAFETY_MISS"
    elif frac < 0.60:
        star_verdict = "STAR_PROXY_NOT_REFERENCE_AGREEING"
    elif frac >= 0.90:
        star_verdict = "STAR_PROXY_AGREES"
    else:
        star_verdict = "STAR_PROXY_PARTIAL"

    pred = bar["registered_predictions"]
    prediction_held = (star_verdict == "STAR_PROXY_PARTIAL"
                       and len(rows) - star_agree == pred["predicted_disagreement_count"])

    out = {
        "schema": "cyp4f2-getrm-concordance-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "cell": "pgx:human:cyp4f2",
        "acceptance_bar": "wiki/cyp4f2_star_acceptance_bar.json",
        "truth_source": bar["truth"]["source"],
        "independence": bar["truth"]["independence"],
        "n_scored": len(rows),
        "n_parse_failures": len(parse_failures),
        "parse_failures": parse_failures,
        "star_name_axis": star_axis,
        "function_axis_control": {
            **func_control, "IS_A_CONTROL_NOT_VALIDATION": True,
            "why": "GeT-RM's own *3 assignment derives from rs2108622, the one SNP this caller reads, so "
                   "a near-perfect number here is definitional and is not evidence.",
        },
        "safety_axis": safety,
        "af_corroboration_axis": {
            "denominator_note": "computed over the FULL 3202-sample 1000G panel, NOT the 64-sample truth "
                                "overlap -- the right denominator for a population-frequency claim and a "
                                "DIFFERENT one from every other number in this artifact",
            "asserted_by_docstring": ASSERTED_AF,
            "tolerance": AF_TOLERANCE,
            "measured": af,
            "checks": af_checks,
            "claim_holds": af_claim_holds,
        },
        "per_reference_allele": by_allele,
        "star_verdict": star_verdict,
        "af_verdict": "AF_CLAIM_HOLDS" if af_claim_holds else "AF_CLAIM_WRONG",
        "registered_prediction_held": prediction_held,
        "registered_predictions": pred,
        "scope_limits": bar["scope_limits_registered_in_advance"],
        "rows": rows,
    }

    stamp = datetime.date.today().isoformat()
    jf = REPO / "wiki" / f"cyp4f2_getrm_concordance_{stamp}.json"
    jf.write_text(json.dumps(out, indent=2), encoding="utf-8")

    lines = [
        f"# CYP4F2 single-SNP surface vs GeT-RM CDC consensus ({stamp})", "",
        f"Scored **{len(rows)}** 1000G samples with CDC multi-lab consensus CYP4F2 diplotypes. Bar and "
        f"predictions frozen beforehand: `wiki/cyp4f2_star_acceptance_bar.json`.", "",
        "This cell carries the SAME *\"never an independent star-concordance number\"* claim as SLCO1B1, "
        "plus an **AF-corroboration** claim cited as a validation leg. Both are tested.", "",
        "## 1. Star-NAME axis", "",
        f"- exact star-name agreement **{star_agree}/{len(rows)} = {frac}**",
        f"- {n_invisible}/{len(rows)} samples carry `*2`, which this caller cannot see (it reads exactly "
        f"one SNP, rs2108622)",
        f"- **verdict: {star_verdict}** · registered prediction (~0.81, 12 disagreements) "
        f"**{'HELD' if prediction_held else 'DID NOT HOLD'}**", "",
        "## 2. FUNCTION axis (a CONTROL, not validation)", "",
        f"- `*3` carriership: sens **{func_control['sensitivity']}** · spec "
        f"**{func_control['specificity']}** · PPV **{func_control['ppv']}** (n={func_control['n']}), "
        f"dosage-exact {func_control['zygosity_exact']}/{func_control['n']}",
        "- definitional: GeT-RM's `*3` assignment derives from the one SNP we read.", "",
        "## 3. SAFETY axis", "",
        f"- truth `*3` carriers: **{n_star3_carriers}**; called Normal Function anyway: "
        f"**{safety['n_called_normal_function_anyway']}**",
        f"- *Scope:* {safety['scope']}", "",
        "## 4. AF-corroboration axis (the second shipped claim)", "",
        f"Computed over the **full 3202-sample panel**, not the 64-sample overlap.", "",
        "| pop | asserted | measured | delta | within +/-0.10 |", "|---|---|---|---|---|",
    ]
    for p, c in af_checks.items():
        lines.append(f"| {p} | {c['asserted']} | {c['measured']} | {c['delta']} | "
                     f"{'yes' if c['within_tolerance'] else '**NO**'} |")
    lines += ["", f"- **verdict: {out['af_verdict']}**", "",
              "All measured superpopulation frequencies:", "",
              "| pop | n | ALT alleles | freq |", "|---|---|---|---|"]
    for p, v in af.items():
        lines.append(f"| {p} | {v['n_samples']} | {v['alt_alleles']} | {v['alt_freq']} |")
    lines += ["", "## Which reference allele drives the naming gap", "",
              "| allele | samples | name agrees |", "|---|---|---|"]
    for a, b in sorted(by_allele.items(), key=lambda kv: -kv[1]["n"]):
        lines.append(f"| `{a}` | {b['n']} | {b['name_agrees']} |")
    lines += ["", "## Scope limits", ""] + [f"- {h}" for h in out["scope_limits"]]
    if parse_failures:
        lines.append(f"- {len(parse_failures)} truth diplotype(s) did not parse and were EXCLUDED.")
    mf = REPO / "wiki" / f"cyp4f2_getrm_concordance_{stamp}.md"
    mf.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"scored={len(rows)}  *3 carriers={n_star3_carriers}  *2 carriers={n_invisible}  "
          f"parse_failures={len(parse_failures)}")
    print(f"STAR-NAME  {star_agree}/{len(rows)} = {frac}  -> {star_verdict} "
          f"(prediction {'HELD' if prediction_held else 'DID NOT HOLD'})")
    print(f"FUNCTION   (control) exact={func_control['zygosity_exact']}/{func_control['n']}")
    print(f"SAFETY     *3 carriers called Normal = {safety['n_called_normal_function_anyway']}")
    for p, c in af_checks.items():
        print(f"AF {p}      asserted={c['asserted']} measured={c['measured']} "
              f"delta={c['delta']} -> {'ok' if c['within_tolerance'] else 'OUT OF BAND'}")
    print(f"AF VERDICT {out['af_verdict']}")
    print(f"wrote {jf.name} + {mf.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
