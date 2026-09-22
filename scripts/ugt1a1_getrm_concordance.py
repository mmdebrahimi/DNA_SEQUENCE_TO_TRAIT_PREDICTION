"""UGT1A1 tag-SNP surface vs GeT-RM CDC consensus on real 1000G genomes -- and a head-to-head
DIRECT TA-repeat arm that tests the catalog's own impossibility claim.

WHY THIS EXISTS
---------------
`dna_decode/pgx/ugt1a1_catalog.py` ships at NEAR_INDEPENDENT while never having been scored against
GeT-RM. Its major reduced-function allele *28 is a promoter TA-dinucleotide repeat, and v0 substitutes
rs887829 (*80) as an LD-tag PROXY, asserting "EUR r^2 ~0.9+" -- a number never measured here. This repo
has a demotion precedent for precisely that pattern: the HLA B*58:01 tag reached sens 0.61 / PPV 0.18 on
1000G truth and was declared unsafe for an SJS/TEN screen.

The acceptance bar was FROZEN BEFORE any result was read: wiki/ugt1a1_tag_acceptance_bar.json.

WHAT IS AND IS NOT THE METRIC
-----------------------------
GeT-RM truth speaks in repeat alleles (*28/*36/*37/*60/*27/*6/*7) the shipped SNP caller structurally
cannot emit, so EXACT-DIPLOTYPE concordance is not well-posed and is deliberately not the bar. The bar is
*28 CARRIER DETECTION (+ zygosity, which is what the CPIC activity score actually turns on).

THE DIRECT ARM
--------------
The catalog states the TA repeat is "unresolvable from a short-read SNP VCF". The 1000G 30x phased panel
represents it at chr2:233760233 as indel alleles: C>CAT (+1 TA = TA7 = *28), C>CATAT (+2 TA = TA8 = *37),
CAT>C (-1 TA = TA5 = *36). That makes the asserted impossibility TESTABLE rather than assumed, so this
script scores a direct repeat call head-to-head against the tag on the same samples and same truth.

A win for the direct arm is PANEL-SPECIFIC and must be reported that way (see the bar's scope_limit).

TRUTH PROVENANCE (exact, so nobody assumes a merge that did not happen)
----------------------------------------------------------------------
The truth TSV is built from the CDC CONSOLIDATED table ALONE:
    uv run python scripts/getrm_consolidated_ingest.py \
        --consolidated data/pgx_getrm/getrm_consolidated_pgx_hla.xlsx
NOT merged with the 137-sample panel table (that xlsx is not on this host). For UGT1A1 specifically the
two agree: the committed coverage artifact records UGT1A1 = 197 and the consolidated-only build yields
the same 197 rows, so this cell's truth set is unaffected by the omission. Do NOT re-run that ingest
WITHOUT --panel137 and commit the result -- it overwrites wiki/pgx_getrm_consolidated_coverage_*.json
with a consolidated-only version, silently dropping its panel137 block.
"""
from __future__ import annotations

import collections
import csv
import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TRUTH_TSV = REPO / "data" / "pgx_getrm" / "getrm_consolidated_truth.tsv"
VCF = REPO / "data" / "pgx_1000g" / "ugt1a1_1000g.vcf"
PED = REPO / "data" / "pgx_1000g" / "1000G_3202_samples.ped"
BAR = REPO / "wiki" / "ugt1a1_tag_acceptance_bar.json"

# --- loci (GRCh38, chr2 plus strand; coords verified in ugt1a1_catalog.py) ---
TAG_POS = 233759924          # rs887829 C>T  == *80, the shipped LD-tag proxy for *28
STAR6_POS = 233760498        # rs4148323 G>A == *6 (direct SNP; like-for-like control)
REPEAT_POS = 233760233       # UGT1A1 promoter TA repeat, represented as indels in this panel

# TA-repeat indel alleles -> star allele. Reference is A(TA)6TAA == *1.
REPEAT_ALLELE_TO_STAR = {
    ("C", "CAT"): "*28",      # +1 TA -> TA7
    ("C", "CATAT"): "*37",    # +2 TA -> TA8
    ("CAT", "C"): "*36",      # -1 TA -> TA5
}


def load_truth() -> dict[str, str]:
    """{coriell_id: raw UGT1A1 diplotype} from the CDC consolidated consensus."""
    out: dict[str, str] = {}
    with TRUTH_TSV.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] == "UGT1A1" and r["diplotype"].strip():
                out[r["coriell_id"]] = r["diplotype"].strip()
    return out


def truth_star_count(diplotype: str, star: str) -> int | None:
    """Count haplotypes carrying `star`. Returns None if the diplotype does not parse as two haplotypes.

    GeT-RM writes compound haplotypes with parentheses/plus, e.g. `(*28 + *60)/(*36 + *60)` or `*28/(*60)`.
    Splitting on the single top-level '/' yields the two haplotypes; each may list several alleles.
    A parse failure is RETURNED AS None and excluded upstream -- never silently defaulted to 0, which
    would score an unparsed row as a true negative.
    """
    parts = diplotype.split("/")
    if len(parts) != 2:
        return None
    n = 0
    for hap in parts:
        alleles = re.findall(r"\*\d+", hap)
        if not alleles:
            return None
        if star in alleles:
            n += 1
    return n


def truth_reduced_count(diplotype: str) -> int | None:
    """Count haplotypes carrying a REDUCED-FUNCTION promoter repeat allele (*28 = TA7 or *37 = TA8).

    POST-HOC framing, added after the frozen bar was scored and labelled as such everywhere it is
    reported. Rationale: the diagnosis showed rs887829-T tags the reduced-function promoter CLASS, and
    the clinical decision (irinotecan dosing) turns on reduced function, not on *28 specifically.
    *36 (TA5) is deliberately EXCLUDED -- a shorter repeat is INCREASED function, not reduced.
    """
    parts = diplotype.split("/")
    if len(parts) != 2:
        return None
    n = 0
    for hap in parts:
        alleles = re.findall(r"\*\d+", hap)
        if not alleles:
            return None
        if "*28" in alleles or "*37" in alleles:
            n += 1
    return n


def load_ped() -> dict[str, str]:
    """{sample: superpopulation}."""
    pop: dict[str, str] = {}
    with PED.open(encoding="utf-8") as fh:
        hdr = fh.readline().split()
        si, pi = hdr.index("SampleID"), hdr.index("Superpopulation")
        for line in fh:
            f = line.split()
            if len(f) > pi:
                pop[f[si]] = f[pi]
    return pop


def _dose(gt: str) -> int | None:
    """ALT dosage from a VCF genotype field. None when uncalled."""
    g = gt.split(":")[0].replace("|", "/")
    if "." in g:
        return None
    return sum(1 for a in g.split("/") if a != "0")


def read_vcf(samples: set[str]) -> tuple[dict[str, dict], list[str]]:
    """Per-sample ALT dosages at the tag, *6, and each TA-repeat indel allele."""
    calls: dict[str, dict] = {s: {} for s in samples}
    order: list[str] = []
    with VCF.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                order = f[9:]
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            key = None
            if pos == TAG_POS and (ref, alt) == ("C", "T"):
                key = "tag80"
            elif pos == STAR6_POS and (ref, alt) == ("G", "A"):
                key = "star6"
            elif pos == REPEAT_POS and (ref, alt) in REPEAT_ALLELE_TO_STAR:
                key = "repeat" + REPEAT_ALLELE_TO_STAR[(ref, alt)]
            if key is None:
                continue
            for name, gt in zip(order, f[9:]):
                if name in calls:
                    d = _dose(gt)
                    if d is not None:
                        calls[name][key] = d
    return calls, order


def confusion(pairs: list[tuple[int, int]]) -> dict:
    """2x2 on CARRIERSHIP (count>0) plus exact-zygosity agreement, from (truth_n, pred_n) pairs."""
    tp = sum(1 for t, p in pairs if t > 0 and p > 0)
    fn = sum(1 for t, p in pairs if t > 0 and p == 0)
    fp = sum(1 for t, p in pairs if t == 0 and p > 0)
    tn = sum(1 for t, p in pairs if t == 0 and p == 0)
    exact = sum(1 for t, p in pairs if t == p)

    def rate(a, b):
        return round(a / (a + b), 4) if (a + b) else None

    return {
        "n": len(pairs), "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sensitivity": rate(tp, fn), "specificity": rate(tn, fp),
        "ppv": rate(tp, fp), "npv": rate(tn, fn),
        "zygosity_exact": exact,
        "zygosity_exact_frac": round(exact / len(pairs), 4) if pairs else None,
    }


def main(argv=None) -> int:
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    truth = load_truth()
    pop = load_ped()

    with VCF.open(encoding="utf-8") as fh:
        vcf_samples = set()
        for line in fh:
            if line.startswith("#CHROM"):
                vcf_samples = set(line.rstrip("\n").split("\t")[9:])
                break

    scored_ids = sorted(set(truth) & vcf_samples)
    if not scored_ids:
        print("ERROR: no overlap between GeT-RM truth and the 1000G panel", file=sys.stderr)
        return 2

    calls, _ = read_vcf(set(scored_ids))

    rows, parse_failures = [], []
    for sid in scored_ids:
        dip = truth[sid]
        t28 = truth_star_count(dip, "*28")
        t6 = truth_star_count(dip, "*6")
        if t28 is None:
            parse_failures.append({"sample": sid, "diplotype": dip})
            continue
        c = calls.get(sid, {})
        direct28 = c.get("repeat*28")
        rows.append({
            "sample": sid, "superpopulation": pop.get(sid, "UNKNOWN"), "truth_diplotype": dip,
            "truth_star28_n": t28, "truth_star6_n": t6,
            "truth_reduced_fn_n": truth_reduced_count(dip),
            "tag_star80_n": c.get("tag80"), "direct_star28_n": direct28,
            "direct_star36_n": c.get("repeat*36"), "direct_star37_n": c.get("repeat*37"),
            "called_star6_n": c.get("star6"),
        })

    n_carriers = sum(1 for r in rows if r["truth_star28_n"] > 0)
    if n_carriers == 0:
        print("REFUSED: truth overlap contains zero *28 carriers -- a tag cannot be tested against "
              "absent positives.", file=sys.stderr)
        return 3

    tag_pairs = [(r["truth_star28_n"], r["tag_star80_n"]) for r in rows if r["tag_star80_n"] is not None]
    dir_pairs = [(r["truth_star28_n"], r["direct_star28_n"]) for r in rows
                 if r["direct_star28_n"] is not None]
    six_pairs = [(r["truth_star6_n"], r["called_star6_n"]) for r in rows
                 if r["truth_star6_n"] is not None and r["called_star6_n"] is not None]

    # POST-HOC arm (NOT in the frozen bar): the clinical axis is reduced FUNCTION, and the diagnosis
    # showed rs887829-T tags the reduced-function promoter class (*28 TA7 + *37 TA8), not *28 alone.
    red_pairs = [(r["truth_reduced_fn_n"], r["tag_star80_n"]) for r in rows
                 if r["truth_reduced_fn_n"] is not None and r["tag_star80_n"] is not None]

    tag = confusion(tag_pairs)
    direct = confusion(dir_pairs) if dir_pairs else {"n": 0, "note": "unresolved -- no direct calls"}
    star6 = confusion(six_pairs) if six_pairs else {"n": 0}

    by_pop = {}
    for p in sorted({r["superpopulation"] for r in rows}):
        sub = [(r["truth_star28_n"], r["tag_star80_n"]) for r in rows
               if r["superpopulation"] == p and r["tag_star80_n"] is not None]
        subd = [(r["truth_star28_n"], r["direct_star28_n"]) for r in rows
                if r["superpopulation"] == p and r["direct_star28_n"] is not None]
        if sub:
            by_pop[p] = {"tag": confusion(sub), "direct": confusion(subd) if subd else {"n": 0}}

    # --- frozen verdicts ---
    sens, ppv = tag["sensitivity"], tag["ppv"]
    failing_pop = [p for p, v in by_pop.items()
                   if v["tag"]["n"] >= 10 and (v["tag"]["sensitivity"] or 0) < 0.80]
    if sens is not None and sens >= 0.90 and (ppv or 0) >= 0.90 and not failing_pop:
        tag_verdict = "TAG_CONFIRMED"
    elif (sens is not None and sens < 0.80) or failing_pop:
        tag_verdict = "TAG_DEMOTE_AND_DISCLOSE"
    else:
        tag_verdict = "TAG_PARTIAL"

    if direct.get("n"):
        dominates = (
            (direct["sensitivity"] or 0) >= (sens or 0)
            and (direct["ppv"] or 0) >= (ppv or 0)
            and (direct["zygosity_exact_frac"] or 0) >= (tag["zygosity_exact_frac"] or 0)
            and ((direct["sensitivity"] or 0) > (sens or 0)
                 or (direct["ppv"] or 0) > (ppv or 0)
                 or (direct["zygosity_exact_frac"] or 0) > (tag["zygosity_exact_frac"] or 0))
        )
        direct_verdict = "DIRECT_STRICTLY_DOMINATES_TAG" if dominates else "DIRECT_DOES_NOT_DOMINATE"
    else:
        direct_verdict = "DIRECT_ARM_UNRESOLVED"

    out = {
        "schema": "ugt1a1-getrm-concordance-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "cell": "pgx:human:ugt1a1",
        "acceptance_bar": str(BAR.relative_to(REPO)).replace("\\", "/"),
        "truth_source": bar["truth"]["source"],
        "independence": bar["truth"]["independence"],
        "n_truth_rows": len(truth),
        "n_scored": len(rows),
        "n_truth_star28_carriers": n_carriers,
        "n_parse_failures": len(parse_failures),
        "parse_failures": parse_failures,
        "tag_arm": tag,
        "direct_repeat_arm": direct,
        "star6_control": star6,
        "post_hoc_reduced_function_arm": {
            "post_hoc": True,
            "not_in_frozen_bar": True,
            "target": "reduced-function promoter carriership (*28 TA7 OR *37 TA8); *36 (TA5) excluded "
                      "as INCREASED function",
            "why": "Diagnosis of the frozen arm's errors showed rs887829-T tags the reduced-function "
                   "promoter CLASS, not *28 specifically. Irinotecan dosing turns on reduced function. "
                   "Reported BESIDE the frozen primary, never in place of it.",
            **confusion(red_pairs),
        },
        "diagnosed_tag_mechanism": (
            "rs887829-T tags the reduced-function promoter haplotype class -- BOTH *28 (TA7) and *37 "
            "(TA8) -- not *28 alone. Confirmed directly: the three *28-specific tag errors (HG01190, "
            "NA19920, NA19239) are haplotypes where the C>CATAT (+2 TA = *37) allele is present. "
            "Against a *28-specific target these read as false positives; against reduced FUNCTION "
            "they are correct calls."),
        "unresolved_discrepancy": {
            "sample": "NA20509",
            "truth": "(*28)/(*28 + *60) -- GeT-RM consensus = 2 copies of *28",
            "observed": "1 copy, by BOTH independent assays (rs887829 tag GT 0|1 AND the C>CAT repeat "
                        "indel GT 0|1, phase-consistent on the same haplotype); the other haplotype is "
                        "reference TA6 at every repeat allele",
            "reading": "Two independent assays of different variant types agree with each other and "
                       "against the reference consensus. Reported as an UNRESOLVED panel-vs-consensus "
                       "discrepancy -- this run does NOT adjudicate which is right, and it is the single "
                       "error in the direct arm.",
        },
        "by_superpopulation": by_pop,
        "tag_verdict": tag_verdict,
        "direct_verdict": direct_verdict,
        "failing_superpopulations": failing_pop,
        "scope_limit": bar["direct_repeat_arm"]["scope_limit_even_if_it_wins"],
        "honest_limits": bar["honest_limits_registered_in_advance"],
        "rows": rows,
    }

    stamp = datetime.date.today().isoformat()
    jf = REPO / "wiki" / f"ugt1a1_getrm_concordance_{stamp}.json"
    jf.write_text(json.dumps(out, indent=2), encoding="utf-8")

    lines = [
        f"# UGT1A1 tag-SNP surface vs GeT-RM CDC consensus ({stamp})", "",
        f"Scored **{len(rows)}** 1000G samples with CDC multi-lab consensus UGT1A1 diplotypes "
        f"({n_carriers} carry `*28`). Bar frozen beforehand: `{jf.parent.name}/{BAR.name}`.", "",
        "## The shipped tag arm (rs887829 / `*80` as a proxy for `*28`)", "",
        f"- sensitivity **{tag['sensitivity']}** · specificity **{tag['specificity']}** · "
        f"PPV **{tag['ppv']}** · NPV **{tag['npv']}** (n={tag['n']})",
        f"- exact `*28` **zygosity** agreement: {tag['zygosity_exact']}/{tag['n']} "
        f"({tag['zygosity_exact_frac']})",
        f"- **verdict: {tag_verdict}**", "",
        "## The direct TA-repeat arm (the catalog's asserted structural wall)", "",
    ]
    if direct.get("n"):
        lines += [
            f"- sensitivity **{direct['sensitivity']}** · specificity **{direct['specificity']}** · "
            f"PPV **{direct['ppv']}** · NPV **{direct['npv']}** (n={direct['n']})",
            f"- exact `*28` **zygosity** agreement: {direct['zygosity_exact']}/{direct['n']} "
            f"({direct['zygosity_exact_frac']})",
            f"- **verdict: {direct_verdict}**", "",
            f"**Scope limit.** {out['scope_limit']}", "",
        ]
    else:
        lines += [f"- **{direct_verdict}** -- no direct repeat calls resolved.", ""]
    lines += ["## By superpopulation (the catalog names non-EUR LD breakdown as an unmeasured residual)",
              "", "| pop | n | tag sens | tag PPV | direct sens |", "|---|---|---|---|---|"]
    for p, v in sorted(by_pop.items()):
        d = v["direct"].get("sensitivity") if v["direct"].get("n") else None
        lines.append(f"| {p} | {v['tag']['n']} | {v['tag']['sensitivity']} | {v['tag']['ppv']} | {d} |")
    rf = out["post_hoc_reduced_function_arm"]
    lines += ["", "## POST-HOC: the clinical axis (reduced FUNCTION, not `*28` specifically)", "",
              "*Not in the frozen bar. Added after diagnosing the frozen arm's errors; reported beside "
              "the frozen primary, never in place of it.*", "",
              f"**Diagnosed mechanism.** {out['diagnosed_tag_mechanism']}", "",
              f"- tag vs reduced-function carriership (`*28` or `*37`): sensitivity **{rf['sensitivity']}** "
              f"· specificity **{rf['specificity']}** · PPV **{rf['ppv']}** (n={rf['n']}), "
              f"zygosity-exact {rf['zygosity_exact']}/{rf['n']}", "",
              "## Unresolved panel-vs-consensus discrepancy", "",
              f"- **{out['unresolved_discrepancy']['sample']}** — truth "
              f"`{out['unresolved_discrepancy']['truth']}`; observed "
              f"{out['unresolved_discrepancy']['observed']}. {out['unresolved_discrepancy']['reading']}",
              "",
              "## Control: `*6` (rs4148323, a direct SNP -- not a tag)", "",
              f"- n={star6.get('n')} · sensitivity {star6.get('sensitivity')} · "
              f"PPV {star6.get('ppv')} · zygosity-exact {star6.get('zygosity_exact_frac')}", "",
              "## Honest limits", ""]
    lines += [f"- {h}" for h in out["honest_limits"]]
    if parse_failures:
        lines += ["", f"- {len(parse_failures)} truth diplotype(s) did not parse and were EXCLUDED "
                      f"(never defaulted to non-carrier)."]
    mf = REPO / "wiki" / f"ugt1a1_getrm_concordance_{stamp}.md"
    mf.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"scored={len(rows)}  *28 carriers={n_carriers}  parse_failures={len(parse_failures)}")
    print(f"TAG    sens={tag['sensitivity']} spec={tag['specificity']} ppv={tag['ppv']} "
          f"zyg={tag['zygosity_exact_frac']}  -> {tag_verdict}")
    if direct.get("n"):
        print(f"DIRECT sens={direct['sensitivity']} spec={direct['specificity']} ppv={direct['ppv']} "
              f"zyg={direct['zygosity_exact_frac']}  -> {direct_verdict}")
    for p, v in sorted(by_pop.items()):
        print(f"   {p:8} n={v['tag']['n']:3} tag_sens={v['tag']['sensitivity']}")
    print(f"wrote {jf.name} + {mf.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
