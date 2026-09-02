"""The control that decides whether PRJNA1322038's 60 S-labelled `rmt` carriers are evidence or artifact.

THE SETUP. The inverted specificity hunt found 60 S-labelled `rmt` carriers -- the counter-examples the
deployed gentamicin v2 rule had never been tested against. Every one traces to a SINGLE BioProject
(PRJNA1322038, University of Queensland): inside it 0 R / 60 S, outside it 146 R / 0 S across 23 other
projects. Perfect separation by submission is gate G2 (study == class) at its maximum.

WHY A CONTROL IS STILL NEEDED. "One project" is not by itself a verdict. Two hypotheses survive and they
have opposite consequences:

  H-artifact  That project's gentamicin labels are systematically odd (encoding, method, or breakpoint),
              so its S calls say nothing about the rule. Then the rule's specificity stays UNTESTED.
  H-real      That project's labels are ordinary, and it genuinely contains rmt carriers that test
              susceptible. Then the deployed rule really does over-call and that is a defect to fix.

THE DISCRIMINATOR is the project's OWN internal behaviour, which needs no external comparison:
  1. Base rate -- what fraction of ALL its gentamicin-labelled isolates are R? A project that is ~100% S
     everywhere is reporting something unusual; one with a normal R/S mix is labelling normally.
  2. The aac(3) control -- `aac(3)` is the CLASSIC gentamicin determinant that the FROZEN rule already
     counted before the v2 rescue existed. If this project labels its aac(3) carriers R, its labels track
     gentamicin resistance correctly and its rmt-carrier S calls are a specific, real signal. If it labels
     aac(3) carriers S too, the whole project's gentamicin column is not behaving like the others.

That second check is the load-bearing one: it uses a determinant nobody disputes as an internal yardstick.

Writes wiki/gentamicin_rmt_project_control.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dna_decode.data.pd_ast import ast_label_for  # noqa: E402
from gentamicin_rmt_specificity_hunt import (  # noqa: E402
    RESCUE, RESCUE_RE, latest_metadata_url, parse_amr_genotypes,
)

AAC3_RE = re.compile(r"^aac\(3\)")
ARMA_RE = re.compile(r"^armA\d*$")


def scan(group: str, project: str, drug: str) -> dict:
    r = urllib.request.urlopen(latest_metadata_url(group), timeout=600)
    cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    gi, ai = idx["AMR_genotypes"], idx["AST_phenotypes"]
    bi = idx.get("bioproject_acc")

    inside = {"all": collections.Counter(), "rmt": collections.Counter(),
              "aac3_no_rmt": collections.Counter(), "arma_no_rmt": collections.Counter(),
              "no_known_gent": collections.Counter()}
    outside = {"all": collections.Counter(), "rmt": collections.Counter(),
               "aac3_no_rmt": collections.Counter()}
    other_drugs_inside = collections.Counter()

    for line in r:
        f = line.decode("utf8", "replace").rstrip("\n").split("\t")
        if len(f) <= max(gi, ai):
            continue
        in_proj = bi is not None and len(f) > bi and f[bi] == project
        label = ast_label_for(f[ai], drug)
        if in_proj:
            # how does this project label OTHER drugs? an all-S project is a different story
            for d in ("ciprofloxacin", "ceftriaxone", "meropenem", "ampicillin"):
                lab = ast_label_for(f[ai], d)
                if lab:
                    other_drugs_inside[f"{d}:{lab}"] += 1
        if label not in ("R", "S", "I"):
            continue
        syms = parse_amr_genotypes(f[gi])
        has_rmt = any(RESCUE_RE.match(s) for s in syms)
        has_aac3 = any(AAC3_RE.match(s) for s in syms)
        has_arma = any(ARMA_RE.match(s) for s in syms)
        tgt = inside if in_proj else outside
        tgt["all"][label] += 1
        if has_rmt:
            tgt["rmt"][label] += 1
        elif has_aac3:
            tgt["aac3_no_rmt"][label] += 1
        elif in_proj and has_arma:
            tgt["arma_no_rmt"][label] += 1
        elif in_proj and not (has_aac3 or has_arma):
            tgt["no_known_gent"][label] += 1

    return {"group": group,
            "inside": {k: dict(v) for k, v in inside.items()},
            "outside": {k: dict(v) for k, v in outside.items()},
            "other_drugs_inside": dict(other_drugs_inside)}


def rate(c: dict) -> str:
    r, s = c.get("R", 0), c.get("S", 0)
    return f"{r}R/{s}S" + (f"  ({r / (r + s):.0%} R)" if r + s else "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", default="PRJNA1322038")
    ap.add_argument("--groups", default="Klebsiella,Pseudomonas_aeruginosa,Escherichia_coli_Shigella")
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "wiki" /
                    "gentamicin_rmt_project_control.json")
    a = ap.parse_args()

    print(f"control on {a.project}; deployed rescue = {RESCUE}\n")
    per_group = []
    for g in [x.strip() for x in a.groups.split(",") if x.strip()]:
        res = scan(g, a.project, a.drug)
        per_group.append(res)
        i, o = res["inside"], res["outside"]
        print(f"{g}")
        print(f"   INSIDE  all={rate(i['all']):22} rmt={rate(i['rmt']):20} "
              f"aac(3),no-rmt={rate(i['aac3_no_rmt'])}")
        print(f"           armA,no-rmt={rate(i['arma_no_rmt']):14} "
              f"no-known-gent-gene={rate(i['no_known_gent'])}")
        print(f"   OUTSIDE all={rate(o['all']):22} rmt={rate(o['rmt']):20} "
              f"aac(3),no-rmt={rate(o['aac3_no_rmt'])}\n", flush=True)

    def tot(side, key):
        c = collections.Counter()
        for g in per_group:
            c.update(g[side][key])
        return dict(c)

    ins_all, ins_rmt = tot("inside", "all"), tot("inside", "rmt")
    ins_aac3, ins_none = tot("inside", "aac3_no_rmt"), tot("inside", "no_known_gent")
    out_all, out_rmt, out_aac3 = tot("outside", "all"), tot("outside", "rmt"), tot("outside", "aac3_no_rmt")

    print("=" * 74)
    print(f"INSIDE  {a.project}:  all {rate(ins_all)} | rmt {rate(ins_rmt)} | "
          f"aac(3) {rate(ins_aac3)} | no known gent gene {rate(ins_none)}")
    print(f"OUTSIDE            :  all {rate(out_all)} | rmt {rate(out_rmt)} | aac(3) {rate(out_aac3)}")

    # The verdict rule, stated before the numbers were seen: if the project labels its aac(3) carriers
    # R at a normal rate while calling rmt carriers S, its labels track gentamicin and the signal is
    # specific. If it calls aac(3) carriers S too, its gentamicin column is not comparable.
    a3r, a3s = ins_aac3.get("R", 0), ins_aac3.get("S", 0)
    aac3_rate = a3r / (a3r + a3s) if (a3r + a3s) else None
    out_a3r, out_a3s = out_aac3.get("R", 0), out_aac3.get("S", 0)
    out_rate = out_a3r / (out_a3r + out_a3s) if (out_a3r + out_a3s) else None

    if aac3_rate is None:
        verdict, why = "INCONCLUSIVE", "the project has no aac(3)-carrying, rmt-free labelled isolate"
    elif aac3_rate < 0.5 and (out_rate or 0) > 0.8:
        verdict = "LABEL_ARTIFACT"
        why = (f"the project calls its aac(3) carriers R only {aac3_rate:.0%} of the time vs "
               f"{out_rate:.0%} elsewhere -- a determinant nobody disputes. Its gentamicin column "
               f"does not behave like the rest of PD, so its rmt S calls cannot test the rule.")
    else:
        verdict = "SPECIFIC_TO_RMT"
        why = (f"the project calls its aac(3) carriers R {aac3_rate:.0%} of the time -- comparable to "
               f"{out_rate if out_rate is None else format(out_rate, '.0%')} elsewhere -- so its labels "
               f"track gentamicin, and S on rmt carriers is a SPECIFIC signal that needs explaining.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "gentamicin-rmt-project-control-v1", "project": a.project, "drug": a.drug,
           "deployed_rescue_pattern": RESCUE, "per_group": per_group,
           "inside_totals": {"all": ins_all, "rmt": ins_rmt, "aac3_no_rmt": ins_aac3,
                             "no_known_gent_gene": ins_none},
           "outside_totals": {"all": out_all, "rmt": out_rmt, "aac3_no_rmt": out_aac3},
           "aac3_R_rate_inside": aac3_rate, "aac3_R_rate_outside": out_rate,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "This tests INTERNAL consistency of one project's gentamicin column against a determinant "
               "nobody disputes. It cannot prove the labels wrong -- only that they do or do not behave "
               "like every other submission.",
               "A LABEL_ARTIFACT verdict leaves the deployed rule's specificity UNTESTED, not vindicated. "
               "It means these 60 isolates cannot serve as the test, not that no counter-example exists.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
