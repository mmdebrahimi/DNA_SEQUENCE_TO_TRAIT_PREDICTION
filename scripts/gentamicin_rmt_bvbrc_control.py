"""The control that decides whether BV-BRC's 67 susceptible `rmt` carriers are evidence or artifact.

SAME QUESTION, SAME YARDSTICK AS THE PD CONTROL -- deliberately, so the two verdicts are comparable.
On NCBI-PD, 60 susceptible carriers all came from one BioProject and were killed as a LABEL_ARTIFACT:
that project called its `aac(3)` carriers resistant 2% of the time against 97% everywhere else, so its
gentamicin column did not track genotype at all.

BV-BRC's susceptible carriers are ALSO concentrated (94% from one study, pmid 36801013), so the same
suspicion applies and the same discriminator settles it. `aac(3)` is the classic gentamicin-modifying
enzyme that the FROZEN rule already counted before the v2 rescue existed, which is what makes it a
yardstick independent of the rescue under scrutiny.

THE DECISION RULE, stated before the numbers are read (identical to the PD control):
  - the study calls its aac(3) carriers R at a rate far BELOW everyone else  -> LABEL_ARTIFACT
  - the study calls its aac(3) carriers R at a comparable rate               -> SPECIFIC_TO_RMT

A SPECIFIC_TO_RMT verdict is the one that costs us something: it means the labels are sound and the
susceptible carriers are real counter-examples to the deployed rule.

Network-only. Writes wiki/gentamicin_rmt_bvbrc_control.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

API = "https://www.bv-brc.org/api"
UA = {"Accept": "application/json",
      "User-Agent": "dna_decode/0.13 (research; genotype-phenotype validation)"}


def q(coll: str, query: str, limit: int = 25000, offset: int = 0) -> list[dict]:
    url = f"{API}/{coll}/?{query}&limit({limit},{offset})"
    raw = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300).read()
    data = json.loads(raw.decode("utf8", "replace"))
    if isinstance(data, dict):                     # BV-BRC wraps outages in HTTP 200
        raise RuntimeError(f"BV-BRC error envelope: {str(data)[:200]}")
    return data


def has(genes: set[str], prefix: str) -> bool:
    return any(g.lower().startswith(prefix) for g in genes)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pmid", default="36801013", help="the dominant study among susceptible carriers")
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / "gentamicin_rmt_bvbrc_control.json")
    a = ap.parse_args()

    sel = ("&select(genome_id,genome_name,resistant_phenotype)")
    inside = q("genome_amr", f"and(eq(antibiotic,{a.drug}),eq(evidence,Laboratory%20Method),"
                             f"eq(pmid,{a.pmid})){sel}")
    ph = {str(r["genome_id"]): r["resistant_phenotype"] for r in inside if r.get("resistant_phenotype")}
    names = {str(r["genome_id"]): r.get("genome_name", "") for r in inside}
    print(f"study pmid {a.pmid}: {len(ph)} measured {a.drug} records", flush=True)

    gids = sorted(ph)
    carr: dict[str, set[str]] = {}
    for i in range(0, len(gids), 150):
        chunk = ",".join(gids[i:i + 150])
        for r in q("sp_gene", f"and(in(genome_id,({chunk})),eq(source,CARD))&select(genome_id,gene)"):
            carr.setdefault(str(r["genome_id"]), set()).add(str(r.get("gene", "")))

    def stratum(pred) -> dict:
        picked = [g for g in gids if pred(carr.get(g, set()))]
        c = collections.Counter(ph[g] for g in picked)
        r, s = c.get("Resistant", 0), c.get("Susceptible", 0)
        return {"R": r, "S": s, "I": c.get("Intermediate", 0), "n": r + s,
                "r_rate": (r / (r + s)) if (r + s) else None}

    inside_tot = {
        "aac3_no_rmt": stratum(lambda g: has(g, "aac(3)") and not has(g, "rmt")),
        "rmt": stratum(lambda g: has(g, "rmt")),
        "no_known_determinant": stratum(
            lambda g: not has(g, "aac(3)") and not has(g, "rmt") and not has(g, "arma")),
    }

    # The comparator: the same aac(3) stratum computed OUTSIDE this study.
    outside = q("genome_amr", f"and(eq(antibiotic,{a.drug}),eq(evidence,Laboratory%20Method),"
                              f"ne(pmid,{a.pmid})){sel}")
    oph = {str(r["genome_id"]): r["resistant_phenotype"] for r in outside if r.get("resistant_phenotype")}
    ogids = sorted(oph)[:6000]                       # bounded: this is a comparator, not a census
    ocarr: dict[str, set[str]] = {}
    for i in range(0, len(ogids), 150):
        chunk = ",".join(ogids[i:i + 150])
        for r in q("sp_gene", f"and(in(genome_id,({chunk})),eq(source,CARD))&select(genome_id,gene)"):
            ocarr.setdefault(str(r["genome_id"]), set()).add(str(r.get("gene", "")))
    oc = collections.Counter(oph[g] for g in ogids
                             if has(ocarr.get(g, set()), "aac(3)") and not has(ocarr.get(g, set()), "rmt"))
    o_r, o_s = oc.get("Resistant", 0), oc.get("Susceptible", 0)
    outside_aac3 = {"R": o_r, "S": o_s, "n": o_r + o_s,
                    "r_rate": (o_r / (o_r + o_s)) if (o_r + o_s) else None,
                    "n_genomes_sampled": len(ogids)}

    for k, v in inside_tot.items():
        rate = f"{v['r_rate']:.0%}" if v["r_rate"] is not None else "n/a"
        print(f"  INSIDE  {k:24} R={v['R']:4d} S={v['S']:4d}  ({rate} R)")
    rate_o = f"{outside_aac3['r_rate']:.0%}" if outside_aac3["r_rate"] is not None else "n/a"
    print(f"  OUTSIDE aac3_no_rmt              R={o_r:4d} S={o_s:4d}  ({rate_o} R), "
          f"sampled {len(ogids)} genomes")

    ai, ao = inside_tot["aac3_no_rmt"]["r_rate"], outside_aac3["r_rate"]
    if ai is None or ao is None:
        verdict, why = "INCONCLUSIVE", "no aac(3)-carrying, rmt-free labelled isolate on one side"
    elif ai < 0.5 and ao > 0.8:
        verdict = "LABEL_ARTIFACT"
        why = (f"the study calls its aac(3) carriers R only {ai:.0%} of the time vs {ao:.0%} elsewhere; "
               "its gentamicin column does not behave like the rest of the archive")
    else:
        verdict = "SPECIFIC_TO_RMT"
        why = (f"the study calls its aac(3) carriers R {ai:.0%} of the time vs {ao:.0%} elsewhere, and "
               f"its isolates with no known determinant R "
               f"{inside_tot['no_known_determinant']['r_rate']:.0%} of the time -- its labels TRACK "
               "genotype, so susceptibility on rmt carriers is a specific signal that needs explaining")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "gentamicin-rmt-bvbrc-control-v1", "archive": "BV-BRC", "pmid": a.pmid,
           "drug": a.drug, "inside": inside_tot, "outside_aac3_no_rmt": outside_aac3,
           "aac3_R_rate_inside": ai, "aac3_R_rate_outside": ao,
           "verdict": verdict, "why": why,
           "honest_limits": [
               "This tests whether ONE study's gentamicin column tracks genotype, using a determinant "
               "nobody disputes. It cannot prove any individual label correct.",
               "The outside comparator is a bounded sample of genomes, not a full census -- it is a "
               "reference rate, not an archive-wide statistic.",
               "Carrier calls are CARD/BLAT via BV-BRC sp_gene, a different tool from our AMRFinder.",
               "A SPECIFIC_TO_RMT verdict makes the susceptible carriers real counter-examples; it does "
               "NOT by itself say the deployed rule is wrong in its own validated organism scope, which "
               "must be read per-organism.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
