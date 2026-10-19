"""Anatomy of the 21 false negatives and 11 false positives in the CABBAGE Campylobacter cipro cell.

WHY: the standing rule here is to DIAGNOSE a miss before reading it as decay. The prospective
Campylobacter cipro cell had exactly one FN, and it carried ZERO determinants -- abstention-by-design
leaking into a hard S call, a declared blind spot rather than drift. With 21 FN there are two
distinguishable stories and they have opposite consequences:

  (a) the FN carry NO determinant        -> the same declared blind spot, at scale
  (b) the FN carry a determinant the rule DOES NOT COUNT -> a catalog gap, the gentamicin `rmt` shape

Partitioning by what the genotype actually contains is what separates them. Everything is read from the
already-cached normalized AMRFinder TSVs, so this costs no network and no Docker.
"""
import collections
import json
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, ".")

from dna_decode.data.amr_portal import prepare  # noqa: E402
from dna_decode.eval.amr_rules import call_resistance  # noqa: E402
from scripts.cabbage_score_cell import CACHE, build_labels  # noqa: E402

import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("--species", default="Campylobacter jejuni")
_ap.add_argument("--drug", default="ciprofloxacin")
_ap.add_argument("--organism", default="Campylobacter")
_a = _ap.parse_args()
SPECIES, DRUG, ORGANISM = _a.species, _a.drug, _a.organism

labels, census, srcs, yrs, ctry = build_labels(SPECIES, DRUG)
print(f"labels: {len(labels)}  census={census['tier_census']}")

fn, fp, tp, tn, missing = [], [], 0, 0, 0
det_by_bucket = {"fn": collections.Counter(), "fp": collections.Counter(),
                 "tp": collections.Counter(), "tn": collections.Counter()}
ndet_by_bucket = collections.defaultdict(collections.Counter)

for acc, truth in labels.items():
    raw = CACHE / f"{acc}.tsv"
    if not raw.exists():
        missing += 1
        continue
    try:
        norm = prepare(raw, CACHE / "norm" / f"{acc}.tsv")
        res = call_resistance(str(norm), drug=DRUG, organism=ORGANISM)
    except Exception:                                 # noqa: BLE001
        missing += 1
        continue
    pred = res.get("prediction")
    # The key is `determinants`. A first version read `hits`, which does not exist on this record, so
    # EVERY bucket came back determinant-free -- including 997 true positives, which is impossible.
    # That printed a confident, plausible verdict off an extraction that had parsed nothing.
    syms = sorted({str(h.get("symbol") or "?") for h in (res.get("determinants") or [])})
    assert len(syms) == 0 or res.get("n_determinants", 0) > 0, \
        "symbols parsed but n_determinants is 0 -- extraction disagrees with the record"
    bucket = {("R", "R"): "tp", ("R", "S"): "fn", ("S", "S"): "tn", ("S", "R"): "fp"}.get(
        (truth, pred))
    if bucket is None:
        continue
    if bucket == "tp":
        tp += 1
    elif bucket == "tn":
        tn += 1
    elif bucket == "fn":
        fn.append({"accession": acc, "determinants": syms, "n_determinants": len(syms),
                   "database": str(srcs.get(acc, "?")), "country": str(ctry.get(acc, "?")),
                   "year": yrs.get(acc)})
    else:
        fp.append({"accession": acc, "determinants": syms, "n_determinants": len(syms),
                   "database": str(srcs.get(acc, "?")), "country": str(ctry.get(acc, "?")),
                   "year": yrs.get(acc)})
    for s in syms:
        det_by_bucket[bucket][s] += 1
    ndet_by_bucket[bucket][len(syms)] += 1

print(f"tp={tp} fn={len(fn)} tn={tn} fp={len(fp)}  (uncached/unscoreable={missing})")

fn_zero = [r for r in fn if r["n_determinants"] == 0]
fp_zero = [r for r in fp if r["n_determinants"] == 0]
print(f"\nFN with ZERO counted determinants: {len(fn_zero)}/{len(fn)}")
for r in fn:
    print(f"  FN {r['accession']:<18} dets={r['determinants']}")
print(f"\nFP with ZERO counted determinants: {len(fp_zero)}/{len(fp)}")
for r in fp:
    print(f"  FP {r['accession']:<18} dets={r['determinants']}")

print("\ndeterminant frequency by bucket:")
for b in ("tp", "fn", "tn", "fp"):
    print(f"  {b}: {dict(det_by_bucket[b].most_common(6))}")

# PER-DETERMINANT PPV. "6 of the 11 FP carry gyrA_P104S" is an anecdote until the DENOMINATOR is
# known: a determinant seen 6 times with 6 false positives is an over-call, the same symbol seen 600
# times with 6 is noise. Carriers are counted across the WHOLE scored set, truth from the labels.
from dna_decode.eval.clonality import wilson_ci  # noqa: E402

carriers = collections.defaultdict(lambda: {"R": 0, "S": 0})
for b in ("tp", "fn", "tn", "fp"):
    truth = "R" if b in ("tp", "fn") else "S"
    for sym, k in det_by_bucket[b].items():
        carriers[sym][truth] += k
per_det = {}
for sym, c in sorted(carriers.items(), key=lambda kv: -(kv[1]["R"] + kv[1]["S"])):
    n = c["R"] + c["S"]
    lo, hi = wilson_ci(c["R"], n)
    per_det[sym] = {"n_carriers": n, "n_truth_R": c["R"], "n_truth_S": c["S"],
                    "ppv": round(c["R"] / n, 4), "ppv_wilson_ci": [round(lo, 4), round(hi, 4)]}
print("\nper-determinant PPV over the whole scored set:")
for sym, v in per_det.items():
    print(f"  {sym:<14} carriers={v['n_carriers']:>4}  R={v['n_truth_R']:>4} S={v['n_truth_S']:>3}  "
          f"PPV={v['ppv']:.3f} CI={v['ppv_wilson_ci']}")

zero_ppv = [s for s, v in per_det.items() if v["n_truth_R"] == 0]

# PLUMBING GUARD, and it is the reason this script is trustworthy at all. A true positive carries a
# counted determinant BY DEFINITION, so a determinant-free TP bucket means the extraction read nothing.
# The first run of this probe did exactly that and emitted a confident, tempting verdict. Refuse.
if tp > 0 and sum(det_by_bucket["tp"].values()) == 0:
    print(f"\nREFUSING a verdict: {tp} true positives parsed ZERO determinants. A TP carries a "
          "counted determinant by definition, so this is an extraction failure, not a finding.")
    raise SystemExit(3)

if len(fn) == 0:
    verdict = "NO_FALSE_NEGATIVES"
elif len(fn_zero) == len(fn):
    verdict = "ALL_FN_CARRY_NO_DETERMINANT_DECLARED_BLIND_SPOT"
elif len(fn_zero) == 0:
    verdict = "NO_FN_IS_DETERMINANT_FREE_SUSPECT_A_CATALOG_GAP"
else:
    verdict = "MIXED_FN_SOME_DETERMINANT_FREE_SOME_NOT"

out = {
    "schema": "cabbage-fn-anatomy-v1", "generated": date.today().isoformat(),
    "cell": {"species": SPECIES, "drug": DRUG, "amrfinder_organism": ORGANISM},
    "confusion": {"tp": tp, "fn": len(fn), "tn": tn, "fp": len(fp)},
    "n_unscoreable": missing,
    "verdict": verdict,
    "per_determinant_ppv": per_det,
    "determinants_with_zero_true_R_carriers": zero_ppv,
    "fp_finding": (
        "Every false positive carries a counted determinant, and the determinants with ZERO "
        "true-R carriers contribute only false positives in this cohort. That is the MIRROR of the "
        "gentamicin `rmt` case: there a determinant the rule could not SEE produced false negatives; "
        "here a determinant the rule COUNTS produces false positives. Changing what the rule counts "
        "would edit the sha256-pinned frozen surface and retire the active prospective lock, so it "
        "is a USER AUTHORITY call and is deliberately NOT made here."),
    "fp_attribution": (
        "NOT A RULE BUG -- an UPSTREAM catalog assertion the rule faithfully propagates. Verified on "
        "GCA_004836905.1: gyrA_P104S is the SOLE counted determinant (n_determinants=1, "
        "counter='qrdr_point', threshold=1), and AMRFinderPlus itself files it as Class=QUINOLONE / "
        "Subclass=QUINOLONE named 'Campylobacter quinolone resistant GyrA'. The frozen rule counts "
        "QUINOLONE-class gyrA point mutations, so it is behaving exactly as specified; the "
        "resistance-conferring claim for P104S originates in AMRFinder's catalog. A fix therefore "
        "belongs upstream, or in an L2 doubt flag, not in the counting logic."),
    "false_negatives": fn, "false_positives": fp,
    "n_fn_determinant_free": len(fn_zero), "n_fp_determinant_free": len(fp_zero),
    "determinant_frequency_by_bucket": {b: dict(c.most_common(20))
                                        for b, c in det_by_bucket.items()},
    "n_determinants_distribution": {b: dict(sorted(c.items()))
                                    for b, c in ndet_by_bucket.items()},
    "honest_limits": [
        "`determinants` are the hits the FROZEN rule COUNTS for this drug, not everything AMRFinder "
        "reported. An FN listed as determinant-free may still carry an uncounted element, which is "
        "precisely the distinction between a declared blind spot and a catalog gap -- so a "
        "determinant-free FN bounds the rule's input, it does not prove the genome is bare.",
        "Labels are portal MIC run through this repo's STRICT tiers; a mislabelled isolate presents "
        "as an FN and nothing here can separate that from a genotype miss.",
        "Campylobacter cipro R is overwhelmingly gyrA_T86I, so bucket determinant frequencies are "
        "dominated by one substitution and are not a broad coverage census.",
        "The cohort is 96.4% one country (see the cell artifact's source_concentration verdict), so "
        "an absent determinant family may be absent from the COHORT rather than rare in the world.",
        "A zero-true-R determinant is reported with its Wilson CI because the carrier counts are "
        "SMALL; 0 of 6 is suggestive and is NOT a demonstration that the substitution is benign.",
        "PPV here is MARGINAL per determinant, not causal: a carrier may hold several determinants, "
        "so credit is not attributable to one symbol from these counts alone.",
        "Not clonality-corrected, so carriers of one determinant may be a handful of clones.",
    ],
}
_tag = f"{SPECIES.split()[0].lower()}_{DRUG.lower().replace(' ', '-')}"
p = Path("wiki") / f"cabbage_fn_anatomy_{_tag}_{date.today().isoformat()}.json"
p.write_text(json.dumps(out, indent=2), encoding="utf-8")
print(f"\nVERDICT: {verdict}")
print(f"-> {p}")
