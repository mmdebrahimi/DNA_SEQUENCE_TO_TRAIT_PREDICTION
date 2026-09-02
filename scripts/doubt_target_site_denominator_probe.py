"""Can the doubt layer's COMPLETENESS screen extend to the target-site catalogs? Measure, don't build.

THE OPEN QUESTION (F-A candidate 2): "extend the completeness screen to the target-site catalogs -- one
vocabulary or two?" The AMR arm's screen ranks determinant FAMILIES the deployed rule cannot represent by
a PURITY signature: a family whose carriers are labelled R and never S is a candidate rule gap. It found
exactly one confirmed gap (`rmtE1`) after a family-wise correction.

WHY MEASURE FIRST. The same question has now twice been settled by a census rather than by building: the
per-cell regime field looked like ceremony until the census showed the column was not constant, and the
NNRTI curation looked obvious until the recovery was measured against the incumbent. The screen's purity
signature needs a specific shape to mean anything:

  1. a NEGATIVE class -- some carriers labelled susceptible, or "zero S carriers" is true by construction
     and carries no information;
  2. a candidate unit the deployed rule CANNOT already represent -- otherwise the flag can never fire on
     anything actionable;
  3. enough carriers per unit that zero-S is surprising rather than small-n.

This probe measures all three on the HIV NNRTI cell -- the target-site cell with a free, independent,
isolate-level wet-lab label (Stanford PhenoSense fold-change) and therefore the BEST case. If the shape
fails on the best case, it fails on the others, and the answer is "two vocabularies", settled cheaply.

Reads the gitignored Stanford dataset; skips cleanly when absent. Writes
wiki/doubt_target_site_denominator_probe.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dna_decode.data.hiv_amr import NNRTI_RT_MAJOR_DRMS  # noqa: E402
from dna_decode.eval.doubt import completeness_surprise, completeness_tier  # noqa: E402

DATA = Path("data/raw/hiv/NNRTI_DataSet.txt")
CUTOFF = 3.0          # the sourced Stanford DRMcv.R clinical cutoff for EFV/NVP/ETR/RPV
SUB = re.compile(r"^([A-Z])(\d+)([A-Za-z*]+)$")
MIN_CARRIERS = 5


def catalogued_positions() -> set[int]:
    return {int("".join(c for c in m if c.isdigit())) for m in NNRTI_RT_MAJOR_DRMS}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="EFV")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "wiki" /
                    "doubt_target_site_denominator_probe.json")
    a = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    path = root / DATA
    if not path.is_file():
        print(f"SKIP: {path} absent (gitignored dataset).")
        return 0

    rows = [ln.rstrip("\n").split("\t") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    hdr, rows = rows[0], rows[1:]
    idx = {c: i for i, c in enumerate(hdr)}
    di = idx.get(a.drug)
    mi = idx.get("CompMutList")
    if di is None or mi is None:
        print(f"SKIP: columns {a.drug}/CompMutList not found in {hdr[:8]}")
        return 0

    cat_pos = catalogued_positions()
    # unit -> Counter({R,S}); the candidate unit is a SUBSTITUTION outside the catalogued positions,
    # i.e. exactly what the deployed mutant-level rule cannot represent.
    per_sub: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    n_used = n_R = n_S = 0
    for r in rows:
        if len(r) <= max(di, mi):
            continue
        raw = r[di].strip()
        if not raw or raw.upper() in {"NA", "NULL", "."}:
            continue
        try:
            fold = float(raw.lstrip("<>=~"))
        except ValueError:
            continue
        label = "R" if fold >= CUTOFF else "S"
        n_used += 1
        n_R += label == "R"
        n_S += label == "S"
        seen = set()
        for tok in r[mi].split(","):
            m = SUB.match(tok.strip())
            if not m:
                continue
            wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
            if wt == mut:            # self-to-self = a mixture marker, not a substitution
                continue
            if pos in cat_pos:       # already representable by the deployed catalog
                continue
            for aa in mut:
                seen.add(f"{wt}{pos}{aa}")
        for s in seen:
            per_sub[s][label] += 1

    base_s = n_S / n_used if n_used else 0.0
    eligible = {s: c for s, c in per_sub.items() if (c["R"] + c["S"]) >= MIN_CARRIERS}
    pure = {s: c for s, c in eligible.items() if c["S"] == 0}
    n_tested = len(eligible)

    strong = []
    for s, c in pure.items():
        p = completeness_surprise(c["R"], c["S"], base_s)
        if completeness_tier(p, n_tested) == "strong":
            strong.append((s, c["R"], p))

    print(f"HIV NNRTI / {a.drug}: {n_used} isolates with a fold value  ({n_R} R, {n_S} S at fold>={CUTOFF})")
    print(f"  base susceptible rate                 : {base_s:.3f}")
    print(f"  candidate units (non-catalogued subs) : {len(per_sub)}")
    print(f"  with >= {MIN_CARRIERS} carriers                    : {n_tested}")
    print(f"  of those, PURE (zero S carriers)      : {len(pure)}")
    print(f"  surviving family-wise correction      : {len(strong)}")
    for s, n, p in sorted(strong, key=lambda x: x[2])[:10]:
        print(f"      {s:10} {n} carriers, all R   p={p:.2e}")

    # The three shape conditions, answered from the numbers rather than asserted.
    has_negative = n_S >= MIN_CARRIERS
    has_units = n_tested >= 10
    informative = len(pure) < n_tested          # if EVERY unit is pure, purity separates nothing
    shape_ok = has_negative and has_units and informative
    verdict = "ONE_VOCABULARY" if shape_ok and strong else (
        "SHAPE_OK_NO_HITS" if shape_ok else "TWO_VOCABULARIES")

    why = {
        "ONE_VOCABULARY": "the purity signature has a real denominator here AND fires, so the AMR screen's "
                          "vocabulary transfers to the target-site arm.",
        "SHAPE_OK_NO_HITS": "the signature is well-formed here (negative class, many units, purity is "
                            "discriminating) but nothing survives correction -- the screen COULD extend; "
                            "on this cell it finds nothing, which is a result, not a blocker.",
        "TWO_VOCABULARIES": "the signature is malformed on the best-case target-site cell, so forcing the "
                            "AMR vocabulary onto this arm would report structure that is not there.",
    }[verdict]

    print(f"\n  negative class present (>= {MIN_CARRIERS} S) : {has_negative}  (n_S={n_S})")
    print(f"  enough candidate units (>= 10)        : {has_units}")
    print(f"  purity is discriminating (not all)    : {informative}  ({len(pure)}/{n_tested} pure)")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "doubt-target-site-denominator-probe-v1", "cell": f"HIV NNRTI / {a.drug}",
           "fold_cutoff": CUTOFF, "n_isolates_scored": n_used, "n_R": n_R, "n_S": n_S,
           "base_susceptible_rate": base_s, "catalogued_positions": sorted(cat_pos),
           "n_candidate_units": len(per_sub), "n_units_with_min_carriers": n_tested,
           "n_pure_units": len(pure), "n_surviving_familywise": len(strong),
           "surviving": [{"substitution": s, "carriers": n, "p": p} for s, n, p in strong],
           "shape": {"negative_class_present": has_negative, "enough_units": has_units,
                     "purity_discriminating": informative},
           "verdict": verdict, "why": why,
           "honest_limits": [
               "Measured on ONE target-site cell (HIV NNRTI) -- deliberately the BEST case, since it is "
               "the only one with a free independent isolate-level wet-lab label. A negative here "
               "generalises downward; a positive would need re-checking per cell.",
               "The R/S split is a threshold on a continuous fold-change at the sourced DRMcv.R cutoff; "
               "the AMR arm's labels are categorical AST. The two are not the same measurement.",
               "In-distribution: catalog and label both trace to Stanford HIVDB.",
               "This probe answers whether the SHAPE supports the screen. It is not itself a curation "
               "recommendation -- data-derived NNRTI curation was measured and DECLINED on 2026-09-01.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
