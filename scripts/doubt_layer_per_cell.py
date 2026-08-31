"""Measure the L2 doubt layer PER CELL. Read-only; exit 0 always (a report, not a gate).

TWO ARMS, NOT ONE NUMBER. The doubt layer carries two signal kinds and they are NOT on a common
scale, so this script refuses to pool them:

  determinant_completeness (AMR cells)   per DETERMINANT FAMILY. "Does the deployed rule fail to
                                         represent a family whose labelled carriers are uniformly
                                         resistant?" There is exactly ONE known gap in this arm, so
                                         its recall is 1/1 -- a single case, never a rate.
  position_novelty (target-site cells)   per ISOLATE, and genuinely measured: median sensitivity
                                         0.604 on the EFV catalog-negative blind spot, lift 4.69.
                                         That number is the INCUMBENT this layer must be read
                                         against; it is reproduced here from its committed artifact,
                                         not recomputed.

Averaging a per-family screen against a per-isolate detector would produce a headline that describes
neither. Every figure below is per cell.

Run: uv run python scripts/doubt_layer_per_cell.py
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
WIKI = ROOT / "wiki"

INCUMBENT = {"source": "wiki/hiv_blindspot_position_novelty_2026-07-11.json",
             "metric": "median flag sensitivity on the EFV catalog-negative blind spot",
             "value": 0.604, "lift": 3.98}

# The one determinant-completeness gap this project has independently confirmed. Used to state the
# AMR arm's recall honestly as 1-of-1 rather than dressing a single case as a rate.
KNOWN_AMR_GAPS = {("gentamicin", "rmt")}


def base_s_rates(census: dict) -> dict:
    """Per-drug susceptible-rate among labelled genomes -- the null the purity test runs against."""
    c: dict[str, Counter] = defaultdict(Counter)
    for rec in (census.get("labels") or {}).values():
        for drug, call in (rec.get("calls") or {}).items():
            c[drug][str(call).upper()] += 1
    out = {}
    for drug, n in c.items():
        r, s = n.get("R", 0), n.get("S", 0)
        if r + s:
            out[drug] = {"n_labelled_r": r, "n_labelled_s": s, "base_s_rate": s / (r + s)}
    return out


def score_amr_arm(screen: dict, rates: dict) -> list[dict]:
    """Tier every uncounted determinant family, per drug. No pooling across drugs."""
    from dna_decode.eval.doubt import STRONG, completeness_signal

    out = []
    for d in screen.get("drugs", []):
        drug = d["drug"]
        rate = rates.get(drug)
        cands = d.get("candidates") or []
        n_tested = len(cands)
        if not rate:
            out.append({"cell_kind": "determinant_completeness", "drug": drug,
                        "status": "no_labels", "n_families_uncounted": n_tested,
                        "note": "no labelled genome carries a call for this drug -- unassessable, "
                                "which is not the same as clean"})
            continue

        sigs = [completeness_signal(c["symbol"], c["subclass"], c["r_carriers"], c["s_carriers"],
                                    rate["base_s_rate"], n_tested) for c in cands]
        strong = [s for s in sigs if s.tier == STRONG]
        weak = [s for s in sigs if s.tier == "weak"]
        # raw signature = step 1's uncorrected heuristic (>=3 R carriers, 0 S). Reported beside the
        # corrected count so the correction's effect is visible rather than asserted.
        raw = [c for c in cands if c["r_carriers"] >= 3 and c["s_carriers"] == 0]
        hit = [s for s in strong
               if any(s.evidence["symbol"].lower().startswith(fam) for dr, fam in KNOWN_AMR_GAPS
                      if dr == drug)]
        out.append({
            "cell_kind": "determinant_completeness", "drug": drug, "status": "scored",
            "n_genomes_scanned": d.get("n_genomes_scanned"),
            "n_determinants_probed": d.get("n_distinct_determinants_probed"),
            "n_counted_by_rule": d.get("n_counted_by_rule"),
            "n_families_uncounted": n_tested,
            "n_raw_signature": len(raw), "n_strong": len(strong), "n_weak": len(weak),
            "base_s_rate": round(rate["base_s_rate"], 4),
            "n_labelled_r": rate["n_labelled_r"], "n_labelled_s": rate["n_labelled_s"],
            "known_gap_recovered": bool(hit),
            "strong": [s.as_dict() for s in strong],
        })
    return out


def amr_arm_specificity(amr: list[dict]) -> dict:
    """How often does the AMR arm fire where no gap has ever been confirmed?

    THE CATEGORIES ARE PREDECLARED AND DELIBERATELY NOT CALLED "CLEAN". A drug with no confirmed
    completeness gap is `unconfirmed`, not gap-free -- both known gaps in this project were invisible
    until an independent label set arrived, so absence of a confirmed gap is absence of evidence. A
    STRONG hit on an `unconfirmed` drug is therefore AMBIGUOUS between a false positive and a new
    discovery, and nothing here can separate them without new labels. Reporting it as a false-positive
    RATE would assert the negative this project has twice been wrong about.

    What the number does bound: the screen's firing rate on the surface where no gap is known. Zero
    hits across that surface is the honest statement -- it is compatible with a specific screen AND
    with undiscovered gaps, and says so.
    """
    confirmed = {d for d, _fam in KNOWN_AMR_GAPS}
    conf_rows = [c for c in amr if c["drug"] in confirmed and c.get("status") == "scored"]
    unconf = [c for c in amr if c["drug"] not in confirmed and c.get("status") == "scored"]
    unassessable = [c["drug"] for c in amr if c.get("status") != "scored"]

    fam_unconf = sum(c.get("n_families_uncounted") or 0 for c in unconf)
    strong_unconf = sum(c.get("n_strong") or 0 for c in unconf)
    raw_unconf = sum(c.get("n_raw_signature") or 0 for c in unconf)
    return {
        "predeclared_categories": {
            "confirmed_gap": sorted(confirmed),
            "unconfirmed": sorted(c["drug"] for c in unconf),
            "unassessable_no_labels": sorted(unassessable),
        },
        "confirmed_gap_drugs": {
            "n_drugs": len(conf_rows),
            "n_families_screened": sum(c.get("n_families_uncounted") or 0 for c in conf_rows),
            "n_strong": sum(c.get("n_strong") or 0 for c in conf_rows),
        },
        "unconfirmed_drugs": {
            "n_drugs": len(unconf),
            "n_families_screened": fam_unconf,
            "n_strong_after_correction": strong_unconf,
            "n_raw_signature_before_correction": raw_unconf,
            "correction_removed": raw_unconf - strong_unconf,
        },
        "interpretation": (
            f"Across {fam_unconf} candidate families on {len(unconf)} drugs with no CONFIRMED "
            f"completeness gap, the corrected screen fires {strong_unconf} time(s); the raw purity "
            f"signature fired {raw_unconf}. This is NOT a false-positive rate -- a hit here would be "
            "ambiguous between a false positive and an undiscovered gap, and no evidence in hand can "
            "separate those."),
    }


def read_target_site_arm() -> list[dict]:
    """Reproduce the position-novelty arm from its committed artifact. Never recomputed here."""
    p = WIKI / "hiv_blindspot_position_novelty_2026-07-11.json"
    if not p.exists():
        return []
    d = json.loads(p.read_text(encoding="utf-8"))
    out = []
    for drug, r in (d.get("per_drug") or {}).items():
        out.append({"cell_kind": "position_novelty", "drug": drug, "status": "scored",
                    "organism": "HIV-1",
                    "sens_on_blindspot": r.get("flag_sens_on_blindspot"),
                    "fp_on_catalog_negative_S": r.get("flag_fp_on_catneg_S"),
                    "lift": r.get("lift"), "powered": r.get("powered"),
                    "n_blindspot_true_R": r.get("n_blindspot_true_R"),
                    "n_catalog_negative": r.get("n_catalog_negative")})
    return sorted(out, key=lambda x: x["drug"])


def main() -> int:
    screen_files = sorted(WIKI.glob("determinant_completeness_screen_*.json"))
    if not screen_files:
        print("no completeness-screen artifact -- run scripts/determinant_completeness_screen.py first")
        return 0
    screen = json.loads(screen_files[-1].read_text(encoding="utf-8"))

    census_p = WIKI / "unscored_genome_label_census.json"
    census = json.loads(census_p.read_text(encoding="utf-8")) if census_p.exists() else {}
    rates = base_s_rates(census)

    amr = score_amr_arm(screen, rates)
    tgt = read_target_site_arm()

    out = {
        "schema": "doubt-layer-per-cell-v1", "generated": date.today().isoformat(),
        "incumbent": INCUMBENT,
        "contract": ("A DOUBT signal qualifies a call and explains itself; it never emits one. The "
                     "two arms measure different objects and are deliberately NOT pooled."),
        "screen_source": screen_files[-1].name,
        "determinant_completeness_arm": amr,
        "position_novelty_arm": tgt,
        "amr_arm_known_gap_recall": {
            "n_known_gaps": len(KNOWN_AMR_GAPS),
            "n_recovered": sum(1 for c in amr if c.get("known_gap_recovered")),
            "note": ("One known gap exists in this arm, so this is a single case, NOT a rate. It "
                     "bounds nothing about gaps that have never been independently confirmed."),
        },
        "amr_arm_specificity": amr_arm_specificity(amr),
    }

    print(f"\nincumbent: {INCUMBENT['metric']} = {INCUMBENT['value']} (lift {INCUMBENT['lift']})\n")
    print("determinant-completeness arm (per drug; families the deployed rule cannot represent)")
    print(f"  {'drug':16} {'uncounted':>9} {'raw-sig':>7} {'STRONG':>6} {'weak':>5}  {'baseS':>6}  known-gap")
    for c in amr:
        if c["status"] != "scored":
            print(f"  {c['drug']:16} {c['n_families_uncounted']:>9}  -- {c['status']}")
            continue
        print(f"  {c['drug']:16} {c['n_families_uncounted']:>9} {c['n_raw_signature']:>7} "
              f"{c['n_strong']:>6} {c['n_weak']:>5}  {c['base_s_rate']:>6.3f}  "
              f"{'RECOVERED' if c['known_gap_recovered'] else '-'}")
        for s in c["strong"]:
            e = s["evidence"]
            print(f"      STRONG  {e['symbol']:14} {e['carriers_labelled_r']}R/"
                  f"{e['carriers_labelled_s']}S  p={e['purity_surprise_p']:.2e}")

    print("\nposition-novelty arm (per drug; per-ISOLATE sensitivity -- NOT the same object as above)")
    print(f"  {'drug':16} {'sens':>6} {'FP':>6} {'lift':>6}  powered")
    for c in tgt:
        print(f"  {c['drug']:16} {c['sens_on_blindspot']:>6.3f} {c['fp_on_catalog_negative_S']:>6.3f} "
              f"{c['lift']:>6.2f}  {c['powered']}")

    sp = out["amr_arm_specificity"]
    u, c = sp["unconfirmed_drugs"], sp["confirmed_gap_drugs"]
    print(f"\nAMR-arm firing rate (predeclared categories, NOT a false-positive rate)")
    print(f"  confirmed-gap drugs   {c['n_drugs']} drug(s), {c['n_families_screened']:>4} families "
          f"-> {c['n_strong']} STRONG")
    print(f"  unconfirmed drugs     {u['n_drugs']} drug(s), {u['n_families_screened']:>4} families "
          f"-> {u['n_strong_after_correction']} STRONG "
          f"(raw signature fired {u['n_raw_signature_before_correction']}; "
          f"correction removed {u['correction_removed']})")
    print(f"  unassessable          {sp['predeclared_categories']['unassessable_no_labels']}")
    print("  A hit on an unconfirmed drug would be AMBIGUOUS between a false positive and a new")
    print("  gap -- both known gaps were invisible until independent labels arrived.")

    dest = WIKI / f"doubt_layer_per_cell_{out['generated']}.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote wiki/{dest.name}")
    print("Per cell by construction. The two arms are NOT pooled -- they measure different objects.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
