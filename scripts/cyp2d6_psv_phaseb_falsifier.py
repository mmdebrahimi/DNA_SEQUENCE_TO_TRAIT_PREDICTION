"""CYP2D6 hybrid-identity Phase-B FALSIFIER — does the PSV signal separate at the FULL N?

Phase A proved the signal on n=1/type. The /brainstorm (R2) flagged that anchor-reproduction is
NECESSARY-BUT-NOT-SUFFICIENT: a caller can match a few anchors while failing on pure-dup / deletion /
low-depth / mixed-copy. This runs the Phase A evidence extractor over the FULL structural panel (all 13
committed hybrids + 25 non-hybrids + the fetchable Cyrius *13 anchor HG01161) and asks, per allele:
  - do the hybrids get their CORRECT directional/dip SIGNAL? (per-allele sensitivity)
  - do the non-hybrids (normal / pure-dup / deletion) stay FLAT? (specificity — the confound test)

This is DIAGNOSTIC (the coarse Phase-A signal, not the eventual Phase-B classifier). It is the go/no-go for
building the classifier: clean separation at full N -> build it; muddy -> depth-only presence stays the ceiling
and identity is documented as long-read-required.
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from cyp2d6_psv_evidence import load_psvs, sample_evidence  # noqa: E402

# expected signal BUCKET (must match cyp2d6_psv_evidence.identity_signal.split(' ')[0] exactly).
_EXPECT = [("*68", "directional_5p_high_3p_low"), ("*13", "directional_5p_low_3p_high"),
           ("*36", "exon9_tip_dip")]


def _expected_signal(truth: str) -> str | None:
    """The signal a sample's truth should produce (None = non-hybrid -> flat). *68 checked before *36
    (a *68 tandem may also carry *4/*10; the *68 breakpoint dominates the profile)."""
    for star, sig in _EXPECT:
        if star in truth:
            return sig
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="CYP2D6 hybrid-identity Phase-B full-N falsifier.")
    ap.add_argument("--pileup-dir", type=Path, required=True)
    ap.add_argument("--panel", type=Path, required=True, help="panel file: 'name kind truth cram' per line")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "cyp2d6_psv_phaseb_falsifier.json")
    args = ap.parse_args(argv)
    psvs = load_psvs()

    rows = []
    for line in args.panel.read_text(encoding="utf-8").splitlines():
        f = line.split()
        if len(f) < 3:
            continue
        name, kind, truth = f[0], f[1], f[2]
        if not (args.pileup_dir / f"{name}.d6.txt").exists():
            continue
        ev = sample_evidence(name, args.pileup_dir, psvs)
        exp = _expected_signal(truth)
        got = ev["identity_signal"]
        got_bucket = got.split(" ")[0]
        exp_hyb = exp is not None
        got_hyb = got_bucket != "flat_nonhybrid"
        correct = (got_bucket == exp) if exp_hyb else (not got_hyb)
        rows.append({"sample": name, "truth": truth, "expected_signal": exp or "flat_nonhybrid",
                     "got_signal": got_bucket, "n_callable": ev["n_callable"],
                     "five_prime_minus_three_prime": ev["five_prime_minus_three_prime"],
                     "exon9_tip_dip": ev["exon9_tip_dip"], "is_hybrid": exp_hyb, "correct": correct})

    # per-allele sensitivity + non-hybrid specificity
    by_allele = {}
    for star, _ in _EXPECT:
        hits = [r for r in rows if star in r["truth"]]
        by_allele[star] = {"n": len(hits), "correct_signal": sum(r["correct"] for r in hits),
                           "rate": round(sum(r["correct"] for r in hits) / len(hits), 3) if hits else None}
    non = [r for r in rows if not r["is_hybrid"]]
    spec = round(sum(r["correct"] for r in non) / len(non), 3) if non else None
    hyb = [r for r in rows if r["is_hybrid"]]
    sens = round(sum(r["correct"] for r in hyb) / len(hyb), 3) if hyb else None

    # verdict: hybrids get their signal (well-powered alleles) AND non-hybrids stay flat (the confound)
    well = [a for a in ("*68", "*36") if by_allele[a]["n"] >= 3]
    well_ok = all((by_allele[a]["rate"] or 0) >= 0.6 for a in well)
    go = bool(well) and well_ok and (spec is not None and spec >= 0.85)
    rep = {
        "schema": "cyp2d6-psv-phaseb-falsifier-v0", "analysis_date": datetime.date.today().isoformat(),
        "n_samples": len(rows), "n_hybrid": len(hyb), "n_nonhybrid": len(non),
        "overall_hybrid_signal_sensitivity": sens, "nonhybrid_specificity": spec,
        "per_allele": by_allele,
        "verdict": "GO_BUILD_CLASSIFIER" if go else "PARTIAL_DEPTH_ONLY_CEILING",
        "interpretation": ("Full-N separation: the well-powered hybrid alleles (*68 n>=4, *36 n>=8) get their "
                           "directional/exon9 signal AND non-hybrids (normal/dup/deletion) stay flat (the "
                           "confound test). GO -> the Phase-B abstaining classifier is justified. PARTIAL -> "
                           "some alleles muddy at scale; keep presence-only as the ceiling for those."),
        "rows": sorted(rows, key=lambda r: (not r["is_hybrid"], r["truth"])),
    }
    args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    print(f"Phase-B falsifier: {rep['verdict']}")
    print(f"  hybrid signal sensitivity {sens} (n={len(hyb)}); non-hybrid specificity {spec} (n={len(non)})")
    for star, d in by_allele.items():
        print(f"  {star}: {d['correct_signal']}/{d['n']} correct-signal ({d['rate']})")
    print(f"[falsifier -> {args.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
