"""Screen a deep-research harvest through BOTH standing screens, and refuse to call anything promising
until it has cleared them.

WHY THIS EXISTS. A literature find is a LEAD, not a finding. This project has a recorded pattern of
candidates that looked promising in prose and fell the moment they were acted on (the `/innovate` sweep:
most survivors were refuted or retracted once measured), plus a sharper one the same week — a run named
"score gLM2-650M vs the curated baseline" as the highest-VOI move when `regime.screen_proposal` puts it in
`natural x molecular x zero_shot = LOSES_TO_CATALOG`, i.e. it had proposed re-running a recorded negative.

So every harvested candidate passes two INDEPENDENT screens, and they answer different questions:

  * `rejection_gates.screen_candidate` (G1-G10) -- does a usable LABEL exist, and is the decoder's rule
    scoreable against a genotype at all? G1 (wet-lab vs tool-derived) and G3 (assay reading vs collection
    context) are JUDGMENTS: they need a human evidence string plus an explicit assertion and REFUSE the
    whole verdict without one. They never default to pass.
  * `regime.screen_proposal` -- is this proposal's (population, endpoint, method) cell a RECORDED NEGATIVE?
    A dataset candidate may legitimately have no regime triple; a METHOD candidate must have one, because
    that is the screen which catches a re-run of a closed direction.

Neither screen is a gate on the user; both are gates on ME calling something promising.

HONEST SCOPE. A `CLEARS` verdict bounds the LABEL question only. It is not a build recommendation:
artifact reachability, regime fit and worth-doing are separate, and PEAR is the standing proof — it clears
every applicable gate and still was not buildable as specified (its data ships as serialized ggplot2
objects). Read-only; writes one artifact; exit 0 unless the harvest itself is malformed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from dna_decode.eval.regime import screen_proposal  # noqa: E402
from dna_decode.eval.rejection_gates import screen_candidate  # noqa: E402

# A locator is the difference between a citation and a rumour. These tokens mark a candidate whose
# locator the harvesting agent could NOT confirm resolves -- they are preserved, never silently dropped,
# because an unverified accession is the fabrication hazard this project guards hardest against.
UNVERIFIED_MARKERS = ("UNVERIFIED", "unverified", "TODO", "?")

PROMISING = "PROMISING"
LEAD_ONLY = "LEAD_ONLY"
REJECTED = "REJECTED"
UNSCREENABLE = "UNSCREENABLE"


def _locator_state(cand: dict) -> tuple[bool, str]:
    """Is this candidate's locator present and claimed-resolvable?"""
    loc = (cand.get("locator") or "").strip()
    if not loc:
        return False, "no locator — a candidate with no DOI/URL/accession cannot be verified at all"
    if any(m in loc for m in UNVERIFIED_MARKERS):
        return False, f"locator flagged unverified by the harvester: {loc!r}"
    return True, loc


def screen_one(cand: dict) -> dict:
    """Run both screens on one candidate and classify it. Never upgrades on missing evidence."""
    out = dict(cand)
    ok_loc, loc_note = _locator_state(cand)
    out["locator_ok"] = ok_loc
    out["locator_note"] = loc_note

    # --- label side (G1-G10) -------------------------------------------------------------------
    lab = cand.get("label_screen")
    if lab:
        res = screen_candidate(cand["name"], lab["intended_layer"], lab.get("evidence", {}))
        out["label_verdict"] = res.verdict
        out["label_reason"] = res.reason
        out["label_gates"] = {g.gate: g.verdict for g in res.gates}
        out["label_gate_reasons"] = {g.gate: g.reason for g in res.gates if g.verdict == "trip"}
    else:
        out["label_verdict"] = None
        out["label_reason"] = "no label_screen supplied — not screened, which is NOT a pass"
        out["label_gates"] = {}
        out["label_gate_reasons"] = {}

    # --- regime side ---------------------------------------------------------------------------
    reg = cand.get("regime_screen")
    if reg:
        r = screen_proposal(reg["population"], reg["endpoint"], reg["method"],
                            curated_catalog_exists=reg.get("curated_catalog_exists", False))
        out["regime_verdict"] = r.verdict
        out["regime_note"] = getattr(r, "evidence", "") or getattr(r, "note", "") or ""
    else:
        out["regime_verdict"] = None
        # A DATASET candidate legitimately has no regime triple; a METHOD candidate without one is
        # exactly the gLM2 hole, so the classifier below treats a method with no regime as UNSCREENABLE.
        out["regime_note"] = "no regime triple supplied"

    out["screen_class"] = _classify(out)
    return out


def _classify(o: dict) -> str:
    """PROMISING requires BOTH screens to have actually run and neither to condemn it."""
    if not o["locator_ok"]:
        return UNSCREENABLE
    if o["label_verdict"] in ("REJECTED", "INCOMPLETE") or o["regime_verdict"] in (
            "CLOSED_NEGATIVE", "LOSES_TO_CATALOG"):
        return REJECTED
    is_method = (o.get("kind") or "").lower() == "method"
    if o["label_verdict"] is None and not is_method:
        return LEAD_ONLY            # a dataset we could not screen for labels yet
    if is_method and o["regime_verdict"] is None:
        return UNSCREENABLE         # a method with no regime triple is the gLM2 hole
    if o["label_verdict"] == "CLEARS" or o["regime_verdict"] in ("WORKS", "OPEN",
                                                                 "REQUIRES_DECONFOUNDING"):
        return PROMISING
    return LEAD_ONLY


def screen_harvest(path: Path) -> dict:
    doc = json.loads(path.read_text(encoding="utf-8"))
    cands = doc.get("candidates", [])
    screened = [screen_one(c) for c in cands]
    counts: dict[str, int] = {}
    for s in screened:
        counts[s["screen_class"]] = counts.get(s["screen_class"], 0) + 1
    by_family: dict[str, dict[str, int]] = {}
    for s in screened:
        fam = s.get("family", "?")
        by_family.setdefault(fam, {})
        k = s["screen_class"]
        by_family[fam][k] = by_family[fam].get(k, 0) + 1
    return {
        "artifact": "deep_research_screen",
        "schema": "deep-research-screen-v1",
        "generated": doc.get("generated"),
        "window": doc.get("window"),
        "n_candidates": len(screened),
        "class_counts": counts,
        "by_family": by_family,
        "candidates": screened,
        "honest_limits": [
            "a CLEARS verdict bounds the LABEL question ONLY -- not artifact reachability, not regime "
            "fit, not worth-doing (PEAR clears every applicable gate and was still not buildable)",
            "G1 and G3 are JUDGMENTS carried on a human evidence string; they refuse without one and "
            "never default to pass",
            "a candidate with no label_screen is LEAD_ONLY, never PROMISING -- unscreened is not passed",
            "a METHOD candidate with no regime triple is UNSCREENABLE, because the regime screen is "
            "precisely what catches a re-run of a recorded negative",
        ],
    }


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--harvest", type=Path,
                    default=REPO / "wiki" / "deep_research_candidates_2026-09-28.json")
    ap.add_argument("--out", type=Path,
                    default=REPO / "wiki" / "deep_research_screen_2026-09-28.json")
    a = ap.parse_args(argv)

    if not a.harvest.exists():
        print(f"no harvest at {a.harvest} — nothing to screen", file=sys.stderr)
        return 2
    rep = screen_harvest(a.harvest)
    a.out.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")

    print(f"screened {rep['n_candidates']} candidates")
    for k in (PROMISING, LEAD_ONLY, REJECTED, UNSCREENABLE):
        if rep["class_counts"].get(k):
            print(f"  {k:12s} {rep['class_counts'][k]}")
    for fam, c in sorted(rep["by_family"].items()):
        print(f"  {fam}: {c}")
    if rep["class_counts"].get(PROMISING):
        print("\nPROMISING (label-screen cleared or regime-open) — still LEADS, not findings:")
        for s in rep["candidates"]:
            if s["screen_class"] == PROMISING:
                print(f"  [{s.get('family')}] {s['name']}")
                print(f"      label={s['label_verdict']} regime={s['regime_verdict']}")
                print(f"      {s.get('locator')}")
    # Display a repo-relative path when the output IS under the repo, else the absolute one. A bare
    # relative_to() raises on any --out outside the tree, which crashed the first validation run against
    # a scratchpad path after all the real work had already succeeded.
    try:
        shown = a.out.relative_to(REPO).as_posix()
    except ValueError:
        shown = a.out.as_posix()
    print(f"\nwrote {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
