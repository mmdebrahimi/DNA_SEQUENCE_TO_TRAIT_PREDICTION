"""The two archives DISAGREE about Klebsiella, and the disagreement is testable.

WHY THIS EXISTS. On 2026-09-03 we published a Klebsiella over-call for the gentamicin `rmt` rescue
(BV-BRC: 58R/64S, PPV 0.475), then found it fails this project's own source-diversity bar (66.4% from
one study). The obvious next question -- what does the OTHER archive say about Klebsiella specifically?
-- had never been asked, because the PD sweep was read as a pooled number and the per-organism split
was sitting unexamined in the committed artifact.

It says the opposite. NCBI-PD holds 53 Klebsiella `rmt` carriers with a gentamicin label and ZERO
susceptible ones, across 12 BioProjects with a largest-source share of 0.26 -- which PASSES the bar
that the BV-BRC Klebsiella cell fails.

THE ASYMMETRY IS PRINCIPLED, NOT CONVENIENT, and that has to be stated plainly because it is the part
a reader should be most suspicious of. Both archives had a dominant source, and both got the SAME
pre-registered `aac(3)` control:
  - PD's PRJNA1322038 FAILED it (calls aac(3) carriers R 2% vs 97% elsewhere) -> excluded as a label
    artifact. All 42 of PD's susceptible Klebsiella carriers live in that project.
  - BV-BRC's pmid 36801013 PASSED it (99% vs 83%) -> kept.
So each archive's counter-examples were judged by the same rule, and they came out differently.

WHAT IS TESTED HERE. If the true Klebsiella PPV were the BV-BRC estimate (0.475), observing 53 carriers
with 0 susceptible is a binomial tail event. That is computed exactly, with no distributional
assumption beyond independence -- which is itself flagged as the main limit, since carriers within a
BioProject are not independent (clonality). The test is reported as an inconsistency measure, NOT as a
proof the rule is safe.

Offline: reads two committed artifacts. Writes wiki/rmt_klebsiella_archive_conflict_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from datetime import date as _date
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from source_diverse_validate import MAX_SOURCE_SHARE  # noqa: E402

PD_ARTIFACT_PROJECT = "PRJNA1322038"   # failed its aac(3) control -> excluded
BVBRC_DOMINANT_PMID = "36801013"       # passed its aac(3) control -> kept


def concentration(records: list[dict], key: str) -> dict:
    """Distinct sources + largest-source share, and whether it clears the project's own bar."""
    if not records:
        return {"n": 0, "n_sources": 0, "largest_share": None, "passes_bar": None}
    c = collections.Counter((r.get(key) or "NO_SOURCE") for r in records)
    top, n = c.most_common(1)[0]
    share = n / len(records)
    return {"n": len(records), "n_sources": len(c), "largest_source": top,
            "largest_share": share, "passes_bar": share <= MAX_SOURCE_SHARE}


def exact_binomial_tail_zero(n: int, p_success: float) -> float:
    """P(0 successes in n trials) where success = 'carrier is susceptible'.

    Computed with Fraction so a tiny probability is not silently flushed to 0.0 by float underflow --
    reporting 0.0 would overstate the result as impossible rather than merely extreme.
    """
    q = Fraction(1) - Fraction(p_success).limit_denominator(10**6)
    return float(q ** n)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pd", type=Path, default=ROOT / "wiki" / "gentamicin_rmt_specificity_hunt.json")
    ap.add_argument("--bvbrc", type=Path, default=ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" /
                    f"rmt_klebsiella_archive_conflict_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    # ---- NCBI-PD, Klebsiella only -------------------------------------------------------------
    pd = json.loads(a.pd.read_text(encoding="utf-8"))
    grp = next((g for g in pd["per_group"] if g["group"] == "Klebsiella"), None)
    if grp is None:
        print("no Klebsiella group in the PD artifact", file=sys.stderr)
        return 2
    pd_R_all = grp.get("rmt_R_records") or []
    pd_S_all = grp.get("rmt_S_records") or []

    def not_artifact(r): return (r.get("bioproject_acc") or "") != PD_ARTIFACT_PROJECT
    pd_R = [r for r in pd_R_all if not_artifact(r)]
    pd_S = [r for r in pd_S_all if not_artifact(r)]
    pd_conc = concentration(pd_R + pd_S, "bioproject_acc")
    pd_ppv = len(pd_R) / (len(pd_R) + len(pd_S)) if (pd_R or pd_S) else None

    # ---- BV-BRC, Klebsiella only --------------------------------------------------------------
    bv = json.loads(a.bvbrc.read_text(encoding="utf-8"))

    def src(h):
        p = h.get("pmid")
        if isinstance(p, list):
            p = ",".join(str(x) for x in p) if p else ""
        return str(p) if p else "NO_PMID"

    kleb = [h for h in bv["all_hits"]
            if "klebsiella" in (h.get("genome_name") or "").lower()
            and h.get("phenotype") in ("Resistant", "Susceptible")]
    bv_R = [h for h in kleb if h["phenotype"] == "Resistant"]
    bv_S = [h for h in kleb if h["phenotype"] == "Susceptible"]
    bv_conc = concentration(kleb, "__src")  # placeholder; recompute with the normalised key below
    c = collections.Counter(src(h) for h in kleb)
    top, n = c.most_common(1)[0]
    bv_conc = {"n": len(kleb), "n_sources": len(c), "largest_source": top,
               "largest_share": n / len(kleb), "passes_bar": (n / len(kleb)) <= MAX_SOURCE_SHARE}
    bv_ppv = len(bv_R) / len(kleb)

    print(f"Klebsiella rmt carriers, per archive (bar = largest-source share <= {MAX_SOURCE_SHARE:.0%})\n")
    print(f"  NCBI-PD  (artifact project excluded): {len(pd_R)}R / {len(pd_S)}S  PPV="
          f"{pd_ppv:.3f}  sources={pd_conc['n_sources']}  largest={pd_conc['largest_share']:.3f}  "
          f"{'PASSES' if pd_conc['passes_bar'] else 'FAILS'}")
    print(f"  BV-BRC   (dominant study kept)      : {len(bv_R)}R / {len(bv_S)}S  PPV="
          f"{bv_ppv:.3f}  sources={bv_conc['n_sources']}  largest={bv_conc['largest_share']:.3f}  "
          f"{'PASSES' if bv_conc['passes_bar'] else 'FAILS'}")

    # ---- is PD's zero consistent with BV-BRC's PPV? --------------------------------------------
    p_sus = 1.0 - bv_ppv
    p_zero = exact_binomial_tail_zero(len(pd_R) + len(pd_S), p_sus) if pd_S == [] else None
    # CLONALITY-CONSERVATIVE version. Carriers inside one BioProject may be near-clonal, so the
    # per-isolate test is optimistic. Collapsing to ONE vote per BioProject is the pessimistic bound:
    # it throws away all within-project replication and asks only whether 12 INDEPENDENT sources all
    # came back clean. Report BOTH; the conservative one is the one to quote.
    n_eff = pd_conc["n_sources"]
    p_zero_eff = exact_binomial_tail_zero(n_eff, p_sus) if pd_S == [] else None
    if p_zero is not None:
        print(f"\n  If the true Klebsiella PPV were BV-BRC's {bv_ppv:.3f}, seeing 0 susceptible among")
        print(f"    {len(pd_R)} PD carriers (per-isolate, assumes independence) : p = {p_zero:.3e}")
        print(f"    {n_eff} PD BioProjects (one vote each, clonality-safe)      : p = {p_zero_eff:.3e}"
              f"   <- quote THIS one")

    conflict = (pd_conc["passes_bar"] and not bv_conc["passes_bar"]
                and pd_ppv is not None and abs(pd_ppv - bv_ppv) > 0.2)
    if conflict:
        verdict = "ARCHIVES_CONFLICT_DIVERSE_SIDE_SHOWS_NO_OVERCALL"
        why = (f"the two archives disagree about Klebsiella, and the one that CLEARS the "
               f"source-diversity bar ({pd_conc['n_sources']} BioProjects, largest "
               f"{pd_conc['largest_share']:.0%}) shows {len(pd_R)}R/{len(pd_S)}S -- no over-call. The "
               f"archive reporting the over-call FAILS that bar ({bv_conc['largest_share']:.0%} from "
               f"one study). The over-call is NOT retracted, but it is now contradicted by the more "
               "source-diverse evidence and rests on a single study.")
    else:
        verdict = "ARCHIVES_CONSISTENT_OR_UNDECIDABLE"
        why = "the archives do not conflict in the way the diversity bar would make decisive."
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "rmt-klebsiella-archive-conflict-v1",
        "question": "do NCBI-PD and BV-BRC agree about the gentamicin rmt rescue in Klebsiella?",
        "bar_imported_from": "scripts/source_diverse_validate.py::MAX_SOURCE_SHARE",
        "bar": MAX_SOURCE_SHARE,
        "ncbi_pd": {"R": len(pd_R), "S": len(pd_S), "ppv": pd_ppv,
                    "excluded_project": PD_ARTIFACT_PROJECT,
                    "S_before_exclusion": len(pd_S_all),
                    "concentration": pd_conc},
        "bvbrc": {"R": len(bv_R), "S": len(bv_S), "ppv": bv_ppv,
                  "dominant_source_kept": BVBRC_DOMINANT_PMID,
                  "concentration": bv_conc},
        "consistency_test": {
            "question": "under BV-BRC's Klebsiella PPV, how likely is PD's zero-susceptible result?",
            "p_susceptible_under_bvbrc": p_sus,
            "n_pd_carriers": len(pd_R) + len(pd_S),
            "p_observing_zero_susceptible_per_isolate": p_zero,
            "per_isolate_caveat": "assumes independence between carriers, which is FALSE under "
                                  "clonality -- carriers within a BioProject may be near-clonal, so "
                                  "this OVERSTATES the evidence and should not be the quoted figure.",
            "n_effective_sources": n_eff,
            "p_observing_zero_susceptible_per_source": p_zero_eff,
            "quote_this": "p_observing_zero_susceptible_per_source -- one vote per BioProject, which "
                          "discards all within-project replication and is the clonality-safe bound.",
        },
        "verdict": verdict, "why": why,
        "why_the_exclusions_are_not_special_pleading": (
            "both archives had a dominant source and both got the SAME pre-registered aac(3) control. "
            f"PD's {PD_ARTIFACT_PROJECT} FAILED it (aac(3) carriers called R 2% vs 97% elsewhere) and "
            f"was excluded; BV-BRC's pmid {BVBRC_DOMINANT_PMID} PASSED it (99% vs 83%) and was kept. "
            "The asymmetry is the control's output, not a choice made after seeing which way it cut."
        ),
        "honest_limits": [
            "PD's 42 susceptible Klebsiella carriers ALL live in the excluded project, so PD's zero is "
            "a zero-after-exclusion. If that exclusion were wrong, PD would read 53R/42S = PPV 0.558 "
            "and would CORROBORATE the over-call rather than contradict it.",
            "Absence of susceptible carriers is weaker than presence of resistant ones; 53R/0S bounds "
            "the over-call but does not prove the rule safe in Klebsiella.",
            "Clonality is not corrected here. The binomial p-value assumes independent carriers and is "
            "therefore optimistic; a lineage-collapsed version would have a smaller effective N.",
            "Different label provenance: PD's AST_phenotypes vs BV-BRC's publication-curated MICs. A "
            "systematic difference between those label sources is not excluded.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
