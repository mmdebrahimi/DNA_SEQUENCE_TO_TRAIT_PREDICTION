"""`dna-tb` -- M. tuberculosis drug resistance from a VCF, via the pinned WHO catalogue.

    dna-tb --vcf isolate.vcf --drug rifampicin
    dna-tb --vcf isolate.vcf --drug isoniazid --json
    dna-tb --list-drugs

WHY THIS EXISTS: the decoder did, and nothing could reach it. `organism_rules/tb_amr.py` has
`score_rif` / `score_inh` / `score_drug`, the WHO catalogue loader verifies three sha256 pins, and
there are validated artifacts for twelve drugs -- but there was NO console script and ZERO
`cell_registry` cells, so the whole arm was invisible to the trust surface and unrunnable by a user.
`DrugCall` even carries `rule_status` / `rule_scope` / `catalogue_commit`, i.e. it was built to be a
contracted cell and never wired up. Same class as the 2026-09-01 HCMV registry gap, one step worse:
HCMV was at least CLI-routable.

SEPARATE CONSOLE SCRIPT, not folded into `dna-amr`, for three reasons: TB lives in the NON-frozen
`organism_rules` package rather than the frozen `DRUG_RULE` engine; its input is a VCF against H37Rv
rather than a genome FASTA or `--observed`; and `dna-flowering` already sets the
organism_rules-gets-its-own-entry-point precedent. This also keeps the five sha256-pinned frozen files
entirely untouched.

TWO EVIDENCE TIERS, AND THEY MUST NOT BE QUOTED AS ONE:
  * rifampicin + isoniazid have an INDEPENDENT number (EBI AMR-Portal, provenance-disjoint from the
    CRyPTIC set the WHO catalogue was partly built on).
  * the other ten have a CRyPTIC IN-DISTRIBUTION baseline only -- the catalogue was built partly from
    CRyPTIC, so scoring on it is a knowledge baseline, not validation.

AND FOR RIF/INH THE HEADLINE IS THE LINEAGE-COLLAPSED NUMBER, not the raw one. TB is heavily clonal;
raw counts one vote per isolate, so 2,845 isolates collapse to ~67 barcode lineages. Raw RIF
0.920/0.955 is clonality-INFLATED; the honest figure is 0.444/0.979.

REFUSES rather than guesses. An absent or pin-failing catalogue exits non-zero with an actionable
message -- it never returns a susceptible call, because "no determinant found" and "no catalogue
loaded" would otherwise be indistinguishable in the output.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

#: Drugs with an INDEPENDENT (provenance-disjoint) number, and that number -- lineage-collapsed,
#: because raw is clonality-inflated. Read from wiki/tb_independent_amr_portal_lineage_collapsed.json.
INDEPENDENT_DRUGS = {
    "rifampicin": {"lineage_sens": 0.444, "lineage_spec": 0.979,
                   "raw_sens": 0.920, "raw_spec": 0.955},
    "isoniazid": {"lineage_sens": 0.321, "lineage_spec": 0.972,
                  "raw_sens": 0.879, "raw_spec": 0.962},
}


def _catalogue_state(cat_dir: Path) -> tuple[bool, str]:
    """(usable, reason). Refuses on absent OR pin-failing -- a drifted catalogue is not the one the
    published numbers were measured against."""
    from dna_decode.data import tb_who_catalogue as cat

    if not cat.catalogue_available(cat_dir):
        return False, (f"WHO TB catalogue not found under {cat_dir}. The 37 MB master file is "
                       f"gitignored (only CHECKSUMS is tracked), so it must be fetched on this host "
                       f"before TB calls are possible.")
    try:
        pins = cat.verify_pins(cat_dir)
    except Exception as exc:
        return False, f"catalogue pin verification failed: {type(exc).__name__}: {exc}"
    bad = sorted(k for k, v in pins.items() if not v)
    if bad:
        return False, (f"catalogue pin MISMATCH on {bad} -- this is not the catalogue the published "
                       f"numbers were measured against, so a call would not be the validated rule.")
    return True, "ok"


def normalize_filter_column(vcf_text: str) -> tuple[str, bool]:
    """(normalized_text, rewrote_any). Rewrites FILTER `.` -> `PASS`, same transform the validated
    AMR-Portal runner applies (`run_tb_independent_amr_portal.py::_pass_mask`).

    WHY THIS IS NEEDED AND WHY IT IS SAFE. `parse_masked_calls` requires `FILTER==PASS`, which is right
    for the CRyPTIC **minos** VCFs it was written against. A minimap2/paftools.js VCF -- which is what
    the 2,845-isolate AMR-Portal arm actually produced -- writes FILTER `.` on every record. Without
    this, a VCF carrying 1,382 genuine `GT=1/1` variant rows parses to ZERO calls and every drug comes
    back S. Measured on a real isolate before writing this, not assumed.

    `.` in VCF means "no filters applied", so promoting it cannot smuggle a FAILING record through: a
    caller that does filter marks failures with a filter NAME (`MIN_DP`), never `.`. Only `.` is
    rewritten; a real PASS/filter vocabulary passes through untouched.

    THE CONSEQUENCE IS DISCLOSED, NOT HIDDEN: `.` also means there is NO caller-applied quality floor,
    so `parse_masked_calls`'s "PASS subsumes the MIN_DP/MIN_FRS floor" reasoning does not hold for such
    a VCF. The caller returns the flag so the record can say so.
    """
    out, rewrote = [], False
    for ln in vcf_text.splitlines():
        if ln.startswith("#"):
            out.append(ln)
            continue
        parts = ln.split("\t")
        if len(parts) > 6 and parts[6] == ".":
            parts[6] = "PASS"
            rewrote = True
            out.append("\t".join(parts))
        else:
            out.append(ln)
    return "\n".join(out), rewrote


def count_data_rows(vcf_text: str) -> int:
    """Data (non-header) rows, used to tell an EMPTY VCF from a VCF we failed to PARSE."""
    return sum(1 for ln in vcf_text.splitlines() if ln.strip() and not ln.startswith("#"))


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="dna-tb",
        description="M. tuberculosis drug resistance from a VCF against H37Rv, using the pinned WHO "
                    "mutation catalogue v2 (2023). Abstains on uncallable determinant positions and "
                    "REFUSES when the catalogue is absent or pin-mismatched.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", type=Path, help="masked VCF called against H37Rv (NC_000962.3)")
    ap.add_argument("--drug", default=None, help="drug name (see --list-drugs)")
    ap.add_argument("--regeno", type=Path, default=None,
                    help="optional regenotyped VCF for CALLABILITY (present-and-failed -> ABSTAIN). "
                         "Absent positions are treated as callable, per the ratified fail-only rule.")
    ap.add_argument("--catalogue-dir", type=Path, default=Path("data/raw/who_tb_catalogue"))
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--list-drugs", action="store_true",
                    help="print supported drugs and which have an INDEPENDENT number")
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = build_parser()
    a = ap.parse_args(argv)

    from dna_decode.data import tb_who_catalogue as cat
    from dna_decode.data.routable_drugs import all_routable_tb_drugs

    # ONE definition, shared with cell_registry's coverage test -- so a drug cannot be routable here
    # while the registry believes it is not and ships it with no evidence contract (the HCMV failure).
    longs = sorted(all_routable_tb_drugs())

    if a.list_drugs:
        usable, reason = _catalogue_state(a.catalogue_dir)
        print(f"catalogue: {'usable' if usable else 'UNUSABLE -- ' + reason}")
        print(f"\n{len(longs)} supported drugs. Evidence tier differs and must not be quoted as one:")
        for d in longs:
            ind = INDEPENDENT_DRUGS.get(d)
            if ind:
                print(f"   {d:<14} INDEPENDENT (AMR-Portal, provenance-disjoint): lineage-collapsed "
                      f"sens {ind['lineage_sens']} / spec {ind['lineage_spec']}")
            else:
                print(f"   {d:<14} CRyPTIC IN-DISTRIBUTION baseline only -- the WHO catalogue was "
                      f"built partly from CRyPTIC, so this is a knowledge baseline, NOT validation")
        print("\nFor rifampicin/isoniazid the LINEAGE-COLLAPSED number is the headline; the raw "
              "per-isolate figure (0.920/0.879 sens) is clonality-INFLATED because ~2,845 isolates "
              "collapse to ~67 barcode lineages.")
        return 0

    if not a.vcf or not a.drug:
        ap.error("--vcf and --drug are both required (or use --list-drugs)")

    drug = a.drug.strip().lower()
    if drug not in cat.DRUG_CATALOGUE_NAME:
        print(f"ERROR: {a.drug!r} is not in the WHO catalogue. Supported: {longs}", file=sys.stderr)
        return 2
    if not a.vcf.exists():
        print(f"ERROR: VCF not found: {a.vcf}", file=sys.stderr)
        return 2

    # REFUSAL comes BEFORE any scoring: a missing catalogue must never surface as a susceptible call.
    usable, reason = _catalogue_state(a.catalogue_dir)
    if not usable:
        print(f"CATALOGUE UNUSABLE: {reason}", file=sys.stderr)
        print("Refusing to score -- 'no determinant found' and 'no catalogue loaded' would be "
              "indistinguishable in the output.", file=sys.stderr)
        return 3

    from dna_decode.organism_rules import tb_amr, tb_vcf

    determinants = cat.load_determinants(drug, a.catalogue_dir)
    raw = a.vcf.read_text(encoding="utf-8", errors="replace")
    n_rows = count_data_rows(raw)
    text, rewrote = normalize_filter_column(raw)
    calls = tb_vcf.parse_masked_calls(text)

    # A VCF WITH DATA ROWS THAT PARSES TO ZERO CALLS IS A PARSE FAILURE, NOT A CLEAN GENOME.
    # Both states would otherwise print `PREDICTION: S`, which is the indistinguishable-failure
    # shape -- and it is exactly how the FILTER-vocabulary bug above presented: 1,382 real variant
    # rows, zero parsed calls, a confident S. Refusing here means a future caller shape that this
    # normalizer does not understand fails LOUDLY instead of reporting every isolate susceptible.
    if n_rows > 0 and not calls:
        print(f"PARSE FAILURE: {a.vcf.name} has {n_rows} data rows but ZERO non-reference calls were "
              f"parsed. That is a parser/format mismatch, not a susceptible genome -- refusing rather "
              f"than reporting S. Check the FILTER and GT columns against "
              f"tb_vcf.parse_masked_calls's contract (FILTER==PASS and a GT allele index >= 1).",
              file=sys.stderr)
        return 4

    regeno = (a.regeno.read_text(encoding="utf-8", errors="replace")
              if a.regeno and a.regeno.exists() else None)
    if regeno is not None:
        regeno, _ = normalize_filter_column(regeno)
    res = tb_amr.score_drug(drug, calls, determinants, regeno_text=regeno)

    ind = INDEPENDENT_DRUGS.get(drug)
    rec = {
        "organism": "Mycobacterium_tuberculosis",
        "drug": drug,
        "prediction": res.prediction,
        # `variant` already carries the gene prefix (`rpoB_p.Ser450Leu`), so joining with `gene` again
        # yields `rpoB_rpoB_p.Ser450Leu`. Use the variant alone, falling back for any that lacks it.
        "matched_determinants": [d.variant if d.variant.startswith(d.gene) else f"{d.gene}_{d.variant}"
                                 for d in res.matched],
        "coverage_scope": res.coverage_scope,
        "n_determinant_positions": res.n_determinant_positions,
        "n_uncallable_positions": res.n_uncallable_positions,
        "callability_assessed": res.callability_assessed,
        "rule_status": res.rule_status,
        "rule_scope": res.rule_scope,
        "input_type": res.input_type,
        "catalogue_commit": res.catalogue_commit,
        "n_nonref_calls": len(calls),
        "n_vcf_data_rows": n_rows,
        "filter_column_normalized": rewrote,
        **({"no_caller_quality_floor": "this VCF's FILTER column was '.' (no filters applied by the "
                                       "caller), so there is NO caller-applied depth/fraction floor "
                                       "behind these calls -- PASS does not subsume a quality gate here"}
           if rewrote else {}),
        "evidence": ({"tier": "independent_provenance_disjoint",
                      "headline_lineage_collapsed": {"sens": ind["lineage_sens"],
                                                     "spec": ind["lineage_spec"]},
                      "raw_per_isolate_CLONALITY_INFLATED": {"sens": ind["raw_sens"],
                                                             "spec": ind["raw_spec"]},
                      "source": "wiki/tb_independent_amr_portal_lineage_collapsed.json"}
                     if ind else
                     {"tier": "cryptic_in_distribution_knowledge_baseline",
                      "note": "the WHO catalogue was built partly from CRyPTIC, so a CRyPTIC-scored "
                              "number is in-distribution, NOT independent validation",
                      "source": f"wiki/tb_{drug[:3]}_cryptic_parquet_baseline_*.json"}),
        "does_not_support": [
            "a clinical treatment decision",
            "drugs outside the WHO catalogue's grade 1/2 determinant set",
            "an S call ruling out resistance by an UNCATALOGUED mechanism -- absence of a catalogued "
            "determinant is not absence of resistance",
        ],
    }
    if a.json:
        print(json.dumps(rec, indent=2))
        return 0

    print(f"ORGANISM   : Mycobacterium tuberculosis")
    print(f"DRUG       : {drug}")
    print(f"PREDICTION : {res.prediction}")
    print(f"determinants matched: {', '.join(rec['matched_determinants']) or '(none)'}")
    print(f"positions  : {res.n_determinant_positions} catalogued, "
          f"{res.n_uncallable_positions} uncallable"
          f"{'' if res.callability_assessed else ' (callability NOT assessed -- no --regeno)'}")
    print(f"calls      : {len(calls)} non-reference from {n_rows} VCF data rows")
    if rewrote:
        print("             NOTE: FILTER was '.' (no caller filters) -- normalized to PASS, so there "
              "is NO caller-applied quality floor behind these calls.")
    # `rule_status` is the library's label for how the RULE was DERIVED (a curated catalogue, hence
    # KNOWLEDGE_BASELINE). The `evidence` line below is about the COHORT it was SCORED on, which for
    # rifampicin/isoniazid is provenance-disjoint. Printing the two bare and adjacent reads as a
    # contradiction -- "KNOWLEDGE_BASELINE" directly above "INDEPENDENT" -- so each is labelled with
    # what it actually describes.
    print(f"rule deriv.: {res.rule_status} / {res.rule_scope} (how the RULE was built: a curated "
          f"catalogue)  catalogue {res.catalogue_commit[:8]}")
    if ind:
        print(f"cohort ev. : INDEPENDENT (AMR-Portal). Headline is LINEAGE-COLLAPSED "
              f"sens {ind['lineage_sens']} / spec {ind['lineage_spec']}; the raw per-isolate "
              f"{ind['raw_sens']}/{ind['raw_spec']} is clonality-INFLATED.")
    else:
        print(f"cohort ev. : CRyPTIC IN-DISTRIBUTION baseline only -- the catalogue was built partly "
              f"from CRyPTIC, so this is a knowledge baseline, NOT independent validation.")
    print("scope      : not a clinical decision; an S call does NOT rule out an uncatalogued "
          "resistance mechanism.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
