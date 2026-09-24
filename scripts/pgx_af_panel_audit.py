"""Audit every asserted PGx allele-frequency claim against the 1000G panel ON DISK.

WHY THIS EXISTS (a real defect, found 2026-09-24)
-------------------------------------------------
`dna_decode/pgx/cyp4f2.py` asserts the *3 allele is "~29% EUR / ~79% EAS", and the registry contract cites
"AF-corroborated" as one of this cell's validation legs. Measured over the 3202-sample panel that ships
beside it, EAS *3 is **0.2128**, not 0.79 -- and 1 - 0.2128 = 0.7872, i.e. the asserted figure matches the
REFERENCE-allele frequency to within 0.3 percentage points. (That is consistent with an allele-polarity
slip; this audit measures frequencies and does NOT claim to establish how the number got there.)

The existing `scripts/pgx_af_corroboration.py` did not catch it, and the reason is structural, not a typo:

  * it checks ONE population per variant (CYP4F2: EUR only -- the number that is CORRECT), while the
    docstrings assert TWO;
  * its AF values are HARDCODED from a 2026-07-07 Ensembl fetch, so by default it compares a hand-entered
    number against a hand-entered band -- literature vs literature, deriving nothing at run time;
  * it never reads the 1000G panel sitting in data/pgx_1000g/, which is a free, local, per-population
    source for exactly these frequencies.

This audit closes that: EVERY asserted population is checked, and every frequency is DERIVED from the local
panel. It is the same "derive, don't assert" discipline `scripts/project_status.py` applies to scope.

REFUSES RATHER THAN PASSES. A claim whose variant is not in any local VCF is reported
UNMEASURABLE_NO_LOCAL_PANEL -- never a pass. If nothing at all was measurable the audit exits non-zero
rather than printing a clean-looking report over zero checks (the vacuous-negative trap).

SCOPE: this bounds whether a SHIPPED FREQUENCY CLAIM matches this panel. It is not a phenotype validation,
and a matching frequency does not make a cell correct.
"""
from __future__ import annotations

import datetime
import json
import sys
from collections import Counter
from pathlib import Path

from scripts.ugt1a1_getrm_concordance import _dose, load_ped

REPO = Path(__file__).resolve().parent.parent
VCF_DIR = REPO / "data" / "pgx_1000g"

TOLERANCE = 0.10          # percentage-point band a "~N%" claim must land within

# Every asserted population-frequency claim on the PGx surface, with WHERE it is asserted so a failure
# points at the file to fix. `asserted` maps superpopulation -> claimed ALT (functional-allele) frequency.
CLAIMS = [
    {"gene": "CYP4F2", "allele": "*3", "rsid": "rs2108622", "chrom": "19", "pos": 15879621,
     "ref": "C", "alt": "T", "vcf": "cyp4f2_1000g.vcf",
     "asserted": {"EUR": 0.29, "EAS": 0.79},
     "asserted_in": "dna_decode/pgx/cyp4f2.py docstring + cell_registry CYP4F2 contract",
     "note": "warfarin dose-up allele"},
    {"gene": "UGT1A1", "allele": "*80", "rsid": "rs887829", "chrom": "2", "pos": 233759924,
     "ref": "C", "alt": "T", "vcf": "ugt1a1_1000g.vcf",
     "asserted": {"EUR": 0.30},
     "asserted_in": "dna_decode/pgx/runner.py UGT1A1 block + pgx_af_corroboration.py",
     "note": "*28 LD-tag; the claim is that EUR AF == the *28 frequency"},
    {"gene": "UGT1A1", "allele": "*6", "rsid": "rs4148323", "chrom": "2", "pos": 233760498,
     "ref": "G", "alt": "A", "vcf": "ugt1a1_1000g.vcf",
     "asserted": {"EAS": 0.14},
     "asserted_in": "scripts/pgx_af_corroboration.py",
     "note": "decreased function, EAS-common"},
    {"gene": "CYP2B6", "allele": "*6-proxy", "rsid": "rs3745274", "chrom": "19", "pos": 41006936,
     "ref": "G", "alt": "T", "vcf": "cyp2b6_1000g.vcf",
     "asserted": {"ALL": 0.320},
     "asserted_in": "dna_decode/pgx/cyp2b6_catalog.py PROVENANCE block",
     "note": "global (not per-superpopulation) claim: '*6+*9 combined global frequency'"},
    # --- no local panel: reported UNMEASURABLE, never a pass ---
    {"gene": "ABCG2", "allele": "141K", "rsid": "rs2231142", "chrom": "4", "pos": 88131171,
     "ref": "G", "alt": "T", "vcf": "abcg2_1000g.vcf",
     "asserted": {"EUR": 0.09, "EAS": 0.29},
     "asserted_in": "dna_decode/pgx/abcg2.py docstring + cell_registry ABCG2 contract",
     "note": "rosuvastatin; NO local VCF -- fetch the region to measure"},
    {"gene": "NUDT15", "allele": "*3", "rsid": "rs116855232", "chrom": "13", "pos": 48037782,
     "ref": "C", "alt": "T", "vcf": "nudt15_1000g.vcf",
     "asserted": {"EAS": 0.10},
     "asserted_in": "dna_decode/pgx/nudt15_catalog.py + runner.py NUDT15 block",
     "note": "thiopurine toxicity; NO local VCF -- fetch the region to measure"},
]


def panel_frequencies(vcf: Path, pos: int, ref: str, alt: str, pop: dict[str, str]) -> dict | None:
    """Per-superpopulation ALT frequency for one variant, derived from the panel. None if absent."""
    alt_n, tot_n = Counter(), Counter()
    found = False
    with vcf.open(encoding="utf-8") as fh:
        order: list[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                order = f[9:]
                continue
            if int(f[1]) != pos or f[3] != ref or f[4] != alt:
                continue
            found = True
            for name, gt in zip(order, f[9:]):
                d = _dose(gt)
                if d is None:
                    continue
                p = pop.get(name)
                if not p:
                    continue
                alt_n[p] += d
                tot_n[p] += 2
            break
    if not found or not tot_n:
        return None
    out = {p: round(alt_n[p] / tot_n[p], 4) for p in sorted(tot_n)}
    out["ALL"] = round(sum(alt_n.values()) / sum(tot_n.values()), 4)
    out["_n_samples"] = sum(tot_n.values()) // 2
    return out


def audit(claims=CLAIMS, tolerance: float = TOLERANCE) -> dict:
    pop = load_ped()
    results = []
    for c in claims:
        vcf = VCF_DIR / c["vcf"]
        measured = panel_frequencies(vcf, c["pos"], c["ref"], c["alt"], pop) if vcf.exists() else None
        if measured is None:
            results.append({**{k: c[k] for k in ("gene", "allele", "rsid", "asserted", "asserted_in",
                                                 "note")},
                            "status": "UNMEASURABLE_NO_LOCAL_PANEL",
                            "reason": f"{c['vcf']} absent from data/pgx_1000g/"
                                      if not vcf.exists() else
                                      f"variant {c['chrom']}:{c['pos']} {c['ref']}>{c['alt']} not in "
                                      f"{c['vcf']}",
                            "checks": {}})
            continue
        checks = {}
        for p, asserted in c["asserted"].items():
            m = measured.get(p)
            complement = None if m is None else round(1 - m, 4)
            checks[p] = {
                "asserted": asserted, "measured": m,
                "delta": None if m is None else round(m - asserted, 4),
                "within_tolerance": m is not None and abs(m - asserted) <= tolerance,
                # a claim that matches 1-ALT far better than ALT is the allele-polarity signature
                "matches_complement_better": (
                    m is not None and abs(asserted - complement) < abs(asserted - m)),
            }
        results.append({**{k: c[k] for k in ("gene", "allele", "rsid", "asserted", "asserted_in", "note")},
                        "status": "OK" if all(v["within_tolerance"] for v in checks.values())
                                  else "CLAIM_OUT_OF_BAND",
                        "measured_all_populations": {k: v for k, v in measured.items()
                                                     if not k.startswith("_")},
                        "n_samples": measured["_n_samples"], "checks": checks})

    measurable = [r for r in results if r["status"] != "UNMEASURABLE_NO_LOCAL_PANEL"]
    bad = [r for r in results if r["status"] == "CLAIM_OUT_OF_BAND"]
    return {
        "schema": "pgx-af-panel-audit-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "method": "every asserted population frequency DERIVED from the local 1000G panel in "
                  "data/pgx_1000g/; a claim with no local panel is UNMEASURABLE, never a pass",
        "tolerance": tolerance,
        "n_claims": len(results),
        "n_measurable": len(measurable),
        "n_unmeasurable": len(results) - len(measurable),
        "n_out_of_band": len(bad),
        "verdict": ("NO_CLAIM_MEASURABLE" if not measurable
                    else "ALL_MEASURABLE_CLAIMS_HOLD" if not bad else "CLAIMS_OUT_OF_BAND"),
        "scope": "bounds whether a shipped FREQUENCY claim matches this panel; NOT a phenotype validation",
        "claims": results,
    }


def render_md(rep: dict) -> str:
    L = [f"# PGx allele-frequency panel audit ({rep['analysis_date']})", "",
         f"Every asserted population frequency on the PGx surface, **derived from the 1000G panel on "
         f"disk** rather than compared against a hand-entered value. Tolerance +/-{rep['tolerance']}.", "",
         f"**{rep['verdict']}** — {rep['n_measurable']} of {rep['n_claims']} claims measurable, "
         f"{rep['n_out_of_band']} out of band, {rep['n_unmeasurable']} unmeasurable (no local panel).", "",
         "| gene | allele | pop | asserted | measured | delta | verdict |", "|---|---|---|---|---|---|---|"]
    for r in rep["claims"]:
        if r["status"] == "UNMEASURABLE_NO_LOCAL_PANEL":
            L.append(f"| {r['gene']} | `{r['allele']}` | — | "
                     f"{', '.join(f'{k} {v}' for k, v in r['asserted'].items())} | — | — | "
                     f"UNMEASURABLE ({r['reason']}) |")
            continue
        for p, c in r["checks"].items():
            flag = "ok" if c["within_tolerance"] else "**OUT OF BAND**"
            if c.get("matches_complement_better"):
                flag += " — matches 1−ALT better than ALT"
            L.append(f"| {r['gene']} | `{r['allele']}` | {p} | {c['asserted']} | {c['measured']} | "
                     f"{c['delta']} | {flag} |")
    L += ["", "## Where each claim is asserted", ""]
    for r in rep["claims"]:
        L.append(f"- **{r['gene']} `{r['allele']}`** ({r['rsid']}) — {r['asserted_in']}")
    L += ["", f"_{rep['scope']}._", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    rep = audit()
    stamp = rep["analysis_date"]
    (REPO / "wiki" / f"pgx_af_panel_audit_{stamp}.json").write_text(json.dumps(rep, indent=2),
                                                                   encoding="utf-8")
    (REPO / "wiki" / f"pgx_af_panel_audit_{stamp}.md").write_text(render_md(rep), encoding="utf-8")
    for r in rep["claims"]:
        if r["status"] == "UNMEASURABLE_NO_LOCAL_PANEL":
            print(f"{r['gene']:8} {r['allele']:10} UNMEASURABLE ({r['reason']})")
            continue
        for p, c in r["checks"].items():
            tag = "ok " if c["within_tolerance"] else "OUT"
            extra = "  <- matches 1-ALT better" if c.get("matches_complement_better") else ""
            print(f"{r['gene']:8} {r['allele']:10} {p:4} asserted={c['asserted']:<6} "
                  f"measured={c['measured']:<7} delta={c['delta']:<8} {tag}{extra}")
    print(f"VERDICT {rep['verdict']}  (measurable {rep['n_measurable']}/{rep['n_claims']}, "
          f"out-of-band {rep['n_out_of_band']})")
    if rep["verdict"] == "NO_CLAIM_MEASURABLE":
        print("REFUSED: no claim was measurable -- a clean report over zero checks is not a result.",
              file=sys.stderr)
        return 3
    return 1 if rep["n_out_of_band"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
