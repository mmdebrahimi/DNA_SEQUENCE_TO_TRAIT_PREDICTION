"""Score the EXPERIMENTAL scorer-local TMP-SMX rule on the in-hand external cohorts (Oxford + Sci234).

TMP-SMX (co-trimoxazole) is a genuinely-new resistance MECHANISM (acquired folate-pathway sul + dfr)
distinct from the 6 frozen decoder cells. Its rule `(>=1 sul) AND (>=1 dfr)` is a shape the frozen
`DRUG_RULE` cannot express, so it lives in the NON-FROZEN overlay (`experimental_drug_rules`) and is
validated ONLY here on the external-validation arm — never the frozen report card / shipped surface.

Pure join of each cohort's OWN genotype + measured co-trimoxazole MIC, scored by the overlay AND rule.
No download / assembly / Docker. Emits ONE top-level `external-validation-v1` artifact PER cohort,
branded EXPERIMENTAL_SCORED / scorer_local, with the 4 genotype strata (sul+dfr / sul-only / dfr-only /
neither). The strata-reproduction gate (sul+dfr is the highest-R stratum AND sul-only R-rate < 0.5 on
BOTH cohorts) is what separates a SCORED cell from INDETERMINATE — the real de-risk that the AND rule
generalizes, not just an aggregate sens/spec.

Binary (MIC>=4 R / <=2 S) is the PRIMARY metric for clean measured MIC; the 4x-margin strict tier
(HIGH_R needs MIC>=16, HIGH_S needs MIC<=0.5) over-excludes co-trimoxazole and is S-degenerate on Oxford
(min MIC 1) — reported but not headline (matches the Oxford/Sci234 binary-primary precedent).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.experimental_drug_rules import (
    DRUG, RULE_SCOPE, RULE_STATUS, cotrimoxazole_tier, tmp_smx_call,
)
from scripts.independent_cohort_validate import _conf
import scripts.oxford_score as oxford_score
import scripts.sci234_score as sci234_score

REGISTRY_ORGANISM = "Escherichia_coli_Shigella"
_STRICT_R, _STRICT_S = "HIGH_R", "HIGH_S"
_RELAXED_R, _RELAXED_S = {"HIGH_R", "DECISIVE_R"}, {"HIGH_S", "DECISIVE_S"}


def _mic(cell) -> float | None:
    s = str(cell)
    for op in ("<=", ">=", ">", "<", "="):
        s = s.replace(op, "")
    try:
        return float(s.strip())
    except ValueError:
        return None


def load_sci234():
    """-> (genotype {key:[gene symbols]}, trisul_mic {key:float})."""
    import openpyxl
    _qrdr, geno_acq, _mics = sci234_score.load_cohort()
    wb = openpyxl.load_workbook(sci234_score.SUPPL, read_only=True)
    it = wb["Supplementary data 2"].iter_rows(values_only=True)
    h = [str(c) for c in next(it)]
    ti = h.index("TRISUL")
    trisul = {}
    for r in it:
        m = _mic(r[ti]) if ti < len(r) else None
        if m is not None:
            trisul[str(r[0]).strip()] = m
    return geno_acq, trisul


def load_oxford():
    """-> (genotype {guuid:[Element symbols]}, cotrim_mic {guuid:float})."""
    base = Path("data/raw/oxford")
    header, by_guuid = oxford_score.group_amrfinder(base / "amrfinder.tsv")
    cols = header.split("\t")
    esym = cols.index("Element symbol")
    geno = {g: [ln.split("\t")[esym].strip() for ln in lines if len(ln.split("\t")) > esym]
            for g, lines in by_guuid.items()}
    cot = {}
    for r in csv.DictReader(open(base / "main_data.csv", encoding="utf-8")):
        up = (r.get("Cotrim_upper") or "").strip()
        if up in ("", "NA"):
            continue
        try:
            cot[r["guuid"].strip()] = 2.0 ** float(up)
        except ValueError:
            continue
    return geno, cot


def _binary(mic: float) -> str | None:
    if mic >= 4:
        return "R"
    if mic <= 2:
        return "S"
    return None


def score_cohort(geno: dict, mics: dict) -> dict:
    strata = {k: {"R": 0, "S": 0} for k in ("sul+dfr", "sul-only", "dfr-only", "neither")}
    strict_pairs, relaxed_pairs, binary_pairs = [], [], []
    for key, mic in mics.items():
        if key not in geno:
            continue
        call = tmp_smx_call(geno[key])
        pred = call["prediction"]
        has_sul, has_dfr = bool(call["matched_sul"]), bool(call["matched_dfr"])
        stratum = ("sul+dfr" if (has_sul and has_dfr) else
                   "dfr-only" if has_dfr else "sul-only" if has_sul else "neither")
        b = _binary(mic)
        if b is not None:
            binary_pairs.append((pred, 1 if b == "R" else 0))
            strata[stratum]["R" if b == "R" else "S"] += 1
        tier = cotrimoxazole_tier([mic])
        if tier in (_STRICT_R, _STRICT_S):
            strict_pairs.append((pred, 1 if tier == _STRICT_R else 0))
        if tier in (_RELAXED_R | _RELAXED_S):
            relaxed_pairs.append((pred, 1 if tier in _RELAXED_R else 0))
    strata_tbl = {}
    for k, v in strata.items():
        n = v["R"] + v["S"]
        strata_tbl[k] = {"n": n, "R": v["R"], "S": v["S"],
                         "r_rate": round(v["R"] / n, 3) if n else None}
    rates = {k: t["r_rate"] for k, t in strata_tbl.items() if t["r_rate"] is not None}
    sd = strata_tbl["sul+dfr"]["r_rate"]
    so = strata_tbl["sul-only"]["r_rate"]
    reproduced = bool(rates and sd is not None and sd == max(rates.values())
                      and (so is None or so < 0.5))
    return {"strict": _conf(strict_pairs), "relaxed": _conf(relaxed_pairs),
            "binary": _conf(binary_pairs), "strata": strata_tbl, "strata_reproduced": reproduced}


def build_artifact(cohort: str, cells: dict, run_id: str, independence: str, leakage: str) -> dict:
    headline = "SCORED" if cells["strata_reproduced"] else "INDETERMINATE"
    return {
        "_schema": "external-validation-v1", "date": _date.today().isoformat(),
        "cohort": cohort, "organism": REGISTRY_ORGANISM, "drug": DRUG, "run_id": run_id,
        "strict": cells["strict"], "relaxed": cells["relaxed"], "binary": cells["binary"],
        "strata": cells["strata"], "strata_reproduced": cells["strata_reproduced"],
        "headline": headline,
        "rule_status": RULE_STATUS, "rule_scope": RULE_SCOPE, "not_in_shipped_surface": True,
        "rule_text": tmp_smx_call([])["rule_text"],
        "primary_metric": ("binary (clean measured MIC; strict HIGH_S needs MIC<=0.5 -> over-excludes "
                           "co-trimoxazole, S-degenerate on Oxford)"),
        "evidence_tier": "external_clinical_experimental",
        "independence_tier": independence,
        "leakage_status": "DISJOINT_PROJECT_LEVEL", "leakage_evidence": leakage,
        "fidelity_caveat": ("acquired-gene AND rule; folP/folA target POINT-mutation TMP-R is invisible "
                            "to it (the sul+dfr-negative R isolates). Genotype callers (Oxford AMRFinder / "
                            "Sci234 ResFinder-style) share the curated acquired-gene vocabulary -> "
                            "independence is of the pipeline, not the caller class."),
        "note": "EXPERIMENTAL scorer-local rule (dna_decode/data/experimental_drug_rules.py); NOT in the "
                "frozen DRUG_RULE / shipped_decoder_surface. Direct supplement/AMRFinder join, no Docker.",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", default=f"tmpsmx{_date.today().isoformat().replace('-', '')}")
    a = ap.parse_args()
    out_paths = []
    for cohort, loader, indep, leak in [
        ("sci234", load_sci234,
         "Sci Rep 2023 234-isolate E. coli (Spain): own ResFinder/PointFinder-style genotype + measured "
         "TRISUL MIC, joined by study Key.",
         "PRJNA854358 deposits 0 NCBI assemblies (reads-only); decoder tuning is 100% GCA/GCF."),
        ("oxford", load_oxford,
         "Oxford ecoli_mic_arg (UK): own AMRFinder genotype + measured Cotrim MIC, joined by guuid.",
         "PRJNA604975 + PRJNA1007570 deposit 0 NCBI assemblies (reads-only); tuning is 100% GCA/GCF."),
    ]:
        geno, mics = loader()
        cells = score_cohort(geno, mics)
        art = build_artifact(cohort, cells, a.run_id, indep, leak)
        b, st = art["binary"], art["strict"]
        print(f"\n{cohort} [{art['headline']}] TMP-SMX: "
              f"binary n={b['n_scored']} acc={b['acc']} sens={b['sens']} spec={b['spec']} "
              f"(R{b['tp']+b['fn']}/S{b['tn']+b['fp']}) | strict n={st['n_scored']} sens={st['sens']} spec={st['spec']}")
        print(f"  strata: " + " | ".join(
            f"{k} n={v['n']} R-rate={v['r_rate']}" for k, v in art["strata"].items()))
        print(f"  strata_reproduced={art['strata_reproduced']}")
        p = Path(f"wiki/external_validation_{cohort}tmpsmx_{a.run_id}_{_date.today().isoformat()}.json")
        p.write_text(json.dumps(art, indent=2), encoding="utf-8")
        out_paths.append(p)
        print(f"  artifact -> {p}")
    print("\nroll up with: python scripts/build_external_validation_report.py "
          f"--run-id {a.run_id} --no-clonality")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
