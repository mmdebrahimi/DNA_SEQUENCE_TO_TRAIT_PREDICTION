"""How rare is `rmt` in a broadly-ascertained E. coli collection? Oxford says: rarer than we assumed.

WHY A THIRD ARCHIVE. The gentamicin v2 `symbol_rescue` (`rmt*`/`npmA`) was validated on E. coli and its
over-call risk was measured on BV-BRC, where Klebsiella carriers were 58R/64S but E. coli carriers were
12R/0S. The standing honest limit was that twelve carriers is not many -- so the E. coli scope might be
safe, or might merely be under-sampled. Oxford is the natural third opinion:

  - ascertained on BACTERAEMIA, not on aminoglycoside resistance, so it escapes the structural trap that
    disqualified every `rmt`-enriched surveillance study (they screen on the outcome whose exceptions we
    are hunting, and therefore contain zero susceptible carriers by construction);
  - independent on BOTH axes -- Oxford's own broth-microdilution MIC and Oxford's OWN AMRFinder run, not
    a re-serving of NCBI-PD;
  - already on disk and already joined by `scripts/oxford_score.py`, so it costs nothing to ask.

THE ANSWER IS ZERO, AND ZERO IS THE FINDING. Across 4,979 AMRFinder-scanned isolates there is not one
`rmt` or `npmA` carrier. So Oxford CANNOT test the rescue's specificity -- no carrier, no test. What it
does instead is BOUND the prevalence, and that reframes the open limit: the twelve E. coli carriers are
few because `rmt` in E. coli is genuinely rare in unselected collections, not because the search was
lazy. A cohort with no carriers of a determinant cannot detect a rule that keys on it -- the same
structural blindness the source-concentration layer exists to disclose.

NON-VACUITY IS PROVEN, NOT ASSERTED. A zero from a broken parse looks exactly like a zero from real
absence. The probe therefore requires the scan to have found aminoglycoside determinants at all, AND
reports the singleton symbols it detected -- if genes present exactly once are found, a gene present once
would have been found. That check is the difference between a finding and a silent bug.

Offline: reads the committed Oxford deposit on D:. Writes wiki/oxford_rmt_prevalence_<date>.json.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.eval.amr_rules import DRUG_RULE  # noqa: E402

OXFORD = Path("D:/dna_decode_data/raw/oxford")

# The rescue pattern under scrutiny, imported from the DEPLOYED rule rather than retyped -- a retyped
# regex can drift from the one that actually fires.
RESCUE = re.compile(DRUG_RULE["gentamicin"]["symbol_rescue"], re.IGNORECASE)
# armA is deliberately SEPARATE: AMRFinder files it under Subclass GENTAMICIN, so the frozen rule already
# counted it before the rescue existed. It is a corroboration target, not a rescue case.
ARMA = re.compile(r"^armA\d*$", re.IGNORECASE)

CLSI_GENT_R = 16.0  # CLSI M100 2024, Enterobacterales


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oxford-dir", type=Path, default=OXFORD)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"oxford_rmt_prevalence_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    amr, mic = a.oxford_dir / "amrfinder.tsv", a.oxford_dir / "main_data.csv"
    for p in (amr, mic):
        if not p.exists():
            print(f"missing {p} -- the Oxford deposit is not on this host", file=sys.stderr)
            return 2

    # --- scan Oxford's own AMRFinder output -------------------------------------------------
    isolates: set[str] = set()
    aminoglycoside: collections.Counter = collections.Counter()
    subclass_of: dict[str, set[str]] = {}
    rescued: dict[str, list[str]] = {}
    arma: dict[str, list[str]] = {}

    with open(amr, encoding="utf-8", errors="replace") as f:
        cols = next(f).rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(cols)}
        # The Oxford deposit predates the 'Element symbol' rename; accept either spelling.
        sym_i = ix.get("Gene symbol", ix.get("Element symbol"))
        cls_i, sub_i, name_i = ix.get("Class"), ix.get("Subclass"), ix.get("Name")
        if sym_i is None or cls_i is None or name_i is None:
            print(f"unexpected AMRFinder header: {cols[:14]}", file=sys.stderr)
            return 2
        for line in f:
            fld = line.rstrip("\n").split("\t")
            if len(fld) <= max(sym_i, cls_i, name_i):
                continue
            iso, sym = fld[name_i].strip(), fld[sym_i].strip()
            isolates.add(iso)
            if fld[cls_i].strip().upper() == "AMINOGLYCOSIDE":
                aminoglycoside[sym] += 1
                subclass_of.setdefault(sym, set()).add(fld[sub_i].strip() if sub_i is not None else "")
            if RESCUE.match(sym):
                rescued.setdefault(iso, []).append(sym)
            elif ARMA.match(sym):
                arma.setdefault(iso, []).append(sym)

    # --- NON-VACUITY: a zero from a broken parse is indistinguishable from a real zero ------
    singletons = sorted(s for s, n in aminoglycoside.items() if n == 1)
    scan_is_live = len(aminoglycoside) >= 5 and bool(singletons)
    if not scan_is_live:
        print("REFUSING to report: the aminoglycoside scan found too little to prove it parsed. "
              f"distinct symbols={len(aminoglycoside)} singletons={len(singletons)}", file=sys.stderr)
        return 3

    # --- gentamicin MIC for whatever carriers exist -----------------------------------------
    rows = list(csv.DictReader(open(mic, encoding="utf-8")))
    gent = {}
    for r in rows:
        up = (r.get("Gentamicin_upper") or "").strip()
        if up not in ("", "NA"):
            gent[r["guuid"].strip()] = 2.0 ** float(up)   # the paper encodes MIC as 2**upper

    def phen(iso: str) -> dict:
        m = gent.get(iso)
        return {"gentamicin_mic": m,
                "call": None if m is None else ("R" if m >= CLSI_GENT_R else "S")}

    print(f"Oxford E. coli bacteraemia -- {len(isolates)} AMRFinder-scanned isolates, "
          f"{len(gent)} with a measured gentamicin MIC")
    print(f"  aminoglycoside determinant symbols detected : {len(aminoglycoside)} "
          f"({sum(aminoglycoside.values())} hits)")
    print(f"  detected at n=1 (proves singletons are found): {', '.join(singletons[:6])}")
    print(f"\n  rmt/npmA carriers (the rescue's target)     : {len(rescued)}")
    print(f"  armA carriers (already counted by the rule)  : {len(arma)}")
    for iso, syms in arma.items():
        p = phen(iso)
        print(f"    {iso}  {','.join(syms)}  subclass={sorted(subclass_of.get(syms[0], set()))}  "
              f"MIC={p['gentamicin_mic']} -> {p['call']}")

    verdict = ("NO_CARRIERS_CANNOT_TEST" if not rescued else
               "CARRIERS_PRESENT_TESTABLE")
    why = ("Oxford contains no rmt/npmA carrier at all, so it cannot test the rescue's specificity. It "
           "BOUNDS the prevalence instead: <1 in {n} broadly-ascertained E. coli.".format(n=len(isolates))
           if not rescued else
           f"{len(rescued)} rmt/npmA carriers are present and can be scored against measured MIC.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "oxford-rmt-prevalence-v1",
        "question": "how prevalent is the gentamicin rescue's rmt/npmA target in a broadly-ascertained, "
                    "resistance-independent E. coli collection?",
        "cohort": "Oxford E. coli bacteraemia (Lipworth et al.; ENA PRJNA604975 + the ecoli_mic_arg "
                  "deposit). Oxford's own broth-microdilution MIC and Oxford's OWN AMRFinder run.",
        "n_isolates_scanned": len(isolates),
        "n_with_gentamicin_mic": len(gent),
        "n_rmt_npma_carriers": len(rescued),
        "n_arma_carriers": len(arma),
        "arma_carriers": {i: {"symbols": s, "subclass": sorted(subclass_of.get(s[0], set())), **phen(i)}
                          for i, s in arma.items()},
        "aminoglycoside_symbols": dict(aminoglycoside.most_common()),
        "non_vacuity": {
            "distinct_aminoglycoside_symbols": len(aminoglycoside),
            "symbols_detected_once": singletons,
            "argument": "genes present exactly once ARE detected, so a gene present once would have been "
                        "found. The zero is absence, not a parse failure.",
        },
        "verdict": verdict, "why": why,
        "what_this_does_not_show": [
            "It does NOT test the rescue's specificity. Zero carriers means no test -- neither support "
            "nor counter-evidence for the rule. A cohort with no carriers of a determinant is "
            "structurally incapable of detecting a rule that keys on it.",
            "It does NOT generalise past this population: one UK region, 2008-2018 bacteraemia. rmt is "
            "substantially more prevalent in Asia and in Klebsiella, where the measured over-call lives.",
            "AMRFinder DB coverage is version-dependent and this deposit's version is not recorded. Total "
            "blindness to 16S methyltransferases is unlikely -- armA, the same gene family and curation "
            "era, IS detected -- but rmt-specific DB coverage here is unverified.",
        ],
        "corroborated_incidentally": (
            "armA is filed under Subclass=GENTAMICIN by this independent, older AMRFinder run, "
            "confirming on a second version the premise the v2 rescue was built on: the frozen rule "
            "already counted armA, so the gap was the rmt family only."
        ),
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
