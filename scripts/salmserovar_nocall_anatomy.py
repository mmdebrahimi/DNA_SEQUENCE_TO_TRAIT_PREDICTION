"""Why does the Salmonella serovar caller abstain? Measure the causes before fixing any of them.

THIS RUN EXISTS BECAUSE THE PREVIOUS DIAGNOSIS WAS WRONG. The 2026-09-04 validation reported that "33
of 59 no-calls have an empty H2" and named phase-2 flagellin handling as the cell's single largest
defect. The COUNT was right. The CAUSAL ATTRIBUTION was not, and a fix built on it would have fixed
almost nothing.

An antigenic formula is O:H1:H2. A formula printed as `4:H?:-` has an empty H2 -- but its H2 is empty
because NOTHING was resolved on the H axes, not because the genome is phase-1-only. Counting trailing
`-` conflates "the phase-2 antigen is genuinely absent" with "the caller failed upstream", and those
demand opposite fixes.

So this partitions every abstention by the FIRST axis that actually failed, and reports the size of the
only bucket a phase-2 fix could reach: formulas whose O:H1 pair is VALID in the White-Kauffmann table
and which therefore fail solely because H2 is missing.

IT ALSO MEASURES THE HEADROOM OF THE OBVIOUS FIX before anyone builds it. The tempting move is: when
H2 is absent, resolve on O:H1 alone if that pair maps to exactly ONE serovar. That is measured here and
the answer is that it recovers nothing -- the H2-blocked formulas are precisely the AMBIGUOUS ones,
which is why they need H2 in the first place. A fix whose headroom is zero should not be written.

Offline: reads the committed W-K-L table + the cached validation results.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

UNRESOLVED_O = {"O?", "", "-"}
UNRESOLVED_H = {"H?", "", "-"}


def load_table(p: Path) -> tuple[set[tuple[str, str]], dict[tuple[str, str], set[str]]]:
    """-> (valid O:H1 pairs, O:H1 -> serovar names)."""
    pairs: set[tuple[str, str]] = set()
    by_pair: dict[tuple[str, str], set[str]] = collections.defaultdict(set)
    with open(p, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:4] != ["O", "H1", "H2", "Serovar"]:
            raise ValueError(f"unexpected serovar table header: {header[:4]}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or not f[3].strip():
                continue
            key = (f[0].strip(), f[1].strip())
            pairs.add(key)
            by_pair[key].add(f[3].strip())
    return pairs, by_pair


def classify(formula: str | None, pairs: set) -> str:
    """The FIRST axis that failed -- not merely whether H2 is empty."""
    if not formula or formula.count(":") != 2:
        return "no_formula_at_all"
    o, h1, h2 = [x.strip() for x in formula.split(":")]
    o_bad = o in UNRESOLVED_O or o.startswith("O?")
    h1_bad = h1 in UNRESOLVED_H
    if o_bad and h1_bad:
        return "both_O_and_H1_unresolved"
    if o_bad:
        return "O_antigen_unresolved"
    if h1_bad:
        return "H1_phase1_flagellin_unresolved"
    if (o, h1) in pairs:
        return "O_H1_valid_only_H2_blocks_it"
    return "O_H1_called_but_pair_absent_from_table"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/results.jsonl"))
    ap.add_argument("--table", type=Path,
                    default=ROOT / "data" / "salmserovar_db" / "serovar_table.tsv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"salmserovar_nocall_anatomy_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.results.exists():
        print(f"validation results absent at {a.results}", file=sys.stderr)
        return 2
    pairs, by_pair = load_table(a.table)
    rows = [json.loads(ln) for ln in a.results.read_text(encoding="utf-8").splitlines() if ln.strip()]
    ok = [r for r in rows if r.get("status") == "ok"]
    nocalls = [r for r in ok if not r.get("pred_serovar")]
    if not nocalls:
        print("no abstentions to analyse -- nothing to diagnose", file=sys.stderr)
        return 3

    causes = collections.Counter(classify(r.get("pred_formula"), pairs) for r in nocalls)

    # HEADROOM of the tempting fix: resolve on O:H1 alone when it is UNIQUE in the table.
    recoverable, would_be_right = 0, 0
    for r in nocalls:
        f = r.get("pred_formula") or ""
        if f.count(":") != 2:
            continue
        o, h1, h2 = [x.strip() for x in f.split(":")]
        if h2 not in ("-", ""):
            continue
        if o in UNRESOLVED_O or h1 in UNRESOLVED_H:
            continue
        cand = by_pair.get((o, h1), set())
        if len(cand) == 1:
            recoverable += 1
            if (r.get("label") or "").lower() == next(iter(cand)).lower():
                would_be_right += 1

    n = len(nocalls)
    print(f"abstentions: {n} of {len(ok)} scored\n")
    for k, v in causes.most_common():
        print(f"  {v:>3}  ({v / n:5.1%})  {k}")
    h2_bucket = causes["O_H1_valid_only_H2_blocks_it"]
    print(f"\n  bucket a PHASE-2 fix could reach : {h2_bucket} ({h2_bucket / n:.1%})")
    print(f"  headroom of the O:H1-unique fallback: {recoverable} recoverable, "
          f"{would_be_right} of them would match the wet-lab label")

    biggest = causes.most_common(1)[0]
    if h2_bucket and biggest[0] == "O_H1_valid_only_H2_blocks_it":
        verdict = "H2_IS_THE_DOMINANT_CAUSE"
    else:
        verdict = "H2_IS_NOT_THE_DOMINANT_CAUSE"
    why = (f"the largest cause is `{biggest[0]}` at {biggest[1]} of {n} ({biggest[1] / n:.1%}); the "
           f"phase-2-reachable bucket is {h2_bucket} ({h2_bucket / n:.1%}). The O:H1-unique fallback "
           f"recovers {recoverable} formulas, so its headroom is "
           f"{'zero -- do not build it' if recoverable == 0 else 'nonzero'}.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {
        "schema": "salmserovar-nocall-anatomy-v1", "date": _date.today().isoformat(),
        "question": "which axis actually causes the serovar caller to abstain?",
        "supersedes": ("the 2026-09-04 claim that phase-2 flagellin handling is this cell's single "
                       "largest defect -- the COUNT of empty-H2 formulas was right, the causal "
                       "attribution was not"),
        "n_scored": len(ok), "n_abstentions": n,
        "causes": dict(causes.most_common()),
        "phase2_reachable_bucket": h2_bucket,
        "oh1_unique_fallback_headroom": {"recoverable": recoverable,
                                         "would_match_label": would_be_right},
        "verdict": verdict, "why": why,
        "fix_priority_by_measured_size": [k for k, _ in causes.most_common()],
        "honest_limits": [
            "One cohort (N=200, reference-lab-filtered) and one antigen DB build. The cause mix would "
            "differ with a different O-antigen allele set.",
            "This partitions by the FIRST failing axis. A formula can be wrong on more than one axis; "
            "the counts are of first-failure, not of independent defects.",
            "It measures WHY the caller abstains, NOT whether a resolved call would be correct -- "
            "abstention and error are different failures and are counted separately throughout.",
            "The O-antigen bucket is the largest but its fix is DB coverage (data engineering on the "
            "wzx/wzy allele set), not a code change -- sized here, not attempted.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
