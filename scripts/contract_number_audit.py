"""Check every number a cell contract CITES against the artifact it cites it FROM.

WHY
---
This repo's most frequently-recurring defect class is not wrong code, it is a claim that DRIFTED from the
artifact it rests on. Measured instances, all found by comparing a claim to its own evidence rather than by
any test: `project_status.py` quoting a figure its own data layer forbids; the report card understating the
validated surface by 4 cells; `--help` advertising two HLA screens the repo had already demoted; ResFinder's
"165 beta-lactamase genes"; and this week CYP4F2's contract citing a 54/54 concordance that was a different
axis computed on a smaller parsed set, plus a "~79% EAS" allele frequency whose real value is 0.2128.

`tests/test_cell_registry.py` already pins that every cited artifact EXISTS. Nothing checks whether the
NUMBERS in the prose are in the artifact. That is the gap this closes.

WHAT IT IS, HONESTLY
--------------------
A TRIAGE FUNNEL, not an oracle -- the same posture as scripts/kaggle_staleness_auditor.py, whose field
precision is about 1 in 6. A cited number can legitimately be absent from its artifact: it may be derived
(a ratio, a sum, a percentage of a raw count), rounded differently, or carried over from a superseded run
that the prose correctly describes as superseded. So a miss is CANDIDATE_DRIFT requiring adjudication, and
this script NEVER edits anything. Adjudication is the product.

Scope: only contracts that name a `wiki/...` artifact, and only decimal numbers (a bare integer in prose is
usually a count whose denominator lives elsewhere, which this cannot check without guessing).
"""
from __future__ import annotations

import datetime
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

PROSE_FIELDS = ("claim", "validation_slice", "label_provenance", "demotion_rule", "claim_status")
# Contracts cite artifacts in TWO styles, and both must resolve or the audit invents drift. The second is
# shell-style brace expansion, used throughout this repo:
#     wiki/cyp4f2_getrm_concordance_2026-09-24.{md,json}
#     wiki/forward_inverse_{roundtrip,sweep,deployable}_2026-07-1{6,7}.{md,json}
# Missing it produced a false CANDIDATE_DRIFT on the inverse cell's "+53.0%", whose source
# (margin_vs_better_baseline = 0.5295) sits in a roundtrip artifact the regex never loaded.
ARTIFACT_RE = re.compile(r"wiki/[A-Za-z0-9_][A-Za-z0-9_./{},-]*")
BRACE_RE = re.compile(r"\{([^{}]*)\}")


def _expand_braces(path: str) -> list[str]:
    """Expand one or more {a,b,c} groups into every concrete path."""
    out = [path]
    while any(BRACE_RE.search(p) for p in out):
        nxt = []
        for p in out:
            m = BRACE_RE.search(p)
            if not m:
                nxt.append(p)
                continue
            for alt in m.group(1).split(","):
                nxt.append(p[:m.start()] + alt.strip() + p[m.end():])
        out = nxt
    return out
# A decimal immediately preceded by a letter/underscore is part of an IDENTIFIER, not a measurement --
# canFam3.1, GCA_000005845.2, v1.2. Capturing those produced a false CANDIDATE_DRIFT on the dog
# morphology cell's "3.1", which is the assembly name canFam3.1.
DECIMAL_RE = re.compile(r"(?<![\w.])(\d{1,4}\.\d{1,4})(?![\d.])")

# Numbers that are not measurements: dates, tool versions, guideline years, p-value exponents.
NON_MEASUREMENT = re.compile(r"^(19|20)\d{2}\.")

# Numbers ADJUDICATED once as legitimately absent from any artifact, with the reason. Deliberately a small
# explicit table keyed by (cell_id, number): a blanket "is it anywhere in the source tree" check was tried
# and is vacuous (it accepted 12/12 randomly-generated numbers), so each exemption is named or it flags.
ADJUDICATED_BENIGN: dict[tuple[str, str], str] = {
    ("typing:Escherichia_coli:serotype", "0.60"):
        "not a measurement: the project's own source-diversity BAR, defined in code as "
        "scripts/source_diverse_validate.MAX_SOURCE_SHARE. A threshold has no artifact to live in.",
}


def _cited_artifacts(blob: str) -> list[str]:
    out, seen = [], set()
    for raw in ARTIFACT_RE.findall(blob):
        for p in _expand_braces(raw.rstrip(".,;:)")):
            p = p.rstrip(".,;:)_")
            # a contract may cite a stem without its extension (…_2026-09-24) -> try both sidecars
            for cand in ([p] if p.endswith((".json", ".md")) else [p + ".json", p + ".md", p]):
                if cand not in seen:
                    seen.add(cand)
                    out.append(cand)
    return out


def _load_artifact_text(rel: str) -> str | None:
    p = REPO / rel
    if not p.exists() or p.is_dir():
        return None
    try:
        raw = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    if p.suffix == ".json":
        try:
            # re-dump so nested numeric values are searchable as plain text
            return raw + "\n" + json.dumps(json.loads(raw))
        except json.JSONDecodeError:
            return raw
    return raw


def _artifact_numbers(text: str) -> list[float]:
    """Every numeric token in the artifact, as floats, for rounding-aware comparison."""
    out = []
    for m in re.findall(r"-?\d+(?:\.\d+)?", text):
        try:
            out.append(float(m))
        except ValueError:
            pass
    return out


def _matches_by_rounding(tok: str, art_nums: list[float]) -> float | None:
    """Is the cited number a ROUNDING of some number in the artifact (or vice versa)?

    Load-bearing: prose legitimately rounds. This audit's first run flagged SLCO1B1's cited 0.398 as drift
    when the artifact holds 0.3977 -- a false positive produced by string matching alone, which would have
    sent a reader to adjudicate a number that is correct. Compare numerically at the cited precision, and
    also allow the percent/fraction pairing.
    """
    try:
        cited = float(tok)
    except ValueError:
        return None
    dec = len(tok.split(".")[1]) if "." in tok else 0
    # Tolerance, NOT round(): round(52.95, 1) is 52.9 because 52.95 has no exact binary representation,
    # so a cited "53.0" would spuriously fail against an artifact's 0.5295. Half a unit at the cited
    # precision is what "rounds to" means, and it is immune to float-repr edge cases.
    tol = 0.5 * (10.0 ** -dec) + 1e-9
    for a in art_nums:
        for cand in (a, a * 100.0, a / 100.0):
            if abs(cand - cited) <= tol:
                return a
    return None


_SOURCE_BLOB: str | None = None


def _package_source_blob() -> str:
    """Concatenated package + scripts source, for tracing a cited THRESHOLD to the code that defines it."""
    global _SOURCE_BLOB
    if _SOURCE_BLOB is None:
        parts = []
        for d in ("dna_decode", "scripts"):
            for p in sorted((REPO / d).rglob("*.py")):
                try:
                    parts.append(p.read_text(encoding="utf-8", errors="replace"))
                except OSError:
                    pass
        _SOURCE_BLOB = "\n".join(parts)
    return _SOURCE_BLOB


def _in_package_source(tok: str) -> bool:
    return any(v in _package_source_blob() for v in _number_variants(tok))


def _number_variants(tok: str) -> list[str]:
    """A cited 0.83 may live in the artifact as 0.8300, 0.83, or 83.0 (percent). Accept those forms."""
    v = {tok, tok.rstrip("0").rstrip(".") if "." in tok else tok}
    try:
        f = float(tok)
    except ValueError:
        return sorted(v)
    for nd in (1, 2, 3, 4):
        v.add(f"{f:.{nd}f}")
    if f < 1:                      # fraction cited, percent stored (or vice versa)
        v.add(f"{f * 100:.1f}")
        v.add(f"{f * 100:.0f}")
        v.add(f"{f * 100:g}")
    if f > 1:
        v.add(f"{f / 100:.4f}".rstrip("0"))
    return sorted(x for x in v if x)


def audit() -> dict:
    sys.path.insert(0, str(REPO))
    from dna_decode.data.cell_registry import cells

    rows = []
    # Cells carrying candidate numbers that cite NO artifact. `scope` below has always said the audit
    # covers "cells whose prose names a wiki/ artifact", so the exclusion is disclosed -- but its SIZE
    # was not, and that is what tells a reader whether the audited set is near-complete or a fraction of
    # the surface. It is a fraction: the skipped set is LARGER than the checked one. A cell citing only
    # its SCRIPT (typing:bacteria:mlst names scripts/mlst_serotype_purity.py, never the wiki artifact
    # holding the numbers) drops out entirely, so its measured purity figures are unverifiable here.
    # Reported, never counted as drift -- an uncited number is unchecked, not wrong.
    unverifiable = []
    for c in cells():
        blob = " ".join(str(getattr(c, f, "") or "") for f in PROSE_FIELDS)
        arts = _cited_artifacts(blob)
        if not arts:
            bare = [t for t in dict.fromkeys(DECIMAL_RE.findall(blob))
                    if not NON_MEASUREMENT.match(t) and (c.cell_id, t) not in ADJUDICATED_BENIGN]
            if bare:
                unverifiable.append({"cell_id": c.cell_id, "n_numbers": len(bare),
                                     "numbers": bare, "why": "prose cites no wiki/ artifact"})
            continue
        texts = {a: _load_artifact_text(a) for a in arts}
        resolved = {a: t for a, t in texts.items() if t is not None}
        if not resolved:
            rows.append({"cell_id": c.cell_id, "status": "NO_ARTIFACT_RESOLVED",
                         "cited": arts, "checked": 0, "found": 0, "candidate_drift": []})
            continue
        haystack = "\n".join(resolved.values())
        nums = [t for t in dict.fromkeys(DECIMAL_RE.findall(blob)) if not NON_MEASUREMENT.match(t)]
        art_nums = _artifact_numbers(haystack)
        found, missing, rounded, in_code = [], [], {}, []
        for tok in nums:
            if any(v in haystack for v in _number_variants(tok)):
                found.append(tok)
                continue
            src = _matches_by_rounding(tok, art_nums)
            if src is not None:
                found.append(tok)
                rounded[tok] = src          # matched only after rounding; recorded for transparency
                continue
            # A cited number may be a THRESHOLD / BAR rather than a measurement, correctly absent from
            # every artifact. Those are handled by the SMALL explicit ADJUDICATED_BENIGN table above --
            # NOT by searching the package source, which was tried and is VACUOUS: variant matching of a
            # 2-4 digit token against a multi-megabyte source blob accepted 17/17 tokens including 12
            # randomly-generated ones, turning real signal into a clean-looking zero.
            key = (c.cell_id, tok)
            if key in ADJUDICATED_BENIGN:
                found.append(tok)
                in_code.append({"number": tok, "why": ADJUDICATED_BENIGN[key]})
            else:
                missing.append(tok)
        rows.append({
            "cell_id": c.cell_id,
            "status": "CANDIDATE_DRIFT" if missing else "ALL_CITED_NUMBERS_PRESENT",
            "cited": sorted(resolved), "unresolved_citations": sorted(set(arts) - set(resolved)),
            "checked": len(nums), "found": len(found), "candidate_drift": missing,
            "matched_only_after_rounding": rounded,
            "adjudicated_benign": in_code,
        })

    checked = sum(r["checked"] for r in rows)
    drift = sum(len(r["candidate_drift"]) for r in rows)
    return {
        "schema": "contract-number-audit-v1",
        "analysis_date": datetime.date.today().isoformat(),
        "posture": "TRIAGE FUNNEL, not an oracle -- a cited number may legitimately be derived, rounded, or "
                   "describe a superseded run. Every hit needs adjudication; this script edits nothing.",
        "scope": "cells whose prose names a wiki/ artifact; decimal numbers only",
        "n_cells_audited": len(rows),
        "n_numbers_checked": checked,
        "n_candidate_drift": drift,
        "hit_rate": round(drift / checked, 4) if checked else None,
        "verdict": "NOTHING_TO_ADJUDICATE" if not drift else "ADJUDICATION_REQUIRED",
        # Coverage, so a clean verdict cannot be read as "every contract number is verified". These do
        # NOT enter n_numbers_checked / n_candidate_drift / verdict -- unchecked is not drift.
        "n_cells_with_numbers_but_no_citation": len(unverifiable),
        "n_numbers_unverifiable": sum(u["n_numbers"] for u in unverifiable),
        "cells_unverifiable": sorted(unverifiable, key=lambda u: u["cell_id"]),
        "cells": rows,
    }


def main(argv=None) -> int:
    rep = audit()
    if rep["n_numbers_checked"] == 0:
        print("REFUSED: no cited numbers were checkable -- a clean report over zero checks is not a result.",
              file=sys.stderr)
        return 3
    stamp = rep["analysis_date"]
    (REPO / "wiki" / f"contract_number_audit_{stamp}.json").write_text(json.dumps(rep, indent=2),
                                                                      encoding="utf-8")
    for r in rep["cells"]:
        flag = "" if r["status"] == "ALL_CITED_NUMBERS_PRESENT" else f"  <- {r['status']}"
        print(f"{r['cell_id']:42} checked={r['checked']:3} found={r['found']:3} "
              f"drift={len(r['candidate_drift'])}{flag}")
        if r["candidate_drift"]:
            print(f"    cited but absent from {', '.join(r['cited'])}: "
                  f"{', '.join(r['candidate_drift'])}")
    print(f"\n{rep['verdict']}  {rep['n_candidate_drift']}/{rep['n_numbers_checked']} numbers need "
          f"adjudication across {rep['n_cells_audited']} cells")
    if rep["cells_unverifiable"]:
        print(f"COVERAGE: {rep['n_numbers_unverifiable']} number(s) in "
              f"{rep['n_cells_with_numbers_but_no_citation']} further cell(s) are UNVERIFIABLE here — "
              f"they cite no wiki/ artifact, so the clean verdict above does not cover them:")
        for u in rep["cells_unverifiable"]:
            print(f"    {u['cell_id']:42} {u['n_numbers']:3} number(s): "
                  f"{', '.join(u['numbers'][:8])}{' …' if u['n_numbers'] > 8 else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
