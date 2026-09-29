"""Does a cited number match the right QUANTITY, or merely the right VALUE?

THE GAP THIS MEASURES. `scripts/contract_number_audit.py` verifies that a number quoted in a cell
contract is PRESENT in the artifact the contract cites, as a whole number (after the 2026-09-28
boundary-guard fix). Presence is necessary and not sufficient: a contract claiming `sens 0.993` is
"verified" by an artifact carrying `spec: 0.993`, because the audit compares VALUES and never asks
which FIELD the value came from. The audit's own residual says so in as many words -- "whole-number
matching is not SEMANTIC" -- and this script is the measurement of how much that costs.

WHAT IT DOES. For every cited number it extracts (a) the quantity LABEL the prose claims the number
is, and (b) every numeric FIELD PATH in the cited artifact whose value equals that number. It then
asks whether any matching path is label-CONSISTENT.

    prose "sens 0.993"  +  artifact metrics.sens = 0.993   -> field_confirmed
    prose "sens 0.993"  +  artifact metrics.spec = 0.993    -> field_mismatch   <-- the defect
    prose "0.993"       (no label)                          -> label_unlabeled  (not measurable)

HONESTY, and it is load-bearing in three places.

1. **A `field_mismatch` is a LEAD, not a proven defect.** Equal values legitimately co-occur -- a
   cohort with spec 1.000 and a normalized share of 1.0, or a sens that genuinely equals a recall
   elsewhere in the same file. Every flag is reported for ADJUDICATION with its competing paths
   attached; this script never edits a contract.
2. **Label extraction is HEURISTIC.** It reads a window around the number for a known quantity word.
   Prose that names the quantity far from the number, or in a table header, is `label_unlabeled` and
   is counted as NOT MEASURED rather than silently passed. The measurable fraction is reported first,
   because a check that can only see a third of its input should say so before it says anything else.
3. **`field_confirmed` is weaker than it sounds.** It says a label-consistent field carries that
   value -- not that the contract's sentence is a faithful summary of the artifact.

DELIBERATELY NOT A GATE. It reports; exit code is 0 unless the run itself failed. The repo's standing
lesson is that a screen which rewrites what it screens becomes its own control case.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from scripts.contract_number_audit import (  # noqa: E402
    DECIMAL_RE,
    _bounded_in,
    NON_MEASUREMENT,
    PROSE_FIELDS,
    _cited_artifacts,
    _load_artifact_text,
)

# Quantity vocabulary: canonical label -> (prose words, field-name tokens).
# Prose words are matched case-insensitively in a window around the number; field tokens are matched
# against the artifact's dotted field path. Kept deliberately SMALL and explicit -- a sprawling
# synonym list would make almost any path "consistent" and quietly turn the check vacuous.
LABEL_VOCAB: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {
    "sens":      (("sens", "sensitivity", "recall", "tpr"),        ("sens", "sensitivity", "recall", "tpr")),
    "spec":      (("spec", "specificity", "tnr"),                  ("spec", "specificity", "tnr")),
    "acc":       (("acc", "accuracy"),                             ("acc", "accuracy")),
    "prec":      (("prec", "precision", "ppv"),                    ("prec", "precision", "ppv")),
    "auc":       (("auc", "auroc", "auprc", "roc"),                ("auc", "auroc", "auprc", "roc")),
    "corr":      (("spearman", "pearson", "rho", "correlation", "corr"),
                  ("spearman", "pearson", "rho", "corr")),
    "mcc":       (("mcc",),                                        ("mcc",)),
    "f1":        (("f1",),                                         ("f1",)),
    "jaccard":   (("jaccard",),                                    ("jaccard",)),
    "balacc":    (("balacc", "balanced accuracy", "bal-acc", "balanced"),
                  ("balacc", "balanced")),
    "share":     (("share", "fraction", "proportion"),             ("share", "fraction", "proportion")),
    "rate":      (("rate",),                                       ("rate",)),
    "pvalue":    (("p-value", "pvalue", "p ="),                    ("p_value", "pvalue", "pval")),
    "vme":       (("vme", "very major error"),                     ("vme",)),
    "identity":  (("identity", "ident"),                           ("identity", "ident", "pct_identity")),
    "coverage":  (("coverage", "cov"),                             ("coverage", "cov")),
    # Added 2026-09-28, driven by an OBSERVED occurrence rather than anticipation: the cyp4f2 cell quotes
    # panel-derived allele frequencies (EUR 0.2773 / EAS 0.2128) that live at `...alt_freq` paths, and with
    # no `freq` quantity `path_is_consistent` refused every binding to them. Note the PROSE side barely
    # matters here -- that prose says "EUR 0.2773" with no quantity word at all, so only a BINDING can reach
    # these numbers; the vocabulary entry exists so the binding's declared quantity can be checked against
    # the path. Grown per observed need on purpose: this vocabulary must stay INDEPENDENT of the artifacts'
    # own field names, because deriving it from them would make every path "consistent" and the mismatch
    # check vacuous.
    "freq":      (("frequency", "freq", "allele frequency"),       ("freq", "frequency", "af")),
    # Deliberately NOT an alias of `acc`. Concordance is agreement between two CALLERS where neither is
    # truth; accuracy is measured against a label. That distinction is the whole basis of this repo's
    # tier system (FAITHFUL_TO_TOOL vs INDEPENDENT_MEASURED), so collapsing them would let a tool-agreement
    # figure verify a prose claim of accuracy -- the exact right-value/wrong-quantity error this check is for.
    # `agreement` is a near-synonym and is EXCLUDED: measured at 2 conversions against 1 FALSE LEAD, and a
    # false lead is strictly worse than an honestly-unlabeled number, so +1 net does not pay for it.
    "concordance": (("concordance",),                              ("concordance",)),
}

# How far to look for a label. Asymmetric on purpose: English writes "sens 0.993" and "0.993 recall",
# but the leading form dominates in these contracts, so the backward window is the wider one.
WINDOW_BEFORE = 44
WINDOW_AFTER = 0  # trailing labels are not used; see label_for

ST_CONFIRMED = "field_confirmed"
ST_MISMATCH = "field_mismatch"
ST_UNLABELED = "label_unlabeled"
ST_NO_FIELD = "no_numeric_field"
ST_ABSENT = "absent_from_artifact"
# The value is present but only in a markdown-cited artifact (or reached via artifact pooling), so
# this check cannot adjudicate its FIELD. Not measurable -- deliberately NOT counted as a mismatch.
ST_MD_ONLY = "md_only_or_pooled"
# --- binding states (added 2026-09-28) -----------------------------------------------------------
# A DECLARED binding resolved: the prose token is present, the cited artifact carries that value at the
# declared field path, and the path is consistent with the declared quantity. This is the only state that
# does NOT rest on the prose-window heuristic.
ST_BOUND = "field_confirmed_by_binding"
# A binding was declared but its token does not appear in ANY prose field of that cell. This is a DEFECT,
# not a pass: it means the declaration has drifted from the prose it is supposed to make checkable, and it
# is exactly how a binding could quietly become a second truth surface.
ST_BOUND_NO_PROSE = "binding_declared_but_prose_token_absent"
# A binding was declared but the cited artifact does not carry that value at the declared path.
ST_BOUND_UNRESOLVED = "binding_declared_but_path_unresolved"
MEASURABLE_BY_BINDING = (ST_BOUND,)

MEASURABLE = (ST_CONFIRMED, ST_MISMATCH)


def numeric_leaves(obj, prefix: str = "") -> list[tuple[str, float]]:
    """Every numeric leaf in a parsed artifact, as (dotted_path, value).

    Booleans are excluded: `bool` is a subclass of `int` in Python, so `"verified": true` would
    otherwise enter the pool as the value 1 and could "confirm" a prose 1.0.
    """
    out: list[tuple[str, float]] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.extend(numeric_leaves(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out.extend(numeric_leaves(v, f"{prefix}[{i}]"))
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float)):
        out.append((prefix, float(obj)))
    return out


def value_matches(token: str, value: float) -> bool:
    """Does `value` equal the prose token at the token's own precision?

    Rounding to the token's decimals is what lets a prose `0.993` match a stored `0.9928` -- the
    contracts quote rounded figures, so an exact-string compare would report a false absence.
    """
    decimals = len(token.split(".")[1]) if "." in token else 0
    try:
        return round(value, decimals) == float(token)
    except (ValueError, OverflowError):
        return False


def label_for(blob: str, start: int, end: int,
              prev_end: int = 0, next_start: int | None = None) -> str | None:
    """The quantity label this number occurrence claims, or None if the prose does not name one.

    TWO RULES, both learned by measuring the first version against real contracts.

    (1) **A LEADING label wins outright.** The naive "closest label" rule mis-assigns systematically:
        in `sens 1.00 / spec 0.939` the word `spec` sits 3 chars AFTER `1.00` while `sens` sits 6
        chars before it, so nearest-wins labelled `1.00` as the specificity. A label that FOLLOWS a
        number almost always belongs to the NEXT number, so `before` is consulted first and `after`
        is only a fallback for the genuine trailing form (`0.833 recall`).

    (2) **Windows are clipped to the INTER-NUMBER span.** Without that, a window can reach across a
        neighbouring number and harvest its label. Scoping to (prev_end, next_start) makes each
        number's label search structurally unable to see a sibling's label.
    """
    lo = max(prev_end, start - WINDOW_BEFORE)
    before = blob[lo:start].lower()
    hi = end + WINDOW_AFTER if next_start is None else min(next_start, end + WINDOW_AFTER)
    after = blob[end:max(end, hi)].lower()

    # WORD-BOUNDARY GUARDED, and this is the whole ballgame. The first version matched prose words by
    # raw substring while the FIELD side (path_is_consistent) was already token-aware -- the same
    # strict-on-one-side/loose-on-the-other asymmetry that caused the parent audit's original defect,
    # reproduced one layer up. Adjudicating the 12 flags it produced, 10 were this bug:
    #     "J-acc-ard"        -> acc    (3 flags)
    #     "alphamis-sens-e"  -> sens   (2 flags)
    #     "corr-oborating"   -> corr
    #     "dis-cov-ery"      -> coverage
    #     "subst-rate"       -> rate
    # Every one of those numbers was in fact sitting at its CORRECT field. A checker must be held to
    # the standard it enforces.
    best: tuple[int, str] | None = None
    for canon, (words, _fields) in LABEL_VOCAB.items():
        for w in words:
            for m in re.finditer(rf"(?<![a-z0-9]){re.escape(w)}(?![a-z0-9])", before):
                dist = len(before) - m.start()
                if best is None or dist < best[0]:
                    best = (dist, canon)
    # NO TRAILING FALLBACK. A label that follows a number nearly always belongs to the NEXT number,
    # and clipping to the inter-number span does not help when the label sits BETWEEN them: prose
    # "9.16 -> 2.44 (AMRFinder 2.46), Jaccard 0.4126" made the trailing rule label 2.46 -- a
    # mean-genes-per-genome figure -- as a Jaccard. Measured across this corpus the trailing form
    # produced at least one wrong label and zero correct ones, so numbers whose quantity is named
    # only AFTER them are reported label_unlabeled (not measured) rather than guessed.
    return best[1] if best else None


def path_is_consistent(path: str, canon: str) -> bool:
    """Is this artifact field path plausibly the claimed quantity?

    Token-aware rather than a bare substring test: `spec` must not be satisfied by `..._specimen`,
    and `acc` must not be satisfied by `accession` -- the second is a real path in these artifacts
    and would have manufactured confirmations for every accuracy figure in the repo.
    """
    tokens = re.split(r"[^a-z0-9]+", path.lower())
    fields = LABEL_VOCAB[canon][1]
    return any(t == f or t.startswith(f + "_") or t.endswith("_" + f) for t in tokens for f in fields)


def audit_semantics() -> dict:
    from dna_decode.data.cell_registry import cells

    rows = []
    counts: dict[str, int] = {}
    for c in cells():
        blob = " ".join(str(getattr(c, f, "") or "") for f in PROSE_FIELDS)
        arts = _cited_artifacts(blob)
        if not arts:
            continue
        # Parse every cited artifact that is JSON; a cited .md has no field structure, which is a
        # LIMIT of this check (reported as no_numeric_field), not a pass.
        leaves: list[tuple[str, str, float]] = []
        md_texts: list[str] = []
        resolved = []
        n_unparseable = 0
        for a in arts:
            if _load_artifact_text(a) is None:
                continue
            resolved.append(a)
            if not a.endswith(".json"):
                t = _load_artifact_text(a)
                if t:
                    md_texts.append(t)
                continue
            # Read the JSON from DISK, not via _load_artifact_text: that helper builds a search
            # HAYSTACK and concatenates the brace-expanded siblings, so every json.loads failed with
            # "Extra data: line N column 1" and the whole run reported 0 of 180 measurable -- a broken
            # check that looks exactly like a clean one. Parse the real file instead.
            p = REPO / a
            try:
                leaves.extend((a, path, v)
                              for path, v in numeric_leaves(json.loads(p.read_text(encoding="utf-8"))))
            except (OSError, json.JSONDecodeError, ValueError):
                n_unparseable += 1
        if not resolved:
            continue

        # Index this cell's DECLARED bindings by their prose token. Empty for every cell that declares
        # none, which is the default -- the heuristic path below is unchanged for them.
        bindings = {b.number_token: b for b in getattr(c, "metric_bindings", ())}
        bound_tokens_seen: set[str] = set()

        seen: set[str] = set()
        ms = list(DECIMAL_RE.finditer(blob))
        for idx, m in enumerate(ms):
            tok = m.group(1)
            if tok in seen or NON_MEASUREMENT.match(tok):
                continue
            seen.add(tok)
            prev_end = ms[idx - 1].end(1) if idx else 0
            next_start = ms[idx + 1].start(1) if idx + 1 < len(ms) else None
            canon = label_for(blob, m.start(1), m.end(1), prev_end, next_start)
            hits = [(a, p, v) for (a, p, v) in leaves if value_matches(tok, v)]
            # A DECLARED binding takes precedence over the prose-window heuristic, because it is the one
            # path that does not have to GUESS the quantity. It is still VERIFIED, never trusted: the
            # artifact must actually carry this value at the declared path, and the path must be consistent
            # with the declared quantity. A binding that fails either check is a defect, not a pass.
            bound = bindings.get(tok)
            if bound is not None:
                canon = bound.quantity
                at_path = [(a, p, v) for (a, p, v) in leaves
                           if p == bound.field_path and a == bound.artifact and value_matches(tok, v)]
                if not at_path:
                    status, consistent, competing = ST_BOUND_UNRESOLVED, [], [
                        f"declared {bound.artifact}:{bound.field_path}"]
                elif not path_is_consistent(bound.field_path, bound.quantity):
                    status, consistent, competing = ST_MISMATCH, [], [
                        f"{bound.artifact}:{bound.field_path} (declared quantity {bound.quantity!r} "
                        f"is not consistent with that path)"]
                else:
                    status = ST_BOUND
                    consistent = [f"{bound.artifact}:{bound.field_path}"]
                    competing = []
            elif canon is None:
                status, consistent, competing = ST_UNLABELED, [], []
            elif not hits:
                status, consistent, competing = ST_NO_FIELD, [], []
            else:
                consistent = [f"{a}:{p}" for (a, p, _v) in hits if path_is_consistent(p, canon)]
                competing = [f"{a}:{p}" for (a, p, _v) in hits if not path_is_consistent(p, canon)]
                if consistent:
                    status = ST_CONFIRMED
                elif any(_bounded_in(tok, md_text) for md_text in md_texts):
                    # ARTIFACT POOLING, measured: a cell cites several artifacts and this check pools
                    # their JSONs, so a number whose OWN artifact is markdown-only gets adjudicated
                    # against unrelated files. `forward`'s 0.49 is cited to esm_at_scale_2026-07-17,
                    # which ships as .md ONLY -- it was flagged a mismatch purely because its value
                    # also occurs at unrelated paths in sibling JSONs. Not measurable, not a defect.
                    status, competing = ST_MD_ONLY, competing
                else:
                    status = ST_MISMATCH
            counts[status] = counts.get(status, 0) + 1
            rows.append({
                "cell_id": getattr(c, "cell_id", "?"),
                "number": tok,
                "claimed_label": canon,
                "status": status,
                "consistent_paths": consistent[:4],
                "competing_paths": competing[:4],
                "from_binding": bound is not None,
            })
            if bound is not None:
                bound_tokens_seen.add(tok)

        # A binding whose token never turned up in this cell's prose is a DEFECT. Without this check a
        # binding could sit in the registry describing a number no reader will ever see -- the audit would
        # then be auditing the declaration instead of the prose, which is the failure the binding design
        # exists to avoid.
        for tok, b in bindings.items():
            if tok in bound_tokens_seen:
                continue
            counts[ST_BOUND_NO_PROSE] = counts.get(ST_BOUND_NO_PROSE, 0) + 1
            rows.append({
                "cell_id": getattr(c, "cell_id", "?"),
                "number": tok,
                "claimed_label": b.quantity,
                "status": ST_BOUND_NO_PROSE,
                "consistent_paths": [],
                "competing_paths": [f"declared {b.artifact}:{b.field_path}"],
                "from_binding": True,
            })

    n_total = len(rows)
    n_by_binding = sum(counts.get(s, 0) for s in MEASURABLE_BY_BINDING)
    n_by_heuristic = sum(counts.get(s, 0) for s in MEASURABLE)
    n_measurable = n_by_binding + n_by_heuristic
    mismatches = [r for r in rows if r["status"] == ST_MISMATCH]
    binding_defects = [r for r in rows
                       if r["status"] in (ST_BOUND_NO_PROSE, ST_BOUND_UNRESOLVED)]
    return {
        "artifact": "contract_number_semantics",
        "schema": "contract-number-semantics-v1",
        "scope": ("cells whose prose cites a wiki/ artifact; only numbers whose prose names a "
                  "quantity are measurable at all"),
        "n_numbers": n_total,
        "n_measurable": n_measurable,
        "measurable_fraction": round(n_measurable / n_total, 4) if n_total else None,
        # THE SPLIT IS LOAD-BEARING, not cosmetic. Without it the headline measurability fraction
        # conflates two different things -- "the prose-window heuristic could not infer the quantity" and
        # "no author declared a binding" -- so the number would move with AUTHORING COVERAGE rather than
        # with audit power, and a rising figure could mean either. Report the denominators apart.
        "measurable_by_binding": n_by_binding,
        "measurable_by_heuristic": n_by_heuristic,
        "n_bindings_declared": sum(len(getattr(c, "metric_bindings", ())) for c in cells()),
        "n_binding_defects": len(binding_defects),
        "binding_defects": binding_defects,
        "status_counts": counts,
        "n_mismatch": len(mismatches),
        "mismatches": mismatches,
        "rows": rows,
        "honest_limits": [
            "a field_mismatch is a LEAD for adjudication, not a proven defect -- equal values "
            "legitimately co-occur across quantities in one artifact",
            "label extraction is a HEURISTIC window scan; prose that names the quantity far from the "
            "number is label_unlabeled and counted as NOT measured, never as passing",
            "field_confirmed means a label-consistent field carries that value -- NOT that the "
            "contract's sentence faithfully summarises the artifact",
            "a cited .md has no field structure, so numbers only in markdown are no_numeric_field",
        ],
    }


def main(argv=None) -> int:
    rep = audit_semantics()
    # REFUSE to present a zero-measurable run as a result. The first version of this script reported
    # "180 numbers | 0 measurable | MISMATCH 0" because a JSON loader bug meant no artifact ever
    # parsed -- arithmetically a clean bill of health, actually a check that ran on nothing. A
    # measurable fraction this low is a plumbing signature, not a finding about the contracts.
    if rep["n_numbers"] and rep["measurable_fraction"] is not None and rep["measurable_fraction"] < 0.05:
        print(f"REFUSING a verdict: only {rep['n_measurable']}/{rep['n_numbers']} numbers were "
              f"measurable ({rep['measurable_fraction']}). That is a plumbing signature "
              f"(artifacts not parsing, or labels not extracting), not a clean result.",
              file=sys.stderr)
        return 3
    out = REPO / "wiki" / "contract_number_semantics_2026-09-28.json"
    out.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    print(f"numbers {rep['n_numbers']} | measurable {rep['n_measurable']} "
          f"({rep['measurable_fraction']}) | MISMATCH {rep['n_mismatch']}")
    for s, n in sorted(rep["status_counts"].items(), key=lambda kv: -kv[1]):
        print(f"  {s:22s} {n}")
    if rep["mismatches"]:
        print("\nLEADS for adjudication (value present, claimed quantity absent):")
        for r in rep["mismatches"]:
            print(f"  {r['cell_id']} :: {r['claimed_label']} {r['number']}")
            for p in r["competing_paths"]:
                print(f"      value also/only at: {p}")
    print(f"\nwrote {out.relative_to(REPO).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
