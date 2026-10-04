"""Would extending LABEL_VOCAB convert more contract numbers, or manufacture false leads?

WHY THIS EXISTS. After one vocabulary entry (`freq`) both unlocked bindings AND converted a number via the
prose heuristic unaided, the conclusion recorded was "vocabulary coverage has leverage; bindings do not".
That generalised from n=1. This probe measures the lever across all remaining `label_unlabeled` numbers
before any entry is added, and its answer was the opposite: **no addable entry at the coverage of that
moment.**

CORRECTED THE SAME DAY, and the correction is the more useful result. "Exhausted" was wrong: citing
pneumoserotype's two artifacts moved 5 numbers into the cited set and `concordance` immediately surfaced as
a clean candidate (2 near / 2 would-confirm / 0 would-mismatch) and was added. **The lever's headroom is a
function of WHICH CELLS ARE CITED, so it re-opens whenever coverage grows** -- it is not a fixed quantity
that can be exhausted once. Re-run this probe after any citation change.

IT IS AN UPPER BOUND ON CONVERSIONS, measured: the probe counts a word as converting when it sits in the
window AND the artifact path matches, but `label_for` applies leading-label-wins plus inter-number clipping,
so at most ONE word wins per number. `concordance` was estimated at 2 conversions and delivered 1
(field_confirmed 16 -> 17). Treat its counts as a ceiling, and the `would_mismatch` column likewise.

THE MEASUREMENT, and why both halves are needed. A vocabulary entry only helps when BOTH sides line up:
the prose word must sit near the number AND the cited artifact must carry that value at a path whose tokens
the word matches. If only the prose side matches, the number moves `label_unlabeled` -> `field_mismatch` —
a FALSE LEAD, which is strictly WORSE than unlabeled, because a lead costs adjudication time and an
unlabeled number honestly says "not measured". So this counts `would_confirm` AND `would_mismatch`.

WHAT IT FOUND (2026-09-28, 141 unlabeled numbers):
  * The words that appear to convert cleanly are mostly NOT QUANTITIES — `genes`, `lactam`, `beta`,
    `aminoglycoside`, `star`, `name`, `computed`, `discovery`. They are domain nouns and structural words
    that happen to sit near a number and also appear in a path. Adding them would inflate the measurability
    figure while measuring nothing — the vacuous-control failure this repo has recorded repeatedly.
  * The genuinely quantity-shaped candidates are STATISTICAL MODIFIERS: `mean`, `max`, `null`. These are the
    more dangerous class, because they produce VACUOUS CONFIRMATIONS: prose "max 0.21" against a path
    `max_abs_r` would confirm the quantity as "max" when the number is a CORRELATION.
  * `concordance` was the one genuine quantity that cleared -- 2 estimated / 1 actual conversion, 0 false
    leads -- and it is deliberately NOT an alias of `acc` (agreement between two callers is not accuracy
    against a label; collapsing them would let a tool-agreement figure verify a claim of accuracy). Its
    near-synonym `agreement` was REFUSED at 2 conversions against 1 false lead.
  * The one real quantity found in the first pass, `purity`, scores 1 convert against 4 false mismatches, because the same
    artifact names some fields specifically (`observed_purity`, `null_mean`) and others generically
    (`additional_statistics[0].observed`).

THE STRUCTURE THAT ACTUALLY EXPLAINS IT — and it corrects the earlier ranking. Vocabulary works where the
artifact field name is SPECIFIC; it manufactures false leads where the field name is GENERIC. A generic
field name is not a vocabulary problem and no entry can fix it — it is precisely what a
`ContractNumberBinding` is for, since a binding names the path explicitly. **So vocabulary and bindings are
COMPLEMENTARY, not ranked: vocabulary for specifically-named fields, bindings for generically-named ones.**

Read-only. Exit 0 always; this is a measurement, not a gate.
"""

from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from scripts.contract_number_audit import DECIMAL_RE, PROSE_FIELDS, _cited_artifacts  # noqa: E402
from scripts.contract_number_semantics import (  # noqa: E402
    LABEL_VOCAB,
    WINDOW_BEFORE,
    numeric_leaves,
    value_matches,
)

# Words that are STATISTICAL MODIFIERS rather than quantities. Naming them explicitly is the point: they
# look like the best candidates on frequency and are the worst to add, because a path like `max_abs_r` or
# `null_mean` would "confirm" the modifier while the underlying quantity (a correlation, a purity) goes
# unchecked. A confirmation that does not identify the quantity is vacuous.
MODIFIERS = frozenset({"mean", "max", "min", "median", "null", "total", "overall", "avg", "average",
                       "observed", "expected", "raw", "net", "delta", "gain", "half",
                       # "exact" qualifies a rate ("exact set match rate") -- same vacuity risk
                       "exact"})

STOP = frozenset("""the a an of in on at to for with and or is was are were be been by from as that this
it its their our we you not no but if then than so such which who whom whose when where how why all any
both each few more most other some only own same too very can will just should now over under again
further once here there also per vs versus while does did do done has have had having one two three four
five six seven eight nine ten first second third new old real full every none because about above after
before below between during into through out up down off""".split())


def probe() -> dict:
    from dna_decode.data.cell_registry import cells

    # Resolve the NEWEST artifact. A fixed filename broke the moment the writer started date-stamping,
    # and pinning a date here would reintroduce exactly the staleness that change removed.
    found = sorted((REPO / "wiki").glob("contract_number_semantics_*.json"))
    if not found:
        return {"error": "run scripts/contract_number_semantics.py first "
                         "(no wiki/contract_number_semantics_*.json)"}
    art = found[-1]
    rep = json.loads(art.read_text(encoding="utf-8"))
    unlabeled = [r for r in rep["rows"] if r["status"] == "label_unlabeled"]

    known = {w for words, _f in LABEL_VOCAB.values() for w in words}
    stats: dict[str, dict] = collections.defaultdict(
        lambda: {"near": 0, "would_confirm": 0, "would_mismatch": 0})
    bycell = {getattr(c, "cell_id", ""): c for c in cells()}

    for r in unlabeled:
        c = bycell.get(r["cell_id"])
        if c is None:
            continue
        blob = " ".join(str(getattr(c, f, "") or "") for f in PROSE_FIELDS)
        m = next((m for m in DECIMAL_RE.finditer(blob) if m.group(1) == r["number"]), None)
        if m is None:
            continue
        before = blob[max(0, m.start(1) - WINDOW_BEFORE):m.start(1)].lower()
        words = {w for w in re.split(r"[^a-z]+", before)
                 if len(w) >= 3 and w not in STOP and w not in known}
        if not words:
            continue
        leaves = []
        for a in _cited_artifacts(blob):
            if not a.endswith(".json"):
                continue
            try:
                leaves += [(a, p, v) for p, v in
                           numeric_leaves(json.loads((REPO / a).read_text(encoding="utf-8")))]
            except (OSError, json.JSONDecodeError, ValueError):
                pass
        hits = [(a, p, v) for (a, p, v) in leaves if value_matches(r["number"], v)]
        for w in words:
            st = stats[w]
            st["near"] += 1
            if not hits:
                continue
            if any(w in re.split(r"[^a-z0-9]+", p.lower()) for (_a, p, _v) in hits):
                st["would_confirm"] += 1
            else:
                st["would_mismatch"] += 1
    return {"n_unlabeled": len(unlabeled), "words": dict(stats)}


def verdict_for(word: str, st: dict) -> str:
    """Ordered so the DANGEROUS cases cannot hide behind a benign-sounding label.

    The first version of this function had a real defect worth keeping in mind: its final `else` branch
    printed "no effect (no artifact hit)" for the confirm-but-fewer-than-mismatch case, so `purity` at
    1 convert / 4 false leads read as harmless. A verdict label that understates a net-negative candidate
    is exactly how a bad vocabulary entry gets added.
    """
    c, mm = st["would_confirm"], st["would_mismatch"]
    if word in MODIFIERS:
        return "REFUSE — statistical modifier, would confirm vacuously"
    if c == 0 and mm == 0:
        return "no effect — value never appears in a cited artifact"
    if c == 0:
        return f"REFUSE — {mm} false leads, 0 conversions"
    if mm == 0:
        return f"candidate — {c} clean conversions (check it is a QUANTITY, not a domain noun)"
    if c > mm:
        return f"marginal — {c} convert vs {mm} false leads"
    return f"REFUSE — net negative: {c} convert vs {mm} false leads"


def main(argv=None) -> int:
    rep = probe()
    if "error" in rep:
        print(rep["error"], file=sys.stderr)
        return 2
    rows = sorted(rep["words"].items(),
                  key=lambda kv: (-kv[1]["would_confirm"], -kv[1]["near"]))
    print(f"unlabeled numbers examined: {rep['n_unlabeled']}\n")
    print(f"{'word':18s} {'near':>5s} {'conv':>5s} {'false':>6s}  verdict")
    addable = []
    for w, st in rows:
        if st["near"] < 2 and st["would_confirm"] == 0:
            continue
        v = verdict_for(w, st)
        print(f"{w:18s} {st['near']:5d} {st['would_confirm']:5d} {st['would_mismatch']:6d}  {v}")
        if v.startswith("candidate"):
            addable.append(w)
    print(f"\ncandidates that convert cleanly: {len(addable)} -> {addable}")
    print("NOTE: 'candidate' is NOT 'add it'. Every one must still be a QUANTITY; the cleanly-converting "
          "words in this corpus are domain nouns (genes / lactam / beta / aminoglycoside) and structural "
          "words (star / name / computed / discovery), and adding those would inflate measurability while "
          "measuring nothing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
