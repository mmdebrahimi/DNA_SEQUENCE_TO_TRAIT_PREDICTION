"""The 10-item ground-truth benchmark for a semantic staleness auditor -- built, not hand-typed.

WHY THIS EXISTS. ~91% of the substantive errors in the 2026-08-25/27 sessions were STALE OR UNVERIFIED
SEMANTIC CLAIMS: a sentence asserting something about an artifact that the artifact no longer supports.
`tests/test_claude_md_citations.py` says in its own docstring that it CANNOT catch this class -- it checks
that a cited FILE EXISTS, not that the file still says what the claim says. The ProSST case is the proof:
`wiki/prosst_lift_2026-07-18.md` existed the whole time; what was stale was the CLAIM ABOUT it.

Before spending GPU on a semantic auditor, it needs a benchmark it can FAIL. This builds one.

THE TWO CLASSES, AND WHERE THEY COME FROM

  POSITIVES (should be flagged STALE) -- 5 real stale claims, each already pinned as a regression test in
  `tests/test_claude_md_citations.py` (4) plus the C. auris one found 2026-08-27. Each is quoted VERBATIM
  as it stood before its correction, with the artifact that refutes it.

  NEGATIVES (must NOT be flagged) -- the 5 hits produced by the MECHANICAL proximity screen that was built,
  run, and deliberately NOT shipped because all 5 were false positives. Those hits were described in prose
  but NEVER PERSISTED, so this module RE-RUNS that screen to regenerate them rather than hand-copying a
  list from a comment. That matters: a hand-typed negative set would be my own recollection of what the
  screen found, and the whole point of this benchmark is to not trust recollection.

THE BAR. The mechanical screen scored 0/5 true positives and 5/5 false positives. Pre-registered pass
condition for any candidate auditor: >=3/5 TP AND <=1/5 FP. Set there because the value is catching real
staleness and the cost of a false positive is one adjudication -- but a guard whose noise exceeds ~1-in-5
gets disabled, which is exactly how the mechanical version died.

Run: uv run python scripts/staleness_benchmark.py [--json]
"""
from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLAUDE_MD = ROOT / "CLAUDE.md"


@dataclass(frozen=True)
class Item:
    """One benchmark item: a claim, its cited artifact, and the ground-truth label."""

    item_id: str
    claim: str                 # the claim text AS IT STOOD (verbatim where it was stale)
    artifact: str              # the repo path the claim is about
    label: str                 # "stale" | "supported"
    why: str                   # the ground truth, in one line
    pinned_by: str = ""        # the regression test that records it, when there is one
    source: str = "curated"    # "curated" | "regenerated"


# ---------------------------------------------------------------------------- POSITIVES (label=stale)
# Each was REAL, shipped, and repeated to readers until caught. Quoted as it stood BEFORE correction.

POSITIVES: tuple[Item, ...] = (
    Item(
        item_id="P1_prosst_deferred",
        claim="the real forward pass is deferred to a Kaggle run",
        artifact="wiki/prosst_lift_2026-07-18.md",
        label="stale",
        why=("it RAN the same day it shipped, LOCALLY on CPU -- the artifact reports Spearman 1.0000 "
             "column reproduction over 56 proteins and a paired hybrid lift of +0.0668 (52/56, p=1.1e-11)"),
        pinned_by="test_the_prosst_forward_pass_is_not_described_as_deferred",
    ),
    Item(
        item_id="P2_tb_pending",
        claim="PENDING DATA RUNS (BLOCKED-gated by design, not incomplete code)",
        artifact="wiki/tb_independent_amr_portal_scores.json",
        label="stale",
        why=("both blockers had been resolved -- the CRyPTIC parquet adapter sidestepped the regeno fetch "
             "and the EBI AMR-Portal cohort (N=2,845) delivered the independent number"),
        pinned_by="test_the_tb_pending_data_runs_header_records_that_both_blockers_resolved",
    ),
    Item(
        item_id="P3_genome_map_browser",
        claim="a visual browser deferred",
        artifact="dna_decode/genome_map/browser.py",
        label="stale",
        why=("the browser SHIPPED 2026-07-11 -- the module exists, and the very next bullet in the same "
             "file said 'GRAPHICAL BROWSER SHIPPED'"),
        pinned_by="test_the_genome_map_browser_is_not_described_as_deferred",
    ),
    Item(
        item_id="P4_bvbrc_census",
        claim="strict-MIC 3-drug feasibility census ... deferred to a fresh session",
        artifact="wiki/bvbrc_strict_mic_4drug_census_2026-05-18.md",
        label="stale",
        why="it RAN 2026-05-18 as a FOUR-drug census -- stale on both the drug count and the status",
        pinned_by="test_the_bvbrc_census_is_not_described_as_a_deferred_3_drug_run",
    ),
    Item(
        item_id="P5_caur_no_free_source",
        claim="no free isolate-level phenotype source exists (fungal)",
        artifact="wiki/ar_bank_caur_powered_result_2026-07-20.md",
        label="stale",
        why=("the CDC AR Isolate Bank produced a POWERED result -- 12 isolates (5R/7S), sens 1.00, and "
             "1.00/1.00 on the mechanism-attributable subset"),
        pinned_by="(found 2026-08-27; see wiki/label_constraint_deep_dive_2026-08-27.md)",
    ),
)

# ---------------------------------------------------- NEGATIVES (label=supported) -- REGENERATED, not typed
# The mechanical screen: "a deferral marker within N chars of an existing repo path is a contradiction".
# It returned 5 hits, ALL false positives, and was not shipped. Re-run it to rebuild the negative set.

_DEFERRAL = re.compile(r"\b(?:deferred|pending|blocked|not yet|still to|TODO)\b", re.I)
_PATH = re.compile(r"`([A-Za-z0-9_./-]+/[A-Za-z0-9_./-]+\.(?:py|md|json|yaml|tsv|fna))`")
_WINDOW = 220


def mechanical_screen(text: str, window: int = _WINDOW) -> list[dict]:
    """THE SCREEN THAT FAILED -- a deferral marker within `window` chars of an existing repo path.

    Reproduced here as the negative-set generator, and as the BASELINE any candidate must beat. It is
    kept verbatim in spirit rather than improved: its 0/5 TP and 5/5 FP is the bar, and quietly making it
    better would erase the comparison.
    """
    hits = []
    for m in _PATH.finditer(text):
        path = m.group(1)
        if not (ROOT / path).exists():
            continue                       # the file-exists guard already covers a missing path
        lo, hi = max(0, m.start() - window), min(len(text), m.end() + window)
        region = text[lo:hi]
        d = _DEFERRAL.search(region)
        if d:
            hits.append({"artifact": path, "marker": d.group(0), "region": " ".join(region.split())})
    return hits


def regenerate_negatives(limit: int = 5) -> list[Item]:
    """Rebuild the negative set by re-running the failed screen over the live CLAUDE.md.

    Every hit is a CONFIRMED false positive by construction: the screen's own verdict was that all of its
    hits were false. Labelling them `supported` is therefore ground truth, not a fresh judgment call.
    """
    if not CLAUDE_MD.exists():
        return []
    hits = mechanical_screen(CLAUDE_MD.read_text(encoding="utf-8", errors="replace"))
    out = []
    for i, h in enumerate(hits[:limit], 1):
        out.append(Item(
            item_id=f"N{i}_mechanical_fp",
            claim=h["region"][:400],
            artifact=h["artifact"],
            label="supported",
            why=(f"mechanical FALSE POSITIVE: the marker {h['marker']!r} and the path sit in the same "
                 f"prose region but refer to different things (often the marker is inside CORRECTION "
                 f"text explaining a fix). The screen that produced this scored 5/5 FP and was not shipped"),
            source="regenerated",
        ))
    return out


# ------------------------------------------------------------------------------------- the benchmark

# Pre-registered, and deliberately low: the baseline to beat is 0/5 TP with 5/5 FP.
PASS_MIN_TP = 3
PASS_MAX_FP = 1


@dataclass
class Score:
    tp: int = 0
    fn: int = 0
    fp: int = 0
    tn: int = 0
    missing: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.tp >= PASS_MIN_TP and self.fp <= PASS_MAX_FP

    def as_dict(self) -> dict:
        return {"tp": self.tp, "fn": self.fn, "fp": self.fp, "tn": self.tn,
                "n_positives": self.tp + self.fn, "n_negatives": self.fp + self.tn,
                "pass_condition": f"tp>={PASS_MIN_TP} and fp<={PASS_MAX_FP}",
                "passed": self.passed, "unanswered": self.missing}


def build() -> list[Item]:
    """The full benchmark: curated positives + regenerated negatives."""
    return list(POSITIVES) + regenerate_negatives()


def score(items: list[Item], verdicts: dict[str, str]) -> Score:
    """Score {item_id: 'stale'|'supported'|'unclear'} against ground truth. PURE.

    `unclear` counts as NOT flagged -- an auditor that hedges on a real stale claim has missed it, and
    hedging on a good claim costs nothing. That asymmetry is deliberate: it prevents a model from passing
    by marking everything `unclear`.
    """
    s = Score()
    for it in items:
        v = verdicts.get(it.item_id)
        if v is None:
            s.missing.append(it.item_id)
            v = "unclear"
        flagged = (v == "stale")
        if it.label == "stale":
            s.tp += flagged
            s.fn += not flagged
        else:
            s.fp += flagged
            s.tn += not flagged
    return s


def baseline_verdicts(items: list[Item]) -> dict[str, str]:
    """What the MECHANICAL screen says about each item -- the bar to beat, computed not asserted.

    It flags exactly what its regex finds: every regenerated negative (by construction) and any positive
    whose text happens to trip it.
    """
    out = {}
    for it in items:
        if it.source == "regenerated":
            out[it.item_id] = "stale"                    # the screen flagged these -- that is why they exist
        else:
            # would the screen catch this real stale claim? only if the claim text itself trips it AND
            # cites an existing path -- which is precisely what it could not do.
            trips = bool(_DEFERRAL.search(it.claim)) and bool(_PATH.search(it.claim))
            out[it.item_id] = "stale" if trips else "supported"
    return out


def main() -> int:
    items = build()
    base = score(items, baseline_verdicts(items))

    if "--json" in sys.argv:
        print(json.dumps({
            "_schema": "staleness-benchmark-v1",
            "pass_min_tp": PASS_MIN_TP, "pass_max_fp": PASS_MAX_FP,
            "mechanical_baseline": base.as_dict(),
            "items": [{"item_id": i.item_id, "claim": i.claim, "artifact": i.artifact,
                       "label": i.label, "why": i.why, "pinned_by": i.pinned_by,
                       "source": i.source} for i in items],
        }, indent=2))
        return 0

    pos = [i for i in items if i.label == "stale"]
    neg = [i for i in items if i.label == "supported"]
    print(f"Staleness benchmark: {len(pos)} positives (curated) + {len(neg)} negatives (regenerated)\n")
    for i in items:
        print(f"  [{i.label:9s}] {i.item_id:24s} {i.artifact}")
    print(f"\nPre-registered pass condition: tp >= {PASS_MIN_TP} AND fp <= {PASS_MAX_FP}")
    print(f"MECHANICAL BASELINE (the bar to beat): {base.as_dict()}")
    if not neg:
        print("\nWARNING: the negative set is EMPTY -- the mechanical screen found no hits to regenerate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
