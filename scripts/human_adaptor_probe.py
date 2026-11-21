"""Measure WHY the full-230 geometry separates: the flanks are constant and the Markov chain misses them.

The margin run's primary geometry reads CHAIN_IS_WEAKER_ON_HUMAN, and stripping 15 bp from each end
collapses it. That before/after drop is suggestive but underdetermined -- it has an obvious rival
explanation: stripping removes 30 of 230 positions (13%), so a weaker discriminator could simply be
a shorter-sequence effect rather than anything about the adaptors.

This probe decides between them by measuring both halves of the mechanism directly:

  1. Are the flanks actually CONSTANT in natural sequence? (asserted from a column description in the
     first draft of the memo; that is exactly the pattern that has to be measured instead)
  2. Does the Markov null actually FAIL to reproduce them?

If the flanks carry ZERO natural variability, the length rival dies: positions that are identical in
every natural sequence cannot contribute to a discriminator's ability to characterise natural
sequence, so removing them cannot weaken it on that account. They can only remove the
constant-versus-not signal that a stochastic generator fails.

Run: uv run python scripts/human_adaptor_probe.py
"""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path

from dna_decode.glm.baselines import MarkovGenerator
from dna_decode.glm.human_natural import ADAPTOR_BP, FULL_BP, fit_eval_by_fold, human_natural_set

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mpra_human"
OUT = ROOT / "wiki" / "human_adaptor_probe_2026-11-21.json"

CELLS = ("HepG2", "K562")
N_NATURAL = 3000
N_GEN = 120
SEED = 7
MARKOV_ORDER = 3


def conservation(seqs: list[str]) -> list[float]:
    """Per-position share of the most common base. 1.0 means perfectly constant."""
    n = len(seqs)
    return [max(Counter(s[i] for s in seqs).values()) / n for i in range(FULL_BP)]


def probe_cell(cell: str) -> dict:
    seqs, _ = human_natural_set(cell, DATA, limit=N_NATURAL)
    five = Counter(s[:ADAPTOR_BP] for s in seqs)
    three = Counter(s[-ADAPTOR_BP:] for s in seqs)
    cons = conservation(seqs)

    fit, _ev, _ = fit_eval_by_fold(cell, DATA, eval_folds={1}, limit_eval=1200)
    gen = MarkovGenerator(order=MARKOV_ORDER, seed=SEED).fit(fit).generate(N_GEN, FULL_BP)
    p5 = five.most_common(1)[0][0]
    p3 = three.most_common(1)[0][0]

    prefix = []
    for s in gen:
        k = 0
        while k < ADAPTOR_BP and s[k] == p5[k]:
            k += 1
        prefix.append(k)

    return {
        "cell": cell,
        "n_natural": len(seqs),
        "natural_5p": {"distinct": len(five), "modal": p5,
                       "modal_share": round(five.most_common(1)[0][1] / len(seqs), 4)},
        "natural_3p": {"distinct": len(three), "modal": p3,
                       "modal_share": round(three.most_common(1)[0][1] / len(seqs), 4)},
        "mean_max_base_share": {
            "flank_5p": round(sum(cons[:ADAPTOR_BP]) / ADAPTOR_BP, 4),
            "interior": round(sum(cons[ADAPTOR_BP:FULL_BP - ADAPTOR_BP])
                              / (FULL_BP - 2 * ADAPTOR_BP), 4),
            "flank_3p": round(sum(cons[FULL_BP - ADAPTOR_BP:]) / ADAPTOR_BP, 4)},
        "markov": {
            "n_generated": len(gen),
            "exact_5p_matches": sum(1 for s in gen if s[:ADAPTOR_BP] == p5),
            "exact_3p_matches": sum(1 for s in gen if s[-ADAPTOR_BP:] == p3),
            "mean_matching_5p_prefix": round(sum(prefix) / len(prefix), 2),
            "max_matching_5p_prefix": max(prefix),
            "of_possible": ADAPTOR_BP},
    }


def main() -> int:
    t0 = time.time()
    rows = [probe_cell(c) for c in CELLS]
    flanks_constant = all(r["natural_5p"]["distinct"] == 1 and r["natural_3p"]["distinct"] == 1
                          for r in rows)
    markov_never = all(r["markov"]["exact_5p_matches"] == 0 and r["markov"]["exact_3p_matches"] == 0
                       for r in rows)
    verdict = ("ADAPTOR_ARTIFACT_MECHANISM_CONFIRMED" if flanks_constant and markov_never
               else "MECHANISM_NOT_CONFIRMED")

    rec = {
        "record": "human-adaptor-probe-v1", "date": "2026-11-21",
        "question": ("Is the full-230 geometry's separation the constant adaptors, or just a "
                     "shorter-sequence effect from removing 30 of 230 positions?"),
        "verdict": verdict,
        "flanks_are_constant": flanks_constant,
        "markov_never_reproduces_them": markov_never,
        "length_rival_status": (
            "REFUTED. The stripped positions carry ZERO natural variability (max-base share 1.000), "
            "so they cannot contribute to a discriminator's ability to characterise natural "
            "sequence; removing them can only remove the constant-versus-not signal that a "
            "stochastic generator fails." if flanks_constant else "not addressed"),
        "scope": ("Explains the FULL-230 rows only. It says nothing about the element-200 rows, "
                  "where the cells genuinely disagree, and it is NOT evidence about generation."),
        "rows": rows,
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    for r in rows:
        print("%-6s natural flanks: 5p %d distinct (share %.4f) | 3p %d distinct (share %.4f)"
              % (r["cell"], r["natural_5p"]["distinct"], r["natural_5p"]["modal_share"],
                 r["natural_3p"]["distinct"], r["natural_3p"]["modal_share"]))
        print("       max-base share: flank5p %.3f  interior %.3f  flank3p %.3f"
              % (r["mean_max_base_share"]["flank_5p"], r["mean_max_base_share"]["interior"],
                 r["mean_max_base_share"]["flank_3p"]))
        m = r["markov"]
        print("       markov exact flank matches: 5p %d/%d  3p %d/%d | mean prefix %.2f of %d"
              % (m["exact_5p_matches"], m["n_generated"], m["exact_3p_matches"], m["n_generated"],
                 m["mean_matching_5p_prefix"], m["of_possible"]))
    print("\n%s -> %s" % (verdict, OUT.relative_to(ROOT).as_posix()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
