"""Does the generator QUOTE its training corpus? The check the falsifier cannot make.

WHY THIS EXISTS. `heldout_auroc` asks whether a discriminator can separate generated sequence from
natural sequence. A generator that emits VERBATIM training sequence is maximally indistinguishable --
it IS natural sequence -- so it scores perfectly and is worthless. The committed falsifier has no term
that can see this, and the leak gate does not help: it withholds the windows the falsifier SCORES
against, so the model cannot quote those, but nothing stops it reciting the other 60,576 training
windows and passing as natural on their strength.

THE MEASUREMENT. For each candidate sequence, the length of the longest substring it shares verbatim
with the training corpus. Reported against two references that bracket the honest range:

* REAL HELD-OUT windows -- the level genuine E. coli promoter sequence reaches by shared biology
  (conserved motifs, repeated elements, operon structure). This is the floor a generator SHOULD sit at.
* UNIFORM random sequence -- the chance level for a 4-letter alphabet against a 9.7 Mbp corpus.

A generator materially ABOVE the real-held-out level is copying, and its distinguishability score is
measuring recall rather than generation.

CREDIT. The failure mode was found by a sibling session on its own tokenization arms, where the two
BEST-scoring arms were the two quoting at double the natural verbatim length. Its artifact scan
(constant positions, shared prefixes, duplicates, prompt leakage) came back clean and missed this
entirely -- a scan built for a constant-flank artifact cannot see a copying artifact. This is an
INDEPENDENT reimplementation, deliberately: two implementations that agree are worth more here than
one shared one, because the quantity is easy to get subtly wrong (see the separator note below).

THREE TRAPS, each closed here rather than left to the caller.

1. **The corpus must be joined with a NON-ACGT separator.** Concatenating windows directly creates
   junction substrings that exist in no window, inflating every match length for every arm at once --
   so it would not even show up as an arm-vs-arm difference.
2. **A match length is not a copy unless it beats the REAL baseline.** Natural sequence shares long
   substrings with other natural sequence; quoting the median real level is not misconduct. Absolute
   thresholds on verbatim length are therefore meaningless without the real-held-out row, and this
   script REFUSES to emit a verdict without one.
3. **A longest-match statistic is driven by its tail.** The median alone hides a few wholesale copies,
   so the fraction at or above a long threshold is reported beside it, and the maximum too.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

BASES = "ACGT"
#: Separator between corpus windows. Must be outside the DNA alphabet so no match can span a junction.
SEP = "\n"
#: "Long" verbatim share, for the tail fraction. 30 bp is well past anything chance produces at this
#: corpus size and is the length at which a shared substring stops being plausibly coincidental.
LONG_BP = 30


def read_sequences(path: Path) -> list[str]:
    """One-per-line or FASTA. Blank lines and `#` comments skipped; sequence upper-cased."""
    text = path.read_text(encoding="utf-8")
    out: list[str] = []
    if text.lstrip().startswith(">"):
        cur: list[str] = []
        for line in text.splitlines():
            if line.startswith(">"):
                if cur:
                    out.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            out.append("".join(cur).upper())
    else:
        for line in text.splitlines():
            s = line.strip().upper()
            if s and not s.startswith("#"):
                out.append(s)
    return [s for s in out if s]


def longest_shared_substring(seq: str, haystack: str) -> int:
    """Length of the longest substring of `seq` occurring verbatim in `haystack`.

    Binary search on length: if a length-L substring occurs, so does every shorter one, so the
    predicate is monotone. The inner scan breaks on the first hit, which makes feasible lengths cheap
    and bounds the cost to the infeasible half of the search.
    """
    n = len(seq)
    if n == 0:
        return 0
    best = 0
    for i in range(n):
        # No remaining start can beat the incumbent, so stop rather than scan the tail.
        if n - i <= best:
            break
        # Only ever test lengths that would IMPROVE on best; grow while they keep hitting.
        L = best + 1
        while i + L <= n and haystack.find(seq[i:i + L]) >= 0:
            best = L
            L += 1
    return best


def profile(name: str, seqs: list[str], haystack: str, *, verbose: bool = True) -> dict:
    t0 = time.time()
    lens = []
    for i, s in enumerate(seqs):
        lens.append(longest_shared_substring(s, haystack))
        if verbose and (i + 1) % 25 == 0:
            print(f"    {name}: {i + 1}/{len(seqs)} ({time.time() - t0:.0f}s)", flush=True)
    lens.sort()
    return {
        "n": len(lens),
        # Carried so the ABSOLUTE verbatim-reproduction rule can fire without a ratio. Modal length,
        # because a ragged set would otherwise make the ceiling depend on which sequence is longest.
        "seq_len": statistics.mode(len(s) for s in seqs) if seqs else None,
        "median_bp": statistics.median(lens) if lens else None,
        "mean_bp": round(statistics.fmean(lens), 2) if lens else None,
        "max_bp": max(lens) if lens else None,
        "frac_at_or_above_long": (round(sum(1 for x in lens if x >= LONG_BP) / len(lens), 4)
                                  if lens else None),
        "long_bp_threshold": LONG_BP,
        "seconds": round(time.time() - t0, 1),
    }


def uniform_sequences(n: int, length: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice(BASES) for _ in range(length)) for _ in range(n)]


def verdict(gen: dict, real: dict | None, uni: dict | None = None) -> tuple[str, str]:
    """Three outcomes, not two.

    The obvious failure is a generator ABOVE the real level -- it is quoting. But a generator at or
    below the UNIFORM level is its own finding: real promoters share long verbatim stretches with a
    promoter corpus (conserved motifs, paralogues, repeated elements), so a generator that shares no
    more than random sequence does has not captured that shared content either. That is not
    misconduct and it is not a pass on realism; it is a separate, reportable fact, and collapsing it
    into "no copying signal" would hide it.
    """
    if real is None or real["median_bp"] is None:
        return ("NO_VERDICT_WITHOUT_REAL_BASELINE",
                "A verbatim match length is only interpretable against the level real held-out "
                "sequence reaches. Absolute thresholds are meaningless; refusing rather than "
                "inventing one.")
    g, r = gen["median_bp"], real["median_bp"]
    gl, rl = gen["frac_at_or_above_long"], real["frac_at_or_above_long"]

    # ABSOLUTE rule FIRST, before any ratio. Measured 2026-11-22: a positive control of 10 verbatim
    # TRAINING windows scored median 150.0 of 150 with 100% long shares -- total memorisation -- and
    # the ratio rule returned ELEVATED_NOT_CONCLUSIVE because 150/90.5 = 1.66 < 2. A ratio is
    # undefined at the extreme precisely because the baseline it divides by is itself sampled, and
    # the worst case must not depend on a baseline draw. A window reproduced end to end is copying
    # whatever the reference says.
    L = gen.get("seq_len")
    if L and g >= 0.95 * L:
        return ("VERBATIM_REPRODUCTION",
                f"median verbatim share {g} bp of a {L} bp window ({gl:.1%} of sequences share "
                f">={LONG_BP} bp). The generator is reproducing whole training windows. No ratio "
                "against the real baseline is needed or meaningful at this level.")
    if g >= 2 * r or (gl is not None and rl is not None and gl >= 0.10 and gl >= 4 * max(rl, 0.001)):
        return ("COPYING_DETECTED",
                f"median verbatim share {g} bp against real held-out {r} bp, and {gl:.1%} of "
                f"sequences share >={LONG_BP} bp against {rl:.1%} for real sequence. The "
                "distinguishability score for this generator is measuring recall, not generation.")
    if g > r:
        return ("ELEVATED_NOT_CONCLUSIVE",
                f"median verbatim share {g} bp exceeds real held-out {r} bp but not by the 2x "
                "margin this script treats as decisive. Report both rows; do not quote the "
                "distinguishability number without them.")
    u = uni["median_bp"] if uni else None
    # NOT `g <= u`. These are integer-valued medians, so a half-base difference would decide which of
    # two QUALITATIVELY different states is reported -- measured 2026-11-22, where generated 15.5 vs
    # uniform 15.0 flipped the label while the substance was "at chance". The honest test is where g
    # sits in the uniform->real RANGE: within its bottom quarter is at-chance in any useful sense.
    if u is not None and r is not None and r > u and (g - u) <= 0.25 * (r - u):
        return ("AT_CHANCE_LEVEL_NO_SHARED_CONTENT",
                f"median verbatim share {g} bp against real held-out {r} bp and uniform chance "
                f"{u} bp. NOT copying -- but also no more verbatim overlap with the corpus than "
                "random sequence has, so the generator has not captured the shared motif/paralogue "
                "content that makes real promoters resemble a promoter corpus. Report both facts.")
    return ("AT_OR_BELOW_REAL_LEVEL",
            f"median verbatim share {g} bp against real held-out {r} bp (uniform chance {u} bp). "
            "No copying signal; the distinguishability score is not explained by recall.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--corpus", required=True, help="training corpus, one window per line or FASTA")
    ap.add_argument("--generated", required=True, help="generated sequences to audit")
    ap.add_argument("--held-out", help="REAL held-out sequences (the baseline; verdict refuses without it)")
    ap.add_argument("--sample-real", type=int, default=120,
                    help="subsample the held-out set to this many (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", help="write a JSON artifact here")
    a = ap.parse_args()

    corpus = read_sequences(Path(a.corpus))
    generated = read_sequences(Path(a.generated))
    if not corpus or not generated:
        print("REFUSING: empty corpus or generated set", file=sys.stderr)
        return 2
    haystack = SEP.join(corpus)
    print(f"corpus {len(corpus)} windows / {len(haystack):,} chars; generated {len(generated)}")

    real: list[str] = []
    if a.held_out:
        real = read_sequences(Path(a.held_out))
        # A held-out set that is actually IN the corpus is not held out. Drop and count.
        inside = sum(1 for s in real if s in corpus_set(corpus))
        if inside:
            print(f"  NOTE: {inside} of {len(real)} 'held-out' sequences are present in the corpus "
                  "verbatim; they are NOT held out and are dropped")
            cs = corpus_set(corpus)
            real = [s for s in real if s not in cs]
        if a.sample_real and len(real) > a.sample_real:
            random.Random(a.seed).shuffle(real)
            real = real[:a.sample_real]

    length = len(generated[0])
    uni = uniform_sequences(min(len(generated), 120), length, a.seed)

    print("profiling (binary search on verbatim length) ...")
    res = {
        "record": "glm-copying-penalty-v1",
        "corpus": {"path": str(a.corpus), "n_windows": len(corpus), "n_chars": len(haystack)},
        "separator": "newline (non-ACGT, so no match spans a window junction)",
        "generated": profile("generated", generated, haystack),
        "uniform_chance": profile("uniform", uni, haystack),
    }
    res["real_held_out"] = profile("real", real, haystack) if real else None
    v, why = verdict(res["generated"], res["real_held_out"], res["uniform_chance"])
    res["verdict"], res["verdict_detail"] = v, why

    print("\n" + "=" * 74)
    print(f"{'set':<22}{'median':>8}{'mean':>8}{'max':>7}{f'  >={LONG_BP}bp':>10}")
    print("-" * 74)
    for key, label in (("real_held_out", "REAL held-out"), ("uniform_chance", "uniform (chance)"),
                       ("generated", "GENERATED")):
        r = res.get(key)
        if not r:
            print(f"{label:<22}{'(absent)':>8}")
            continue
        print(f"{label:<22}{r['median_bp']:>8}{r['mean_bp']:>8}{r['max_bp']:>7}"
              f"{r['frac_at_or_above_long']:>9.1%}")
    print("=" * 74)
    print(f"VERDICT: {v}\n  {why}")

    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=2), encoding="utf-8")
        print(f"\n-> {a.out}")
    return 0


_CS: set[str] | None = None


def corpus_set(corpus: list[str]) -> set[str]:
    global _CS
    if _CS is None:
        _CS = set(corpus)
    return _CS


if __name__ == "__main__":
    raise SystemExit(main())
