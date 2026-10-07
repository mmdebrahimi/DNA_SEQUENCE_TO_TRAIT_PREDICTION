"""Does a generator produce sequence that is actually DNA-like — and does it beat a dumb null?

The design doc's original bar was "constraint-pass rate above a 3-mer Markov baseline". The prior-art scan
showed that is far too weak: two independent 2026 audits separate generated from natural sequence at
**93% sens / 97.9% spec** and at **AUROC 0.97 (eukaryote) / 0.82 (prokaryote)** against two different
published genome models, with long-range organisation named as the failure. So passing hard constraints
says almost nothing, and the real question is **can a classifier tell your output from real DNA**.

Three numbers, and the first two exist to make the third trustworthy:

1. **`control_self` — natural vs natural.** Split the real corpus in half and discriminate the halves.
   They are the same distribution, so the honest answer is 0.5. Anything higher is the instrument reacting
   to the SPLIT rather than to the sequence, and it bounds how much signal is noise.
2. **`control_uniform` — natural vs i.i.d. uniform.** Obviously not DNA, so this must come out near 1.0.
   If it does not, the discriminator has no power and its verdict on a real generator is void.
3. **`distinguishability` — natural vs the candidate.** LOWER IS BETTER: 0.5 means indistinguishable from
   real DNA, 1.0 means trivially detectable.

**The tolerance is DERIVED, not asserted.** `control_self` is re-run over several seeds and the MAXIMUM
observed value becomes the noise ceiling; "distinguishable" means exceeding that ceiling. This mirrors the
standing practice in this repo of clearing a null's maximum rather than its p95, and it means no threshold
in this module is a number someone picked.

**`verdict` REFUSES rather than reporting** when either control fails (`INSTRUMENT_INVALID` /
`INSTRUMENT_UNDERPOWERED`). A measurement whose own instrument is broken must not emit a number — a
discriminator AUROC of 0.52 reads as "indistinguishable, excellent" and is equally consistent with "the
classifier learned nothing because the features were empty".

**Honest note on the hard constraints.** `dna_decode/constraints/` is CDS-shaped — reading frame, internal
stop codons, codon reachability. A promoter is NOT a coding sequence, so those constraints are
**INAPPLICABLE here**, and this module says so (`constraint_scope`) instead of reporting a vacuous 100%
pass rate. What IS meaningful for regulatory sequence is checked: alphabet validity and exact length. A
filter that cannot fail is not a control, which this repo has already had to retract twice.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product

# Verdicts
BEATS_NULL = "BEATS_MARKOV_NULL"
NO_BETTER = "NO_BETTER_THAN_MARKOV_NULL"
WORSE = "WORSE_THAN_MARKOV_NULL"
INSTRUMENT_INVALID = "INSTRUMENT_INVALID_SELF_CONTROL_FAILED"
INSTRUMENT_UNDERPOWERED = "INSTRUMENT_UNDERPOWERED_NO_SEPARATING_POWER"
INSUFFICIENT = "INSUFFICIENT_SEQUENCES"

VERDICTS = frozenset({BEATS_NULL, NO_BETTER, WORSE, INSTRUMENT_INVALID,
                      INSTRUMENT_UNDERPOWERED, INSUFFICIENT})

#: Minimum per-class sequences for a held-out discriminator to mean anything.
MIN_PER_CLASS = 40
#: `control_uniform` must reach this or the instrument has no power. Uniform random vs real DNA is a
#: gross difference; anything below this means the features or the split are broken, not that uniform
#: noise is subtly DNA-like.
MIN_POWER_AUROC = 0.90


def kmer_features(seqs: list[str], k: int = 4) -> list[list[float]]:
    """PURE: normalized k-mer frequency vectors. No sklearn, no numpy needed to build them."""
    if k < 1:
        raise ValueError("k must be >= 1")
    vocab = {"".join(p): i for i, p in enumerate(product("ACGT", repeat=k))}
    rows = []
    for s in seqs:
        s = s.upper()
        v = [0.0] * len(vocab)
        tot = 0
        for i in range(len(s) - k + 1):
            j = vocab.get(s[i:i + k])
            if j is not None:
                v[j] += 1.0
                tot += 1
        if tot:
            v = [x / tot for x in v]
        rows.append(v)
    return rows


def positional_features(seqs: list[str], *, bins: int = 10, k: int = 1) -> list[list[float]]:
    """PURE: per-POSITION-BIN composition — features that know WHERE a motif sits, not just that it exists.

    Why this exists, and it was measured rather than assumed. A k-mer histogram is TRANSLATION-INVARIANT,
    and so is a Markov chain: both describe local composition with no notion of position. Measured on real
    *E. coli* promoters, an order-3 Markov null reached distinguishability 0.5747 against a derived noise
    floor of 0.5447 — i.e. a trivial null very nearly SATURATES a k-mer discriminator, leaving ~0.03 AUROC
    of headroom in which to rank a real generator. That is a ceiling effect in the instrument, not a
    property of the generator.

    A real promoter has ARCHITECTURE: the -10 and -35 boxes sit at characteristic distances from the
    start codon. Splitting the window into bins and taking each bin's composition makes that visible,
    which is the cheap end of the "long-range organisation" axis the published audits name as the thing
    generated sequence fails to reproduce.

    Returns `bins * 4**k` features: each bin's normalized k-mer frequencies.
    """
    if bins < 1:
        raise ValueError("bins must be >= 1")
    if k < 1:
        raise ValueError("k must be >= 1")
    vocab = {"".join(p): i for i, p in enumerate(product("ACGT", repeat=k))}
    width = len(vocab)
    rows = []
    for s in seqs:
        s = s.upper()
        n = len(s)
        v = [0.0] * (bins * width)
        for b in range(bins):
            lo = (n * b) // bins
            hi = (n * (b + 1)) // bins
            seg = s[lo:hi]
            tot = 0
            base = b * width
            for i in range(len(seg) - k + 1):
                j = vocab.get(seg[i:i + k])
                if j is not None:
                    v[base + j] += 1.0
                    tot += 1
            if tot:
                for j in range(width):
                    v[base + j] /= tot
        rows.append(v)
    return rows


def build_features(seqs: list[str], *, mode: str = "kmer", k: int = 4, bins: int = 10
                   ) -> list[list[float]]:
    """Feature builder. `mode` is 'kmer' (translation-invariant), 'positional', or 'both'."""
    if mode == "kmer":
        return kmer_features(seqs, k=k)
    if mode == "positional":
        return positional_features(seqs, bins=bins, k=min(k, 2))
    if mode == "both":
        a = kmer_features(seqs, k=k)
        b = positional_features(seqs, bins=bins, k=min(k, 2))
        return [x + y for x, y in zip(a, b)]
    raise ValueError(f"mode must be 'kmer', 'positional' or 'both', got {mode!r}")


def heldout_auroc(a: list[str], b: list[str], *, k: int = 4, seed: int = 0,
                  mode: str = "kmer", bins: int = 10) -> float | None:
    """AUROC of a logistic-regression discriminator between two sequence sets, on a HELD-OUT half.

    In-sample separation is meaningless here — with 4^k features and a few hundred sequences a classifier
    can memorise. So the data is split 50/50 by a seeded shuffle, fitted on one half and scored on the
    other. Returns None when either class is too small to support the split.
    """
    import random as _random

    if len(a) < 4 or len(b) < 4:
        return None
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import roc_auc_score
        from sklearn.preprocessing import StandardScaler
    except ImportError:  # pragma: no cover - sklearn is a project dep
        return None

    X = build_features(list(a) + list(b), mode=mode, k=k, bins=bins)
    y = [0] * len(a) + [1] * len(b)
    idx = list(range(len(y)))
    _random.Random(seed).shuffle(idx)
    cut = len(idx) // 2
    tr, te = idx[:cut], idx[cut:]
    ytr = [y[i] for i in tr]
    yte = [y[i] for i in te]
    if len(set(ytr)) < 2 or len(set(yte)) < 2:
        return None
    Xa = np.asarray(X, dtype=float)
    sc = StandardScaler().fit(Xa[tr])
    clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Xa[tr]), ytr)
    p = clf.predict_proba(sc.transform(Xa[te]))[:, 1]
    return float(roc_auc_score(yte, p))


def self_control_ceiling(natural: list[str], *, k: int = 4, seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
                         mode: str = "kmer", bins: int = 10) -> tuple[float | None, list[float]]:
    """DERIVE the noise ceiling: split the natural corpus against ITSELF over several seeds.

    Returns `(max_observed, all_values)`. The maximum — not the mean — is the ceiling a candidate must
    exceed to count as distinguishable, so a lucky split cannot manufacture a finding.

    **Re-derived per feature mode.** A richer feature set has more capacity to fit split noise, so its
    ceiling is higher; reusing a k-mer ceiling to judge a positional measurement would understate the
    noise and manufacture a finding.
    """
    import random as _random

    vals: list[float] = []
    for s in seeds:
        shuf = list(natural)
        _random.Random(1000 + s).shuffle(shuf)
        half = len(shuf) // 2
        auc = heldout_auroc(shuf[:half], shuf[half:2 * half], k=k, seed=s, mode=mode, bins=bins)
        if auc is not None:
            vals.append(auc)
    if not vals:
        return None, []
    return max(vals), vals


def alphabet_and_length_ok(seqs: list[str], length: int) -> dict:
    """The constraint checks that are actually MEANINGFUL for non-coding sequence.

    The CDS-shaped universal constraints (frame, internal stop, codon reachability) do NOT apply to a
    promoter, so they are deliberately not run and `constraint_scope` records why.
    """
    ok_alpha = sum(1 for s in seqs if s and not (set(s.upper()) - set("ACGT")))
    ok_len = sum(1 for s in seqs if len(s) == length)
    n = len(seqs) or 1
    return {
        "n": len(seqs),
        "alphabet_valid_rate": round(ok_alpha / n, 4),
        "exact_length_rate": round(ok_len / n, 4),
        "constraint_scope": (
            "dna_decode/constraints is CDS-shaped (frame / internal stop / codon reachability) and is "
            "INAPPLICABLE to non-coding regulatory sequence, so it is not run here and no vacuous "
            "100%-pass is reported. Alphabet + exact length are the checks that do apply."
        ),
    }


@dataclass
class FalsifierResult:
    """One generator's result. `verdict` is REFUSAL-capable by design."""

    candidate_name: str
    n_natural: int
    n_candidate: int
    k: int
    feature_mode: str = "kmer"
    bins: int | None = None
    distinguishability: float | None = None
    null_distinguishability: float | None = None
    control_self_max: float | None = None
    control_self_values: list = field(default_factory=list)
    control_uniform: float | None = None
    gc_natural: float | None = None
    gc_candidate: float | None = None
    checks: dict = field(default_factory=dict)
    verdict: str = INSUFFICIENT
    reason: str = ""

    def as_dict(self) -> dict:
        return {
            "record": "glm-falsifier-v1",
            "candidate": self.candidate_name,
            "n_natural": self.n_natural,
            "n_candidate": self.n_candidate,
            "kmer_k": self.k,
            "feature_mode": self.feature_mode,
            "positional_bins": self.bins,
            "distinguishability_auroc_LOWER_IS_BETTER": self.distinguishability,
            "markov_null_distinguishability": self.null_distinguishability,
            "control_self_natural_vs_natural_max": self.control_self_max,
            "control_self_values": self.control_self_values,
            "control_uniform_vs_natural": self.control_uniform,
            "gc_natural": self.gc_natural,
            "gc_candidate": self.gc_candidate,
            "checks": self.checks,
            "verdict": self.verdict,
            "reason": self.reason,
        }


def run_falsifier(
    natural: list[str],
    candidate: list[str],
    *,
    candidate_name: str,
    null_sequences: list[str] | None = None,
    uniform_sequences: list[str] | None = None,
    k: int = 4,
    length: int | None = None,
    seed: int = 7,
    mode: str = "kmer",
    bins: int = 10,
) -> FalsifierResult:
    """Score one candidate generator's output against the natural corpus, controls first.

    `null_sequences` is the Markov null's output; omit it only when scoring the null itself.
    `mode` selects the feature set: 'kmer' (translation-invariant), 'positional' (architecture-aware)
    or 'both'. Every control is re-derived under the SAME mode.
    """
    from dna_decode.glm.corpus import gc_fraction

    length = length if length is not None else (len(candidate[0]) if candidate else 0)
    res = FalsifierResult(candidate_name=candidate_name, n_natural=len(natural),
                          n_candidate=len(candidate), k=k)
    res.feature_mode = mode
    res.bins = bins if mode in ("positional", "both") else None
    res.checks = alphabet_and_length_ok(candidate, length)
    if natural:
        res.gc_natural = round(sum(gc_fraction(s) for s in natural) / len(natural), 4)
    if candidate:
        res.gc_candidate = round(sum(gc_fraction(s) for s in candidate) / len(candidate), 4)

    if len(natural) < MIN_PER_CLASS or len(candidate) < MIN_PER_CLASS:
        res.verdict = INSUFFICIENT
        res.reason = (f"need >= {MIN_PER_CLASS} per class for a held-out discriminator; "
                      f"have natural={len(natural)}, candidate={len(candidate)}")
        return res

    # --- controls FIRST: the instrument is validated before it is used -------------------------------
    ceiling, vals = self_control_ceiling(natural, k=k, mode=mode, bins=bins)
    res.control_self_max = None if ceiling is None else round(ceiling, 4)
    res.control_self_values = [round(v, 4) for v in vals]
    if uniform_sequences:
        u = heldout_auroc(natural, uniform_sequences, k=k, seed=seed, mode=mode, bins=bins)
        res.control_uniform = None if u is None else round(u, 4)

    if ceiling is None:
        res.verdict = INSTRUMENT_INVALID
        res.reason = "natural-vs-natural control could not be computed; no noise ceiling exists"
        return res
    if res.control_uniform is not None and res.control_uniform < MIN_POWER_AUROC:
        res.verdict = INSTRUMENT_UNDERPOWERED
        res.reason = (f"uniform-vs-natural control reached only {res.control_uniform} "
                      f"(< {MIN_POWER_AUROC}); the discriminator cannot separate i.i.d. noise from real "
                      "DNA, so any verdict on a real generator would be unsupported")
        return res

    # --- the measurement ----------------------------------------------------------------------------
    d = heldout_auroc(natural, candidate, k=k, seed=seed, mode=mode, bins=bins)
    res.distinguishability = None if d is None else round(d, 4)
    if null_sequences and len(null_sequences) >= MIN_PER_CLASS:
        nd = heldout_auroc(natural, null_sequences, k=k, seed=seed, mode=mode, bins=bins)
        res.null_distinguishability = None if nd is None else round(nd, 4)

    if res.distinguishability is None:
        res.verdict = INSUFFICIENT
        res.reason = "discriminator could not be fitted on the candidate set"
        return res

    if res.null_distinguishability is None:
        # Scoring the null itself, or no null supplied: report against the derived ceiling only.
        res.verdict = NO_BETTER
        res.reason = (f"no Markov null supplied for comparison; candidate distinguishability "
                      f"{res.distinguishability} against a derived noise ceiling of {res.control_self_max}")
        return res

    # The null's own spread is unknown, so the comparison tolerance reuses the DERIVED noise ceiling's
    # distance above chance -- how much apparent separation this instrument produces on identical
    # distributions -- rather than an asserted epsilon.
    tol = max(0.0, res.control_self_max - 0.5)
    if res.distinguishability < res.null_distinguishability - tol:
        res.verdict = BEATS_NULL
    elif res.distinguishability > res.null_distinguishability + tol:
        res.verdict = WORSE
    else:
        res.verdict = NO_BETTER
    res.reason = (f"candidate {res.distinguishability} vs Markov null {res.null_distinguishability}; "
                  f"tolerance {round(tol, 4)} derived from the self-control ceiling "
                  f"{res.control_self_max}. LOWER distinguishability is better.")
    return res
