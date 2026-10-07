"""The EXPRESSION ORACLE — supervised sequence → measured promoter strength.

This is the scoring head the generative loop has been missing, and it is the one component of
AlphaGenome's design that genuinely transfers to bacteria: **supervised sequence-to-function, learned on
real measured tracks, with a rich representation and a deliberately trivial head.** (AlphaGenome itself is
human/mouse and its heads predict splicing, histone marks and chromatin contacts — quantities a bacterium
does not have. Its *logic* transfers; the model does not. Its licence separately forbids training on its
outputs.)

**Substrate: GEO GSE108535** (Urtecho, Tripp, Insigne, Kim & Kosuri, *Biochemistry* 2019;58(11):1539–1551),
a fully-crossed σ70 promoter grid measured by genomically-integrated MPRA in *E. coli*. Verified by download
rather than from a summary: **10,898 sequence→expression pairs, every sequence exactly 150 bp, ACGT-only,
10,898 distinct, expression spanning 615x.**

**TWO TRAPS, both measured, both encoded here rather than left in prose.**

1. **The join silently loses 65% of the data.** The expression table writes UP-element names with HYPHENS
   (`gourse-326fold-up`) and the barcode map writes them with UNDERSCORES (`gourse_136fold_up`). A naive
   join on `name` returns **3,821** pairs and looks perfectly healthy; normalising `-` to `_` returns all
   **10,898**. `_norm_key` is that fix, and `load_pairs` asserts a minimum join rate so the failure can
   never go quiet again.
2. **The files are SPACE-delimited, not tab.** A `split('\t')` parses zero rows, and a zero-row parse that
   is then averaged reads as "no signal" rather than "no data".

**The split is the whole experiment.** This repo has measured three times over that **random splits inflate
results ~3x** (GB1 1-vs-rest 0.32 -> 0.92 random; one-hot 0.579 -> 0.027 across positional splits; R^2 0.90
-> 0.23 stratified). A designed grid makes the honest alternative available and cheap: the table ships the
ELEMENT IDENTITIES (`UP_element`, `Minus35`, `Spacer`, `Minus10`, `Background`), so entire elements can be
held out. **`leave_element_out` is the bacterial analogue of Enformer's held-out-chromosome split** — it
asks whether the model learned promoter grammar or memorised this grid's cells. A random split cannot tell
those apart.
"""
from __future__ import annotations

import gzip
import statistics
from dataclasses import dataclass, field
from pathlib import Path

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108535/suppl"
EXPR_FILE = "GSE108535_sigma70_variant_data.txt.gz"
MAP_FILE = "GSE108535_barcode_mapping.txt.gz"

#: Where the cached MPRA files live. D:, per the standing rule that heavy files never land on C:.
DEFAULT_CACHE = Path("D:/dna_decode_cache/mpra")

#: Element columns in the expression table, in table order. These ARE the split axes.
ELEMENT_COLUMNS = ("UP_element", "Minus35", "Spacer", "Minus10", "Background")

#: Measured on the real file. A join below this means the key formats diverged again.
MIN_JOIN_RATE = 0.95

SEQ_LENGTH = 150


class ExpressionDataError(RuntimeError):
    """The MPRA substrate could not be loaded in a trustworthy state."""


#: Mutation notation (`35T->A`) also contains a hyphen. It is protected during normalisation.
_ARROW = "->"
_ARROW_SENTINEL = "\x00ARROW\x00"


def _norm_key(k: str) -> str:
    """Normalise a variant name for joining. THE hard-won fix — see the module docstring.

    The expression table and the barcode map disagree on hyphen vs underscore inside the UP-element name
    (`gourse-326fold-up` vs `gourse_136fold_up`), which costs 7,077 of 10,898 pairs if unhandled.

    **The mutation arrow is PROTECTED.** A blanket `replace("-", "_")` also rewrites `35T->A` to `35T_>A`.
    That still joins — it is applied to both sides — but it is lossy, and two genuinely different variants
    could in principle collide once the arrow stops being distinguishable from a separator. Protecting the
    arrow keeps the key injective over the parts that carry meaning.
    """
    return k.replace(_ARROW, _ARROW_SENTINEL).replace("-", "_").replace(_ARROW_SENTINEL, _ARROW)


@dataclass
class Pair:
    """One measured promoter."""

    name: str
    seq: str
    expression: float
    elements: dict = field(default_factory=dict)


@dataclass
class LoadReport:
    """What the load actually produced. Reported so a degraded substrate is visible, never inferred."""

    n_expression_rows: int = 0
    n_sequences: int = 0
    n_joined: int = 0
    join_rate: float = 0.0
    n_wrong_length: int = 0
    n_non_acgt: int = 0
    naive_join: int = 0

    def as_dict(self) -> dict:
        return {
            "n_expression_rows": self.n_expression_rows,
            "n_sequences_available": self.n_sequences,
            "n_joined": self.n_joined,
            "join_rate": round(self.join_rate, 4),
            "n_wrong_length_dropped": self.n_wrong_length,
            "n_non_acgt_dropped": self.n_non_acgt,
            "naive_join_would_have_given": self.naive_join,
            "note": ("naive_join is what a join without hyphen/underscore normalisation returns; the gap "
                     "is the silent 65% data loss this loader exists to prevent"),
        }


def ensure_cached(cache_dir: Path | str = DEFAULT_CACHE, *, allow_download: bool = True) -> Path:
    """Make sure both GSE108535 files are present locally; download them if not. Returns the directory."""
    d = Path(cache_dir)
    d.mkdir(parents=True, exist_ok=True)
    missing = [f for f in (EXPR_FILE, MAP_FILE) if not (d / f).exists()]
    if missing and not allow_download:
        raise ExpressionDataError(f"missing {missing} in {d} and downloading is disabled")
    for f in missing:
        import urllib.request
        url = f"{GEO_BASE}/{f}"
        tmp = d / (f + ".part")
        urllib.request.urlretrieve(url, tmp)       # noqa: S310 - fixed NCBI FTP-over-HTTPS URL
        tmp.replace(d / f)
    return d


def _rows(path: Path):
    """Yield whitespace-split fields. SPACE-delimited, not tab — see the module docstring."""
    with gzip.open(path, "rt", errors="replace") as fh:
        header = next(fh).split()
        for line in fh:
            parts = line.split()
            if parts:
                yield header, parts


def load_pairs(cache_dir: Path | str = DEFAULT_CACHE, *, allow_download: bool = True,
               ) -> tuple[list[Pair], LoadReport]:
    """Load (sequence, measured expression, element labels) triples. Fails loudly on a degraded join."""
    d = ensure_cached(cache_dir, allow_download=allow_download)
    rep = LoadReport()

    expr: dict[str, tuple[float, dict]] = {}
    for header, p in _rows(d / EXPR_FILE):
        if len(p) < 7:
            continue
        try:
            val = float(p[6])
        except ValueError:
            continue
        expr[p[0]] = (val, dict(zip(ELEMENT_COLUMNS, p[1:6])))
    rep.n_expression_rows = len(expr)

    seqs: dict[str, str] = {}
    for _h, p in _rows(d / MAP_FILE):
        if len(p) >= 3:
            seqs.setdefault(p[2], p[1])
    rep.n_sequences = len(seqs)

    rep.naive_join = sum(1 for n in expr if n in seqs)

    by_norm: dict[str, str] = {}
    for n, s in seqs.items():
        by_norm.setdefault(_norm_key(n), s)

    out: list[Pair] = []
    for name, (val, elems) in expr.items():
        s = by_norm.get(_norm_key(name))
        if s is None:
            continue
        s = s.upper()
        if len(s) != SEQ_LENGTH:
            rep.n_wrong_length += 1
            continue
        if set(s) - set("ACGT"):
            rep.n_non_acgt += 1
            continue
        out.append(Pair(name=name, seq=s, expression=val, elements=elems))

    rep.n_joined = len(out)
    rep.join_rate = len(out) / rep.n_expression_rows if rep.n_expression_rows else 0.0
    if rep.join_rate < MIN_JOIN_RATE:
        raise ExpressionDataError(
            f"join rate {rep.join_rate:.3f} < {MIN_JOIN_RATE}: only {rep.n_joined} of "
            f"{rep.n_expression_rows} expression rows matched a sequence. The key formats have diverged "
            f"again (a naive join gives {rep.naive_join}); fix _norm_key rather than lowering this floor."
        )
    return out, rep


# ---------------------------------------------------------------------------------------------------
# splits -- the part that decides whether the number means anything
# ---------------------------------------------------------------------------------------------------
def random_split(pairs: list[Pair], *, frac: float = 0.25, seed: int = 0):
    """A random held-out split. **Reported for comparison, NOT as the headline.**

    Measured three times in this repo: random splits inflate results ~3x versus a biologically-motivated
    split, because near-identical sequences land on both sides.
    """
    import random as _r
    idx = list(range(len(pairs)))
    _r.Random(seed).shuffle(idx)
    cut = int(len(idx) * (1 - frac))
    return [pairs[i] for i in idx[:cut]], [pairs[i] for i in idx[cut:]]


def leave_element_out(pairs: list[Pair], column: str, *, held: set[str] | None = None,
                      frac: float = 0.25, seed: int = 0):
    """Hold out ENTIRE element identities — the honest generalisation test.

    The bacterial analogue of a held-out-chromosome split: the test set contains only promoters built from
    element variants the model never saw, so it measures learned grammar rather than grid interpolation.
    Disjointness is BY CONSTRUCTION — membership is a pure function of the element label, so no sequence
    can appear on both sides.
    """
    if column not in ELEMENT_COLUMNS:
        raise ValueError(f"column must be one of {ELEMENT_COLUMNS}, got {column!r}")
    values = sorted({p.elements.get(column, "") for p in pairs})
    if len(values) < 2:
        raise ValueError(f"column {column!r} has <2 distinct values; cannot split on it")
    if held is None:
        import random as _r
        vs = list(values)
        _r.Random(seed).shuffle(vs)
        k = max(1, int(round(len(vs) * frac)))
        held = set(vs[:k])
    train = [p for p in pairs if p.elements.get(column, "") not in held]
    test = [p for p in pairs if p.elements.get(column, "") in held]
    return train, test, sorted(held)


# ---------------------------------------------------------------------------------------------------
# features
# ---------------------------------------------------------------------------------------------------
def gc_feature(seqs: list[str]) -> list[list[float]]:
    """GC content alone. The crudest possible baseline — and this repo has measured that crude baselines
    win more often than not, so it is mandatory rather than decorative."""
    return [[(s.count("G") + s.count("C")) / len(s)] for s in seqs]


def kmer_feature(seqs: list[str], k: int = 4) -> list[list[float]]:
    """Normalised k-mer frequencies. Reuses the falsifier's definition so the two cannot drift apart."""
    from dna_decode.glm.falsifier import kmer_features
    return kmer_features(seqs, k=k)


def positional_onehot(seqs: list[str]) -> list[list[float]]:
    """Flat per-position one-hot. Position-aware and linear-model-friendly.

    This is the representation that matters for a promoter: the -10 and -35 boxes sit at characteristic
    offsets, and a translation-invariant k-mer histogram cannot see where anything is.
    """
    idx = {b: i for i, b in enumerate("ACGT")}
    out = []
    for s in seqs:
        v = [0.0] * (len(s) * 4)
        for i, b in enumerate(s):
            j = idx.get(b)
            if j is not None:
                v[i * 4 + j] = 1.0
        out.append(v)
    return out


FEATURES = {
    "gc": gc_feature,
    "kmer4": lambda s: kmer_feature(s, k=4),
    "kmer6": lambda s: kmer_feature(s, k=6),
    "onehot": positional_onehot,
}


# ---------------------------------------------------------------------------------------------------
# the oracle
# ---------------------------------------------------------------------------------------------------
def spearman(a: list[float], b: list[float]) -> float:
    """Spearman rho via scipy, which handles ties correctly. Ties matter: this assay has repeated values."""
    from scipy.stats import spearmanr
    r = spearmanr(a, b).statistic
    return float(r) if r == r else 0.0


def fit_predict_ridge(train: list[Pair], test: list[Pair], *, feature: str = "onehot",
                      alpha: float = 1.0, log_target: bool = True,
                      return_pred: bool = False) -> dict:
    """Ridge regression on `feature`. The TRIVIAL HEAD half of AlphaGenome's design principle.

    `log_target` is on by default and is not cosmetic: expression spans 615x, so a linear model on the raw
    scale is dominated by the few brightest promoters. Spearman is rank-based and therefore invariant to
    the monotone transform, but the FIT is not.

    `return_pred` adds `pred` + `measured` to the result. OFF by default on purpose: callers that dump this
    dict straight into a JSON artifact would otherwise bloat it with thousands of floats. Downstream
    consumers (the selection measurement) need the predictions and must come through THIS function rather
    than re-fitting, so the two cannot drift apart.
    """
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    if feature not in FEATURES:
        raise ValueError(f"unknown feature {feature!r}; have {sorted(FEATURES)}")
    fn = FEATURES[feature]
    Xtr = np.asarray(fn([p.seq for p in train]), dtype=float)
    Xte = np.asarray(fn([p.seq for p in test]), dtype=float)
    ytr = np.asarray([p.expression for p in train], dtype=float)
    yte = [p.expression for p in test]
    if log_target:
        ytr = np.log1p(ytr)
    sc = StandardScaler().fit(Xtr)
    model = Ridge(alpha=alpha).fit(sc.transform(Xtr), ytr)
    pred = model.predict(sc.transform(Xte))
    out = {
        "feature": feature,
        "n_features": int(Xtr.shape[1]),
        "n_train": len(train),
        "n_test": len(test),
        "spearman": round(spearman(list(pred), yte), 4),
        "alpha": alpha,
        "log_target": log_target,
    }
    if return_pred:
        out["pred"] = [float(v) for v in pred]
        out["measured"] = [float(v) for v in yte]
    return out
