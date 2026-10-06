"""Join essentiality LABELS to a genome ANNOTATION, with a gated and composition-reported join rate.

WHY THIS IS THE INTEGRITY CRUX OF THE LADDER. The conserved-core decoder reads PRODUCT TEXT, so a rung
only scores if its label ids join to that organism's annotation. The FBA arm's "walls" for S. aureus and
P. aeruginosa are MODEL-GENE-JOIN walls and genuinely do not apply here -- but the wall does not vanish,
it MOVES to "label ids must join the exact annotation table", which is this module.

TWO RECORDED FAILURES SHAPE IT:
  * `gene_id` vs `gene_symbol` is this repo's most expensive identifier trap: a cross-strain join on the
    strain-unique id gave 0% vocabulary overlap and a confident AUROC of 0.000. A low join rate must be a
    WALL, not a quieter number.
  * A join rate alone cannot see BIAS. 70% of labels joining is fine; 70% joining because the joined
    subset is the conserved core and the unjoined remainder is the lineage-specific tail is not -- and
    the rate is identical in both cases. So composition is reported, not just the rate.

MIN_JOIN_RATE is an ASSERTED bar (0.70), flagged as asserted. Below it a rung emits
WALL_JOIN_RATE_BELOW_FLOOR and carries NO coverage key at all -- absence, not a None sitting next to a
number that looks computed.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass, field

# ASSERTED, not derived. Stated so a reader is never misled into treating it as measured.
MIN_JOIN_RATE = 0.70

WALL_JOIN_RATE_BELOW_FLOOR = "WALL_JOIN_RATE_BELOW_FLOOR"
WALL_ANNOTATION_UNREACHABLE = "WALL_ANNOTATION_UNREACHABLE"

ANNOTATION_KINDS = ("ncbi_feature_table", "sgd_features", "ncbi_gene_info")


class AnnotationError(RuntimeError):
    """Annotation could not be loaded or has the wrong shape."""


def _open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, "rt", encoding="utf-8", errors="replace")


def iter_feature_table(path):
    """CDS rows of an NCBI feature table -> (symbol, product name).

    PROMOTED VERBATIM from dna_decode/essentiality/cli.py::_iter_feature_table so there is exactly one
    definition; the CLI now imports this. Behaviour is byte-identical by construction, pinned by test.
    """
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            isym = header.index("symbol")
            ina = header.index("name")
        except ValueError:
            return
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(isym, ina) and p[0] == "CDS" and p[isym]:
                yield p[isym], p[ina]


def _iter_gene_info(path):
    """NCBI gene_info -> (Symbol, description) for protein-coding genes, keyed later by GeneID."""
    with _open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        try:
            gid, sym, dsc, ty = (h.index("GeneID"), h.index("Symbol"),
                                 h.index("description"), h.index("type_of_gene"))
        except ValueError as e:
            raise AnnotationError("gene_info missing a required column: %s" % e) from e
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(gid, sym, dsc, ty) and p[ty] == "protein-coding":
                yield p[gid], "%s %s" % (p[sym], p[dsc])


def _iter_sgd_features(path):
    """SGD_features.tab -> (systematic ORF name, description).

    Column layout per SGD's own README: [3]=feature type, [5]=systematic name, [15]=description.
    Only ORF rows are kept; anything else cannot carry a deletion-collection label.
    """
    with _open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 16:
                continue
            if p[3].strip().upper() != "ORF":
                continue
            orf = p[5].strip()
            if orf:
                yield orf, "%s %s" % (p[4].strip(), p[15].strip())


def load_annotation(path, kind: str) -> dict[str, str]:
    """identifier -> product/description text. Raises rather than returning an empty dict silently."""
    if kind not in ANNOTATION_KINDS:
        raise AnnotationError("unknown annotation kind %r (known: %s)" % (kind, ANNOTATION_KINDS))
    it = {"ncbi_feature_table": iter_feature_table,
          "ncbi_gene_info": _iter_gene_info,
          "sgd_features": _iter_sgd_features}[kind]
    out: dict[str, str] = {}
    for key, text in it(path):
        out.setdefault(key, text)
    if not out:
        raise AnnotationError(
            "annotation %s (%s) yielded ZERO entries -- an empty annotation joins nothing and would "
            "produce a clean-looking 0.0 coverage on no data" % (path, kind))
    return out


# ---------------------------------------------------------------------------------------------------
# Identifier normalisation -- EXPLICIT per key space. Never a generic fuzzy fallback: a fuzzy match is
# how a strain-unique id silently joins to the wrong gene.
# ---------------------------------------------------------------------------------------------------
def normalise_key(key: str, key_space: str) -> str:
    k = (key or "").strip()
    if not k:
        return ""
    if key_space == "gene_symbol":
        return k                                   # case IS meaningful for gene symbols
    if key_space == "systematic_orf":
        return k.upper()                           # YAL001C vs yal001c
    if key_space in ("sausa300_locus_tag", "pa14_locus_tag"):
        # SAUSA300_0001 / SAUSA300_1 -> zero-pad the numeric tail to 4, uppercase the prefix
        up = k.upper().replace("-", "_")
        if "_" in up:
            pre, _, tail = up.rpartition("_")
            if tail.isdigit():
                return "%s_%s" % (pre, tail.zfill(4))
        return up
    if key_space == "entrez_gene_id":
        return k.lstrip("0") or k                  # numeric id; leading zeros are not meaningful
    raise AnnotationError("no normalisation rule for key space %r -- add one explicitly rather than "
                          "falling back to a fuzzy match" % key_space)


@dataclass
class JoinResult:
    n_labels: int
    n_joined: int
    join_rate: float
    n_unjoined: int
    unjoined_sample: list = field(default_factory=list)
    rows: list = field(default_factory=list)              # (key, text, label) for joined labels
    composition: dict = field(default_factory=dict)
    wall: str | None = None

    def as_dict(self) -> dict:
        """The record a rung carries. NOTE: no coverage key here -- coverage is computed by the caller
        ONLY when wall is None, so a walled rung cannot carry a coverage number at all."""
        d = {"n_labels": self.n_labels, "n_joined": self.n_joined,
             "join_rate": round(self.join_rate, 4), "n_unjoined": self.n_unjoined,
             "unjoined_sample": self.unjoined_sample[:8],
             "joined_vs_unjoined_composition": self.composition,
             "min_join_rate_asserted": MIN_JOIN_RATE}
        if self.wall:
            d["wall"] = self.wall
        return d


def _composition(joined_labels, unjoined_labels) -> dict:
    """Label-class makeup of the joined vs unjoined sets.

    THE BIAS THE RATE CANNOT SEE. If the joined subset is enriched for essentials relative to what did
    not join, the rung's coverage is measured on a non-representative slice -- and the join RATE is
    identical whether or not that is true.
    """
    def frac(xs):
        return (sum(1 for x in xs if x) / len(xs)) if xs else None
    je, ue = frac(joined_labels), frac(unjoined_labels)
    out = {"joined_n": len(joined_labels), "unjoined_n": len(unjoined_labels),
           "joined_essential_fraction": None if je is None else round(je, 4),
           "unjoined_essential_fraction": None if ue is None else round(ue, 4)}
    out["essential_enrichment_in_joined"] = (None if (je is None or ue is None)
                                             else round(je - ue, 4))
    return out


def join_labels(annotation: dict[str, str], essential_ids, nonessential_ids, *,
                key_space: str, min_join_rate: float = MIN_JOIN_RATE) -> JoinResult:
    """Join a label set to an annotation; gate on the join rate; always report composition."""
    norm_ann = {}
    for k, v in annotation.items():
        norm_ann.setdefault(normalise_key(k, key_space), v)

    labels = [(i, True) for i in essential_ids] + [(i, False) for i in nonessential_ids]
    if not labels:
        raise AnnotationError("join_labels called with an empty label set")

    rows, joined_lbl, unjoined_lbl, unjoined_keys = [], [], [], []
    for raw, is_ess in labels:
        nk = normalise_key(raw, key_space)
        text = norm_ann.get(nk)
        if text is None:
            unjoined_lbl.append(is_ess)
            unjoined_keys.append(raw)
        else:
            rows.append((raw, text, is_ess))
            joined_lbl.append(is_ess)

    rate = len(rows) / len(labels)
    res = JoinResult(n_labels=len(labels), n_joined=len(rows), join_rate=rate,
                     n_unjoined=len(unjoined_lbl), unjoined_sample=unjoined_keys[:8],
                     rows=rows, composition=_composition(joined_lbl, unjoined_lbl))
    # INCLUSIVE at the floor, matching the documented bar.
    if rate < min_join_rate:
        res.wall = WALL_JOIN_RATE_BELOW_FLOOR
        res.rows = []        # a walled rung must not hand scorable rows downstream
    return res
