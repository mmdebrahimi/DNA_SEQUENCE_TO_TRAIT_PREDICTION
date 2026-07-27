"""GEUVADIS data loaders + the population-structure de-confounding split.

Light path only: the quantified RPKM matrix (462 x 23722) + the E-GEUV-1 sample
annotation (ancestry category = population). Genotypes are sliced per-gene cis-window
elsewhere. Pure stdlib + numpy; no network.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ExprMatrix:
    genes: list[str]          # TargetID (ENSG...)
    chrom: list[str]
    coord: list[int]          # gene coordinate (GRCh37)
    samples: list[str]        # HG.../NA... ids, column order
    # values[g][s] RPKM; kept as list-of-lists to avoid a hard numpy import here
    values: list[list[float]]


def load_expr(path: str | Path) -> ExprMatrix:
    op = gzip.open if str(path).endswith(".gz") else open
    genes, chrom, coord, values = [], [], [], []
    with op(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[4:]                       # cols 0-3 = TargetID,Gene_Symbol,Chr,Coord
        for line in fh:
            p = line.rstrip("\n").split("\t")
            genes.append(p[0]); chrom.append(p[2]); coord.append(int(p[3]))
            values.append([float(x) for x in p[4:]])
    return ExprMatrix(genes, chrom, coord, samples, values)


def parse_sample_population(sdrf_path: str | Path) -> dict[str, str]:
    """{sample_id -> population} from the E-GEUV-1 SDRF ancestry-category column."""
    rows = Path(sdrf_path).read_text(encoding="utf-8", errors="replace").splitlines()
    header = rows[0].split("\t")

    def col(sub: str) -> int:
        for i, h in enumerate(header):
            if sub.lower() in h.lower():
                return i
        raise KeyError(sub)

    i_ind, i_anc = col("Characteristics[individual]"), col("Characteristics[ancestry category]")
    out: dict[str, str] = {}
    for r in rows[1:]:
        p = r.split("\t")
        if len(p) > max(i_ind, i_anc):
            out[p[i_ind].strip()] = p[i_anc].strip()
    return out


# ancestry-category label -> 1000G population code (SDRF uses short labels)
_POP_CANON = {
    "utah": "CEU", "finnish": "FIN", "british": "GBR", "tuscan": "TSI", "yoruba": "YRI",
    # long-form fallbacks (other SDRF vintages)
    "utah residents (ceph) with northern and western european ancestry": "CEU",
    "finnish in finland": "FIN", "british in england and scotland": "GBR",
    "toscani in italia": "TSI", "yoruba in ibadan, nigeria": "YRI",
}


def canon_pop(label: str) -> str:
    return _POP_CANON.get(label.strip().lower(), label.strip())
