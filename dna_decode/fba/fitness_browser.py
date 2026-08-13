"""Conditional essentiality from the Fitness Browser RB-TnSeq compendium — 25 carbon sources, not 4 media.

`conditional_essentiality.py` scores the medium-dependent essentiality switch against the Orth 2011 screen:
**1,075 genes x 4 media, 68 conditionally essential, 268 gene x condition cells**. That was the best
substrate available offline. It is no longer the best available.

The Fitness Browser (Price & Deutschbauer, figshare 25236931, CC BY 4.0 — 48 organisms / 7,552 experiments
/ 27.4M gene-fitness measurements) carries the same organism (*E. coli* BW25113, orgId `Keio`) across
**168 experiments, 68 of them carbon-source**. Joining is direct: `Gene.sysName` is the b-number, verified
2026-08-12 at 1,075/1,075 against the Orth gold standard and 1,515/1,516 against iML1515.

**What this is and is not.** RB-TnSeq measures a CONTINUOUS fitness value per gene per experiment (plus a
t-statistic), so the GENE axis is genuinely two-sided. It does NOT make the growth/no-growth carbon
benchmark two-sided — these assays only ran on sources the organism grows on. Two different axes; only
the first is fixed here.

**The essentiality threshold is inherited, not invented:** `fit < -2` is the same cutoff the shipped FBA
cell's Keio validation already uses (Bernstein 2023 method, recorded in the `fba` evidence contract as
"fitness<-2 = essential-on-glucose"). Reusing it keeps this comparable to the existing number instead of
introducing a second, incompatible convention.
"""
from __future__ import annotations

import sqlite3
import statistics
from collections import defaultdict
from pathlib import Path

from dna_decode.fba.carbon_growth import build_exchange_name_index, match_carbon_exchange
from dna_decode.fba.conditional_essentiality import GeneRecord

DEFAULT_DB = Path("D:/dna_decode_cache/fitness_browser/feba.db")

# Inherited from the shipped fba cell's Keio validation (Bernstein 2023). NOT a fresh invention.
ESSENTIAL_FITNESS = -2.0
ORG_ID = "Keio"


def open_db(path: str | Path | None = None) -> sqlite3.Connection:
    """Read-only connection. Raises FileNotFoundError with the provenance if the 7 GB db is absent."""
    p = Path(path) if path else DEFAULT_DB
    if not p.exists():
        raise FileNotFoundError(
            f"Fitness Browser db not found at {p}. It is NOT committed (7.4 GB). Fetch feba.db.gz from "
            "figshare 10.6084/m9.figshare.25236931 (CC BY 4.0) and gunzip it. See "
            "wiki/label_acquisition_candidates_2026-08-12.md")
    return sqlite3.connect(f"file:{p.as_posix()}?mode=ro", uri=True)


def carbon_conditions(conn: sqlite3.Connection, model) -> dict[str, str]:
    """{carbon-source label -> iML1515 exchange id} for every carbon source this model can represent.

    An unmappable label (racemic mixture, casamino acids) is DROPPED rather than guessed — see
    `carbon_growth.match_carbon_exchange`. Measured 2026-08-12: 25 of 28 map.
    """
    idx = build_exchange_name_index(model)
    out: dict[str, str] = {}
    for (cond,) in conn.execute(
            "SELECT DISTINCT condition_1 FROM Experiment WHERE orgId=? AND expGroup='carbon source'",
            (ORG_ID,)):
        ex = match_carbon_exchange(cond, idx)
        if ex:
            out[cond] = ex
    return out


def load_records(conn: sqlite3.Connection, conditions: dict[str, str],
                 gene_filter: set[str] | None = None,
                 threshold: float = ESSENTIAL_FITNESS,
                 min_abs_t: float | None = None,
                 min_abs_t_mode: str = "per_cell") -> list[GeneRecord]:
    """Experimental conditional essentiality per (gene, carbon source), as GeneRecords.

    Replicate experiments for the same carbon source are averaged (62 experiments over 25 sources, so
    ~2.5 replicates each). A gene is kept only if it has a value in EVERY condition — a partial row would
    make the switch pattern incomparable across genes.

    `paper_fba` is left empty: unlike the Orth substrate there is no published FBA column here, so this
    substrate has NO built-in reproduction gate. That is a real difference and is recorded in the artifact.

    **The t-statistic.** `GeneFitness` carries a per-measurement `t` beside `fit` (verified on the live
    db: columns are `['orgId','locusId','expName','fit','t']`). It was never selected by this loader — a
    published memo claimed it was read-but-unused, which was false in the "read" half.

    `min_abs_t` applies a confidence bar, and **`min_abs_t_mode` is load-bearing — the wrong mode
    destroys the very phenomenon being measured** (measured 2026-08-13, see
    `wiki/fba_label_threshold_sweep_2026-08-13.md`):

    - **`per_cell` (DEFAULT, correct):** a cell is called essential only if the fitness clears the bar
      AND that measurement is confident. A low-|t| cell simply is not claimed as essential; the gene
      stays. This is the only mode that answers "are the essential CALLS trustworthy?".
    - **`all_conditions` (the v1 shape, kept for the record):** drop the gene unless `abs(mean t) >=
      min_abs_t` in EVERY condition. **This is ANTI-SELECTIVE for conditional essentiality** — a
      switcher is confidently essential in one condition and confidently NEUTRAL (t ≈ 0) in the others,
      so it fails the bar by construction. Measured: 15 of 15 grid settings collapsed to 100% constant
      predictions and ZERO commitments, because the filter had removed every switcher.

    `min_abs_t=None` is the default and changes nothing — the inherited `fit < -2` cutoff still defines
    every shipped number.
    """
    agg: dict[tuple[str, str], list[float]] = defaultdict(list)
    agg_t: dict[tuple[str, str], list[float]] = defaultdict(list)
    exp_to_cond = {}
    for name, cond in conn.execute(
            "SELECT expName, condition_1 FROM Experiment WHERE orgId=? AND expGroup='carbon source'",
            (ORG_ID,)):
        if cond in conditions:
            exp_to_cond[name] = cond

    q = ("SELECT g.sysName, f.expName, f.fit, f.t FROM GeneFitness f "
         "JOIN Gene g ON g.orgId=f.orgId AND g.locusId=f.locusId WHERE f.orgId=?")
    for sysname, expname, fit, tstat in conn.execute(q, (ORG_ID,)):
        cond = exp_to_cond.get(expname)
        if cond is None or not sysname:
            continue
        if gene_filter is not None and sysname not in gene_filter:
            continue
        agg[(sysname, cond)].append(fit)
        if tstat is not None:
            agg_t[(sysname, cond)].append(tstat)

    by_gene: dict[str, dict[str, float]] = defaultdict(dict)
    for (gene, cond), vals in agg.items():
        by_gene[gene][cond] = statistics.mean(vals)
    t_by_gene: dict[str, dict[str, float]] = defaultdict(dict)
    for (gene, cond), vals in agg_t.items():
        t_by_gene[gene][cond] = statistics.mean(vals)

    keys = sorted(conditions)
    records: list[GeneRecord] = []
    for gene, per_cond in by_gene.items():
        if len(per_cond) != len(keys):
            continue                                  # incomplete row -> not comparable
        per_t = t_by_gene.get(gene, {})
        if min_abs_t is not None and min_abs_t_mode == "all_conditions":
            # ANTI-SELECTIVE for switchers -- kept only so the sweep can demonstrate it. See docstring.
            if len(per_t) != len(keys) or any(abs(per_t[c]) < min_abs_t for c in keys):
                continue
        if min_abs_t is not None and min_abs_t_mode == "per_cell":
            # Fail-closed per CELL: an unconfident measurement cannot support an ESSENTIAL claim, but
            # it is no reason to discard the gene (its other conditions may be perfectly measured).
            exp = {c: (per_cond[c] < threshold and abs(per_t.get(c, 0.0)) >= min_abs_t) for c in keys}
        else:
            exp = {c: (per_cond[c] < threshold) for c in keys}
        records.append(GeneRecord(gene_id=gene, gene=gene, experimental=exp, paper_fba={},
                                  conditionally_essential=any(exp.values()) and not all(exp.values())))
    return records


def mean_fitness_matrix(conn: sqlite3.Connection, conditions: dict[str, str],
                        genes: set[str]) -> dict[str, dict[str, float]]:
    """{condition: {gene: mean fitness}} — the CONTINUOUS values, for threshold-free scoring."""
    agg: dict[tuple[str, str], list[float]] = defaultdict(list)
    exp_to_cond = {}
    for name, cond in conn.execute(
            "SELECT expName, condition_1 FROM Experiment WHERE orgId=? AND expGroup='carbon source'",
            (ORG_ID,)):
        if cond in conditions:
            exp_to_cond[name] = cond
    q = ("SELECT g.sysName, f.expName, f.fit, f.t FROM GeneFitness f "
         "JOIN Gene g ON g.orgId=f.orgId AND g.locusId=f.locusId WHERE f.orgId=?")
    for sysname, expname, fit, _t in conn.execute(q, (ORG_ID,)):
        cond = exp_to_cond.get(expname)
        if cond is not None and sysname in genes:
            agg[(sysname, cond)].append(fit)
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for (gene, cond), vals in agg.items():
        out[cond][gene] = statistics.mean(vals)
    return dict(out)


def mean_t_matrix(conn: sqlite3.Connection, conditions: dict[str, str],
                  genes: set[str]) -> dict[str, dict[str, float]]:
    """{condition: {gene: mean t-statistic}} — mirrors `mean_fitness_matrix` on the confidence axis.

    A sweep needs the t DISTRIBUTION before it can pick a bar; this makes it inspectable without
    re-querying the 7 GB table through `load_records`.
    """
    agg: dict[tuple[str, str], list[float]] = defaultdict(list)
    exp_to_cond = {}
    for name, cond in conn.execute(
            "SELECT expName, condition_1 FROM Experiment WHERE orgId=? AND expGroup='carbon source'",
            (ORG_ID,)):
        if cond in conditions:
            exp_to_cond[name] = cond
    q = ("SELECT g.sysName, f.expName, f.t FROM GeneFitness f "
         "JOIN Gene g ON g.orgId=f.orgId AND g.locusId=f.locusId WHERE f.orgId=?")
    for sysname, expname, tstat in conn.execute(q, (ORG_ID,)):
        cond = exp_to_cond.get(expname)
        if cond is not None and sysname in genes and tstat is not None:
            agg[(sysname, cond)].append(tstat)
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for (gene, cond), vals in agg.items():
        out[cond][gene] = statistics.mean(vals)
    return dict(out)


def apply_carbon_condition(model, exchange: str, uptake: float = 10.0,
                           all_carbon: tuple[str, ...] = ()) -> None:
    """Set `exchange` as the sole carbon source, IN PLACE.

    Every other candidate carbon exchange is closed first, so a residual glucose uptake can never leak in
    and make every condition silently score as glucose (the failure mode `_ALL_CARBON` guards in the
    4-media path).
    """
    have = {r.id for r in model.exchanges}
    if exchange not in have:
        raise KeyError(f"model {model.id} has no exchange {exchange!r}")
    medium = dict(model.medium)
    for ex in set(all_carbon) | {"EX_glc__D_e"}:
        medium.pop(ex, None)
    medium[exchange] = uptake
    model.medium = medium
