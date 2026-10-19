"""Re-derive every number REPORTED for the tile-rung selection/hybrid plan, from raw per-seed data.

    uv run python scripts/verify_tile_arm_claims.py     # exit 0 = every claim reproduces

Run on user request (2026-10-14) to double-check the reported findings. 52 checks; each states the
value CLAIMED in the session report and the value recomputed here from `per_seed` arrays, the
prediction tables on D:, and the committed artifacts.

**A PRECISION HAZARD THIS SCRIPT FOUND IN ITSELF, now fixed rather than documented as a wart.** The
hybrid artifact stores `per_seed` rounded to 4 decimals, so a paired median recomputed from it can
disagree with the artifact's own full-precision value in the 4th decimal. It fired on exactly one row:
`conv` vs `max_tx` gave **-0.0766** from the rounded array against the reported **-0.0767**. The
full-precision median is **-0.076666** -- 0.00002 from the 0.07665 boundary -- so the REPORT was right
and the first version of this audit was the less precise instrument.

C6 therefore recomputes that row from the PREDICTION VECTORS, and **C6b keeps the rounded computation
visible beside it** so the hazard is recorded rather than buried: a future reader can see that a
paired median taken off `per_seed` is good to ~1e-4, not better. The other two paired rows (C2, C8)
matched under rounding only because they do not sit near a boundary -- that was luck, not rigour.

**What this does NOT check.** It verifies that the reported numbers are faithful to the artifacts and
that they reproduce. It does NOT re-validate the science: whether the asserted 0.50 headroom bar is the
right bar, whether peak tiles are a representative pool, or whether the scripts that produced the
artifacts are themselves correct. The checks that DO bear on correctness are the external ones -- the
four published medians reproducing (A2), the designed-grid anchor (B13/B14), and the k=1 degeneracy
(B6) -- because those were fixed by someone else before this plan existed.

Deliberately recomputes medians, paired margins and means from `per_seed` / `arms[...]["per_seed"]`
arrays rather than reading `median` / `_paired_*` keys, because trusting a summary field would make
this check a tautology over the same arithmetic it is supposed to audit.
"""
from __future__ import annotations

import gzip
import json
import statistics as st
import sys
from pathlib import Path

REPO = Path(r"C:\Users\Farshad\PythonProjects\dna_decode")
PRED = Path("D:/dna_decode_cache/promcalc/tile_predictions")

rows: list[tuple[str, str, str, bool]] = []


def chk(name: str, claimed, got, ok: bool) -> None:
    rows.append((name, str(claimed), str(got), ok))


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def median_of(xs):
    return st.median(xs)


# ---------------------------------------------------------------- A. prediction table
tabs = []
for p in sorted(PRED.glob("seed-*.json.gz")):
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        tabs.append(json.load(fh))

chk("A1 seed tables", 10, len(tabs), len(tabs) == 10)

PUB = {"gc": 0.3112, "conv": 0.3234, "max_tx": 0.3968, "decomposition": 0.3925}
for sc, pub in PUB.items():
    vals = [t["spearman_recomputed"][sc] for t in tabs]
    m = median_of(vals)
    chk(f"A2 {sc} median == published {pub}", pub, round(m, 5), close(m, pub, 3e-5))

chk("A3 all four within 3e-5", "yes",
    "yes" if all(close(median_of([t["spearman_recomputed"][s] for t in tabs]), v, 3e-5)
                 for s, v in PUB.items()) else "no",
    all(close(median_of([t["spearman_recomputed"][s] for t in tabs]), v, 3e-5) for s, v in PUB.items()))

# every table must carry all four scorers (the merge claim)
allfour = all(set(t["predictions"]) == set(PUB) for t in tabs)
chk("A4 every seed has all 4 scorers", True, allfour, allfour)
nosmoke = all(not t.get("smoke") for t in tabs)
chk("A5 no smoke table", True, nosmoke, nosmoke)

man = json.loads((PRED / "manifest.json").read_text(encoding="utf-8"))
chk("A6 n_tiles", 44106, man["n_tiles"], man["n_tiles"] == 44106)
chk("A7 manifest scorers (post-merge)", "4", str(len(man["scorers"])), len(man["scorers"]) == 4)

# ---------------------------------------------------------------- B. selection
sel = json.loads((REPO / "wiki/glm_tile_selection_2026-10-14.json").read_text(encoding="utf-8"))
chk("B1 verdict", "SELECTION_IS_THE_BOTTLENECK", sel["verdict"],
    sel["verdict"] == "SELECTION_IS_THE_BOTTLENECK")


def cell(seed_rec, sc, k):
    for c in seed_rec["cells"]:
        if c["scorer"] == sc and c["k"] == k:
            return c
    return None


h10 = [cell(s, "max_tx", 10)["headroom_captured"] for s in sel["per_seed"]]
chk("B2 max_tx k=10 median", 0.4807, round(median_of(h10), 4), close(median_of(h10), 0.4807, 5e-5))
chk("B3 max_tx k=10 mean", 0.4847, round(st.fmean(h10), 4), close(st.fmean(h10), 0.4847, 5e-5))
n_pass = sum(1 for x in h10 if x >= 0.50)
chk("B4 seeds >= 0.50 bar", "2 of 10", f"{n_pass} of {len(h10)}", n_pass == 2 and len(h10) == 10)

KTREND = {2: 0.3880, 5: 0.4468, 10: 0.4847, 20: 0.5366}
for k, claimed in KTREND.items():
    vals = [cell(s, "max_tx", k)["headroom_captured"] for s in sel["per_seed"]]
    got = st.fmean(vals)
    chk(f"B5 max_tx k={k} mean", claimed, round(got, 4), close(got, claimed, 5e-5))

k1_none = all(cell(s, sc, 1)["headroom_captured"] is None
              for s in sel["per_seed"] for sc in sel["scorers"])
chk("B6 every k=1 headroom is None", True, k1_none, k1_none)

nb = [s["permutation_null_by_k"]["10"] for s in sel["per_seed"]]
nmean = st.fmean([n["mean"] for n in nb])
nsd = st.fmean([n["stdev"] for n in nb])
chk("B7 null k=10 mean-of-means", -0.0035, round(nmean, 4), close(nmean, -0.0035, 5e-5))
chk("B8 null k=10 mean stdev", 0.0485, round(nsd, 4), close(nsd, 0.0485, 5e-5))

dec10 = [cell(s, "decomposition", 10)["headroom_captured"] for s in sel["per_seed"]]
chk("B9 decomposition k=10 median", 0.4722, round(median_of(dec10), 4),
    close(median_of(dec10), 0.4722, 5e-5))
n_dec_wins = sum(1 for a, b in zip(dec10, h10) if a > b)
chk("B10 decomposition beats max_tx on", "3 of 10", f"{n_dec_wins} of 10", n_dec_wins == 3)

gc10 = st.fmean([cell(s, "gc", 10)["headroom_captured"] for s in sel["per_seed"]])
cv10 = st.fmean([cell(s, "conv", 10)["headroom_captured"] for s in sel["per_seed"]])
chk("B11 gc k=10 mean", 0.2948, round(gc10, 4), close(gc10, 0.2948, 5e-5))
chk("B12 conv k=10 mean", 0.3634, round(cv10, 4), close(cv10, 0.3634, 5e-5))

# ---------------------------------------------------------------- B13. the designed-grid anchor
grid = json.loads((REPO / "wiki/glm_oracle_selection_2026-10-07.json").read_text(encoding="utf-8"))


def gcell(ax, sc, k):
    for c in ax["cells"]:
        if c["scorer"] == sc and c["k"] == k:
            return c
    return None


oh = [gcell(a, "onehot", 10)["headroom_captured"] for a in grid["per_axis"]]
gg = [gcell(a, "gc", 10)["headroom_captured"] for a in grid["per_axis"]]
chk("B13 grid onehot cross-AXIS mean", 0.7937, round(st.fmean(oh), 4), close(st.fmean(oh), 0.7937, 5e-5))
chk("B14 grid gc cross-AXIS mean", 0.1328, round(st.fmean(gg), 4), close(st.fmean(gg), 0.1328, 5e-5))
chk("B15 grid verdict", "ORACLE_READY_AS_SCORING_HEAD", grid["verdict"],
    grid["verdict"] == "ORACLE_READY_AS_SCORING_HEAD")
chk("B16 grid axes", 5, len(grid["per_axis"]), len(grid["per_axis"]) == 5)

# ---------------------------------------------------------------- C. hybrid
hy = json.loads((REPO / "wiki/glm_tile_hybrid_arm_2026-10-14.json").read_text(encoding="utf-8"))
A = hy["arms"]
chk("C1 verdict", "HYBRID_ADDS", hy["verdict"], hy["verdict"] == "HYBRID_ADDS")

hyb = A["hybrid_conv_maxtx"]["per_seed"]
mtx = A["max_tx"]["per_seed"]
cnv = A["conv"]["per_seed"]
ctl = A["control_shuffled_conv"]["per_seed"]

# RECOMPUTE the paired margin from per-seed values rather than reading _paired_*
d_h = [a - b for a, b in zip(hyb, mtx)]
chk("C2 paired hybrid-vs-max_tx median", 0.0339, round(median_of(d_h), 4),
    close(median_of(d_h), 0.0339, 5e-5))
chk("C3 hybrid wins", "10 of 10", f"{sum(1 for d in d_h if d > 0)} of {len(d_h)}",
    sum(1 for d in d_h if d > 0) == 10)
chk("C4 hybrid median spearman", 0.4254, round(median_of(hyb), 4), close(median_of(hyb), 0.4254, 5e-5))
chk("C5 max_tx median (hybrid artifact)", 0.3968, round(median_of(mtx), 4),
    close(median_of(mtx), 0.3968, 5e-5))

# C6 recomputes from the FULL-PRECISION prediction vectors, not from the 4dp `per_seed` array.
# Why: this row sits 0.00002 from the 0.07665 rounding boundary, so the rounded inputs give -0.0766
# against the artifact's (correct) -0.0767. Fixing the instrument beats documenting its wart -- but the
# rounded value is computed too and printed beside it, so the hazard stays visible rather than hidden.
from dna_decode.glm.genomewide import spearman  # noqa: E402  (needs the venv; only this row does)

_conv_full, _mtx_full = [], []
for _t in tabs:
    _m = _t["measured"]
    _conv_full.append(spearman(_t["predictions"]["conv"], _m))
    _mtx_full.append(spearman(_t["predictions"]["max_tx"], _m))
d_c = [a - b for a, b in zip(_conv_full, _mtx_full)]
d_c_rounded = [a - b for a, b in zip(cnv, mtx)]
chk("C6 conv-alone vs max_tx median (full precision)", -0.0767, round(median_of(d_c), 4),
    close(median_of(d_c), -0.0767, 5e-5))
chk("C6b same from 4dp per_seed (shows the hazard)", "-0.0766 (NOT the claim)",
    round(median_of(d_c_rounded), 4), close(median_of(d_c_rounded), -0.0766, 5e-5))
chk("C7 conv-alone wins", "0 of 10", f"{sum(1 for d in d_c if d > 0)} of {len(d_c)}",
    sum(1 for d in d_c if d > 0) == 0)

d_ctl = [a - b for a, b in zip(ctl, mtx)]
chk("C8 control vs max_tx median", -0.1232, round(median_of(d_ctl), 4),
    close(median_of(d_ctl), -0.1232, 5e-5))
chk("C9 control wins", "0 of 10", f"{sum(1 for d in d_ctl if d > 0)} of {len(d_ctl)}",
    sum(1 for d in d_ctl if d > 0) == 0)

ci = A["_paired_hybrid_vs_comparator"]["ci95"]
chk("C10 CI95 excludes zero", "[0.0238, 0.0409]", str([round(x, 4) for x in ci]),
    ci[0] > 0 and close(ci[0], 0.0238, 5e-5) and close(ci[1], 0.0409, 5e-5))

o = hy["orthogonality"]
chk("C11 orthogonality median", 0.3862, o["median"], close(o["median"], 0.3862, 5e-5))
chk("C12 orthogonality min", 0.359, o["min"], close(o["min"], 0.359, 5e-4))
chk("C13 orthogonality max", 0.479, o["max"], close(o["max"], 0.479, 5e-4))
chk("C14 orthogonality recomputed median", o["median"], round(median_of(o["per_seed"]), 4),
    close(median_of(o["per_seed"]), o["median"], 5e-5))

f_mtx, f_hyb = median_of(mtx) / 0.941, median_of(hyb) / 0.941
chk("C15 ceiling fraction max_tx", 0.422, round(f_mtx, 3), close(f_mtx, 0.422, 5e-4))
chk("C16 ceiling fraction hybrid", 0.452, round(f_hyb, 3), close(f_hyb, 0.452, 5e-4))

d_h3 = [a - b for a, b in zip(A["hybrid_conv_maxtx_gc"]["per_seed"], mtx)]
chk("C17 hybrid+gc paired median", 0.0302, round(median_of(d_h3), 4),
    close(median_of(d_h3), 0.0302, 5e-5))
chk("C18 3-way is WORSE than 2-way", True, median_of(d_h3) < median_of(d_h),
    median_of(d_h3) < median_of(d_h))

# ---------------------------------------------------------------- D. ledger
CRIT = ("ceiling: best achievable prediction of bacterial promoter activity from sequence, "
        "reported against the replicate-derived 0.941 ceiling")
led = json.loads((REPO / "soraya_runs/arms-ledger.json").read_text(encoding="utf-8"))
arms = led["criteria"][CRIT]
meas = {n: a["metric"] for n, a in arms.items() if a["status"] == "measured"}
plan = [n for n, a in arms.items() if a["status"] == "planned"]
chk("D1 total arms", 10, len(arms), len(arms) == 10)
chk("D2 measured / planned", "6 / 4", f"{len(meas)} / {len(plan)}", len(meas) == 6 and len(plan) == 4)
HY = "promcalc-plus-conv-rank-average-hybrid-150bp"
chk("D3 hybrid metric on disk", 0.4254, arms[HY]["metric"], close(arms[HY]["metric"], 0.4254, 1e-9))
chk("D4 hybrid status", "measured", arms[HY]["status"], arms[HY]["status"] == "measured")
best = max(meas, key=lambda k: meas[k])
chk("D5 hybrid is best measured", HY, best, best == HY)

print(f"{'check':46s} {'claimed':>30s} {'recomputed':>30s}  verdict")
print("-" * 118)
bad = 0
for n, c, g, ok in rows:
    if not ok:
        bad += 1
    print(f"{n:46s} {c:>30s} {g:>30s}  {'ok' if ok else 'MISMATCH'}")
print("-" * 118)
print(f"{len(rows)} checks, {bad} MISMATCH")
sys.exit(1 if bad else 0)
