"""EP-4 v0.1 ExPEC-recall hardening gate (mission ep4-v01-expec-recall).

CROSS-AXIS rule of record (user round-2 commitment 2026-06-03): an ExPEC rescue requires BOTH
extraintestinal axes — >=1 iron-acquisition gene AND >=1 capsule/serum gene (each >=0.80). This
structurally excludes capsule-only genomes (the lone-traT JSPG) that the earlier flat-K=1 burden
over-rescued. The committed decision: "a clean 0.833 beats an overfit 0.85" — so cohort ExPEC recall
is CAPPED at 0.833 (10/12) = `blocked:strategy-cap`, NOT pushed to the in-sample 0.917, and the two
un-rescued strains are exactly JSPG (capsule-only, structurally excluded) and JSLG (empty / no
support, unrescuable without weakening the strong-marker rule).

Endpoints asserted on the 24-genome H4 cohort:
  1. ExPEC recall == 10/12 (the cross-axis strategy-cap; JSPG + JSLG correctly NOT rescued)
  2. confident-supported-call precision == 1.0  (the precision invariant — unchanged)
  3. EPEC recall == 1.0 with the DEC-module (LEE) gate above the ExPEC support branch
  4. rescued ExPEC calls are LOW_CONFIDENCE only
  5. lowest non-ExPEC support-axis margin reported (overfit canary): JSPG is single-axis

Runs fully OFFLINE from two committed artifacts (no ENA fetch, no genome read, no sklearn):
  - data/pathotype_cov_cache/<gid>_v1.json      cluster -> [best_gene, coverage]
  - data/pathotype_pergene_cache/<gid>_v1.json  support gene_prefix -> coverage

Labels: Horesh-2021 F1 roster (ExPEC=Salipante isolation-site, EPEC=Hazen DECA-curated), selected
identically to scripts/build_pergene_support_cache.py. Runnable via pytest OR standalone.
"""
import csv
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pathotype.resolve import resolve_call
from dna_decode.pathotype.expec_score import (
    support_gene_count, meets_cross_axis_support, iron_capsule_counts,
)
from dna_decode.pathotype.markers import EXPEC_SUPPORT_GENE_PREFIXES

F1 = REPO / "data/external/horesh2021_F1_genome_metadata.csv"
COV = REPO / "data/pathotype_cov_cache"
PERGENE = REPO / "data/pathotype_pergene_cache"
CONFIDENT_COV, PARTIAL_COV = 0.80, 0.65
N_PER_CLASS = 12
WGS_MASTER = re.compile(r"^[A-Z]{4}\d{8}(\.fa)?$", re.I)

EXPEC_OK = {"UPEC_COMPATIBLE", "ExPEC_COMPATIBLE"}
EPEC_OK = {"tEPEC_COMPATIBLE", "aEPEC_COMPATIBLE"}
# Cross-axis strategy-cap (committed): 10/12. NOT a >=0.85 bar — JSPG (capsule-only) is
# structurally excluded and JSLG (empty) is unrescuable; a clean 0.833 beats an overfit 0.917.
EXPEC_RECALL_CAP = 10 / 12


def _roster():
    """(gid, label) for the 24-genome cohort; label 0=ExPEC, 1=EPEC. Mirrors the build script."""
    rows = list(csv.DictReader(open(F1, encoding="utf-8")))
    clean = [r for r in rows if "(predicted)" not in r["Pathotype"]
             and r["Pathotype"] not in ("Not determined", "")]
    ok = lambda r: bool(WGS_MASTER.match(r["Assembly_name"].strip()))
    expec = [r for r in clean if r["Pathotype"].strip().startswith("ExPEC")
             and r["Source"].startswith("Salipante") and ok(r)][:N_PER_CLASS]
    epec = [r for r in clean if r["Pathotype"].strip().startswith("EPEC")
            and r["Source"].startswith("Hazen") and ok(r)][:N_PER_CLASS]
    return [(r["ID"], 0) for r in expec] + [(r["ID"], 1) for r in epec]


def _profile_and_support(gid):
    """Build (cluster_profile, partial_clusters, support_gene_count, cross_axis) from caches."""
    cov = json.loads((COV / f"{gid}_v1.json").read_text())
    profile, partial = {}, set()
    for cluster, (_gene, c) in cov.items():
        if c >= CONFIDENT_COV:
            profile[cluster] = True
        elif c >= PARTIAL_COV:
            partial.add(cluster)
    pf = PERGENE / f"{gid}_v1.json"
    pergene = json.loads(pf.read_text()) if pf.exists() else {}
    sgc = support_gene_count(pergene)
    cross = meets_cross_axis_support(pergene)
    return profile, frozenset(partial), sgc, cross


def _resolve_cohort():
    rows = []
    for gid, y in _roster():
        if not (COV / f"{gid}_v1.json").exists():
            continue
        profile, partial, sgc, cross = _profile_and_support(gid)
        call = resolve_call(profile, partial_clusters=partial,
                            support_gene_count=sgc, cross_axis_support=cross)
        rows.append({"gid": gid, "y": y, "call": call["primary"], "tier": call["confidence_tier"]})
    return rows


def test_expec_recall_equals_cross_axis_cap():
    """Cross-axis caps cohort ExPEC recall at 10/12; JSPG + JSLG are the two correctly-unrescued."""
    rows = _resolve_cohort()
    expec = [r for r in rows if r["y"] == 0]
    assert len(expec) == 12, f"expected 12 ExPEC genomes, got {len(expec)}"
    recall = sum(1 for r in expec if r["call"] in EXPEC_OK) / len(expec)
    missed = sorted(r["gid"][:4] for r in expec if r["call"] not in EXPEC_OK)
    assert recall == EXPEC_RECALL_CAP, f"ExPEC recall {recall:.3f} != cap {EXPEC_RECALL_CAP:.3f} (missed {missed})"
    assert missed == ["JSLG", "JSPG"], f"unexpected missed set {missed} (want JSLG empty + JSPG capsule-only)"


def test_jspg_capsule_only_structurally_excluded():
    """JSPG carries capsule-only support (lone traT, 0 iron) -> cross-axis must NOT rescue it."""
    pf = PERGENE / "JSPG00000000_v1.json"
    pergene = json.loads(pf.read_text())
    iron, caps = iron_capsule_counts(pergene)
    assert iron == 0 and caps >= 1, f"JSPG expected iron=0,caps>=1; got iron={iron},caps={caps}"
    assert not meets_cross_axis_support(pergene), "JSPG (capsule-only) must fail cross-axis"


def test_confident_supported_precision_is_1():
    """Among CONFIDENT supported calls, every one must match its label (precision invariant)."""
    rows = _resolve_cohort()
    conf = [r for r in rows if r["tier"] == "CONFIDENT" and r["call"] in (EXPEC_OK | EPEC_OK)]
    assert conf, "expected >=1 confident supported call to test precision against"
    correct = sum(1 for r in conf
                  if (r["y"] == 0 and r["call"] in EXPEC_OK) or (r["y"] == 1 and r["call"] in EPEC_OK))
    precision = correct / len(conf)
    assert precision == 1.0, f"confident-supported precision regressed to {precision:.3f}"


def test_rescued_calls_are_low_confidence_not_confident():
    """The per-gene rescue must NEVER produce a CONFIDENT ExPEC call (structural precision guard)."""
    rows = _resolve_cohort()
    for r in rows:
        if r["call"] == "ExPEC_COMPATIBLE":
            assert r["tier"] == "LOW_CONFIDENCE", f"{r['gid']} ExPEC call tier={r['tier']} (must be LOW)"


def test_epec_recall_preserved_and_lee_gate_above_expec():
    """EPEC strains stay tEPEC/aEPEC (LEE-gated above the ExPEC support branch). Q3 regression:
    an LEE + FULL cross-axis support synthetic profile resolves to EPEC, never ExPEC."""
    rows = _resolve_cohort()
    epec = [r for r in rows if r["y"] == 1]
    epec_recall = sum(1 for r in epec if r["call"] in EPEC_OK) / len(epec)
    assert epec_recall == 1.0, f"EPEC recall regressed to {epec_recall:.3f}"
    # Q3 synthetic: LEE + full cross-axis support (both axes) must still gate to EPEC, not ExPEC.
    r = resolve_call({"LEE": True}, support_gene_count=len(EXPEC_SUPPORT_GENE_PREFIXES),
                     cross_axis_support=True)
    assert r["primary"] in ("tEPEC_COMPATIBLE", "aEPEC_COMPATIBLE")
    assert r["primary"] != "ExPEC_COMPATIBLE"


def test_single_strong_no_support_still_ambiguous():
    """JSLG-shape (1 strong adhesin, no cross-axis support) stays AMBIGUOUS — not rescued."""
    r = resolve_call({"P_FIMBRIAE": True}, support_gene_count=0, cross_axis_support=False)
    assert r["primary"] == "AMBIGUOUS"


def test_lowest_nonexpec_support_margin_reported():
    """Overfit canary (committed): report the closest any UN-rescued ExPEC came to cross-axis.
    JSPG is single-axis (iron=0) — the margin is structural (a whole axis missing), not a tuned
    threshold, confirming the cap is principled not knife-edge."""
    margins = {}
    for gid, y in _roster():
        if y != 0:
            continue
        pf = PERGENE / f"{gid}_v1.json"
        if not pf.exists():
            continue
        iron, caps = iron_capsule_counts(json.loads(pf.read_text()))
        if not (iron >= 1 and caps >= 1):  # un-rescued
            margins[gid[:4]] = (iron, caps)
    # JSPG present as single-axis; JSLG present as zero-axis. No genome sits at a 1-coverage-point edge.
    assert "JSPG" in margins and margins["JSPG"][0] == 0, f"JSPG margin {margins.get('JSPG')}"
    assert "JSLG" in margins and margins["JSLG"] == (0, 0), f"JSLG margin {margins.get('JSLG')}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    raise SystemExit(1 if failed else 0)
