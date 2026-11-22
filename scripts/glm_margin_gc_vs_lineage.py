"""Is the winnable margin driven by GC content or by LINEAGE? The new corpus can finally separate them.

    uv run python scripts/glm_margin_gc_vs_lineage.py [--spec PATH] [--limit-species N]

**The limit this closes, which I declared this morning and could not then close.** Measuring the
winnable margin across the four cached species gave a margin-versus-GC correlation of **-0.941**, and it
had to be reported as uninterpretable: the two low-GC genomes were both *Campylobacter* and the three
high-GC ones were all *Enterobacteriaceae*. GC and phylogeny were perfectly collinear, so the effective
sample size was **2 groups, not 5 points**, and reading it as a GC effect would have repeated the
lineage confound that has already closed four arms in this project. The artifact said what would be
needed to separate them: "substrates that cross them -- a high-GC Campylobacter or a low-GC
enterobacterium -- which this set does not contain."

The expanded corpus does contain them. Sixteen species across roughly twelve phyla spanning GC 0.272 to
0.720 means GC and lineage are no longer the same variable, so the question becomes answerable rather
than confounded.

**Collinearity is TESTED before anything is concluded from it, and that order is the point.** The script
first measures how strongly GC and phylum predict each other in whatever corpus it is handed. If they
are still collinear it says so and refuses to report a GC effect -- because a decoupled corpus is the
precondition for the whole measurement, not an assumption to carry. The four-species case would fail
that test, which is the non-vacuity check: the guard must reject the data that motivated it.

**Two readings are computed and kept apart.**
  * The **raw** correlation of margin against GC, which is what the earlier artifact reported.
  * A **within-phylum** reading: the correlation computed after removing each phylum's mean margin and
    mean GC. If GC matters for its own sake, the relationship survives when lineage can no longer
    explain it. If it was lineage all along, the within-phylum correlation collapses toward zero. This
    is the same move the project's own FBA work used -- a pooled ranking cannot answer a within-group
    question.

Honest scope: this measures the INSTRUMENT's margin against substrate composition. It says nothing about
any generator, and a correlation over a handful of phyla is a description rather than a law.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "refseq_cache"
REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
DEFAULT_SPEC = Path("C:/Users/Laptop Portable/pythonprojects/dna-generator-research/handoffs"
                    "/bacterial-expansion-species-spec.json")
OUT = ROOT / "wiki" / "glm_margin_gc_vs_lineage_2026-11-21.json"

# Held identical to the committed kernel so the margin is the same quantity measured before.
LENGTH = 150
LIMIT = 1200
N_GENERATED = 120
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "positional", "both")
SEEDS = tuple(range(10))

MIN_SPECIES_FOR_A_CLAIM = 8
MIN_PHYLA_FOR_A_CLAIM = 5
COLLINEARITY_CEILING = 0.70      # |corr(GC, phylum-mean GC)| above this means still confounded


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def resolve(values: list[float]) -> dict:
    n = len(values)
    if n < 3:
        return {"n": n, "resolvable_above_zero": None}
    mean = statistics.fmean(values)
    sd = statistics.stdev(values)
    se = sd / math.sqrt(n)
    crit = 2.262 if n == 10 else 1.96
    return {"n": n, "mean": round(mean, 4), "sd": round(sd, 4),
            "ci95": [round(mean - crit * se, 4), round(mean + crit * se, 4)],
            "resolvable_above_zero": bool(mean - crit * se > 0)}


def margin_for_genome(ref_dir: Path) -> dict | None:
    """Null, perfect-generator floor and margin per mode, at the committed geometry.

    Identical construction to `glm_margin_across_organisms` -- spares filtered by SEQUENCE because the
    corpus holds duplicate upstream windows and the discriminator only ever sees sequence.
    """
    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves, gc_fraction
    from dna_decode.glm.falsifier import heldout_auroc

    per_mode: dict[str, list[float]] = {m: [] for m in MODES}
    gcs: list[float] = []
    for s in SEEDS:
        wins, _ = build_natural_set(ref_dir, length=LENGTH, limit=None, seed=s, kind=KIND)
        if len(wins) < LIMIT + N_GENERATED:
            return None
        fit, ev = fit_eval_halves(wins[:LIMIT])
        used = set(fit) | set(ev)
        spare = []
        for w in wins[LIMIT:]:
            if len(spare) >= N_GENERATED:
                break
            if w.seq in used:
                continue
            spare.append(w.seq)
            used.add(w.seq)
        if len(spare) < N_GENERATED:
            return None
        gen = MarkovGenerator(order=MARKOV_ORDER, seed=s).fit(fit).generate(N_GENERATED, LENGTH)
        gcs.append(sum(gc_fraction(x) for x in ev) / len(ev))
        for mode in MODES:
            kw = {"k": KMER_K, "seed": s, "mode": mode}
            if mode in ("positional", "both"):
                kw["bins"] = POSITIONAL_BINS
            null = heldout_auroc(ev, gen, **kw)
            floor = heldout_auroc(ev, spare, **kw)
            if null is not None and floor is not None:
                per_mode[mode].append(null - floor)
    return {"gc": round(statistics.fmean(gcs), 4),
            "margins": {m: resolve(v) for m, v in per_mode.items()}}


def collinearity(gcs: list[float], phyla: list[str]) -> dict:
    """How much of GC is explained by phylum. The precondition, tested rather than assumed.

    Implemented as the correlation between each genome's GC and its PHYLUM MEAN GC: if phylum predicts
    GC almost perfectly that correlation approaches 1 and no GC effect can be separated from lineage.
    """
    means: dict[str, float] = {}
    for p in set(phyla):
        vals = [g for g, q in zip(gcs, phyla) if q == p]
        means[p] = statistics.fmean(vals)
    predicted = [means[p] for p in phyla]
    if len({round(v, 6) for v in predicted}) < 2 or len({round(v, 6) for v in gcs}) < 2:
        return {"corr_gc_vs_phylum_mean": None, "decoupled": False,
                "reason": "too little variation to assess"}
    r = statistics.correlation(gcs, predicted)
    return {
        "corr_gc_vs_phylum_mean": round(r, 3),
        "n_phyla": len(means),
        "ceiling": COLLINEARITY_CEILING,
        "decoupled": bool(abs(r) < COLLINEARITY_CEILING),
        "reading": ("a correlation near 1 means phylum predicts GC and the two cannot be separated; "
                    "below the ceiling they vary independently enough to ask the question"),
    }


def within_group(values: list[float], gcs: list[float], groups: list[str]) -> dict:
    """Correlation of margin against GC after removing each group's means.

    A pooled correlation cannot answer a within-group question -- the same discipline the FBA arm had to
    adopt when a pooled AUROC turned out to be reachable with zero switch signal.
    """
    gm_v: dict[str, float] = {}
    gm_g: dict[str, float] = {}
    for p in set(groups):
        vs = [v for v, q in zip(values, groups) if q == p]
        gs = [g for g, q in zip(gcs, groups) if q == p]
        gm_v[p] = statistics.fmean(vs)
        gm_g[p] = statistics.fmean(gs)
    dv = [v - gm_v[p] for v, p in zip(values, groups)]
    dg = [g - gm_g[p] for g, p in zip(gcs, groups)]
    usable = sum(1 for p in set(groups) if sum(1 for q in groups if q == p) > 1)
    if len({round(x, 9) for x in dg}) < 2 or len({round(x, 9) for x in dv}) < 2:
        return {"corr": None, "n_groups_with_replicates": usable,
                "reason": ("after centring there is no variation left -- every group has a single "
                           "member, so a within-group correlation is undefined. This is the expected "
                           "outcome for one-genome-per-phylum sampling and is NOT a null result.")}
    return {"corr": round(statistics.correlation(dg, dv), 3),
            "n_groups_with_replicates": usable,
            "reading": ("if GC matters for its own sake the relationship survives centring; if it was "
                        "lineage all along it collapses toward zero")}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spec", default=str(DEFAULT_SPEC))
    p.add_argument("--limit-species", type=int, default=0,
                   help="measure only the first N species (deterministic), for a quick pass")
    a = p.parse_args(argv)
    require_sklearn()
    t0 = time.time()

    spec_path = Path(a.spec)
    entries = []
    if spec_path.exists():
        rec = json.loads(spec_path.read_text(encoding="utf-8"))
        for g in rec.get("genomes", []):
            acc = g.get("accession")
            if not acc:
                continue
            d = CACHE / acc
            if not ((d / "genome.fna").exists() and (d / "annotations.gff3").exists()):
                continue
            entries.append({"accession": acc, "species": g.get("species"),
                            "phylum": (g.get("lineage_key") or ";").split(";")[0],
                            "lineage_key": g.get("lineage_key"),
                            "gc_spec": g.get("gc_approx"), "dir": d})
    if not entries:
        print(f"REFUSE: no spec genomes are present in {CACHE}. Fetch them first with "
              "scripts/glm_fetch_corpus_genomes.py; this measurement needs the expanded corpus, "
              "because the four-species cache is exactly the confounded case it exists to resolve.")
        return 2
    if a.limit_species:
        entries = entries[:a.limit_species]

    print(f"measuring {len(entries)} species at the committed geometry "
          f"(length={LENGTH} limit={LIMIT} n={N_GENERATED} order={MARKOV_ORDER})")
    rows = []
    for e in entries:
        got = margin_for_genome(e["dir"])
        if got is None:
            print(f"  {e['accession']:20s} {str(e['species'])[:30]:30s} SKIPPED (too few windows)")
            continue
        row = {**{k: v for k, v in e.items() if k != "dir"},
               "gc_measured": got["gc"],
               "margins": got["margins"]}
        rows.append(row)
        print(f"  {e['accession']:20s} {str(e['species'])[:30]:30s} GC {got['gc']:.3f}  "
              + "  ".join(f"{m}={got['margins'][m]['mean']:+.4f}" for m in MODES))

    if len(rows) < 3:
        print("REFUSE: fewer than three species measured; no correlation is meaningful.")
        return 3

    gcs = [r["gc_measured"] for r in rows]
    phyla = [r["phylum"] or "unknown" for r in rows]
    coll = collinearity(gcs, phyla)

    analysis = {}
    for mode in MODES:
        vals = [r["margins"][mode]["mean"] for r in rows]
        raw = (round(statistics.correlation(gcs, vals), 3)
               if len({round(v, 9) for v in vals}) > 1 else None)
        analysis[mode] = {
            "raw_corr_margin_vs_gc": raw,
            "within_phylum": within_group(vals, gcs, phyla),
            "margin_range": [round(min(vals), 4), round(max(vals), 4)],
            "n_resolvable": sum(1 for r in rows if r["margins"][mode]["resolvable_above_zero"]),
            "n_species": len(rows),
        }

    enough = (len(rows) >= MIN_SPECIES_FOR_A_CLAIM
              and len(set(phyla)) >= MIN_PHYLA_FOR_A_CLAIM)
    if not coll["decoupled"]:
        verdict = "STILL_CONFOUNDED__NO_GC_CLAIM_POSSIBLE"
    elif not enough:
        verdict = "DECOUPLED_BUT_UNDERPOWERED"
    else:
        verdict = "DECOUPLED_AND_POWERED__GC_READING_IS_INTERPRETABLE"

    out = {
        "record": "glm-margin-gc-vs-lineage-v1",
        "date": date.today().isoformat(),
        "question": ("Is the winnable margin driven by GC content or by lineage? The four-species cache "
                     "could not separate them; this corpus can."),
        "what_this_closes": (
            "A declared limit: the cross-organism artifact reported margin-versus-GC at r = -0.941 and "
            "had to call it uninterpretable, because the two low-GC genomes were both Campylobacter and "
            "the three high-GC ones all Enterobacteriaceae -- effective n of 2 groups, not 5 points. It "
            "named what would be needed and said the set did not contain it."),
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "markov_order": MARKOV_ORDER,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "feature_modes": list(MODES), "seeds": list(SEEDS),
                     "note": "identical to the committed kernel, so the margin is the same quantity"},
        "n_species": len(rows), "n_phyla": len(set(phyla)),
        "gc_range": [round(min(gcs), 4), round(max(gcs), 4)],
        "collinearity_precondition": coll,
        "rows": rows,
        "analysis": analysis,
        "verdict": verdict,
        "honest_limits": [
            "Collinearity is tested BEFORE any GC reading is given, and the four-species case fails "
            "that test -- which is the non-vacuity check for the guard: it must reject the data that "
            "motivated building it.",
            "One genome per species, so a within-PHYLUM correlation is undefined wherever a phylum "
            "contributes a single genome. That is reported as undefined rather than as zero, because "
            "an undefined statistic read as a null result is how a sampling design becomes a finding.",
            "A correlation over a handful of phyla is a description of this corpus, not a law. GC is "
            "also only one axis of composition; it does not separate GC from codon usage, repeat "
            "content or gene density, all of which travel with it.",
            "This bounds the INSTRUMENT's margin against substrate composition and says nothing about "
            "any generator.",
            "Draft versus complete assemblies are mixed across the corpus, and contig fragmentation "
            "removes windows at breaks.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print()
    print(f"collinearity: corr(GC, phylum-mean GC) = {coll['corr_gc_vs_phylum_mean']} "
          f"over {coll.get('n_phyla')} phyla -> decoupled={coll['decoupled']}")
    for mode in MODES:
        v = analysis[mode]
        w = v["within_phylum"]
        print(f"  [{mode:10s}] raw r={v['raw_corr_margin_vs_gc']}  within-phylum r={w.get('corr')} "
              f"(groups with replicates {w['n_groups_with_replicates']})  "
              f"margin {v['margin_range']}  resolvable {v['n_resolvable']}/{v['n_species']}")
    print(f"\n{verdict}")
    print(f"{out['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
