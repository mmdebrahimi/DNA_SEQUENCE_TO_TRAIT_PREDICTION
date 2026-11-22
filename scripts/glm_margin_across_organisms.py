"""Does the winnable region survive a change of ORGANISM, at the committed geometry?

    uv run python scripts/glm_margin_across_organisms.py      # read-only; writes one artifact

**The limit this closes.** Every seed-role and margin result in this arm carries "one organism, one
geometry, one null" as an honest limit, and the geometry half has now been swept
(`glm_margin_sample_size_sweep`). The organism half was never tested, even though it is the limit most
likely to matter: the margin is a property of a (substrate, geometry, feature-mode) triple, and nothing
said whether *E. coli* K-12 is representative or peculiar.

**It turns out not to need any download.** `data/refseq_cache/` holds 459 assemblies that carry both a
`genome.fna` and an `annotations.gff3` -- the AMR-cohort genomes -- and they are not all one species:
400 *Campylobacter jejuni*, 26 *C. coli*, 17 *Klebsiella pneumoniae*, 16 *E. coli*. That spans a GC range
from **0.266 to 0.506**, which is the substrate property most likely to move a composition-based
discriminator, and all four support the committed geometry (1200 windows plus 120 spares) with room to
spare.

**Two controls are built into the substrate list rather than argued for.**
  * The MG1655 REFERENCE and an *E. coli* DRAFT assembly are both included. Same species, different
    assembly quality and a different sampling frame (a complete reference genome against a clinical
    isolate's contigs). If their margins agree, assembly quality is not driving the result; if they
    disagree, no cross-ORGANISM claim can be read from the rest of the table.
  * *C. jejuni* and *C. coli* are congeners. Two members of one genus agreeing with each other while
    differing from the *Enterobacteriaceae* is a different and much stronger pattern than one odd
    genome.

**Representative selection is a stated rule, applied before any number was seen:** the
lexicographically smallest accession among the assemblies of that organism that carry both files. It is
not a choice made to favour an outcome, and with 400 *C. jejuni* available the rule matters -- picking
the "best-looking" genome from 400 would be a selection effect large enough to manufacture any result.

Geometry is held IDENTICAL to the committed kernel across every substrate, so the only thing varying is
the organism. Free: pure Python plus sklearn, all data already on disk, no network, no GPU, no D: volume.
"""

from __future__ import annotations

import json
import math
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ_REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
CACHE_DIR = ROOT / "data" / "refseq_cache"
OUT = ROOT / "wiki" / "glm_margin_across_organisms_2026-11-21.json"

# The committed kernel geometry, held identical across every substrate.
LENGTH = 150
LIMIT = 1200
N_GENERATED = 120
KIND = "promoter"
KMER_K = 4
POSITIONAL_BINS = 6
MARKOV_ORDER = 3
MODES = ("kmer", "both")
SEEDS = tuple(range(10))

ORGANISMS = ("Campylobacter jejuni", "Campylobacter coli",
             "Klebsiella pneumoniae", "Escherichia coli")


def require_sklearn() -> None:
    import numpy  # noqa: F401
    import sklearn  # noqa: F401


def pick_representatives() -> dict[str, Path]:
    """Lexicographically smallest accession per organism, among those carrying BOTH required files.

    The rule is fixed in advance and applied mechanically. With 400 *C. jejuni* assemblies on disk, any
    rule that involved looking at results first could produce whatever margin was wanted, so the
    selection is deliberately blind to every downstream number.
    """
    found: dict[str, Path] = {}
    if not CACHE_DIR.exists():
        return found
    for d in sorted(CACHE_DIR.iterdir()):
        if not d.is_dir():
            continue
        if not ((d / "genome.fna").exists() and (d / "annotations.gff3").exists()):
            continue
        try:
            header = (d / "genome.fna").open(encoding="utf-8", errors="replace").readline()
        except OSError:
            continue
        for name in ORGANISMS:
            if name in header and name not in found:
                found[name] = d
    return found


def substrate_cell(ref_dir: Path, seed: int) -> dict:
    """One (substrate, seed) cell at the committed geometry: null, perfect-generator floor, margin.

    Spares are filtered by SEQUENCE, not index: a corpus can hold duplicate upstream windows (operon
    structure, duplicated genes, repeated elements) and the discriminator sees only sequence, so an
    index-disjoint spare that is sequence-identical to an eval window would put the same feature vector
    in both classes and make the floor read cleaner than it is.
    """
    from dna_decode.glm.baselines import MarkovGenerator
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves, gc_fraction
    from dna_decode.glm.falsifier import heldout_auroc

    wins, stats = build_natural_set(ref_dir, length=LENGTH, limit=None, seed=seed, kind=KIND)
    if len(wins) < LIMIT + N_GENERATED:
        return {"seed": seed, "status": "INSUFFICIENT_WINDOWS", "n_windows": len(wins)}
    committed = wins[:LIMIT]
    fit, ev = fit_eval_halves(committed)
    used = set(fit) | set(ev)
    spare, skipped = [], 0
    for w in wins[LIMIT:]:
        if len(spare) >= N_GENERATED:
            break
        if w.seq in used:
            skipped += 1
            continue
        spare.append(w.seq)
        used.add(w.seq)
    if len(spare) < N_GENERATED:
        return {"seed": seed, "status": "INSUFFICIENT_DISJOINT_SPARES", "n_spare": len(spare)}

    gen = MarkovGenerator(order=MARKOV_ORDER, seed=seed).fit(fit).generate(N_GENERATED, LENGTH)
    out = {"seed": seed, "status": "ok",
           "n_windows_available": len(wins),
           "n_fit": len(fit), "n_eval": len(ev),
           "n_spare_skipped_as_duplicates": skipped,
           "gc_natural": round(sum(gc_fraction(s) for s in ev) / len(ev), 4),
           "gc_generated": round(sum(gc_fraction(s) for s in gen) / len(gen), 4),
           "n_cds_total": stats.n_cds_total,
           "n_overlapping_cds": stats.n_overlapping_cds}
    for mode in MODES:
        kw = {"k": KMER_K, "seed": seed, "mode": mode}
        if mode == "both":
            kw["bins"] = POSITIONAL_BINS
        null = heldout_auroc(ev, gen, **kw)
        floor = heldout_auroc(ev, spare, **kw)
        out[mode] = {
            "null": round(null, 4) if null is not None else None,
            "floor": round(floor, 4) if floor is not None else None,
            "margin": round(null - floor, 4) if (null is not None and floor is not None) else None,
        }
    return out


def resolve(margins: list[float]) -> dict:
    n = len(margins)
    if n < 3:
        return {"n": n, "resolvable_above_zero": None}
    mean = statistics.fmean(margins)
    sd = statistics.stdev(margins)
    se = sd / math.sqrt(n)
    crit = 2.262 if n == 10 else 1.96
    lo, hi = mean - crit * se, mean + crit * se
    return {"n": n, "mean": round(mean, 4), "sd": round(sd, 4), "stderr": round(se, 4),
            "ci95": [round(lo, 4), round(hi, 4)],
            "t": round(mean / se, 3) if se else None,
            "resolvable_above_zero": bool(lo > 0),
            "n_positive": sum(1 for v in margins if v > 0)}


def main() -> int:
    require_sklearn()
    t0 = time.time()

    substrates: dict[str, Path] = {}
    if (REFSEQ_REF / "genome.fna").exists():
        substrates["Escherichia coli K-12 MG1655 (REFERENCE)"] = REFSEQ_REF
    reps = pick_representatives()
    for name, d in sorted(reps.items()):
        substrates[f"{name} ({d.name}, draft)"] = d
    if len(substrates) < 2:
        print("REFUSE: fewer than two substrates available; a cross-organism claim needs at least two.")
        return 2

    print(f"geometry held identical across substrates: length={LENGTH} limit={LIMIT} "
          f"n={N_GENERATED} order={MARKOV_ORDER} k={KMER_K} bins={POSITIONAL_BINS}")
    print(f"substrates ({len(substrates)}):")
    for name, d in substrates.items():
        print(f"  {name}  <- {d}")
    print()

    rows, by_substrate = [], {}
    for name, d in substrates.items():
        cells = []
        for s in SEEDS:
            c = substrate_cell(d, s)
            c["substrate"] = name
            cells.append(c)
            rows.append(c)
        ok = [c for c in cells if c["status"] == "ok"]
        if not ok:
            by_substrate[name] = {"status": cells[0]["status"]}
            print(f"  {name}: {cells[0]['status']}")
            continue
        entry = {
            "status": "ok",
            "accession": d.name,
            "n_seeds_ok": len(ok),
            "n_windows_available": ok[0]["n_windows_available"],
            "gc_natural_mean": round(statistics.fmean([c["gc_natural"] for c in ok]), 4),
            "gc_generated_mean": round(statistics.fmean([c["gc_generated"] for c in ok]), 4),
            "n_cds_total": ok[0]["n_cds_total"],
            "overlap_fraction": round(ok[0]["n_overlapping_cds"] / max(1, ok[0]["n_cds_total"]), 3),
        }
        for mode in MODES:
            nulls = [c[mode]["null"] for c in ok if c[mode]["null"] is not None]
            floors = [c[mode]["floor"] for c in ok if c[mode]["floor"] is not None]
            margins = [c[mode]["margin"] for c in ok if c[mode]["margin"] is not None]
            entry[mode] = {
                "null_mean": round(statistics.fmean(nulls), 4),
                "null_sd": round(statistics.stdev(nulls), 4) if len(nulls) > 1 else None,
                "floor_mean": round(statistics.fmean(floors), 4),
                "floor_sd": round(statistics.stdev(floors), 4) if len(floors) > 1 else None,
                "margin": resolve(margins),
            }
        by_substrate[name] = entry
        for mode in MODES:
            e = entry[mode]
            print(f"  {name[:42]:42s} [{mode:5s}] GC {entry['gc_natural_mean']:.3f} "
                  f"null {e['null_mean']:.4f} floor {e['floor_mean']:.4f} -> margin "
                  f"{e['margin']['mean']:+.4f} CI {e['margin']['ci95']} "
                  f"resolvable={e['margin']['resolvable_above_zero']}")

    good = {k: v for k, v in by_substrate.items() if v.get("status") == "ok"}

    # ---- the controls, read mechanically ----
    ecoli_keys = [k for k in good if "Escherichia coli" in k]
    campy_keys = [k for k in good if "Campylobacter" in k]
    controls = {}
    for mode in MODES:
        same_species = None
        if len(ecoli_keys) >= 2:
            vals = [good[k][mode]["margin"]["mean"] for k in ecoli_keys]
            # Overlapping intervals is the honest test, not equality of two point estimates.
            cis = [good[k][mode]["margin"]["ci95"] for k in ecoli_keys]
            overlap = max(c[0] for c in cis) <= min(c[1] for c in cis)
            same_species = {"keys": ecoli_keys, "margins": vals,
                            "spread": round(max(vals) - min(vals), 4),
                            "intervals_overlap": bool(overlap)}
        congeners = None
        if len(campy_keys) >= 2:
            vals = [good[k][mode]["margin"]["mean"] for k in campy_keys]
            cis = [good[k][mode]["margin"]["ci95"] for k in campy_keys]
            congeners = {"keys": campy_keys, "margins": vals,
                         "spread": round(max(vals) - min(vals), 4),
                         "intervals_overlap": bool(
                             max(c[0] for c in cis) <= min(c[1] for c in cis))}
        margins = [v[mode]["margin"]["mean"] for v in good.values()]
        gcs = [v["gc_natural_mean"] for v in good.values()]
        # GC is CONFOUNDED WITH PHYLOGENY here and the confound is reported as a grouping, not a
        # caveat sentence. The low-GC substrates are both Campylobacter and the high-GC ones are all
        # Enterobacteriaceae, so a correlation over five points is really a comparison of two
        # lineage groups and its effective n is 2. Quoting r = -0.94 as a GC effect would be the
        # lineage confound that has already closed four arms in this project.
        lineage = {}
        for k, v in good.items():
            fam = ("Campylobacteraceae" if "Campylobacter" in k else "Enterobacteriaceae")
            lineage.setdefault(fam, []).append(v[mode]["margin"]["mean"])
        controls[mode] = {
            "margin_vs_gc_is_confounded_with_lineage": {
                "by_family": {f: {"n": len(v), "margins": v,
                                  "mean": round(statistics.fmean(v), 4)}
                              for f, v in sorted(lineage.items())},
                "effective_n_for_a_gc_claim": len(lineage),
                "reading": ("the low-GC substrates are both Campylobacter and the high-GC substrates "
                            "are all Enterobacteriaceae, so the GC correlation below is a two-group "
                            "comparison dressed as five points. It is reported as a description of "
                            "this table and must NOT be quoted as a GC effect."),
            },
            "same_species_reference_vs_draft": same_species,
            "congeners_agree": congeners,
            "n_substrates": len(good),
            "n_resolvable": sum(1 for v in good.values()
                                if v[mode]["margin"]["resolvable_above_zero"]),
            "margin_min": round(min(margins), 4),
            "margin_max": round(max(margins), 4),
            "margin_range": round(max(margins) - min(margins), 4),
            "gc_range": [round(min(gcs), 4), round(max(gcs), 4)],
            "margin_vs_gc_correlation": (round(statistics.correlation(gcs, margins), 3)
                                         if len(good) > 2 else None),
        }

    def verdict(mode: str) -> str:
        c = controls[mode]
        if c["n_resolvable"] == c["n_substrates"]:
            return "WINNABLE_REGION_AT_EVERY_ORGANISM_TESTED"
        if c["n_resolvable"] == 0:
            return "NO_RESOLVABLE_WINNABLE_REGION_AT_ANY_ORGANISM"
        return "WINNABLE_REGION_IS_ORGANISM_DEPENDENT"

    rec = {
        "record": "glm-margin-across-organisms-v1",
        "date": date.today().isoformat(),
        "question": ("Does the winnable region measured on E. coli K-12 survive a change of organism, "
                     "at the committed geometry?"),
        "why_it_was_open": ("Every result in this arm carried 'one organism' as an honest limit. The "
                            "geometry half of that limit was swept; the organism half was never tested, "
                            "although the margin is explicitly a property of a (substrate, geometry, "
                            "feature-mode) triple."),
        "geometry": {"window_length_bp": LENGTH, "limit_natural": LIMIT,
                     "n_generated_per_candidate": N_GENERATED, "natural_set_arg": KIND,
                     "kmer_k": KMER_K, "positional_bins": POSITIONAL_BINS,
                     "markov_order": MARKOV_ORDER, "feature_modes": list(MODES),
                     "seeds": list(SEEDS),
                     "note": "held IDENTICAL across substrates so the organism is the only variable"},
        "representative_selection_rule": (
            "lexicographically smallest accession among assemblies of that organism carrying both "
            "genome.fna and annotations.gff3. Fixed in advance and applied blind to every downstream "
            "number, because with 400 C. jejuni assemblies on disk a results-aware rule could produce "
            "any margin wanted."),
        "substrates": {k: str(v) for k, v in substrates.items()},
        "rows": rows,
        "by_substrate": by_substrate,
        "controls": controls,
        "verdict_by_mode": {m: verdict(m) for m in MODES},
        "honest_limits": [
            "ONE representative genome per organism. Within-organism variation is not measured, so an "
            "organism-level claim rests on a single assembly for that organism -- except for E. coli, "
            "which has two (reference and draft) and is the only substrate with any internal replicate.",
            "The non-reference genomes are DRAFT assemblies from the AMR cohorts: contig-level, a "
            "different annotation pipeline, and a clinical-isolate sampling frame rather than a type "
            "strain. The reference-against-draft E. coli pair is included specifically so that "
            "confound is visible rather than assumed away.",
            "An upstream window is a PROXY for a promoter, not a mapped transcription start site, and "
            "that proxy is weaker in compact AT-rich genomes where overlap with the preceding CDS is "
            "common. The overlap fraction is reported per substrate for exactly that reason.",
            "Ten seeds per substrate, so each interval is itself imprecise.",
            "GC content is CONFOUNDED WITH PHYLOGENY in this substrate set: the two low-GC genomes are "
            "both Campylobacter and the three high-GC genomes are all Enterobacteriaceae. The "
            "margin-against-GC correlation is therefore a two-group comparison with an effective n of "
            "2, not a five-point regression, and is reported only as a description of this table. "
            "Reading it as a GC effect would repeat the lineage confound that has already closed four "
            "arms in this project. Separating GC from lineage needs substrates that cross them -- a "
            "high-GC Campylobacter or a low-GC enterobacterium -- which this set does not contain.",
            "The Markov null is REFITTED per substrate, which is correct (a chain fitted on E. coli is "
            "not the right null for Campylobacter) but means the null's strength varies with the "
            "substrate as well as the discriminator's task.",
            "A winnable region is necessary, NOT sufficient: these features see LOCAL composition, so "
            "clearing them would not establish functional realism in any organism.",
            "One geometry, one null order, one window length.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    for mode in MODES:
        c = controls[mode]
        print(f"[{mode}] {verdict(mode)}")
        print(f"      resolvable at {c['n_resolvable']}/{c['n_substrates']} substrates | margin range "
              f"{c['margin_min']:+.4f}..{c['margin_max']:+.4f} | GC range {c['gc_range']}")
        lc = c["margin_vs_gc_is_confounded_with_lineage"]
        fams = "  ".join(f"{f.split('aceae')[0]}aceae n={v['n']} mean {v['mean']:+.4f}"
                         for f, v in lc["by_family"].items())
        print(f"      margin-vs-GC r {c['margin_vs_gc_correlation']} but CONFOUNDED with lineage "
              f"(effective n={lc['effective_n_for_a_gc_claim']}): {fams}")
        if c["same_species_reference_vs_draft"]:
            s = c["same_species_reference_vs_draft"]
            print(f"      CONTROL same-species reference vs draft: margins {s['margins']} "
                  f"(spread {s['spread']}, intervals overlap {s['intervals_overlap']})")
        if c["congeners_agree"]:
            g = c["congeners_agree"]
            print(f"      CONTROL congeners: margins {g['margins']} "
                  f"(spread {g['spread']}, intervals overlap {g['intervals_overlap']})")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
