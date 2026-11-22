"""How much DISTINCT sequence do the 459 cached genomes actually hold? Measure it, don't estimate it.

    uv run python scripts/glm_corpus_distinct_content.py      # read-only; writes one artifact

**Why this exists.** A compute offer for generator pretraining was sized against the raw byte total of
`data/refseq_cache/` -- 0.92 Gbp across 459 annotated assemblies. That number is nearly useless as a
measure of training data, because 400 of the 459 are one species: the directory was built for AMR cohort
work, where 400 strains of *Campylobacter jejuni* is the POINT. A first pass estimated the distinct
content by taking one representative genome per species and got **~14.7 Mbp**, which was reported as an
UPPER-BOUND ESTIMATE with its two failure directions stated: it assumes genomes within a species are
fully redundant (they are not -- bacterial pangenomes carry real accessory content) and that genomes
across species share nothing (they do share some). This script replaces the estimate with a measurement.

**The method, and why this one.** Exact distinct-k-mer counting over 0.92 Gbp would need tens of
gigabytes of memory, so this uses **FracMinHash** -- the scaled-sketch technique sourmash is built on.
Every canonical k-mer is hashed and kept only when `hash % SCALE == 0`, so the sketch is a uniform
1-in-SCALE sample of the k-mer set and `distinct ~= |sketch| * SCALE` is unbiased. Memory falls by
SCALE while the estimate stays good, and because the sample is driven by the HASH rather than by
position, the same k-mer is always either kept or dropped in every genome -- which is exactly what makes
sketches comparable across genomes and lets a union be taken by unioning the sketches.

**What it reports is a RAREFACTION CURVE, which is the thing that actually answers the question.**
Cumulative distinct content as genomes are added, per species: if the curve flattens, the remaining
genomes of that species add almost nothing and the one-representative estimate was close; if it keeps
climbing, the pangenome's accessory content is real training signal and the estimate was too
conservative. The slope at the last genome measured bounds what the un-measured remainder can add, so
no extrapolation is needed and none is made.

**Three things that would silently break this and are guarded.** Canonical k-mers (the minimum of a
k-mer and its reverse complement) are required, or a genome's two strands double-count. The hash must be
DETERMINISTIC across processes -- Python's `hash()` on strings is randomised per process unless
`PYTHONHASHSEED` is set, so a 2-bit packing plus an explicit integer mixer is used instead. And k-mers
spanning a non-ACGT character must be dropped rather than coerced, which matters here because these are
draft assemblies with real ambiguity codes.

Free: numpy only, all data already on disk, no network, no GPU, no D: volume.
"""

from __future__ import annotations

import json
import statistics
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = ROOT / "data" / "refseq_cache"
REFSEQ_REF = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
OUT = ROOT / "wiki" / "glm_corpus_distinct_content_2026-11-21.json"

K = 21                  # the standard genome-comparison k; long enough to be near-unique in a genome
SCALE = 200             # keep 1 k-mer in 200; memory falls 200x, the estimate stays unbiased
MAX_PER_SPECIES = 12    # deterministic cap so the run stays minutes, not hours
ORGANISMS = ("Campylobacter jejuni", "Campylobacter coli",
             "Klebsiella pneumoniae", "Escherichia coli")

# Chinchilla-style compute-optimal ratio, used only to turn a corpus size into a model size.
TOKENS_PER_PARAM = 20


def _mix64(x):
    """splitmix64 finalizer on a numpy uint64 array. Deterministic, fast, well-distributed.

    Used instead of `hash()` because CPython randomises string hashing per process, which would make
    every sketch incomparable between runs and silently so -- the sketch sizes would look fine.
    """
    import numpy as np

    z = (x + np.uint64(0x9E3779B97F4A7C15)).astype(np.uint64)
    z = ((z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)).astype(np.uint64)
    z = ((z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)).astype(np.uint64)
    return (z ^ (z >> np.uint64(31))).astype(np.uint64)


def _segments(seq: str):
    """Split on any non-ACGT run, so no k-mer ever spans an ambiguity code.

    Draft assemblies really do carry N runs between contigs; coercing them to a base would invent
    sequence, and keeping them would make spurious k-mers that inflate the distinct count.
    """
    cur = []
    for ch in seq.upper():
        if ch in "ACGT":
            cur.append(ch)
        else:
            if len(cur) >= K:
                yield "".join(cur)
            cur = []
    if len(cur) >= K:
        yield "".join(cur)


_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def sketch_sequence(seq: str, *, k: int = K, scale: int = SCALE) -> set:
    """FracMinHash sketch of one sequence's CANONICAL k-mers.

    Vectorised with numpy: the k-mer integers are built by `k` shift-and-or passes over the whole
    2-bit-coded array at once, and the reverse complement is built the same way from the complemented,
    reversed codes. Canonicalising by elementwise minimum is what makes the sketch strand-independent.
    """
    import numpy as np

    keep = set()
    thresh = np.uint64(scale)
    for seg in _segments(seq):
        n = len(seg) - k + 1
        if n <= 0:
            continue
        codes = np.frombuffer(seg.encode("ascii"), dtype=np.uint8)
        lut = np.zeros(256, dtype=np.uint64)
        for ch, v in _CODE.items():
            lut[ord(ch)] = v
        c = lut[codes]
        fwd = np.zeros(n, dtype=np.uint64)
        for j in range(k):
            fwd = ((fwd << np.uint64(2)) | c[j:j + n]).astype(np.uint64)
        comp = (np.uint64(3) - c).astype(np.uint64)
        rev = np.zeros(n, dtype=np.uint64)
        for j in range(k):
            idx = k - 1 - j
            rev = ((rev << np.uint64(2)) | comp[idx:idx + n]).astype(np.uint64)
        canon = np.minimum(fwd, rev)
        h = _mix64(canon)
        sel = canon[(h % thresh) == np.uint64(0)]
        keep.update(int(v) for v in sel)
    return keep


def self_check() -> dict:
    """Non-vacuity: on random sequence the estimator must recover a KNOWN distinct count.

    A random 200 kb stretch has essentially no repeated 21-mers, so the true distinct canonical count is
    within a hair of the number of positions. If the scaled estimate did not land near that, the
    packing, the canonicalisation or the sampling would be wrong -- and every downstream number would
    still look plausible, which is exactly why this runs before anything else.
    """
    import random as _random

    rng = _random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(200_000))
    true_positions = len(seq) - K + 1
    est = len(sketch_sequence(seq)) * SCALE
    return {"random_sequence_bp": len(seq),
            "true_distinct_kmers_approx": true_positions,
            "scaled_estimate": est,
            "ratio": round(est / true_positions, 3),
            "passed": bool(0.9 <= est / true_positions <= 1.1)}


def genomes_by_organism() -> dict[str, list[Path]]:
    """Deterministic, lexicographic, blind to any downstream number."""
    found: dict[str, list[Path]] = {}
    if not CACHE_DIR.exists():
        return found
    for d in sorted(CACHE_DIR.iterdir()):
        if not d.is_dir():
            continue
        g, a = d / "genome.fna", d / "annotations.gff3"
        if not (g.exists() and a.exists()):
            continue
        try:
            header = g.open(encoding="utf-8", errors="replace").readline()
        except OSError:
            continue
        for name in ORGANISMS:
            if name in header:
                found.setdefault(name, []).append(d)
                break
    return found


def load_concat(ref_dir: Path) -> str:
    from dna_decode.glm.corpus import load_fasta

    return "".join(load_fasta(ref_dir / "genome.fna").values())


def main() -> int:
    try:
        import numpy  # noqa: F401
    except ImportError:
        print("REFUSE: numpy is required")
        return 2
    t0 = time.time()

    chk = self_check()
    print(f"self-check: {chk['scaled_estimate']:,} estimated vs {chk['true_distinct_kmers_approx']:,} "
          f"true on random sequence -> ratio {chk['ratio']} ({'PASS' if chk['passed'] else 'FAIL'})")
    if not chk["passed"]:
        print("REFUSE: the sketch estimator does not recover a known distinct count.")
        return 3

    by_org = genomes_by_organism()
    if not by_org:
        print(f"REFUSE: no annotated assemblies under {CACHE_DIR}")
        return 2
    print(f"k={K} scale={SCALE} cap={MAX_PER_SPECIES} per species")
    print(f"available: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(by_org.items())))
    print()

    per_species, union_all = {}, set()
    for org in sorted(by_org):
        chosen = by_org[org][:MAX_PER_SPECIES]
        cum, curve, singles = set(), [], []
        for i, d in enumerate(chosen, 1):
            s = sketch_sequence(load_concat(d))
            singles.append(len(s) * SCALE)
            before = len(cum)
            cum |= s
            curve.append({"n_genomes": i, "accession": d.name,
                          "genome_distinct_bp_est": len(s) * SCALE,
                          "cumulative_distinct_bp_est": len(cum) * SCALE,
                          "added_bp_est": (len(cum) - before) * SCALE})
        union_all |= cum
        first = curve[0]["cumulative_distinct_bp_est"]
        last = curve[-1]["cumulative_distinct_bp_est"]
        marginal = curve[-1]["added_bp_est"]
        per_species[org] = {
            "n_available": len(by_org[org]),
            "n_measured": len(chosen),
            "mean_single_genome_distinct_bp_est": round(statistics.fmean(singles)),
            "cumulative_distinct_bp_est": last,
            "growth_over_one_genome": round(last / first, 2) if first else None,
            "marginal_bp_from_last_genome": marginal,
            "marginal_as_fraction_of_total": round(marginal / last, 4) if last else None,
            "saturating": bool(last and marginal / last < 0.02),
            "curve": curve,
        }
        print(f"  {org:24s} n={len(chosen):2d}/{len(by_org[org]):3d}  "
              f"1 genome {first/1e6:5.2f} Mbp -> {len(chosen)} genomes {last/1e6:6.2f} Mbp "
              f"({per_species[org]['growth_over_one_genome']}x)   last genome added "
              f"{marginal/1e6:.3f} Mbp ({per_species[org]['marginal_as_fraction_of_total']*100:.2f}%)")

    # the reference genome, as a cross-check against the committed substrate
    ref = None
    if (REFSEQ_REF / "genome.fna").exists():
        rs = sketch_sequence(load_concat(REFSEQ_REF))
        shared = len(rs & union_all) * SCALE
        ref = {"accession": REFSEQ_REF.name,
               "distinct_bp_est": len(rs) * SCALE,
               "shared_with_cached_union_bp_est": shared,
               "fraction_of_reference_already_in_union": round(
                   shared / max(1, len(rs) * SCALE), 4)}
        union_all |= rs

    total = len(union_all) * SCALE
    sum_species = sum(v["cumulative_distinct_bp_est"] for v in per_species.values())
    supported = total / TOKENS_PER_PARAM

    rec = {
        "record": "glm-corpus-distinct-content-v1",
        "date": date.today().isoformat(),
        "question": ("How much DISTINCT sequence do the cached assemblies hold, as opposed to how many "
                     "bytes, and what model size does that support?"),
        "why": ("A pretrain was sized against a 0.92 Gbp raw byte total over a directory where 400 of "
                "459 assemblies are one species. An earlier one-representative-per-species estimate gave "
                "~14.7 Mbp and was flagged as an upper bound with both failure directions named. This "
                "measures it instead."),
        "method": {
            "sketch": "FracMinHash (scaled MinHash), the sourmash construction",
            "k": K, "scale": SCALE,
            "canonical": "min(kmer, reverse_complement(kmer)), so strand does not double-count",
            "hash": "2-bit packing plus a splitmix64 finalizer -- deterministic across processes, "
                    "unlike CPython's randomised string hash",
            "ambiguity": "sequence is split on any non-ACGT run; no k-mer spans an ambiguity code",
            "estimator": "distinct ~= |sketch| * scale, unbiased because the sample is hash-driven",
            "cap": f"{MAX_PER_SPECIES} genomes per species, lexicographic, chosen blind",
        },
        "self_check": chk,
        "per_species": per_species,
        "reference_cross_check": ref,
        "union": {
            "distinct_bp_est": total,
            "sum_of_per_species_bp_est": sum_species,
            "cross_species_overlap_bp_est": sum_species - total + (
                ref["distinct_bp_est"] if ref else 0),
            "note": ("the union is smaller than the sum of species because species share some k-mers; "
                     "the difference is that shared content and is reported rather than assumed away"),
        },
        "model_size_implied": {
            "tokens_per_param_assumed": TOKENS_PER_PARAM,
            "supported_params": round(supported),
            "supported_params_millions": round(supported / 1e6, 3),
            "shortfall_for": {f"{p/1e6:.0f}M": round(p * TOKENS_PER_PARAM / total, 1)
                              for p in (7e6, 20e6, 50e6)},
        },
        "honest_limits": [
            f"{MAX_PER_SPECIES} genomes per species, not all 459. The rarefaction curve's marginal "
            "addition at the last measured genome BOUNDS what the remainder can add, and no "
            "extrapolation beyond the measured points is made.",
            "A scaled sketch is an estimator. At this scale each per-genome figure carries a few "
            "percent of sampling error, which is far below the effects being discussed but is not zero. "
            "The self-check pins the estimator against a known distinct count on random sequence.",
            "Distinct k-mer content is a measure of SEQUENCE novelty, not of learnable regulatory "
            "signal. Two genomes can be k-mer-distinct and carry the same grammar, so this is an upper "
            "bound on useful training content and not a measure of it.",
            "The 20-tokens-per-parameter ratio is a compute-optimal heuristic from language-model "
            "scaling work, imported here and NOT validated for DNA. It is used to make the shortfall "
            "legible, and a different ratio moves the supported model size proportionally.",
            "k=21 is the conventional genome-comparison choice. A shorter k would find more sharing "
            "between species and a longer one less, so the cross-species overlap figure is k-dependent.",
            "These are draft assemblies of clinical isolates; contig fragmentation removes some k-mers "
            "that span breaks, so the distinct count is a slight UNDER-estimate per genome.",
        ],
        "wall_seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print()
    if ref:
        print(f"  reference {ref['accession']}: {ref['distinct_bp_est']/1e6:.2f} Mbp distinct, "
              f"{ref['fraction_of_reference_already_in_union']*100:.1f}% already in the cached union")
    print(f"\nUNION distinct: {total/1e6:.2f} Mbp "
          f"(sum of species {sum_species/1e6:.2f} Mbp)")
    print(f"supports ~{supported/1e6:.3f}M params at {TOKENS_PER_PARAM} tokens/param")
    for name, mult in rec["model_size_implied"]["shortfall_for"].items():
        print(f"   {name:>4s} params -> {mult}x short")
    print(f"\n{rec['wall_seconds']}s -> {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
