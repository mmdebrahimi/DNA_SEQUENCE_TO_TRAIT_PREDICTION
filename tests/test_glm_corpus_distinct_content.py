"""The distinct-content sketch must be strand-blind, deterministic, and recover a KNOWN count.

This probe exists to replace an estimate with a measurement, so the estimator itself is the thing most
worth testing. Three failure modes would all produce plausible numbers:

1. **Non-canonical k-mers.** If a k-mer and its reverse complement are counted separately, every genome
   double-counts and the corpus looks twice as rich as it is. The decisive test is that a sequence and
   its own reverse complement produce the SAME sketch.
2. **A non-deterministic hash.** CPython randomises string hashing per process, so a sketch built with
   `hash()` would be incomparable between runs while the sketch SIZES stayed perfectly reasonable --
   unions across genomes would silently stop sharing k-mers and the corpus would look far more diverse
   than it is. The 2-bit packing plus explicit mixer exists for that reason and is pinned here.
3. **Coerced ambiguity codes.** These are draft assemblies with real N runs. Treating an N as a base
   invents sequence; keeping it makes spurious k-mers. Either inflates the distinct count.

The self-check inside the script is the non-vacuity guard for the whole thing: on random sequence the
true distinct 21-mer count is essentially the number of positions, so a scaled estimate that misses it
means the packing, the canonicalisation or the sampling is wrong.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_corpus_distinct_content.py"
OUT = ROOT / "wiki" / "glm_corpus_distinct_content_2026-11-21.json"
CACHE_DIR = ROOT / "data" / "refseq_cache"

pytest.importorskip("numpy")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_distinct", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ------------------------------------------------------------------ the hash

def test_the_mixer_is_deterministic_across_calls(mod):
    """A per-process hash would make sketches incomparable between runs, silently."""
    import numpy as np

    x = np.arange(1000, dtype=np.uint64)
    a = mod._mix64(x)
    b = mod._mix64(x)
    assert np.array_equal(a, b)


def test_the_mixer_actually_mixes(mod):
    """Adjacent inputs must not land in adjacent buckets, or the 1-in-scale sample is positional."""
    import numpy as np

    x = np.arange(100_000, dtype=np.uint64)
    h = mod._mix64(x)
    kept = int(((h % np.uint64(mod.SCALE)) == np.uint64(0)).sum())
    # a uniform mixer keeps about 1/SCALE of them; a broken one keeps 0, all, or a periodic pattern
    expected = 100_000 / mod.SCALE
    assert 0.6 * expected < kept < 1.6 * expected, (kept, expected)


def test_the_mixer_does_not_collapse_to_the_identity(mod):
    import numpy as np

    x = np.arange(64, dtype=np.uint64)
    assert not np.array_equal(mod._mix64(x), x)


# ------------------------------------------------------------------ canonical k-mers

def test_a_sequence_and_its_reverse_complement_give_the_SAME_sketch(mod):
    """THE test for canonicalisation. Without it every genome double-counts its two strands."""
    import random

    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    a = mod.sketch_sequence(seq)
    b = mod.sketch_sequence(_revcomp(seq))
    assert a == b, (len(a), len(b), len(a ^ b))
    assert len(a) > 20, "sketch too small for this to be a meaningful comparison"


def test_two_unrelated_sequences_do_NOT_give_the_same_sketch(mod):
    """Control for the test above: equality must mean something."""
    import random

    rng = random.Random(12)
    a = mod.sketch_sequence("".join(rng.choice("ACGT") for _ in range(20_000)))
    b = mod.sketch_sequence("".join(rng.choice("ACGT") for _ in range(20_000)))
    assert a != b
    assert len(a & b) < 0.05 * max(len(a), len(b)), "unrelated random sequence should barely overlap"


def test_an_identical_repeat_does_not_inflate_the_distinct_count(mod):
    """Distinct content, not length: concatenating a sequence with itself adds almost nothing."""
    import random

    rng = random.Random(13)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    once = mod.sketch_sequence(seq)
    twice = mod.sketch_sequence(seq + seq)
    # the junction contributes a handful of genuinely new k-mers, nothing more
    assert len(twice) <= len(once) + 5, (len(once), len(twice))


# ------------------------------------------------------------------ ambiguity codes

def test_segments_split_on_non_acgt_and_drop_short_runs(mod):
    k = mod.K
    seq = "A" * (k + 5) + "N" * 10 + "C" * (k + 7) + "N" + "G" * 3
    segs = list(mod._segments(seq))
    assert len(segs) == 2, segs
    assert set(segs[0]) == {"A"} and set(segs[1]) == {"C"}
    assert all(len(s) >= k for s in segs)


def test_no_kmer_spans_an_ambiguity_code(mod):
    """An N-separated pair must sketch to exactly the union of the two halves sketched apart."""
    import random

    rng = random.Random(14)
    a = "".join(rng.choice("ACGT") for _ in range(5_000))
    b = "".join(rng.choice("ACGT") for _ in range(5_000))
    joined = mod.sketch_sequence(a + "N" * 5 + b)
    apart = mod.sketch_sequence(a) | mod.sketch_sequence(b)
    assert joined == apart


def test_lowercase_is_accepted_not_dropped(mod):
    """Soft-masked FASTA is common; dropping lowercase would silently shrink every genome."""
    import random

    rng = random.Random(15)
    seq = "".join(rng.choice("ACGT") for _ in range(8_000))
    assert mod.sketch_sequence(seq.lower()) == mod.sketch_sequence(seq)


# ------------------------------------------------------------------ the estimator

def test_the_self_check_recovers_a_known_distinct_count(mod):
    """Non-vacuity for the whole pipeline, run against sequence whose answer is known in advance."""
    chk = mod.self_check()
    assert chk["passed"] is True, chk
    assert 0.9 <= chk["ratio"] <= 1.1, chk


def test_a_coarser_scale_gives_a_proportionally_smaller_sketch(mod):
    """If sketch size did not track 1/scale, the `|sketch| * scale` estimator would be wrong."""
    import random

    rng = random.Random(16)
    seq = "".join(rng.choice("ACGT") for _ in range(60_000))
    fine = len(mod.sketch_sequence(seq, scale=50))
    coarse = len(mod.sketch_sequence(seq, scale=500))
    assert fine > coarse
    assert 4 < fine / max(1, coarse) < 25, (fine, coarse)


# ------------------------------------------------------------------ selection

@pytest.mark.skipif(not CACHE_DIR.exists(), reason="refseq_cache absent on this host")
def test_genome_selection_is_lexicographic_and_stable(mod):
    a = mod.genomes_by_organism()
    b = mod.genomes_by_organism()
    assert {k: [p.name for p in v] for k, v in a.items()} == {
        k: [p.name for p in v] for k, v in b.items()}
    for org, paths in a.items():
        names = [p.name for p in paths]
        assert names == sorted(names), org


@pytest.mark.skipif(not CACHE_DIR.exists(), reason="refseq_cache absent on this host")
def test_every_selected_genome_carries_both_required_files(mod):
    for org, paths in mod.genomes_by_organism().items():
        for p in paths:
            assert (p / "genome.fna").exists() and (p / "annotations.gff3").exists(), (org, p.name)


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_artifact_records_a_passing_self_check():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["self_check"]["passed"] is True


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_rarefaction_curve_is_cumulative_and_non_decreasing():
    """A union cannot shrink as genomes are added; if it did, the sketches are not comparable."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for org, v in rec["per_species"].items():
        cum = [c["cumulative_distinct_bp_est"] for c in v["curve"]]
        assert cum == sorted(cum), org
        assert all(c["added_bp_est"] >= 0 for c in v["curve"]), org


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_union_is_no_larger_than_the_sum_of_species():
    """Species share some k-mers, so the union must be smaller; a union LARGER than the sum would
    mean the per-species numbers and the union were computed on different things."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    u = rec["union"]
    assert u["distinct_bp_est"] <= u["sum_of_per_species_bp_est"] + (
        rec["reference_cross_check"]["distinct_bp_est"] if rec.get("reference_cross_check") else 0)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_no_extrapolation_beyond_the_measured_genomes():
    """The cap is real and the artifact must say the remainder is BOUNDED, not estimated."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "BOUNDS" in joined or "bounds" in joined
    assert "no extrapolation" in joined.lower()
    for org, v in rec["per_species"].items():
        assert v["n_measured"] <= v["n_available"], org


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_imported_scaling_heuristic_is_flagged_as_unvalidated_for_dna():
    """20 tokens/param comes from language-model scaling work and is NOT established for DNA.
    The supported-model-size figure moves proportionally with it, so the import must be visible."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "NOT validated for DNA" in joined
    assert rec["model_size_implied"]["tokens_per_param_assumed"] == 20


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_distinct_kmers_are_not_claimed_to_measure_learnable_signal():
    """Two genomes can be k-mer-distinct and carry the same regulatory grammar. This is an upper
    bound on useful training content, and over-reading it is the obvious way to misuse the number."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "upper bound" in joined
    assert "not a measure of it" in joined or "not of learnable" in joined
