"""The cross-organism probe must select representatives blind, and must not let GC stand in for lineage.

Two things make this probe easy to get wrong, and both are pinned.

1. **Selection.** There are 400 *Campylobacter jejuni* assemblies on disk. Any representative rule that
   could see a result would let the probe produce whatever margin was wanted, so the rule is
   lexicographically-smallest-accession, fixed in advance, and tested here directly.
2. **The GC reading.** The low-GC substrates are both *Campylobacter* and the high-GC substrates are all
   *Enterobacteriaceae*, so a margin-against-GC correlation over five points is a two-group comparison
   with an effective n of 2. This repo has already closed four arms to exactly that confound, so the
   artifact must carry the lineage grouping and the effective n, not just the correlation.

The two built-in controls are also pinned, because they are what licenses reading the table as a
cross-ORGANISM result at all: reference-against-draft *E. coli* separates assembly quality from organism,
and the two congeners agreeing separates "a genus differs" from "one genome is odd".
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_margin_across_organisms.py"
OUT = ROOT / "wiki" / "glm_margin_across_organisms_2026-11-21.json"
CACHE_DIR = ROOT / "data" / "refseq_cache"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_across", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ geometry is the controlled variable

def test_the_geometry_is_the_committed_one_and_is_shared(mod):
    """The organism must be the ONLY thing that varies, or the table measures two things at once."""
    assert mod.LENGTH == 150
    assert mod.LIMIT == 1200
    assert mod.N_GENERATED == 120
    assert mod.KMER_K == 4
    assert mod.POSITIONAL_BINS == 6
    assert mod.MARKOV_ORDER == 3


def test_the_seed_set_is_fixed_for_every_substrate(mod):
    assert tuple(mod.SEEDS) == tuple(range(10))


# ------------------------------------------------------------------ the selection rule

def test_the_representative_rule_is_lexicographic_and_blind(mod):
    """Stated in the source and applied mechanically; with 400 candidates this is load-bearing."""
    src = SCRIPT.read_text(encoding="utf-8")
    assert "lexicographically smallest" in src
    assert "blind" in src


@pytest.mark.skipif(not CACHE_DIR.exists(), reason="refseq_cache absent on this host")
def test_the_rule_actually_returns_the_smallest_accession_per_organism(mod):
    """Non-vacuity: re-derive the choice independently and require the same answer.

    A rule that is documented but not followed is worse than no rule, because it reads as a guarantee.
    """
    reps = mod.pick_representatives()
    if not reps:
        pytest.skip("no annotated assemblies found")
    for organism, chosen in reps.items():
        candidates = []
        for d in sorted(CACHE_DIR.iterdir()):
            if not d.is_dir():
                continue
            if not ((d / "genome.fna").exists() and (d / "annotations.gff3").exists()):
                continue
            head = (d / "genome.fna").open(encoding="utf-8", errors="replace").readline()
            if organism in head:
                candidates.append(d.name)
        assert candidates, organism
        assert chosen.name == sorted(candidates)[0], (organism, chosen.name, sorted(candidates)[:3])


@pytest.mark.skipif(not CACHE_DIR.exists(), reason="refseq_cache absent on this host")
def test_at_most_one_representative_per_organism(mod):
    reps = mod.pick_representatives()
    assert len(set(reps.values())) == len(reps)


# ------------------------------------------------------------------ the resolve helper

def test_a_straddling_margin_is_not_resolvable(mod):
    out = mod.resolve([0.05, -0.03, 0.04, -0.01, 0.06, 0.00, -0.04, 0.05, 0.02, 0.07])
    assert out["resolvable_above_zero"] is False


def test_a_tight_positive_margin_resolves(mod):
    out = mod.resolve([0.14, 0.15, 0.13, 0.16, 0.14, 0.15, 0.13, 0.15, 0.14, 0.16])
    assert out["resolvable_above_zero"] is True and out["ci95"][0] > 0


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_more_than_one_organism_was_actually_scored():
    """A cross-organism claim from one organism would be the thing this probe exists to avoid."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    ok = {k: v for k, v in rec["by_substrate"].items() if v.get("status") == "ok"}
    assert len(ok) >= 3, sorted(ok)
    genera = {k.split()[0] for k in ok}
    assert len(genera) >= 2, genera


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_substrates_span_a_real_gc_range():
    """If every substrate had the same composition the probe would not be testing substrate at all."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    lo, hi = rec["controls"]["both"]["gc_range"]
    assert hi - lo > 0.15, (lo, hi)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_gc_reading_carries_its_lineage_confound_and_effective_n():
    """THE load-bearing honesty test. r = -0.94 over five points is two lineage groups."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        lc = rec["controls"][mode]["margin_vs_gc_is_confounded_with_lineage"]
        assert lc["effective_n_for_a_gc_claim"] == 2
        assert set(lc["by_family"]) == {"Campylobacteraceae", "Enterobacteriaceae"}
        assert "must NOT be quoted as a GC effect" in lc["reading"]
    joined = " ".join(rec["honest_limits"])
    assert "CONFOUNDED WITH PHYLOGENY" in joined
    assert "four arms" in joined, "the repo's own precedent for this confound must be named"


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_same_species_control_separates_assembly_quality_from_organism():
    """Reference against draft. If these disagreed, no cross-organism claim could be read at all."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        c = rec["controls"][mode]["same_species_reference_vs_draft"]
        assert c is not None, mode
        assert len(c["keys"]) == 2
        assert c["intervals_overlap"] is True, (mode, c)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_congener_control_is_present_and_agrees():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        g = rec["controls"][mode]["congeners_agree"]
        assert g is not None and g["intervals_overlap"] is True, (mode, g)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_both_mode_generalises_and_kmer_does_not_resolve_everywhere():
    """The finding, pinned: `both` holds at every organism; kmer is positive everywhere but
    under-powered at the committed n, which matches the sample-size sweep from a different direction."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["verdict_by_mode"]["both"] == "WINNABLE_REGION_AT_EVERY_ORGANISM_TESTED"
    ok = {k: v for k, v in rec["by_substrate"].items() if v.get("status") == "ok"}
    # every kmer margin is POSITIVE even where it does not resolve -- that is the shape that makes
    # "unresolved, not absent" the honest reading rather than a hedge.
    assert all(v["kmer"]["margin"]["mean"] > 0 for v in ok.values())
    assert all(v["both"]["margin"]["resolvable_above_zero"] for v in ok.values())


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_floor_sits_near_chance_at_every_organism():
    """The instrument check travels with the substrate: real-against-real is 0.5 in any genome."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    ok = {k: v for k, v in rec["by_substrate"].items() if v.get("status") == "ok"}
    for name, v in ok.items():
        for mode in ("kmer", "both"):
            assert abs(v[mode]["floor_mean"] - 0.5) < 0.04, (name, mode, v[mode]["floor_mean"])


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_draft_assembly_caveats_are_declared():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "DRAFT" in joined and "clinical-isolate" in joined
    assert "ONE representative genome per organism" in joined
