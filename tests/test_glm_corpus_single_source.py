"""The corpus builder must have exactly ONE definition, and it must still reproduce the committed ceiling.

Step 1 of the seed-role plan moved `build_natural`'s body out of `scripts/glm_generator_falsifier.py` into
`dna_decode.glm.corpus.build_natural_set` so the seed-role probes read the same corpus the runner scores
on. The hazard being avoided is the one `geometry_from_args` was created for: two hand-written descriptions
of one thing, which diverged on a single field and were only caught afterwards. A second copy written for
the probes would have rebuilt it one layer down.

So these tests do two different jobs and both are load-bearing:
  1. the delegation is PROVED (identical window sequences), not assumed from reading the code;
  2. the moved body still reproduces the committed seed-7 self-control ceilings, which is the number every
     margin in this arm is measured against.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"

# The committed noise floors, from wiki/glm_generator_falsifier_REAL_2026-10-07.json.
COMMITTED_CEILING_KMER = 0.5447
COMMITTED_CEILING_BOTH = 0.5865
TOL = 5e-4   # the artifact rounds to 4 decimal places

requires_reference = pytest.mark.skipif(
    not (REFSEQ / "genome.fna").exists() or not (REFSEQ / "annotations.gff3").exists(),
    reason="local MG1655 reference absent (data/ is gitignored); see the runner's docstring to fetch it",
)


@requires_reference
def test_the_delegation_is_proved_not_assumed():
    """The runner's `build_natural` and `corpus.build_natural_set` must return the SAME windows.

    Reading the one-line delegation is not evidence: an argument could be dropped or reordered and the
    runner would still work, just on a different corpus. Comparing the sequences is what rules that out.
    """
    import importlib.util

    from dna_decode.glm.corpus import build_natural_set

    spec = importlib.util.spec_from_file_location(
        "_glm_falsifier", ROOT / "scripts" / "glm_generator_falsifier.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(spec and runner or runner)

    via_runner, st_a = runner.build_natural(150, 400, 7, kind="promoter")
    via_module, st_b = build_natural_set(REFSEQ, length=150, limit=400, seed=7, kind="promoter")

    assert [w.seq for w in via_runner] == [w.seq for w in via_module]
    assert len(via_runner) == 400
    assert st_a.n_cds_total == st_b.n_cds_total


@requires_reference
def test_the_moved_body_still_reproduces_the_committed_seed7_ceilings():
    """The eval half at the frozen geometry must still give 0.5447 (kmer) and 0.5865 (both).

    This is the gate the whole seed-role arm rests on. Every per-seed margin is `null - ceiling`, so a
    ceiling that moved would silently move ten margins at once.
    """
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves
    from dna_decode.glm.falsifier import self_control_ceiling

    wins, _ = build_natural_set(REFSEQ, length=150, limit=1200, seed=7, kind="promoter")
    _fit, ev = fit_eval_halves(wins)

    # `self_control_ceiling` returns (max_observed, all_values) -- the MAX is the ceiling, deliberately,
    # so that a lucky split cannot manufacture a finding.
    kmer_max, kmer_vals = self_control_ceiling(ev, k=4, mode="kmer")
    both_max, both_vals = self_control_ceiling(ev, k=4, mode="both", bins=6)

    assert kmer_max == pytest.approx(COMMITTED_CEILING_KMER, abs=TOL), (kmer_max, kmer_vals)
    assert both_max == pytest.approx(COMMITTED_CEILING_BOTH, abs=TOL), (both_max, both_vals)
    assert len(kmer_vals) == 5 and len(both_vals) == 5, "the 5-seed self-control did not run in full"


@requires_reference
def test_the_split_is_the_runners_split_not_a_re_derivation():
    """`fit_eval_halves` must be `len // 2`, fit first, eval second -- the runner's exact rule.

    An off-by-one here changes which windows every ceiling is computed on, so the halving is pinned rather
    than left to each probe to re-derive.
    """
    from dna_decode.glm.corpus import build_natural_set, fit_eval_halves

    wins, _ = build_natural_set(REFSEQ, length=150, limit=401, seed=7, kind="promoter")
    fit, ev = fit_eval_halves(wins)

    assert len(fit) == 200 and len(ev) == 201, (len(fit), len(ev))
    assert fit == [w.seq for w in wins[:200]]
    assert ev == [w.seq for w in wins[200:]]


def test_the_body_lives_in_exactly_one_place():
    """The runner must NOT carry its own copy of the build.

    Pinned on code shape rather than behaviour because the defect this step removes is DUPLICATION, which
    two identical-behaviour copies would pass.
    """
    runner_src = (ROOT / "scripts" / "glm_generator_falsifier.py").read_text(encoding="utf-8")
    body_start = runner_src.index("def build_natural(")
    body_end = runner_src.index("def geometry_from_args(")
    body = runner_src[body_start:body_end]

    assert "build_natural_set" in body, "the runner no longer delegates"
    for rebuilt in ("random_genomic_windows(", "extract_upstream_windows(", "parse_gff3("):
        assert rebuilt not in body, f"the runner rebuilt {rebuilt} instead of delegating"
