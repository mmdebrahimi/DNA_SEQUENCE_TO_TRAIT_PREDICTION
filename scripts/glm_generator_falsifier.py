"""Score GLM sequence generators against real E. coli regulatory sequence. Read-only; writes an artifact.

    uv run python scripts/glm_generator_falsifier.py                  # baselines only, no download
    uv run python scripts/glm_generator_falsifier.py --model auto     # + the real HF generator

**The split that keeps this honest.** The natural corpus is divided into a FIT half and an EVAL half. The
Markov null is fitted on FIT only and every discrimination is against EVAL only, so no generator is ever
compared against sequences it was trained on. Without that, the null would have memorised its own
comparison set and its distinguishability would be flattered.

Needs nothing but a LOCALLY CACHED MG1655 reference at `data/cache/refseq/GCF_000005845.2/{genome.fna,
annotations.gff3}`. No network, no GPU, no labels, once it is there.

**That reference is NOT tracked in git** (`data/` is gitignored; `git ls-files data/cache/refseq/` is
empty), so a fresh clone does not have it and this script REFUSES rather than running on nothing. Fetch it
with the repo's own downloader -- the argument is the PARENT directory, which appends the accession folder:

    from dna_decode.data.refseq import download_genome
    download_genome("GCF_000005845.2", "data/cache/refseq")

This docstring said "the committed MG1655 reference" until 2026-10-07, which is the `--help`-over-claim
pattern: it cost a cross-machine handoff a round trip, because a reader on a fresh clone has every reason
to believe a doc that says the input is already there.
"""
from __future__ import annotations

import argparse
import json
import random
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFSEQ = ROOT / "data" / "cache" / "refseq" / "GCF_000005845.2"
WIKI = ROOT / "wiki"


def build_natural(length: int, limit: int | None, seed: int, kind: str = "promoter"):
    """Thin delegation to `dna_decode.glm.corpus.build_natural_set` -- the ONE definition.

    The body moved there so the seed-role probes read the same corpus this runner scores on, rather than
    growing a second description of it. This wrapper stays because it is this script's established entry
    point and its signature is what `main` and the tests call.
    """
    from dna_decode.glm.corpus import build_natural_set

    return build_natural_set(REFSEQ, length=length, limit=limit, seed=seed, kind=kind)


def geometry_from_args(a) -> dict:
    """The run's FROZEN geometry — the single expression both consumers read.

    There used to be two: the `--compare-to` projection and the artifact dict, hand-written in
    parallel. They diverged on exactly one field. The artifact recorded `positional_bins` as
    `a.bins` unconditionally while the projection recorded `None` unless the mode was
    positional/both, so a `--modes kmer` run compared against its OWN kmer-only baseline was
    REFUSED with `{'positional_bins': (10, None)}` — a false mismatch on a comparable run. It
    could not bite the `kmer,both` runs in flight at the time (both branches yield `a.bins`),
    which is why it survived review of the code that introduced it.

    Changing the one literal would have fixed that instance and left the mechanism — two parallel
    descriptions of one run — fully intact. One expression is what makes a recurrence impossible,
    so the artifact spreads this dict rather than restating its fields.

    Keys are exactly `knob_contract.COMPARISON_GEOMETRY_FIELDS`; a test pins that correspondence,
    because a field named in the tuple but missing here would silently compare as `<absent>`.
    """
    return {
        "window_length_bp": a.length,
        "kmer_k": a.k,
        "positional_bins": a.bins,
        "markov_order": a.markov_order,
        "n_generated_per_candidate": a.n,
        "natural_set_arg": a.natural_set,
        "feature_modes": a.modes,
        "conditional": bool(a.conditional),
        "context_bp": a.context if a.conditional else None,
    }


def make_knob_probe(gen, prompts, length: int):
    """The `generate_at(temperature)` callable `verify_knob_live` drives.

    Module-level, and deliberately so: as a closure inside `main` the choice of generation path
    could only be pinned by reading the source, and this is the one decision in the verification
    that was WRONG. It called the unconditional `generate` while the sweep runs
    `generate_conditional`, which proves the knob moves output on a path no grid point takes.

    One real prompt suffices. Both paths reseed to `self.seed` on entry (`generate.py:222`,
    `:266`), so the only difference between the probe's two calls is temperature.

    The unconditional fallback is unreachable today -- `validate_comparison_mode` leaves
    `promoter + conditional` as the single valid pairing -- but it is kept so the helper is not
    silently wrong if that contract ever widens.
    """
    def _at(t: float):
        gen.temperature = t
        if prompts:
            return gen.generate_conditional(prompts[:1], length)
        return gen.generate(1, length, progress_every=0)
    return _at


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--length", type=int, default=150, help="upstream window length (bp)")
    p.add_argument("--n", type=int, default=400, help="sequences to generate per candidate")
    p.add_argument("--limit-natural", type=int, default=2000)
    p.add_argument("--k", type=int, default=4, help="k-mer size for the discriminator")
    p.add_argument("--markov-order", type=int, default=3)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--modes", default="kmer,positional",
                   help="comma-separated feature modes: kmer | positional | both")
    p.add_argument("--bins", type=int, default=10, help="positional bins")
    p.add_argument("--natural-set", default="promoter",
                   help="promoter | genomic | both -- which natural reference to score against")
    p.add_argument("--model", default=None,
                   help="HF model id, or 'auto' for the default GENERator prokaryote checkpoint")
    p.add_argument("--device", default="auto", help="auto | cpu | cuda")
    p.add_argument("--conditional", action="store_true",
                   help="prompt with real upstream context and score ONLY the continuation; required "
                        "for a meaningful comparison against the ALIGNED promoter set")
    p.add_argument("--context", type=int, default=300, help="prompt length (bp) for --conditional")
    p.add_argument("--cache-dir", default="D:/hf_cache", help="weights cache (keep off C:)")
    # The sampling seam. `dna_decode.glm.generate.HFGenerator` has accepted `temperature` and
    # `top_k` since it was written, but this CLI never passed them, so every recorded run is "at
    # default sampling" BY OMISSION, not by choice. G-D's unblock trigger is a sampling sweep
    # reporting first (so sampling and weights stay separable levers), and that sweep was
    # unrunnable from here without these two flags.
    p.add_argument("--temperature", type=float, default=1.0,
                   help="sampling temperature passed to HFGenerator (1.0 = the recorded default)")
    p.add_argument("--top-k", type=int, default=0,
                   help="top-k truncation passed to HFGenerator (0 = disabled, the recorded default)")
    p.add_argument("--verify-knob", action="store_true",
                   help="BEFORE sweeping, prove the sampling knob moves the OUTPUT (one extra "
                        "model-backed 2-value diff). Required to claim a sweep VARIED sampling; "
                        "without it the artifact records knob_verified_end_to_end=false.")
    p.add_argument("--compare-to", default=None,
                   help="path to a recorded artifact this run will be compared against. Refuses "
                        "BEFORE any work if the frozen geometry differs (the swept `sampling` and "
                        "`device` are deliberately exempt). Added after --bins was left at its "
                        "default 10 against a baseline recorded at 6, voiding a 25-min run.")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
    from dna_decode.glm.falsifier import run_falsifier
    from dna_decode.glm.generate import device_report

    if not (REFSEQ / "genome.fna").exists():
        print(f"REFUSE: reference not found at {REFSEQ}")
        return 2

    # Refuse a meaningless (reference, conditioning) pairing BEFORE any corpus build or
    # model load -- milliseconds, versus ~30 min of generation for an unusable number.
    from dna_decode.glm.knob_contract import validate_comparison_mode

    validate_comparison_mode(a.natural_set, bool(a.conditional))

    # And refuse a run whose FROZEN geometry drifted from the baseline it cites. Same
    # ordering argument: this costs a file read, not a generation pass.
    if a.compare_to:
        from dna_decode.glm.knob_contract import assert_geometry_comparable

        recorded = json.loads(Path(a.compare_to).read_text(encoding="utf-8"))
        assert_geometry_comparable(recorded, geometry_from_args(a))
        print(f"geometry: comparable to {a.compare_to}")

    kinds = ["promoter", "genomic"] if a.natural_set == "both" else [a.natural_set]
    naturals = {}
    for kind in kinds:
        wins, st = build_natural(a.length, a.limit_natural, a.seed, kind=kind)
        if len(wins) < 200:
            print(f"REFUSE: {kind} corpus too small ({len(wins)}) to split into fit/eval halves")
            return 2
        from dna_decode.glm.corpus import fit_eval_halves

        h = len(wins) // 2
        fit, ev = fit_eval_halves(wins)
        naturals[kind] = {"fit": fit, "eval": ev, "stats": st,
                          "eval_windows": wins[h:]}
        print(f"natural[{kind}]: {len(fit) + len(ev)} windows of {a.length}bp -> fit {h} / eval {len(ev)}"
              + (f"  ({st.n_cds_total} CDS, {st.n_overlapping_cds} overlapping)"
                 if kind == "promoter" else ""))

    dev = device_report()
    print(f"device: {dev['device_name'] or 'cpu'}  dtype={dev['recommended_dtype']}"
          + (f"  [{dev['notes'][0]}]" if dev["notes"] else ""))

    uni = UniformGenerator(seed=a.seed).generate(a.n, a.length)

    hf_entry = None
    if a.model:
        model_id = None if a.model == "auto" else a.model
        from dna_decode.glm.generate import DEFAULT_MODEL, HFGenerator
        gen = HFGenerator(model_id or DEFAULT_MODEL, cache_dir=a.cache_dir, seed=a.seed,
                          device=a.device, temperature=a.temperature, top_k=a.top_k)
        from dna_decode.glm.knob_contract import requested_non_default, verify_knob_live

        knob_record = {"requested_non_default": requested_non_default(a.temperature, a.top_k),
                       "knob_wired": bool(gen._sampling_processors()),
                       "knob_verified_end_to_end": False}

        # Prompts are built HERE, before the knob check, so the verification can run on the
        # CONDITIONAL path -- the one the sweep actually uses. It used to call the unconditional
        # `generate`, which proves the knob moves output on a path no grid point takes. Both
        # paths do install the warpers (generate.py:232 and :277), so the risk was low; but the
        # whole subject of this episode is a knob that LOOKED wired and was inert, and a proxy
        # for the path under test is exactly the shape of evidence that failed before.
        # Built once and reused by the real generation below, so the genome is loaded once.
        prompts = None
        if a.conditional:
            # Prompts come from the EVAL windows, so the generated sequence sits at the SAME
            # aligned offsets as the reference it is scored against. Only the continuation is
            # scored.
            from dna_decode.glm.corpus import load_fasta as _lf, window_context
            _g = _lf(REFSEQ / "genome.fna")
            prompts = []
            for w in naturals[kinds[0]]["eval_windows"]:
                c = window_context(_g, w, a.context)
                if c:
                    prompts.append(c)
                if len(prompts) >= a.n:
                    break
            print(f"conditional: {len(prompts)} real prompts of {a.context} bp "
                  f"(continuation ONLY is scored; no prompt leakage)")

        if a.verify_knob and knob_record["requested_non_default"]:
            # The check the 2026-10-19 sweep lacked: three grid points produced three
            # identical numbers because the knob was inert, and the invariance was reported
            # as a finding. Refuse to sweep through an unverified knob.
            print(f"loading {gen.model_id} to verify the sampling knob moves the OUTPUT "
                  f"before sweeping ...", flush=True)
            # try/finally: an exception inside the probe must not leave the generator sitting at
            # the probe's temperature, which would silently mis-sample the run that follows.
            try:
                v = verify_knob_live(make_knob_probe(gen, prompts, a.length))
            finally:
                gen.temperature = a.temperature
            knob_record.update(verification=v, knob_verified_end_to_end=bool(v["live"]),
                               verified_on=("conditional" if prompts else "unconditional"))
            print(f"  knob verification: {v['verdict']} "
                  f"({v['identical']}/{v['n']} identical)", flush=True)
            if not v["live"]:
                print("REFUSING to sweep: the swept variable does not move the output, so "
                      "every grid point would measure the same sampling and any difference "
                      "reported would be noise.", file=sys.stderr)
                return 5
        try:
            # The knob gate above may already have loaded it; saying "loading" twice implies
            # a second load that does not happen.
            if getattr(gen, "_model", None) is None:
                print(f"loading {gen.model_id} (dtype {gen.dtype}) ...")
            if a.conditional:
                # Prompts were built above the knob check (so the verification could use this
                # same path) and are reused here -- the genome is loaded ONCE.
                seqs = gen.generate_conditional(prompts, a.length)
            else:
                seqs = gen.generate(a.n, a.length)
            if gen.fallback_reason:
                print(f"  NOTE: {gen.fallback_reason}")
            hf_entry = (gen.name, seqs)
        except Exception as e:  # noqa: BLE001 - a load/HW failure must not void the baseline result
            print(f"HF generator UNAVAILABLE ({type(e).__name__}: {str(e)[:160]})")
            print("  -> baselines still scored; the artifact records the generator as unavailable")
            hf_entry = (f"hf:{(model_id or DEFAULT_MODEL).split('/')[-1]}", [])

    results = []
    for kind in kinds:
      fit_set, eval_set = naturals[kind]["fit"], naturals[kind]["eval"]
      # The null is re-fitted PER NATURAL SET: a Markov chain fitted on promoters is not the right null
      # for genomic windows (their GC alone differs 0.4528 vs 0.5086).
      markov = MarkovGenerator(order=a.markov_order, seed=a.seed).fit(fit_set)
      null_seqs = markov.generate(a.n, a.length)
      candidates = [(markov.name, null_seqs), ("uniform", uni)]
      if hf_entry:
          candidates.append(hf_entry)
      for mode in a.modes.split(","):
        mode = mode.strip()
        print(f"\n{'=' * 78}\nNATURAL={kind}   FEATURE MODE: {mode}"
              + (f" (bins={a.bins})" if mode in ("positional", "both") else "")
              + f"\n{'=' * 78}")
        for name, seqs in candidates:
            r = run_falsifier(
                eval_set, seqs, candidate_name=name,
                null_sequences=None if name == markov.name else null_seqs,
                uniform_sequences=uni, k=a.k, length=a.length, seed=a.seed,
                mode=mode, bins=a.bins,
            )
            d = r.as_dict(); d["natural_set"] = kind; results.append(d)
            print(f"\n[{name}]  verdict={r.verdict}")
            print(f"   distinguishability (lower=better): {r.distinguishability}"
                  f"   null: {r.null_distinguishability}")
            print(f"   controls: self-max {r.control_self_max} {r.control_self_values}"
                  f"   uniform {r.control_uniform}")
            print(f"   GC natural {r.gc_natural} vs candidate {r.gc_candidate}")
            print(f"   {r.reason}")

    today = date.today().isoformat()
    art = {
        "record": "glm-generator-falsifier-v1",
        "date": today,
        "reference": "GCF_000005845.2 (E. coli K-12 MG1655), locally cached (NOT tracked in git)",
        # The nine frozen-geometry fields come from the SAME expression the --compare-to guard
        # reads, so the two descriptions of one run cannot drift apart again (see
        # `geometry_from_args`). Everything below is per-run output, not geometry.
        **geometry_from_args(a),
        # PROVENANCE: this artifact recorded no seed AT ALL until 2026-11-21, so a committed result
        # could not be traced to the draw that produced it. Additive and deliberately NOT a geometry
        # field -- no entry of COMPARISON_GEOMETRY_FIELDS is a seed, so `--compare-to` is unaffected.
        #
        # ONE seed drives all five roles below.
        #
        # An earlier version of this comment cited `all-coupled sd 0.0140 against generation-only
        # 0.0280` and attributed the gap to the Markov null being fitted on the fit half of the very
        # corpus draw it is scored against. BOTH halves of that are withdrawn. The figures came from a
        # probe running at n_generated=400 while this runner's committed geometry is n=120, so they
        # described a geometry no artifact uses; and the fit-set mechanism could not have been measured
        # by that probe at all, because every one of its arms drew fit and eval from a single corpus
        # call, leaving the matching constant across arms and therefore unable to explain a difference
        # between them. The arm built to test it directly returns an sd ratio of 1.129 with a bootstrap
        # interval of [0.891, 1.445], which contains 1.
        #
        # What the across-seed spread actually IS: the sampling noise of a held-out AUROC at this
        # test-set size and class balance. The Mann-Whitney null sd for about 300 natural against 60
        # generated in the held-out half is 0.0409; the committed grid's own ten nulls spread 0.0410,
        # and the all-coupled arm measured at 240 replicates gives 0.0435 (kmer) and 0.0448 (both),
        # both intervals containing the benchmark. So no seed ROLE carries the spread and splitting the
        # seed has nothing to remove -- the only lever on the error bar is test-set size.
        #
        # On pairing, stated at the strength it was measured and no higher: across all ten
        # checkpoint/null pairs in the committed grid, a resolvable paired-delta benefit appears in 2 of
        # 4 model-vs-NULL pairs in `both` mode (median ratio 0.684) and 0 of 4 in kmer, against 1 of 6
        # checkpoint-vs-checkpoint pairs (median ratio 0.928). The strongest single case is zero-shot
        # minus null in `both`: paired sd 0.0275 against an unpaired 0.0643, correlation +0.832 with CI
        # [0.424, 0.959]. That is suggestive that pairing helps more against the null than between
        # checkpoints, and it is NOT established as a clean rule -- four and six overlapping pairs drawn
        # from five shared series cannot carry one. See wiki/glm_grid_reconcile_2026-11-21.json,
        # wiki/glm_spread_vs_sampling_noise_2026-11-21.json and
        # wiki/glm_seed_role_variance_v2_2026-11-21.json; the v1 artifacts carry supersession stamps.
        "seeds": {"master": a.seed,
                  "drives": ["corpus shuffle + limit truncation", "generator RNG", "Markov null RNG",
                             "uniform null RNG", "discriminator fit/eval split"],
                  "coupled": True,
                  "note": "the across-seed spread is the SAMPLING NOISE of a held-out AUROC at this "
                          "test-set size and class balance (Mann-Whitney benchmark 0.0409; grid 0.0410; "
                          "all-coupled arm 0.0435/0.0448 at 240 replicates), so no seed role carries it "
                          "and the only lever on the error bar is test-set size. Pairing shows a "
                          "resolvable benefit in 2 of 4 model-vs-null pairs in mode both (median ratio "
                          "0.684) against 1 of 6 checkpoint pairs (0.928) -- suggestive, NOT a "
                          "established rule, and not to be quoted as a general benefit of pairing.",
                  "withdrawn_claim": "an earlier artifact attributed the coupling's effect to the null "
                                     "being fitted on the fit half of the corpus draw it is scored "
                                     "against; the arm that tests that directly returns an interval "
                                     "containing 1. See the v1 supersession stamps."},
        "natural_sets": {k: naturals[k]["stats"].as_dict() for k in kinds},
        "fit_eval_split": {"n_fit": len(naturals[kinds[0]]["fit"]),
                           "n_eval": len(naturals[kinds[0]]["eval"]),
                           "note": "Markov null fitted on FIT only; all discrimination on EVAL only, so "
                                   "no generator is scored against sequences it was trained on -- up to "
                                   "duplicate upstream windows, which are real: the E. coli promoter "
                                   "corpus holds 7 sequences repeated across 24 of 4340 windows "
                                   "(operons, duplicated genes, repeated elements), so at the committed "
                                   "geometry 0 to 3 of 600 eval windows also appear in the fit half. "
                                   "Measured, not assumed; too small to move a score and stated because "
                                   "the unqualified claim is not exactly true"},
        "device": dev,
        "knob_contract": (knob_record if a.model else
                          {"requested_non_default": False, "knob_wired": False,
                           "knob_verified_end_to_end": False,
                           "note": "no HF generator in this run"}),
        "sampling": {"temperature": a.temperature, "top_k": a.top_k,
                     "is_recorded_default": (a.temperature == 1.0 and a.top_k == 0),
                     "note": "stamped so a sweep row is attributable; an unstamped artifact cannot "
                             "be distinguished from the default-sampling run already on record"},
        # `conditional` / `context_bp` are supplied by the `geometry_from_args(a)` spread above.
        # They used to be restated here too; a later duplicate key silently WINS over the spread,
        # so leaving them would have re-created the second expression this step exists to remove --
        # and it would have done so invisibly, since both expressions agreed on the value.
        "results": results,
        "honest_limits": [
            "An upstream window is a PROXY for a promoter, not a mapped transcription start site.",
            "One organism (E. coli K-12). Says nothing about transfer to another host.",
            "A k-mer discriminator measures LOCAL composition; the published audits that separate "
            "generated from natural sequence name LONG-RANGE organisation as the failure, which k-mer "
            "features cannot see. So a low distinguishability here is necessary, NOT sufficient.",
            "dna_decode/constraints is CDS-shaped and INAPPLICABLE to non-coding sequence; it is not run "
            "and no vacuous pass-rate is reported.",
        ],
    }
    out = Path(a.out) if a.out else WIKI / f"glm_generator_falsifier_{today}.json"
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
