"""Pins the build-time knobs for the `finetuned-promoteratlas-encoder-150bp` Kaggle kernel.

WHY THESE EXIST. A Kaggle SCRIPT kernel receives NO environment variables. The metadata
`scripts/kaggle_push_poll.py` writes carries `id`, `title`, `code_file`, `language`, `kernel_type`,
`is_private`, `enable_gpu`, `enable_internet`, `machine_shape` and an optional `dataset_sources` --
there is no env field anywhere in it. So the kernel's original `os.environ.get("FT_SMOKE_ONLY")`,
`FT_WALL_HOURS` and `FT_CKPT_IN` reads could only ever return their own defaults on Kaggle, which
left smoke-first sizing and cross-session chaining UNREACHABLE while looking configurable. The
builder now bakes those values into the emitted module constants.

The defect class is the dangerous one: an unreachable knob fails by PASSING. The kernel runs, prints
a projection, and does the non-smoke thing -- there is no error to notice. So these tests assert on
the EMITTED ARTIFACT, never on the parser, and every positive assertion is paired with a negative
one so it cannot be satisfied by both states at once.

Pure static/structural checks: no torch, no network, no GPU, no Kaggle, no D: drive.
"""
from __future__ import annotations

import ast
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "kaggle" / "build_promoteratlas_finetune_kernel.py"
PUSHER = ROOT / "scripts" / "kaggle_push_poll.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def builder():
    return _load(BUILDER, "_bpfk")


@pytest.fixture(scope="module")
def payload(builder):
    """The gzipped source payload, built ONCE -- it is ~190 KB and identical across flag choices."""
    return builder.build_payload()


def _emit(builder, payload, **kw) -> str:
    sources_b64, mods, bar_b64 = payload
    return builder.emit(sources_b64, mods, bar_b64, **kw)


# ---------------------------------------------------------------------------------------------------
# the root finding: there is nowhere to put an environment variable
# ---------------------------------------------------------------------------------------------------
def test_the_push_metadata_has_no_environment_field():
    """The reason baking is necessary at all. If this ever fails, re-examine the baking.

    Asserted against the pusher's source rather than remembered, because the whole of Step 2 rests
    on it. A kernel-metadata schema that gained an env field would make the env reads reachable and
    this design re-decidable.
    """
    src = PUSHER.read_text(encoding="utf-8")
    tree = ast.parse(src)
    # Find the metadata dict literal and read its keys off the AST, not off a regex over prose.
    keys: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict):
            got = {k.value for k in node.keys
                   if isinstance(k, ast.Constant) and isinstance(k.value, str)}
            if "kernel_type" in got and "code_file" in got:
                keys |= got
    assert keys, "could not locate the kernel-metadata dict in the pusher"
    assert "kernel_type" in keys and "enable_gpu" in keys, keys
    forbidden = {k for k in keys if "env" in k.lower() or "variable" in k.lower()}
    assert not forbidden, (
        "the kernel metadata appears to carry an environment field %r. If Kaggle now accepts "
        "environment variables for script kernels, the baked knobs in "
        "build_promoteratlas_finetune_kernel.py are no longer the only reachable mechanism and "
        "that decision should be revisited rather than inherited." % forbidden)


# ---------------------------------------------------------------------------------------------------
# --smoke
# ---------------------------------------------------------------------------------------------------
def test_smoke_build_bakes_a_truthy_smoke_only(builder, payload):
    body = _emit(builder, payload, smoke=True)
    assert "BAKED_SMOKE_ONLY = True" in body


def test_non_smoke_build_does_not_bake_a_truthy_smoke_only(builder, payload):
    """Non-vacuity for the test above: the two builds must differ, or neither assertion means
    anything. `None` (not `False`) is what keeps the env fallback live for a local dry-run."""
    body = _emit(builder, payload, smoke=False)
    assert "BAKED_SMOKE_ONLY = None" in body
    assert "BAKED_SMOKE_ONLY = True" not in body


def test_the_baked_literal_wins_over_the_env_read(builder, payload):
    """Both must be present, in that order of precedence. A baked knob that an unsettable env var
    could override would be theatre; an env read deleted outright would break local dry-runs."""
    body = _emit(builder, payload, smoke=True)
    assert "bool(BAKED_SMOKE_ONLY) if BAKED_SMOKE_ONLY is not None" in body
    assert 'os.environ.get("FT_SMOKE_ONLY", "") == "1"' in body
    baked_at = body.index("BAKED_SMOKE_ONLY = True")
    used_at = body.index("SMOKE_ONLY = (bool(BAKED_SMOKE_ONLY)")
    assert baked_at < used_at, "the baked constant must be defined before it is read"


# ---------------------------------------------------------------------------------------------------
# --wall-hours
# ---------------------------------------------------------------------------------------------------
def test_wall_hours_is_baked_as_a_literal(builder, payload):
    body = _emit(builder, payload, wall_hours=3.5)
    assert "BAKED_WALL_HOURS = 3.5" in body


def test_unset_wall_hours_leaves_the_derived_default(builder, payload):
    """Non-vacuity, and the derivation must travel with the number.

    10.5 h is 12.0 h of session cap minus 1.5 h for the clone, the GEO download and the artifact
    write. A bare 10.5 in the emitted kernel would be a number nobody downstream could check.
    """
    body = _emit(builder, payload)
    assert "BAKED_WALL_HOURS = None" in body
    assert "BAKED_WALL_HOURS = 3.5" not in body
    assert 'os.environ.get("FT_WALL_HOURS", "10.5")' in body
    assert "12.0 - 10.5 = 1.5 h" in body, "the 10.5 h ceiling must ship with its derivation"


def test_the_ceiling_stays_under_the_session_cap(builder, payload):
    """A baked ceiling at or above the 12 h cap would guarantee the artifact is never flushed."""
    body = _emit(builder, payload, wall_hours=11.9)
    assert "BAKED_WALL_HOURS = 11.9" in body
    assert "BOUNDED_BY_WALL_CLOCK" in body, (
        "the kernel must still carry the bounded-vs-complete distinction; a wall-clock stop is "
        "gradeable but is never COMPLETE")


# ---------------------------------------------------------------------------------------------------
# --ckpt-dataset
# ---------------------------------------------------------------------------------------------------
def test_ckpt_dataset_bakes_the_kaggle_input_mount_path(builder, payload):
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    assert "BAKED_CKPT_IN = '/kaggle/input/ft-ckpt-round1'" in body


def test_unset_ckpt_dataset_leaves_the_bare_input_root(builder, payload):
    body = _emit(builder, payload)
    assert "BAKED_CKPT_IN = None" in body
    assert "/kaggle/input/ft-ckpt-round1" not in body
    assert 'os.environ.get("FT_CKPT_IN", "/kaggle/input")' in body


# ---------------------------------------------------------------------------------------------------
# the emitted artifact stays valid, and the closure is flag-independent
# ---------------------------------------------------------------------------------------------------
#: All four knobs, every subset. This was three flags / 8 combinations until --horizon-epochs
#: landed; the old matrix still PASSED afterwards because no build in it baked a horizon, which made
#: it blind to the knob it needed to cover. Enumerated as a power set so adding a fifth flag cannot
#: leave a stale hand-written list behind.
_FLAG_VALUES = {
    "smoke": True,
    "wall_hours": 3.5,
    "ckpt_dataset": "ft-ckpt-round1",
    "horizon_epochs": 38,
}
FLAG_COMBOS = [
    {k: _FLAG_VALUES[k] for k in subset}
    for subset in __import__("itertools").chain.from_iterable(
        __import__("itertools").combinations(sorted(_FLAG_VALUES), r)
        for r in range(len(_FLAG_VALUES) + 1))
]


@pytest.mark.parametrize("kw", FLAG_COMBOS, ids=lambda k: "-".join(sorted(k)) or "none")
def test_emitted_body_parses_under_every_flag_combination(builder, payload, kw):
    """The builder validates its ARTIFACT, not itself -- a collapsed escape has shipped twice. Each
    new format field is a new way to break the emitted body, so every combination is parsed."""
    ast.parse(_emit(builder, payload, **kw))


def _embedded_of(body: str) -> dict:
    """Read the emitted kernel's EMBEDDED map off the AST, not off a regex over the text."""
    embedded = None
    for node in ast.walk(ast.parse(body)):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "EMBEDDED" for t in node.targets):
            embedded = ast.literal_eval(node.value)
    assert embedded, "EMBEDDED not found in the emitted kernel"
    return embedded


@pytest.mark.parametrize("kw", FLAG_COMBOS, ids=lambda k: "-".join(sorted(k)) or "none")
def test_the_embed_set_is_the_payloads_own_map(builder, payload, kw):
    """`emit`'s contract is to render the map it was HANDED.

    The reference here is deliberately `build_payload()`'s map and NOT `closure()`: the payload adds
    the three `__init__.py` package files via `_package_inits`, so comparing against the bare closure
    understates the embed set by 3 and fails for a reason that has nothing to do with the flags.
    """
    embedded = _embedded_of(_emit(builder, payload, **kw))
    expected = {m: p.relative_to(builder.REPO).as_posix()
                for m, p in sorted(payload[1].items())}
    assert embedded == expected
    assert "scripts.finetune_arm_run" in embedded
    assert "scripts.positional_arm_score" in embedded, "transitive import dropped from the closure"
    # The package inits are part of the contract: without them the unpacked payload is not importable.
    assert "scripts" in embedded and embedded["scripts"].endswith("__init__.py")


def test_the_embed_set_is_IDENTICAL_across_every_flag_combination(builder, payload):
    """The invariant the parametrized test above cannot state on its own.

    Each case checks its own build against a fixed reference, which would still pass if every build
    were wrong in the same way. This compares the builds to EACH OTHER: the flags configure RUNTIME
    behaviour and must never move which modules ship, because that would change what the kernel can
    import for a reason no flag is for.
    """
    maps = {("-".join(sorted(kw)) or "none"): _embedded_of(_emit(builder, payload, **kw))
            for kw in FLAG_COMBOS}
    first_id, first = next(iter(maps.items()))
    for cid, m in maps.items():
        assert m == first, "flag combo %r moved the embed set relative to %r" % (cid, first_id)
    assert len(maps) == len(FLAG_COMBOS) >= 8, "the combination matrix shrank"


def test_only_the_knob_constants_differ_between_a_smoke_and_a_full_build(builder, payload):
    """The sharpest form of the claim: --smoke must change the three baked lines and NOTHING else.

    A flag that also perturbed the payload, the bar or the phase order would make the smoke's
    measured rate a measurement of a different program than the full run it is used to size.
    """
    a = _emit(builder, payload).splitlines()
    b = _emit(builder, payload, smoke=True, wall_hours=3.5,
              ckpt_dataset="ft-ckpt-round1", horizon_epochs=38).splitlines()
    assert len(a) == len(b), "a flag changed the kernel's line count"
    differing = [x for x, y in zip(a, b) if x != y]
    assert all(d.startswith("BAKED_") for d in differing), differing
    # Four knobs, so four lines. This read 3 before --horizon-epochs existed and still passed,
    # because neither build baked a horizon -- the count is the whole point of the assertion, so it
    # has to track the number of knobs rather than be left at a number that happens to hold.
    assert len(differing) == 4, differing
    assert {d.split(" =")[0] for d in differing} == {
        "BAKED_SMOKE_ONLY", "BAKED_WALL_HOURS", "BAKED_CKPT_IN", "BAKED_HORIZON_EPOCHS"}, differing


# ---------------------------------------------------------------------------------------------------
# the guards the flags must not weaken
# ---------------------------------------------------------------------------------------------------
def test_the_smoke_never_writes_into_the_real_checkpoint_directory(builder, payload):
    """The recorded hazard. The resume contract reads any existing per-seed checkpoint and never
    recomputes it, so one smoke epoch written into the real dir becomes the real run's trajectory
    prefix -- silently, and in a way the result would not disclose."""
    body = _emit(builder, payload, smoke=True)
    assert 'SMOKE_CKPT = WORK / "smoke_checkpoints_THROWAWAY"' in body
    assert 'CKPT = WORK / "finetune_checkpoints"' in body
    assert "ckpt_dir=SMOKE_CKPT" in body, "the smoke epoch must train into the throwaway dir"
    assert "REFUSED: the smoke run touched the REAL checkpoint dir" in body, (
        "the before/after guard on the real checkpoint dir must survive the flag work")


def test_a_smoke_build_scores_nothing_and_writes_no_run_artifact(builder, payload):
    body = _emit(builder, payload, smoke=True)
    assert '"scored": False' in body
    i_smoke = body.index("if SMOKE_ONLY:")
    i_run = body.index("from scripts.finetune_arm_run import run")
    assert i_smoke < i_run, "the smoke branch must return before the scoring run is imported"


def test_the_smoke_artifact_reports_the_device_it_measured_on(builder, payload):
    """A seconds-per-epoch figure is a property of the hardware that produced it. The T4 pin can
    still yield a P100, so the rate must never be recorded without the device beside it."""
    body = _emit(builder, payload, smoke=True)
    assert '"device": dev' in body
    assert "t4_pin_honoured" in body
    assert "rate_is_for_this_device_only" in body
    assert '"fits_one_session"' in body


def test_the_kernel_still_calls_run_not_main_and_refuses_a_non_kaggle_bar(builder, payload):
    """Two inherited invariants. `main()` holds the local route gate and correctly refuses a
    non-local route, so the kernel calls `run()`; and the kernel refuses a bar that is not routed to
    kaggle_t4, so a stale bar cannot be run on the wrong hardware."""
    body = _emit(builder, payload)
    assert "from scripts.finetune_arm_run import run" in body
    assert "from scripts.finetune_arm_run import main" not in body
    assert 'REFUSED: embedded bar routes to' in body


def test_the_embedded_bar_is_the_frozen_bar_byte_for_byte(builder, payload):
    """The bar travels INSIDE the kernel, so a drifted copy would be scored against silently."""
    import base64
    import gzip
    _, _, bar_b64 = payload
    shipped = gzip.decompress(base64.b64decode(bar_b64))
    on_disk = (builder.REPO / builder.BAR_REL).read_bytes()
    assert shipped == on_disk
    bar = json.loads(shipped.decode("utf-8"))
    assert bar["route"]["chosen"] == "kaggle_t4"
    assert bar["horizon"]["n_seeds"] == 10
    assert bar["frozen_hyperparameters"]["n_epochs_budget"] == 7


def test_every_bar_path_the_kernel_dereferences_EXISTS_in_the_committed_bar(builder, payload):
    """The check that would have saved a Kaggle session, and the reason it did not exist.

    The first real run died 42 s in with `KeyError: 'comparator'`, because
    `measure_seconds_per_epoch` read `bar["comparator"]["split_frac"]` while the frozen bar keeps
    that value at `inherited_not_chosen.comparator.split_frac`. Every local check was STATIC --
    `--self-test` only parses the emitted body, these tests assert on its text, and the payload
    verifier deliberately skipped the timing function. So the one function nobody executed was the
    one that was wrong.

    This resolves each `bar[...]` subscript chain in the emitted kernel against the REAL committed
    bar, which makes a shape mismatch a local failure instead of a remote one.

    SCOPE, stated because it bounds the guard: this reads the EMITTED KERNEL BODY only. Paths
    dereferenced inside the embedded modules (`finetune_arm_run.run()` does its own) are invisible
    to it. That is still exactly where this bug lived.
    """
    body = _emit(builder, payload)
    bar = json.loads((builder.REPO / builder.BAR_REL).read_text(encoding="utf-8"))

    chains: set[tuple[str, ...]] = set()
    for node in ast.walk(ast.parse(body)):
        if not isinstance(node, ast.Subscript):
            continue
        keys: list[str] = []
        cur = node
        while isinstance(cur, ast.Subscript):
            if not (isinstance(cur.slice, ast.Constant) and isinstance(cur.slice.value, str)):
                keys = []
                break
            keys.append(cur.slice.value)
            cur = cur.value
        if keys and isinstance(cur, ast.Name) and cur.id == "bar":
            chains.add(tuple(reversed(keys)))

    assert chains, "found no bar[...] dereferences in the emitted kernel; the scan is vacuous"

    bad = []
    for chain in sorted(chains):
        node = bar
        for key in chain:
            if not isinstance(node, dict) or key not in node:
                bad.append(".".join(chain))
                break
            node = node[key]
    assert not bad, (
        "the emitted kernel dereferences bar path(s) that do not exist in %s: %s"
        % (builder.BAR_REL, bad))


def test_the_bar_shape_guard_runs_before_the_clone(builder, payload):
    """Order is the whole value. Checking the bar AFTER the clone and cohort build costs a session."""
    body = _emit(builder, payload)
    assert "def require_bar_paths(bar):" in body
    i_guard = body.index("require_bar_paths(bar)\n")
    i_clone = body.index('sh(["git", "clone"')
    assert i_guard < i_clone, "the bar-shape guard must run before the 25 s clone"
    assert "inherited_not_chosen" in body, (
        "the required-paths list must name the comparator path the first run got wrong")


def test_the_smoke_uses_the_same_split_frac_as_the_real_run(builder, payload):
    """A smoke that times a DIFFERENT split is not measuring the run it is used to size."""
    body = _emit(builder, payload)
    assert 'inh = bar["inherited_not_chosen"]' in body
    assert 'split_frac=float(inh["comparator"]["split_frac"])' in body
    assert 'split_frac=bar["comparator"]["split_frac"]' not in body, "the shallow path came back"


def test_the_smoke_resolves_a_builder_and_never_passes_a_bare_None(builder, payload):
    """The second Kaggle failure: `TypeError: 'NoneType' object is not callable`.

    `train_one_seed` calls `build()` unconditionally and has no default; only `run()` resolves one.
    The smoke calls `train_one_seed` directly, so passing `build=None` there is always fatal.
    """
    body = _emit(builder, payload, smoke=True)

    # Read the two call sites off the AST. A body-wide "build=None" check would be WRONG: the
    # `run(...)` call legitimately passes None, because run() resolves the builder itself. The
    # invariant is the ASYMMETRY -- None is fine for run(), fatal for train_one_seed.
    builds = {}
    for node in ast.walk(ast.parse(body)):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        if node.func.id not in ("train_one_seed", "run"):
            continue
        for kw in node.keywords:
            if kw.arg == "build":
                builds[node.func.id] = kw.value

    assert "train_one_seed" in builds, "no train_one_seed(build=...) call found; the scan is vacuous"
    ts = builds["train_one_seed"]
    assert not (isinstance(ts, ast.Constant) and ts.value is None), (
        "the smoke passes build=None into train_one_seed, which calls build() unconditionally and "
        "has no default -- this is the exact TypeError the second Kaggle run died on")
    assert isinstance(ts, ast.Name) and ts.id == "_build_trainable", ast.dump(ts)

    assert "run" in builds, "no run(build=...) call found"
    r = builds["run"]
    assert isinstance(r, ast.Constant) and r.value is None, (
        "run(build=None) is CORRECT and deliberate -- run() resolves the builder itself. If this "
        "changed, re-derive whether the smoke should still resolve its own.")

    assert "from scripts.finetune_arm_probe import _build_trainable" in body


def test_the_smoke_uses_THE_SAME_builder_the_real_run_defaults_to(builder, payload):
    """Pins the pair, not just the name.

    The smoke's measured seconds-per-epoch is what the real run's horizon is sized from, so if the two
    trained different models the measurement would not describe the run it is used to size. This reads
    the builder `run()` falls back to out of the RUNNER'S OWN SOURCE rather than restating it here, so
    the two cannot drift apart silently.
    """
    runner = (builder.REPO / "scripts" / "finetune_arm_run.py").read_text(encoding="utf-8")
    assert "from scripts.finetune_arm_probe import _build_trainable" in runner, (
        "run()'s fallback builder import moved; re-derive what the smoke should use")
    body = _emit(builder, payload, smoke=True)
    assert "from scripts.finetune_arm_probe import _build_trainable" in body, (
        "the smoke resolves a DIFFERENT builder than run() falls back to")


def test_the_builder_module_is_in_the_embedded_payload(builder, payload):
    """A resolved import is only resolvable if its module actually shipped."""
    embedded = _embedded_of(_emit(builder, payload, smoke=True))
    assert "scripts.finetune_arm_probe" in embedded, (
        "the smoke imports _build_trainable from scripts.finetune_arm_probe; that module must be in "
        "the AST closure or the import dies on Kaggle")


def test_self_test_still_passes(builder):
    assert builder.self_test() == 0


# ===================================================================================================
# --horizon-epochs: the bake, and the thing the first two Kaggle failures were really about
# ===================================================================================================
def _run_call_keywords(body: str) -> list[dict[str, str]]:
    """Keywords of every call to a BARE name `run` in the emitted kernel, as {kwarg: source}.

    Restricted to `ast.Name` deliberately. The kernel does `from scripts.finetune_arm_run import run`
    and then `run(...)`, while it ALSO calls `subprocess.run` via its `sh` helper -- an
    `ast.Attribute` whose `.attr` is likewise "run". Matching on the attribute name would pull every
    shell invocation into this check and fail for a reason that has nothing to do with the horizon.
    """
    calls = []
    for node in ast.walk(ast.parse(body)):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "run":
            calls.append({k.arg: ast.unparse(k.value) for k in node.keywords if k.arg})
    return calls


def _horizon_reaches_run(body: str) -> tuple[bool, str]:
    """Shared by the live test and its non-vacuity twin, so they cannot drift apart."""
    calls = _run_call_keywords(body)
    if not calls:
        return False, "no call to a bare run() was found in the emitted kernel"
    bad = [c for c in calls if c.get("horizon_epochs") != "HORIZON_EPOCHS"]
    if bad:
        return False, (
            "run() is called without horizon_epochs=HORIZON_EPOCHS (got %r). A baked constant that "
            "never reaches its consumer is EXACTLY the defect class of the first two Kaggle "
            "failures: both were a knob or a value that existed and was never read." % (bad,))
    return True, "run() receives horizon_epochs=HORIZON_EPOCHS at %d call site(s)" % len(calls)


def test_horizon_epochs_is_baked_as_a_literal(builder, payload):
    body = _emit(builder, payload, horizon_epochs=38)
    assert "BAKED_HORIZON_EPOCHS = 38" in body


def test_unset_horizon_leaves_planning_in_charge_and_the_refusal_live(builder, payload):
    """The default must still REFUSE, because that refusal is correct for an unchained run.

    `run()` plans from the ceiling only when `horizon_epochs is None`, and at the measured T4 rate
    `plan_run` refuses every session-sized ceiling. An un-flagged build has to keep reaching that
    refusal; silently acquiring a horizon would turn a correct refusal into a run nobody chose.
    """
    body = _emit(builder, payload)
    assert "BAKED_HORIZON_EPOCHS = None" in body
    assert "BAKED_HORIZON_EPOCHS = 38" not in body
    assert 'os.environ.get("FT_HORIZON_EPOCHS")' in body, (
        "the env read must stay live beneath the bake so a LOCAL dry-run can still override")


def test_the_baked_horizon_REACHES_THE_RUN_CALL(builder, payload):
    """The load-bearing test of Step 1, and the one the plan singled out.

    Baking a constant is not the same as wiring it. The first Kaggle failure read a real value at the
    wrong dict depth; the second passed a real parameter that was never defaulted. Both shipped a
    kernel whose constants looked right. So this asserts the keyword arrives AT THE CALL, not that
    the literal appears somewhere in 200 KB of text.
    """
    ok, why = _horizon_reaches_run(_emit(builder, payload, horizon_epochs=38))
    assert ok, why


def test_the_horizon_pass_through_check_is_NON_VACUOUS(builder, payload):
    """Reintroduce the defect in a throwaway copy; the check must fail and NAME the construct."""
    body = _emit(builder, payload, horizon_epochs=38)
    broken = body.replace("              horizon_epochs=HORIZON_EPOCHS,\n", "", 1)
    assert broken != body, "the defect injection did not apply; this test proves nothing as written"
    ok, why = _horizon_reaches_run(broken)
    assert not ok, "the check PASSED on a kernel that never hands run() a horizon: %s" % why
    assert "horizon_epochs=HORIZON_EPOCHS" in why


def test_the_pass_through_is_UNCONDITIONAL_and_only_the_constant_moves(builder, payload):
    """A deliberate divergence from the plan's own test strategy, recorded rather than hidden.

    The plan expected a horizon build and an un-flagged build to differ in TWO places: the baked
    constant and the `run(` keyword. They differ in ONE. `horizon_epochs=HORIZON_EPOCHS` is passed
    unconditionally and `HORIZON_EPOCHS` resolves to None when nothing was baked -- so there is no
    conditional call path to get wrong, and `run()`'s own `if horizon_epochs is None` branch does the
    switching. That is strictly stronger than what the plan asked for: the keyword cannot be present
    in one build and missing in the other, because it is present in both.
    """
    a = _emit(builder, payload).splitlines()
    b = _emit(builder, payload, horizon_epochs=38).splitlines()
    assert len(a) == len(b), "--horizon-epochs changed the kernel's line count"
    differing = [(x, y) for x, y in zip(a, b) if x != y]
    assert len(differing) == 1, differing
    assert differing[0] == ("BAKED_HORIZON_EPOCHS = None", "BAKED_HORIZON_EPOCHS = 38"), differing
    # And the keyword is in BOTH, which is the part that makes one differing line sufficient.
    for body in ("\n".join(a), "\n".join(b)):
        ok, why = _horizon_reaches_run(body)
        assert ok, why


def test_the_run_artifact_records_the_planning_bypass(builder, payload):
    """A reader who sees the bypass and not the note would reasonably conclude the bar was violated.

    `plan_run`'s refusal text says the run "must not start". Supplying the horizon walks past that
    sentence, so the artifact has to say so in its own words -- and under its OWN key, never folded
    into a field the scorer reads.
    """
    body = _emit(builder, payload, horizon_epochs=38)
    assert '"horizon_divergence"' in body
    assert "finetune-arm-horizon-divergence-v1" in body
    assert "this run must not start" in body, (
        "the divergence block must quote plan_run's own refusal rather than paraphrase it")
    assert "INDETERMINATE_BUDGET_EXHAUSTED" in body, (
        "the block must say that a stop INSIDE the budget phase is not gradeable")
    # Namespace separation: the note is added to a COPY beside rec, not merged into it.
    assert 'payload = dict(rec)' in body
    assert 'payload["horizon_divergence"]' in body


# ===================================================================================================
# Step 2: the zero-copy refusal in seed_checkpoints
# ===================================================================================================
def _seed_refusal_state(body: str) -> tuple[bool, str]:
    """Is the zero-copy refusal present AND fenced behind a baked mount? Shared with its twin."""
    fn = None
    for node in ast.walk(ast.parse(body)):
        if isinstance(node, ast.FunctionDef) and node.name == "seed_checkpoints":
            fn = node
    if fn is None:
        return False, "seed_checkpoints() not found in the emitted kernel"

    # The ZERO-COPY refusal specifically, not merely any SystemExit in the function. There is a
    # second, unrelated refusal here (duplicate basenames under the mount), and an earlier version
    # of this check accepted either -- which let a non-vacuity injection neuter the wrong one and
    # still pass. Require the raise to sit under a guard on `n == 0`.
    raises = []
    for node in ast.walk(fn):
        if not (isinstance(node, ast.If) and "n == 0" in ast.unparse(node.test)):
            continue
        raises += [r for r in ast.walk(node)
                   if isinstance(r, ast.Raise) and r.exc is not None
                   and ast.unparse(r.exc).startswith("SystemExit")]
    if not raises:
        return False, (
            "seed_checkpoints() carries no SystemExit refusal guarded by `n == 0`, so a baked mount "
            "that seeds ZERO files would let the session recompute every seed while looking "
            "entirely normal")

    fence = [n for n in fn.body
             if isinstance(n, ast.If) and "BAKED_CKPT_IN is None" in ast.unparse(n.test)]
    if not fence:
        return False, (
            "the refusal is not fenced behind `BAKED_CKPT_IN is None`, so the DEFAULT "
            "/kaggle/input path would refuse too -- and there a zero is the CORRECT answer")
    if not any(isinstance(n, ast.Return) for n in ast.walk(fence[0])):
        return False, (
            "the `BAKED_CKPT_IN is None` branch does not return, so it does not actually fence the "
            "refusal off from the default path")
    return True, "refusal present, and fenced behind the baked-mount branch"


def test_a_baked_mount_that_seeds_nothing_REFUSES(builder, payload):
    ok, why = _seed_refusal_state(_emit(builder, payload, ckpt_dataset="ft-ckpt-round1"))
    assert ok, why


def test_the_default_input_path_keeps_its_CORRECT_silence(builder, payload):
    """Zero means opposite things in the two cases, and only one of them is a fault.

    With nothing attached, CKPT_IN is the bare /kaggle/input and zero copied files is exactly right.
    The structure must therefore be the same in both builds -- the fence is a runtime branch on the
    baked constant, not a build-time difference -- which is also why this check runs on the UNFLAGGED
    body and still expects the fence to be found.
    """
    ok, why = _seed_refusal_state(_emit(builder, payload))
    assert ok, why
    body = _emit(builder, payload)
    assert "BAKED_CKPT_IN = None" in body
    # The silent branch is the one that returns before the unconditional print.
    assert "Nothing was attached and nothing was asked for" in body


def test_the_refusal_check_is_NON_VACUOUS_when_the_raise_is_removed(builder, payload):
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    # Anchored on the refusal's OWN first line. A bare "        raise SystemExit(" occurs SIX
    # times in the emitted kernel -- including the deeper duplicate-basename raise in this very
    # function -- so replace() hit the wrong one; and `if n == 0:` is not adjacent to the raise
    # (a `found = []` listing sits between them), so the two-line pair does not occur either.
    # Both were tried and both were wrong, which is why this anchors on unique refusal text.
    target = '        raise SystemExit(\n            \"REFUSING to run: --ckpt-dataset baked'
    assert body.count(target) == 1, "the zero-copy refusal is no longer uniquely anchorable"
    broken = body.replace(
        target, '        _ = (\n            \"REFUSING to run: --ckpt-dataset baked', 1)
    assert broken != body, "the defect injection did not apply; this test proves nothing"
    ok, why = _seed_refusal_state(broken)
    assert not ok, "the check PASSED on a seed_checkpoints with no zero-copy refusal: %s" % why
    assert "guarded by `n == 0`" in why


def test_the_refusal_check_is_NON_VACUOUS_when_the_fence_is_removed(builder, payload):
    """Removing the fence makes the DEFAULT path refuse -- a different bug, equally real."""
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    broken = body.replace("    if BAKED_CKPT_IN is None:", "    if False:", 1)
    assert broken != body, "the defect injection did not apply; this test proves nothing"
    ok, why = _seed_refusal_state(broken)
    assert not ok, "the check PASSED with the baked-mount fence removed: %s" % why
    assert "BAKED_CKPT_IN is None" in why


def test_the_refusal_names_the_mount_and_both_causes_of_a_zero(builder, payload):
    """A refusal a reader cannot act on is only half a guard.

    Two causes remain after the two mount findings -- the dataset was never attached, or it is
    genuinely empty -- and they need OPPOSITE fixes. The message has to pose the distinction, name
    the pattern it matched on, and LIST what is actually present, or the reader cannot tell which.
    """
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    msgs = [m for m in _systemexit_messages(body) if "RECOMPUTES EVERY SEED" in m]
    assert len(msgs) == 1, "expected exactly one zero-copy refusal, got %d" % len(msgs)
    m = msgs[0]
    assert "never attached" in m, "the refusal must pose the not-attached cause"
    assert "genuinely empty" in m, "the refusal must pose the empty-dataset cause"
    assert "smoke_*" in m, (
        "the refusal must say smoke output is excluded BY DESIGN, or a reader with a smoke-only "
        "search root chases the wrong cause")
    assert "Present under" in m, "the refusal must list what IS there"
    # The pattern is interpolated, not restated -- a restated pattern would drift from the real one.
    assert "CKPT_NAME_RE.pattern" in body


def test_the_early_return_that_skipped_the_guard_is_GONE(builder, payload):
    """The regression this is really guarding.

    The original `if not CKPT_IN.is_dir(): return 0` ran BEFORE the count, so a baked mount that
    never materialised -- the single most likely form of this failure -- took an early return and
    skipped the refusal entirely. The directory test still exists, but now as a condition on the copy
    loop rather than a return.
    """
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    fn_src = body[body.index("def seed_checkpoints()"):]
    fn_src = fn_src[:fn_src.index("\ndef ")]
    assert "if not CKPT_IN.is_dir():\n        return 0" not in fn_src, (
        "the early return is back; a baked mount that does not exist would skip the refusal")
    assert "if CKPT_IN.is_dir():" in fn_src, "the directory test itself must still be there"


# ===================================================================================================
# Step 5's measured correction: match on the FILENAME, because the mount flattens
# ===================================================================================================
def _systemexit_messages(body: str) -> list[str]:
    """Every `raise SystemExit(...)` message in the emitted kernel, as JOINED text.

    Substring checks against source are fragile for a reason that bit three tests at once: these
    messages are wrapped across source lines, so `"... or it is genuinely " "empty. ..."` contains
    no contiguous "genuinely empty". Python concatenates adjacent string literals at parse time, so
    the AST already holds the joined string -- reading it here makes the assertions independent of
    where the line breaks happen to fall.
    """
    out = []
    for node in ast.walk(ast.parse(body)):
        if not (isinstance(node, ast.Raise) and node.exc is not None):
            continue
        if not (isinstance(node.exc, ast.Call) and isinstance(node.exc.func, ast.Name)
                and node.exc.func.id == "SystemExit"):
            continue
        for arg in node.exc.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                out.append(arg.value)
            elif isinstance(arg, ast.BinOp) and isinstance(arg.left, ast.Constant) \
                    and isinstance(arg.left.value, str):
                out.append(arg.left.value)
    return out


def _runner_ckpt_name_re():
    """The ONE definition, read off the runner's AST so no torch import is ever required."""
    src = (ROOT / "scripts" / "finetune_arm_run.py").read_text(encoding="utf-8")
    for node in ast.walk(ast.parse(src)):
        if (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "CKPT_NAME_RE" for t in node.targets)):
            import re as _re
            return _re.compile(ast.literal_eval(node.value.args[0]))
    raise AssertionError("CKPT_NAME_RE not found in scripts/finetune_arm_run.py")


def _seed_match_rule(body: str) -> tuple[bool, str]:
    """Does the kernel select checkpoints by NAME rather than by parent directory? Shared w/ twin."""
    if 'p.parent.name == "finetune_checkpoints"' in body:
        return False, (
            "the kernel still matches on `p.parent.name == \"finetune_checkpoints\"`. Kaggle "
            "extracts an uploaded directory archive and FLATTENS it, so the parent under the mount "
            "is the dataset SLUG and this rule seeds zero files from a correct dataset.")
    if "CKPT_NAME_RE.match(p.name)" not in body:
        return False, "the kernel does not select checkpoint files by CKPT_NAME_RE.match(p.name)"
    return True, "checkpoints are selected by filename"


def test_checkpoints_are_matched_by_FILENAME_not_by_parent_directory(builder, payload):
    """The measured correction. See wiki/ckpt_dataset_roundtrip_2026-11-15.json.

    Uploading `finetune_checkpoints/seed_99_curve.json` as a Kaggle dataset stores exactly
    `seed_99_curve.json` at the dataset root -- the directory level does not survive. A parent-name
    rule therefore cannot ever seed a chained session.
    """
    ok, why = _seed_match_rule(_emit(builder, payload, ckpt_dataset="ft-ckpt-round1"))
    assert ok, why


def test_the_filename_rule_check_is_NON_VACUOUS(builder, payload):
    """Put the old parent-name rule back in a throwaway copy; the check must name it."""
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    good_rule = ("if not p.is_file() or not CKPT_NAME_RE.match(p.name):"
                 "\n                continue")
    old_rule = ('if not (p.is_file() and p.parent.name == "finetune_checkpoints"):'
                "\n                continue")
    broken = body.replace(good_rule, old_rule, 1)
    assert broken != body, "the defect injection did not apply; this test proves nothing"
    ok, why = _seed_match_rule(broken)
    assert not ok, "the check PASSED on a kernel using the falsified parent-name rule: %s" % why
    assert "FLATTENS" in why


def test_the_pattern_accepts_real_curves_and_REFUSES_smoke_output(builder, payload):
    """The safety property a filename rule must STATE, which the old parent-name rule had by luck.

    The Kaggle rate measurement trains a throwaway epoch at seed 0 with the real hyperparameters,
    split and builder, so `smoke_0_curve.json` is indistinguishable in SHAPE from
    `seed_0_curve.json`. `advance()` READS whatever curve already reaches the requested depth, so
    accepting smoke output would make a discarded measurement epoch the arm's trajectory prefix.
    The pattern is EVALUATED here, not described -- and read from its single source.
    """
    rx = _runner_ckpt_name_re()
    for good in ("seed_0.pt", "seed_0_curve.json", "seed_12.pt", "null_2_curve.json", "null_0.pt"):
        assert rx.match(good), "the pattern rejects a real checkpoint name %r" % good
    for bad in ("smoke_0.pt", "smoke_0_curve.json", "seed_0.pt.tmp", "seed_0_curve.json.tmp",
                "seed_.pt", "notes.txt", "finetune_checkpoints.zip", "seed_0_curve.json.bak",
                "Seed_0.pt", "seed_0.PT"):
        assert not rx.match(bad), "the pattern accepts %r, which it must not" % bad


def test_duplicate_basenames_under_the_mount_REFUSE(builder, payload):
    """The copy is flat, so two same-named files would resolve to one curve for that seed.

    This matters more now than when it was written: the search root falls back to /kaggle/input and
    recurses, so it can legitimately see several attached datasets at once.
    """
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    msgs = [m for m in _systemexit_messages(body) if "both named" in m]
    assert len(msgs) == 1, "expected exactly one duplicate-basename refusal"
    assert "FLAT" in msgs[0], "the refusal must say WHY a duplicate is fatal"
    assert "taken[p.name]" in body, "the refusal must report WHICH two paths collided"


def test_a_missing_baked_mount_FALLS_BACK_to_the_input_root(builder, payload):
    """Measured twice on 2026-11-15, and the second finding replaced this test's original subject.

    An earlier version pinned a "list the parent" message, which existed only because the old code
    could search nowhere but the baked path. Probe 2 then showed the baked path is not where Kaggle
    mounts anything: a dataset attached via `dataset_sources`, reporting status=ready, did not appear
    at /kaggle/input/<slug> and /kaggle/input held only `datasets/`. So the kernel stopped predicting
    the path -- the baked value is a HINT and the fallback search root is /kaggle/input, recursed.
    That SUBSUMES listing the parent, and it is the stronger property, so it is what gets pinned.

    The fallback must widen WHERE we look and never WHETHER a zero is a fault.
    """
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    assert "INPUT_ROOT = Path(" in body, "the kernel must define an input-root constant"
    assert "elif BAKED_CKPT_IN is not None and INPUT_ROOT.is_dir():" in body, (
        "the fallback must be conditioned on a dataset having been baked; falling back with no "
        "dataset attached would search the whole input tree for no reason")
    assert "does not exist; searching" in body, (
        "the fallback must ANNOUNCE itself -- a silent change of search root is how a later reader "
        "mis-attributes a seeded count")
    # The refusal still fires on zero: the fallback is not an escape from the guard.
    assert [m for m in _systemexit_messages(body) if "RECOMPUTES EVERY SEED" in m], (
        "the zero refusal disappeared along with the parent listing")
    # And the no-dataset case must still be silent, which is what keeps the default correct.
    plain = _emit(builder, payload)
    assert "Nothing was attached and nothing was asked for" in plain


def test_the_checkpoint_name_pattern_has_exactly_ONE_definition(builder, payload):
    """The duplication guard, and the failure it prevents is a silent one.

    Step 7's uploader must stage exactly the files the kernel will accept. Two copies of that
    selection rule on opposite ends of a round trip fail by PASSING: the uploader stages a file the
    kernel ignores, seeding reports zero, and the refusal blames the mount -- which is the wrong
    diagnosis, and precisely the diagnosis that cost two probe pushes to correct.
    """
    runner_src = (ROOT / "scripts" / "finetune_arm_run.py").read_text(encoding="utf-8")
    assert runner_src.count("CKPT_NAME_RE = re.compile(") == 1, (
        "scripts/finetune_arm_run.py must define the pattern exactly once -- it owns the names via "
        "_ckpt_path and _curve_path")
    builder_src = BUILDER.read_text(encoding="utf-8")
    assert "CKPT_NAME_RE = re.compile(" not in builder_src, (
        "the kernel template defines its own copy of the pattern again; it must import the one in "
        "scripts/finetune_arm_run.py")
    body = _emit(builder, payload, ckpt_dataset="ft-ckpt-round1")
    assert "from scripts.finetune_arm_run import CKPT_NAME_RE" in body, (
        "the emitted kernel must import the single definition rather than carry a copy")

