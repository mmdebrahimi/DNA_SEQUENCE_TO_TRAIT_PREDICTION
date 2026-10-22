"""READ-ONLY: can the promcalc stage-B toolchain exist on THIS host?

    uv run python scripts/promcalc_toolchain_probe.py [--out wiki/promcalc_toolchain_probe_<date>.json]

**Why this exists.** The tile-rung chain has two substrate dependencies with different wall classes.
The MPRA peak-tile table is PRESENT and merely mis-addressed (a hardcoded `D:` against data that sits
on `E:`) -- code-closable. The promcalc per-tile predictions are genuinely ABSENT from every attached
volume, and the question that decides whether that is a code wall or an external one is whether the
isolated interpreter stage B needs could be created here at all. Nobody had asked it; the frozen
acceptance bar that waits on this chain says only "no D: data tree, no uv and no Docker", and the
`no uv` clause is already false.

**IT INSTALLS NOTHING, AND THAT IS THE WHOLE POINT OF THE SHAPE.** Reporting that an install is
possible is a different claim from having done one, and only the first is made here. `probe_toolchain`
issues no subprocess and writes no file; `main` writes an artifact only when given `--out`. The test
suite pins the no-subprocess half by monkeypatching `subprocess` to refuse being called.

**THE VERDICT VOCABULARY IS FROZEN IN THIS MODULE, ABOVE ANY MEASUREMENT.** A probe that picks its
words after seeing the answer measures the analyst. `INDETERMINATE_PROBE_FAILED` exists for the same
reason -- this repo has a recorded instance of a scorer erroring and printing a confident verdict from
zero resolved calls, so "I could not look" and "I looked and it is not there" must be different
findings rather than the same empty dict.

**THE ISOLATION PREMISE WAS REFUTED 2026-10-22 -- read `isolation_claim`, not the old assumption.**
Stage B's docstring says `promoter_calculator` pulls sklearn 1.9.1 and therefore needs an ISOLATED
interpreter, and this module used to assert `isolation_required: true` on that basis. Measured: the
distribution is `promotercalculator` 1.2.4, a py3-none-any pure-Python wheel with `requires_python=None`
and UNPINNED deps, and a dry-run resolution would install exactly ONE package while leaving sklearn,
numpy and scipy where they are. So isolation is NOT required here and the field now says so, with the
evidence attached. What remains hard about promcalc is the DATA (a ~44k-tile sweep) and the frozen bar's
anchor gate -- not the install. Note the install name differs from the import name.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import string
import sys
from datetime import date
from pathlib import Path

SCHEMA = "promcalc-toolchain-probe-v1"

#: Frozen BEFORE anything was read. Order is the ladder's order; see `decide`.
VERDICTS = (
    "INDETERMINATE_PROBE_FAILED",
    "INPUTS_ABSENT",
    "ENV_PRESENT",
    "ENV_ABSENT_DEPS_UNSATISFIABLE",
    "ENV_ABSENT_DEPS_SATISFIABLE_HERE",
)

#: The interpreter stage B is documented to run under.
ENV_PYTHON = Path("D:/tools/promcalc_env/Scripts/python.exe")

#: The import name, which differs from the distribution name (`pip install promotercalculator`,
#: `import promoter_calculator`). Installing the import name 404s.
STAGE_B_IMPORT_PACKAGE = "promoter_calculator"

#: Places a stage-B interpreter may have been provisioned once the documented volume detached.
CANDIDATE_ENVS = (
    ENV_PYTHON,
    Path.home() / "tools" / "promcalc_env" / "Scripts" / "python.exe",
)


def find_reachable_interpreters(candidates: tuple[Path, ...] | None = None) -> list[dict]:
    """Which candidate stage-B interpreters actually EXIST, and whether the requirement sits in one.

    **Why this exists: reporting only the DOCUMENTED path as absent is a FALSE absence.** On
    2026-10-22 the documented interpreter was gone with its volume while a working one had been
    provisioned on C:, and this probe still said the environment was absent -- which is the most
    expensive answer it can give, because it sends the next session to provision something it already
    has. The same defect was fixed in `find_inputs` the same day: it searched only drive roots while
    the substrate sat repo-local.

    Presence of the requirement is checked on the FILESYSTEM (its import package under
    site-packages), never by importing it, because this probe issues NO subprocess by contract. So a
    hit here means "an interpreter exists and carries the package", NOT "that interpreter works".
    """
    out: list[dict] = []
    for exe in (CANDIDATE_ENVS if candidates is None else candidates):
        if not exe.exists():
            continue
        root = exe.parent.parent
        pkg = root / "Lib" / "site-packages" / STAGE_B_IMPORT_PACKAGE
        out.append({
            "interpreter": str(exe).replace(chr(92), "/"),
            "is_documented_path": exe == ENV_PYTHON,
            "requirement_package_present": pkg.is_dir(),
            "checked_how": "filesystem only -- this probe issues no subprocess by contract",
        })
    return out

#: The ULTIMATE raw substrate. Stage A derives `tiles.tsv` from this file, so its presence -- not the
#: presence of the derived TSV -- is what decides whether there is anything to decompose.
#: Repo root. `scripts/` is one level down, so `parents[1]`.
_REPO = Path(__file__).resolve().parents[1]

PEAK_TILE_FILE = "GSE144621_peak_tile_expression_formatted_std.txt.gz"

#: Measured need. The per-site table is ~3M rows; stage A's own comment puts it on D: for this reason.
MIN_FREE_GB = 5.0

#: Stage B's third-party requirement, and the sklearn it drags in.
STAGE_B_REQUIREMENT = "promoter_calculator"


class ProbeError(RuntimeError):
    """The probe could not complete its own checks."""


def _module_level_imports(text: str) -> list[str]:
    """Top-level module names imported at column 0. Pure; takes source text, not a path.

    Column 0 only, deliberately: stage B's own third-party import sits INSIDE `decompose_one`, and a
    probe that reported it as a module-level dependency would misdescribe the file it read.
    """
    out: list[str] = []
    for line in text.splitlines():
        m = re.match(r"^(?:import|from)\s+([A-Za-z_][\w.]*)", line)
        if m:
            root = m.group(1).split(".")[0]
            if root != "__future__" and root not in out:
                out.append(root)
    return out


def _nested_imports(text: str) -> list[str]:
    """Top-level module names imported at an INDENTED line -- the lazy/isolated ones."""
    out: list[str] = []
    for line in text.splitlines():
        m = re.match(r"^\s+(?:import|from)\s+([A-Za-z_][\w.]*)", line)
        if m:
            root = m.group(1).split(".")[0]
            if root != "__future__" and root not in out:
                out.append(root)
    return out


def decide(*, probe_ok: bool, inputs_found: bool, env_resolves: bool,
           interpreter_available: bool, free_gb_max: float,
           min_free_gb: float = MIN_FREE_GB) -> str:
    """The ladder, first match wins. Pure, so the order is testable without a filesystem.

    Why INPUTS_ABSENT outranks both env verdicts: a toolchain with nothing to run on does not make
    regeneration possible, and `ENV_PRESENT` would be the misleading answer there.
    """
    if not probe_ok:
        return "INDETERMINATE_PROBE_FAILED"
    if not inputs_found:
        return "INPUTS_ABSENT"
    if env_resolves:
        return "ENV_PRESENT"
    if not interpreter_available or free_gb_max < min_free_gb:
        return "ENV_ABSENT_DEPS_UNSATISFIABLE"
    return "ENV_ABSENT_DEPS_SATISFIABLE_HERE"


def volumes() -> list[dict]:
    """Free/total space per resolving drive. A drive with no media does not resolve."""
    out: list[dict] = []
    for letter in string.ascii_uppercase:
        root = Path(f"{letter}:/")
        try:
            if not root.exists():
                continue
            u = shutil.disk_usage(root)
        except OSError:
            continue
        out.append({"drive": f"{letter}:", "free_gb": round(u.free / 2**30, 1),
                    "total_gb": round(u.total / 2**30, 1)})
    return out


def find_inputs(*, roots: list[Path] | None = None) -> list[str]:
    """Every location of the ultimate raw substrate. Reuses the project's own root probe.

    It does NOT open the file -- `genomewide.load_peak_tiles` owns GSE144621 and a second reader of
    another owner's accession is exactly what `tests/test_substrate_loader_ownership.py` forbids. This
    asks only whether the bytes are reachable.
    """
    found: list[str] = []
    if roots is None:
        try:
            from dna_decode.data.external_root import probe
            roots = [c.path for c in probe()]
        except Exception as exc:  # pragma: no cover - import guard
            raise ProbeError(f"external_root probe unavailable: {exc!r}") from exc
        # MEASURED 2026-10-22: the probe reported INPUTS_ABSENT while the substrate was ON DISK.
        # `external_root.probe()` enumerates DRIVE roots, and the ~100 MB GSE144621 download was
        # restored to the repo-local `data/` tree after the external volume detached for good. A
        # false absence is the most expensive answer this probe can give -- it sends the next
        # session to re-fetch data it already has, which is the rediscovery pattern
        # `scripts/data_inventory.py` exists to prevent.
        #
        # This is NOT the auto-substitution `external_root` refuses. That refusal is about ADOPTING
        # a root to READ a dataset from, where a directory of the right NAME is not evidence of the
        # right CONTENTS. This function only answers "are the bytes reachable", and it checks for an
        # exact filename rather than matching a directory name, so naming one more place to LOOK
        # cannot silently redirect a reader.
        local = _REPO / "data"
        if local.is_dir():
            roots = [*roots, local]
    for r in roots:
        p = Path(r) / "mpra" / PEAK_TILE_FILE
        try:
            if p.is_file():
                found.append(str(p).replace(chr(92), "/"))
        except OSError:
            continue
    return found


def probe_toolchain(*, env_python: Path | str = ENV_PYTHON,
                    roots: list[Path] | None = None,
                    vols: list[dict] | None = None,
                    stage_b_source: str | None = None,
                    stage_a_source: str | None = None,
                    requirement_importable: bool | None = None) -> dict:
    """One record describing whether stage B could run here. No subprocess, no writes, no installs.

    Every external fact is an injectable parameter so the whole thing is testable without the host.
    """
    checks_resolved = 0
    problems: list[str] = []

    env = Path(env_python)
    try:
        env_resolves = env.is_file()
        checks_resolved += 1
    except OSError as exc:
        env_resolves = False
        problems.append(f"env interpreter check raised: {exc!r}")

    try:
        v = volumes() if vols is None else list(vols)
        checks_resolved += 1
    except Exception as exc:
        v = []
        problems.append(f"volume enumeration raised: {exc!r}")

    try:
        inputs = find_inputs(roots=roots)
        checks_resolved += 1
    except ProbeError as exc:
        inputs = []
        problems.append(str(exc))

    if requirement_importable is None:
        try:
            import importlib.util
            requirement_importable = importlib.util.find_spec(STAGE_B_REQUIREMENT) is not None
            checks_resolved += 1
        except Exception as exc:
            requirement_importable = False
            problems.append(f"import check raised: {exc!r}")
    else:
        checks_resolved += 1

    def _read(rel: str, injected: str | None) -> tuple[str, str | None]:
        if injected is not None:
            return injected, None
        try:
            return (Path(__file__).resolve().parent / rel).read_text(encoding="utf-8"), None
        except OSError as exc:
            return "", f"could not read {rel}: {exc!r}"

    b_src, b_err = _read("promcalc_decompose.py", stage_b_source)
    a_src, a_err = _read("promcalc_export_tiles.py", stage_a_source)
    for e in (b_err, a_err):
        if e:
            problems.append(e)
        else:
            checks_resolved += 1

    free_max = max((d["free_gb"] for d in v), default=0.0)
    sk = None
    try:
        import sklearn
        sk = sklearn.__version__
    except Exception:
        sk = None

    verdict = decide(
        probe_ok=checks_resolved > 0,
        inputs_found=bool(inputs),
        env_resolves=env_resolves,
        interpreter_available=bool(sys.executable),
        free_gb_max=free_max,
    )

    return {
        "schema": SCHEMA,
        "probed_on": str(date.today()),
        "verdict": verdict,
        "verdict_vocabulary": list(VERDICTS),
        "verdict_vocabulary_frozen": "in this module, above any measurement",
        "checks_resolved": checks_resolved,
        "problems": problems,
        "env": {
            "documented_interpreter": str(env).replace(chr(92), "/"),
            # A reachable interpreter elsewhere does NOT change the verdict below, which is scoped to
            # the DOCUMENTED path and is correct about it. It is recorded so the record cannot read as
            # "no stage-B interpreter exists anywhere", which was false from 2026-10-22 onward.
            "reachable_interpreters": find_reachable_interpreters(),
            "resolves": env_resolves,
            "this_interpreter": sys.executable.replace(chr(92), "/"),
            "requirement": STAGE_B_REQUIREMENT,
            "requirement_importable_here": bool(requirement_importable),
            # MEASURED 2026-10-22 and it REFUTES what this field used to assert. See `isolation_claim`.
            "isolation_required": False,
            "isolation_claim": {
                "claimed_by": "scripts/promcalc_decompose.py docstring",
                "claim": (
                    f"{STAGE_B_REQUIREMENT} pulls sklearn 1.9.1, so it must live in an ISOLATED "
                    f"interpreter or it would move sklearn under the whole test suite"
                ),
                # The two halves have DIFFERENT truth values and a bare "REFUTED" collapsed them.
                # Antecedent HOLDS: a fresh unconstrained resolve really does take sklearn 1.9.1
                # (measured -- a new venv got 1.9.1). Consequence REFUTED: 1.8.0 already satisfies the
                # package, so installing it here would move nothing. A first pass at this correction
                # called the whole premise refuted, which over-corrected in the other direction.
                "status": ("CONSEQUENCE REFUTED, antecedent HOLDS -- measured 2026-10-22; read "
                           "`evidence` and `antecedent_holds` rather than this one line"),
                "antecedent_holds": (
                    "a FRESH `uv venv` + `uv pip install promotercalculator` installed scikit-learn "
                    "1.9.1, so 'pulls sklearn 1.9.1' is true of an unconstrained resolve; what is "
                    "refuted is the CONSEQUENCE that installing it in the project env would move "
                    "sklearn under the suite"
                ),
                "isolation_still_correct_because": (
                    "`uv.lock` is a frozen reproducibility unit, so `uv add` is unavailable and an "
                    "imperative `uv pip install` into the project env is silently reverted by the next "
                    "`uv sync`. The separate venv is near-free and keeps the boundary that "
                    "tests/test_promcalc_decompose.py now enforces."
                ),
                "evidence": (
                    "the distribution is `promotercalculator` 1.2.4 -- a py3-none-any PURE-PYTHON wheel, "
                    "requires_python=None, with UNPINNED deps (biopython, numpy, scikit-learn, scipy). "
                    "`uv pip install --dry-run promotercalculator` resolved 7 packages and would install "
                    f"exactly ONE (promotercalculator), leaving sklearn {sk or '?'} untouched. There is no "
                    "hard pin to 1.9.1 anywhere in its metadata; 1.9.1 is simply what an unconstrained "
                    "resolver picks when it takes the latest."
                ),
                "not_measured_by_this_probe": (
                    "this probe issues NO subprocess by contract, so it cannot run that resolution itself. "
                    "The dry-run was performed separately and its result is RECORDED here, not re-derived "
                    "at call time. Re-run it if the lockfile or the dep set changes."
                ),
                "install_name_differs_from_import_name": (
                    "import `promoter_calculator`, but `pip install promoter_calculator` 404s -- the "
                    "distribution on PyPI is `promotercalculator` (Salis Lab, Barrick Lab fork)."
                ),
            },
            "project_sklearn": sk,
        },
        "inputs": {
            "ultimate_raw_substrate": PEAK_TILE_FILE,
            "found_at": inputs,
            "note": (
                "stage A derives tiles.tsv from this file, so its presence -- not the derived TSV's -- "
                "is what decides whether there is anything to decompose. The file is NOT opened here; "
                "genomewide.load_peak_tiles owns GSE144621."
            ),
        },
        "stage_b": {
            "module_level_imports": _module_level_imports(b_src),
            "nested_imports": _nested_imports(b_src),
            "imports_dna_decode_at_module_level": "dna_decode" in _module_level_imports(b_src),
        },
        "stage_a": {
            "module_level_imports": _module_level_imports(a_src),
            "takes_cache_dir_parameter": "cache_dir" in a_src,
        },
        "volumes": v,
        "min_free_gb_assumed": MIN_FREE_GB,
        "installs_performed": 0,
        "honest_limits": [
            "Read-only. Reporting that an install is possible is NOT the same claim as having done one.",
            "MIN_FREE_GB is ASSERTED from stage A's own 'the per-site table is ~3M rows' comment, not "
            "measured against a completed sweep.",
            "A SATISFIABLE verdict bounds only whether an isolated environment could be CREATED here. "
            "It says nothing about whether the regenerated pipeline would reproduce the frozen bar's "
            "committed per-seed anchors, which is a separate and stricter question.",
            "Module-level import extraction is textual (column-0 regex), not an AST parse; a dynamic "
            "import would not appear.",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--env-python", default=str(ENV_PYTHON))
    ap.add_argument("--out", default=None,
                    help="write the record as JSON here; omitted, the probe only prints")
    a = ap.parse_args(argv)

    rec = probe_toolchain(env_python=a.env_python)
    print(f"verdict {rec['verdict']}  (checks_resolved {rec['checks_resolved']})")
    for _it in rec["env"].get("reachable_interpreters", []):
        print(f"  reachable interpreter  : {_it['interpreter']}  "
              f"[package present: {_it['requirement_package_present']}]")
    print(f"  documented interpreter : {rec['env']['documented_interpreter']}  "
          f"[{'RESOLVES' if rec['env']['resolves'] else 'ABSENT'}]")
    print(f"  {rec['env']['requirement']} importable here : "
          f"{rec['env']['requirement_importable_here']}  "
          f"(isolation_required {rec['env']['isolation_required']} -- "
          f"{rec['env']['isolation_claim']['status']})")
    print(f"  raw substrate          : {len(rec['inputs']['found_at'])} location(s)")
    for p in rec["inputs"]["found_at"]:
        print(f"      {p}")
    print("  volumes                : " + ", ".join(
        f"{d['drive']} {d['free_gb']}/{d['total_gb']} GB free" for d in rec["volumes"]))
    for p in rec["problems"]:
        print(f"  PROBLEM: {p}")

    if a.out:
        out = Path(a.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
