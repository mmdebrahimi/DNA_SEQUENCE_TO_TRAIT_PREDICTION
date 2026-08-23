"""Tests for the unified `dna-decode` dispatcher (dna_decode/cli.py).

Pins: subcommand delegation to the right decoder main(), verbatim argv pass-through, `list` + unknown-trait
handling. Delegation is checked by monkeypatching each decoder main so the test is pure (no Docker / no
AMRFinder / no DB). Runnable via pytest OR standalone.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import dna_decode.cli as uni  # noqa: E402


def test_list_runs(capsys=None):
    rc = uni.main(["list"])
    assert rc == 0


def test_no_args_prints_help():
    assert uni.main([]) == 0


def test_unknown_trait_errors():
    try:
        uni.main(["notatrait"])
    except SystemExit as e:           # argparse .error() raises SystemExit(2)
        assert e.code == 2
        return
    raise AssertionError("expected SystemExit on unknown trait")


def test_amr_delegation_passes_argv(monkeypatch=None):
    captured = {}

    def fake_amr_main(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.amr.cli as amrcli
    orig = amrcli.main
    amrcli.main = fake_amr_main
    try:
        rc = uni.main(["amr", "--drug", "ciprofloxacin", "--amrfinder-run", "X"])
    finally:
        amrcli.main = orig
    assert rc == 0
    assert captured["argv"] == ["--drug", "ciprofloxacin", "--amrfinder-run", "X"]


def test_pathotype_delegation_passes_argv():
    captured = {}

    def fake_patho_main(argv):
        captured["argv"] = argv
        return 0

    import dna_decode.pathotype.cli as pcli
    orig = pcli.main
    pcli.main = fake_patho_main
    try:
        rc = uni.main(["pathotype", "assembly.fna", "--sample-id", "S1"])
    finally:
        pcli.main = orig
    assert rc == 0
    assert captured["argv"] == ["assembly.fna", "--sample-id", "S1"]


def test_traits_registry_matches_console_entries():
    # TRAITS = the deterministic DECODERS (each a console entry). ANALYSES compose them and are kept
    # OUT of TRAITS so this decoder-registry contract is stable.
    #
    # This pin is deliberately hand-maintained: adding a decoder must be a CONSCIOUS act, touching the
    # console entry + TRAITS + this pin + a cell_registry evidence contract. It is not a formality --
    # forward/pigment/flowering (added 2026-07-16) reached the CLI while this pin and the registry both
    # still said they did not exist, and the two guards are what surfaced it.
    assert set(uni.TRAITS) == {"amr", "pathotype", "plasmid", "serotype", "resfinder", "pointfinder",
                               "disinfinder", "mlst", "ktype", "salmserovar", "pneumoserotype", "pgx",
                               "forward", "pigment", "flowering", "inverse", "phage", "kleb", "essentiality",
                               "metabolic", "coatcolor", "morphology", "horsecolor", "catcolor", "plumage",
                               "rabbitcolor", "mousecolor", "cattlecolor", "pigcolor", "sheepcolor",
                               "goatcolor", "alpacacolor", "guineapigcolor", "foxcolor", "donkeycolor",
                               "buffalocolor", "pigeoncolor", "camelcolor", "minkcolor", "roedeercolor",
                               "fba", "motility",
                               # clinvar + hla (2026-08-23): they had console entries and cell_registry
                               # contracts but NO TRAITS row, so they were unreachable via `dna-decode`
                               # and invisible to `dna-decode list` -- while their sibling `pgx` (same
                               # human/VCF shape) was routable. Routing parity, not a new decoder.
                               "clinvar", "hla"}
    # "decode" (added 2026-07-23) is the input-aware ROUTER analysis -- handled inline in cli.py (no
    # delegate module), so it is an ANALYSES entry but not a console script. Conscious addition.
    assert set(uni.ANALYSES) == {"decode", "concordance", "profile", "coloc"}
    assert not (set(uni.TRAITS) & set(uni.ANALYSES))   # disjoint namespaces


def test_every_console_script_is_routable_from_the_unified_cli():
    """The pin above is one-directional: it asserts TRAITS contains exactly a HAND-LISTED set. It cannot
    notice a NEW console script that ships without a TRAITS row -- which is precisely how `clinvar` and
    `hla` stayed unreachable from `dna-decode` (and invisible to `dna-decode list`) while both had console
    entries AND cell_registry contracts. This closes that direction: pyproject is the source of truth for
    what ships, and everything that ships must be reachable through the unified CLI.

    Name mapping mirrors the dispatcher: `dna-pneumo-serotype` -> trait `pneumoserotype` (hyphens dropped).
    """
    import tomllib
    root = Path(__file__).resolve().parent.parent
    scripts = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]["scripts"]
    routable = set(uni.TRAITS) | set(uni.ANALYSES)
    unroutable = sorted(name for name in scripts
                        if name != "dna-decode"
                        and name.removeprefix("dna-").replace("-", "") not in routable)
    assert not unroutable, (
        f"console scripts with no `dna-decode <trait>` route: {unroutable}. Add a TRAITS entry + dispatch "
        f"in dna_decode/cli.py (and the pin above), or the decoder ships unreachable from the unified CLI.")


def test_every_console_script_actually_imports():
    """A console entry whose target module fails to import (or lacks its main) is a broken install-time
    promise that no other test exercises -- nothing imports every decoder in one place."""
    import importlib
    import tomllib
    root = Path(__file__).resolve().parent.parent
    scripts = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]["scripts"]
    broken = []
    for name, spec in sorted(scripts.items()):
        mod, _, fn = spec.partition(":")
        try:
            if not hasattr(importlib.import_module(mod), fn):
                broken.append(f"{name}: {spec} (module lacks {fn}())")
        except Exception as e:  # noqa: BLE001 -- any import failure is the defect
            broken.append(f"{name}: {spec} ({type(e).__name__}: {e})")
    assert not broken, "broken console entries: " + "; ".join(broken)


def test_every_trait_has_an_evidence_contract():
    """The two registries must AGREE: a trait routable from the unified CLI needs a trust-surface contract.

    Without this, the pin above and the cell_registry coverage guard could drift apart -- each satisfied
    while a decoder ships invisibly through the gap between them.

    (clinvar/hla WERE console entries but not dna-decode TRAITS, and this docstring used to say so. That
    gap was closed 2026-08-23: both already had cell_registry contracts whose `route` is dna-clinvar /
    dna-hla, so they satisfied this guard the whole time -- what they lacked was a TRAITS row, which made
    them unreachable from the unified CLI and absent from `dna-decode list`.)
    """
    from dna_decode.data.cell_registry import cells

    contracted = {c.target for c in cells()} | {c.route.removeprefix("dna-") for c in cells()}
    missing = {t for t in uni.TRAITS if t not in contracted}
    assert not missing, f"CLI traits with no cell_registry evidence contract: {sorted(missing)}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
