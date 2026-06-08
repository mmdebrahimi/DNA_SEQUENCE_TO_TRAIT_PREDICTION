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
    # TRAITS = the 5 deterministic DECODERS (each a console entry). ANALYSES compose them and are kept
    # OUT of TRAITS so this decoder-registry contract is stable.
    assert set(uni.TRAITS) == {"amr", "pathotype", "plasmid", "serotype", "resfinder"}
    assert set(uni.ANALYSES) == {"concordance", "profile", "coloc"}
    assert not (set(uni.TRAITS) & set(uni.ANALYSES))   # disjoint namespaces


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
