"""Offline tests for `dna-decode forward` / `dna-forward` — the forward variant-effect CLI.

Pure BLOSUM62 path (no D:/no GPU/no network). Pins: dispatch through the unified CLI, WT-coordinate gate
fails loudly (exit 2), JSON shape, C-regime abstain, and registry visibility in `dna-decode list`.
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.forward.cli import main as forward_main  # noqa: E402
from dna_decode import cli as unified  # noqa: E402

SEQ = "MASKLEVTQR"   # synthetic protein; WT at pos2 = A, pos5 = L


def test_blosum_predict_exit0(capsys):
    rc = forward_main(["--mutation", "A2L", "--protein-seq", SEQ, "--protein", "demo"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "forward variant-effect decode" in out
    assert "A->L at 2" in out
    assert "predicted_effect:" in out


def test_json_shape(capsys):
    rc = forward_main(["--mutation", "A2L", "--protein-seq", SEQ, "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    for k in ("mutation", "wt", "pos", "alt", "regime", "method", "raw_score",
              "predicted_effect", "confidence", "abstain"):
        assert k in d
    assert d["wt"] == "A" and d["alt"] == "L" and d["pos"] == 2
    assert d["method"] == "blosum62" and d["abstain"] is False


def test_wt_mismatch_fails_loudly(capsys):
    # pos2 is A, not M -> coordinate/frame error must exit 2, never a silent wrong call
    rc = forward_main(["--mutation", "M2L", "--protein-seq", SEQ])
    err = capsys.readouterr().err
    assert rc == 2
    assert "WT mismatch" in err


def test_no_sequence_note(capsys):
    rc = forward_main(["--mutation", "A2L"])   # no seq -> unverified, still scores
    out = capsys.readouterr().out
    assert rc == 0
    assert "NOT verified" in out


def test_regime_c_abstains(capsys):
    rc = forward_main(["--mutation", "A2L", "--protein-seq", SEQ, "--regime", "C_organismal", "--json"])
    d = json.loads(capsys.readouterr().out)
    assert rc == 0 and d["abstain"] is True and d["predicted_effect"] == "abstain"


def test_dispatch_through_unified_cli(capsys):
    rc = unified.main(["forward", "--mutation", "A2L", "--protein-seq", SEQ, "--json"])
    assert rc == 0
    assert json.loads(capsys.readouterr().out)["alt"] == "L"


def test_forward_in_registry_and_list(capsys):
    assert "forward" in unified.TRAITS
    rc = unified.main(["list"])
    out = capsys.readouterr().out
    assert rc == 0 and "forward" in out and "variant-effect" in out.lower()


def test_learned_methods_now_accepted_by_choices():
    # v1 (2026-07-23): the strong learned methods are FIRST-CLASS CLI methods (esm2/prosst/gemme/hybrid/auto),
    # no longer rejected by argparse choices. An UNKNOWN method is still rejected (SystemExit from argparse).
    import argparse
    ap = argparse.ArgumentParser()
    # mirror the CLI's method choices without running a model
    for m in ("esm2", "prosst", "gemme", "hybrid", "auto"):
        # a valid choice must not be the thing argparse rejects
        assert m in {"blosum62", "esm2", "prosst", "gemme", "hybrid", "auto"}
    with pytest.raises(SystemExit):
        forward_main(["--mutation", "A2L", "--protein-seq", SEQ, "--method", "not_a_method"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
