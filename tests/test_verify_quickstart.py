"""Pins the verified quickstart: the documented wheel-only / cached / offline decoder paths all work with
NO torch/transformers (the [ml] extra), NO Docker, NO network. If a refactor breaks the quickstart, this
fails. (The script itself is the proof; this just runs it under pytest.)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.verify_quickstart import STEPS, main  # noqa: E402


def test_quickstart_all_steps_pass():
    assert main([]) == 0


def test_quickstart_covers_each_engine_class():
    labels = " ".join(label for label, *_ in STEPS).lower()
    for token in ("hiv", "fungal", "sars", "profile", "offline"):
        assert token in labels, f"quickstart missing a {token} step"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
