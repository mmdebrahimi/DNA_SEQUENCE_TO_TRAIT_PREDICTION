"""Manual wheel-gate: build the wheel and PROVE it ships the trust-surface report cards (the packaging gate).

This is the artifact-boundary check the in-process quickstart cannot do. It builds the wheel via `uv build`,
inspects its contents, and asserts the 4 cards trust_surface loads are present at dna_decode/report_cards/.
Optionally (--fresh-env) it installs the wheel into a throwaway venv OUTSIDE the repo and runs the console
script, proving the badges load from the PACKAGED cards (not repo wiki/).

Run: `uv run python scripts/verify_wheel_ships_cards.py`            (build + contents assert; fast)
     `uv run python scripts/verify_wheel_ships_cards.py --fresh-env` (+ install into a temp venv; heavier)

NOT a pytest test (it builds + optionally installs). Exit 0 = the wheel ships the cards.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CARDS = ["amr_portal_independent_report_card.json", "decoder_validation_report_card.json",
         "hiv_decoder_report_card.json", "tb_report_card.json"]


def _latest_wheel() -> Path | None:
    whls = sorted((REPO / "dist").glob("*.whl"), key=lambda p: p.stat().st_mtime)
    return whls[-1] if whls else None


def build_wheel() -> Path:
    subprocess.run(["uv", "build", "--wheel"], cwd=REPO, check=True, capture_output=True, text=True)
    whl = _latest_wheel()
    if whl is None:
        raise SystemExit("ERROR: no wheel produced by `uv build --wheel`")
    return whl


def assert_cards_in_wheel(whl: Path) -> list[str]:
    names = zipfile.ZipFile(whl).namelist()
    shipped = sorted(n for n in names if "dna_decode/report_cards/" in n and n.endswith(".json"))
    missing = [c for c in CARDS if not any(n.endswith(c) for n in shipped)]
    print(f"wheel: {whl.name}")
    for s in shipped:
        print(f"  ships: {s}")
    if missing:
        raise SystemExit(f"FAIL: wheel is MISSING report cards: {missing}")
    print(f"PASS: all {len(CARDS)} trust cards ship in the wheel")
    return shipped


def fresh_env_smoke(whl: Path) -> None:
    """Install the wheel into a throwaway venv OUTSIDE the repo; prove the badge loads from packaged cards."""
    with tempfile.TemporaryDirectory() as td:
        venv = Path(td) / "venv"
        subprocess.run(["uv", "venv", str(venv)], check=True, capture_output=True, text=True)
        py = venv / ("Scripts/python.exe" if sys.platform == "win32" else "bin/python")
        subprocess.run(["uv", "pip", "install", "--python", str(py), str(whl)],
                       check=True, capture_output=True, text=True)
        # run from the temp dir (NOT the repo) so `import dna_decode` resolves to the INSTALLED wheel
        check = (
            "import dna_decode.data.trust_surface as t;"
            "assert t._PKG_CARDS.exists(), 'packaged cards missing in install';"
            "b=t.lookup_trust('efavirenz');"
            "assert b['tier']=='INDEPENDENT_WETLAB', b;"
            "print('PASS: installed wheel serves trust badges from PACKAGED cards (repo wiki/ not used)')"
        )
        r = subprocess.run([str(py), "-c", check], cwd=td, capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr.strip())
        if r.returncode != 0:
            raise SystemExit("FAIL: fresh-env artifact-boundary check failed")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fresh-env", action="store_true", help="also install into a temp venv + run the badge")
    a = ap.parse_args(argv)
    whl = build_wheel()
    assert_cards_in_wheel(whl)
    if a.fresh_env:
        fresh_env_smoke(whl)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
