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
import tomllib
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
def _force_included_cards() -> list[str]:
    """What the wheel is SUPPOSED to ship -- DERIVED from pyproject, never hand-listed.

    This list used to be four hand-written names while pyproject force-included five, so this
    script -- the artifact-boundary PROOF the pytest suite defers to -- printed
    "PASS: all 4 trust cards ship" without ever checking the HCMV card. A proof that omits
    one of the things it proves is worse than no proof, because it reads as coverage.
    """
    d = tomllib.load(open(REPO / "pyproject.toml", "rb"))
    fi = d["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    return sorted(
        src.split("/", 1)[1]
        for src, dest in fi.items()
        if src.startswith("wiki/") and dest.startswith("dna_decode/report_cards/")
    )


CARDS = _force_included_cards()


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
    if not CARDS:
        raise SystemExit("FAIL: derived no cards from pyproject -- the force-include table moved")
    shipped_names = {n.rsplit("/", 1)[-1] for n in shipped}
    missing = [c for c in CARDS if c not in shipped_names]
    extra = sorted(shipped_names - set(CARDS))
    print(f"wheel: {whl.name}")
    for s in shipped:
        print(f"  ships: {s}")
    if missing:
        raise SystemExit(f"FAIL: wheel is MISSING report cards: {missing}")
    # Asserted BOTH ways. HONEST SCOPE: this branch is NOT currently reachable and was not proven to
    # fire -- CARDS is derived FROM force-include, so `extra` needs the backend to ship something
    # unlisted, and a planted stray json in dna_decode/report_cards/ was MEASURED not to reach the
    # wheel (hatch ships only the force-included files). Kept as a one-line guard against a build-
    # backend behaviour change, not presented as an active control.
    if extra:
        raise SystemExit(f"FAIL: wheel ships report cards NOT in the force-include table: {extra}")
    print(f"PASS: all {len(CARDS)} force-included trust cards ship in the wheel (and nothing else)")
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
