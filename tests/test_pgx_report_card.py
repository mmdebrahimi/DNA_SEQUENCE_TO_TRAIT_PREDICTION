"""Pins the standing PGx trust-surface report card (read-only roll-up; exit 0 always)."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_WRITTEN_ARTIFACTS = (
    REPO / "wiki" / "pgx_report_card.json",
    REPO / "wiki" / "pgx_report_card.md",
)

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def test_report_card_builds_and_covers_the_three_genes():
    from build_pgx_report_card import main
    assert main() == 0
    import json
    rep = json.loads((REPO / "wiki" / "pgx_report_card.json").read_text(encoding="utf-8"))
    genes = {c["gene"] for c in rep["cells"]}
    assert {"CYP2C19", "CYP2C9", "VKORC1"} <= genes
    assert rep["schema"] == "pgx-report-card-v0"
    # no fabricated aggregate headline
    assert "aggregate" not in rep and "overall_concordance" not in rep


def test_report_card_reflects_getrm_numbers_when_present():
    from build_pgx_report_card import main
    main()
    import json
    rep = json.loads((REPO / "wiki" / "pgx_report_card.json").read_text(encoding="utf-8"))
    by = {c["gene"]: c for c in rep["cells"]}
    # if the GeT-RM sidecars are present, the card must reflect them (not fabricate)
    if by["CYP2C19"]["getrm"]:
        assert by["CYP2C19"]["getrm"] == "72/72"
    if by["CYP2C9"]["getrm"]:
        assert by["CYP2C9"]["getrm"] == "73/73"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))


# --------------------------------------------------------------------------- #
# Added 2026-10-20: these tests call the builder's real `main()`, which WRITES the
# committed `wiki/` artifacts. The only field that moves is `analysis_date`, so every
# full-suite run left 2 files dirty with a meaningless date bump -- noise in every
# session's `git status`, and a standing invitation to commit a no-op change. (It was
# manually restored twice on 2026-10-20 before anyone fixed the cause.)
#
# Redirecting the builder's output constant is NOT an option: it also READS its input
# sidecars from the same directory. So the test keeps exercising the real write path and
# simply cleans up after itself. Safe for the currency test, which reads the COMMITTED
# blob via `git show HEAD:...` rather than the working file.
@pytest.fixture(autouse=True)
def _restore_committed_artifacts():
    before = {p: (p.read_bytes() if p.exists() else None) for p in _WRITTEN_ARTIFACTS}
    try:
        yield
    finally:
        for p, blob in before.items():
            if blob is None:
                if p.exists():
                    p.unlink()
            elif p.read_bytes() != blob:
                p.write_bytes(blob)
