"""Shared pytest fixtures."""
from pathlib import Path

import pytest


@pytest.fixture
def fixtures_dir() -> Path:
    """Path to tests/fixtures/. Populated by Step 15 (smoke pipeline + fixtures)."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def project_root() -> Path:
    """Path to the project root."""
    return Path(__file__).parent.parent
