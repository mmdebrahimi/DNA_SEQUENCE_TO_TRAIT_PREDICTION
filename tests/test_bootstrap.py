"""Step 1 acceptance tests: package imports + config loads.

Verifies the bootstrap is sound before Wave 1 agents start adding modules.
"""
from pathlib import Path

import pytest
import yaml


def test_package_imports_cleanly():
    """All 5 subpackages + the top-level package import without side-effects."""
    import dna_decode
    import dna_decode.data
    import dna_decode.eval
    import dna_decode.interp
    import dna_decode.models
    import dna_decode.viz

    assert dna_decode.__version__ == "0.0.1"


def test_datasources_yaml_loads(project_root: Path):
    """config/datasources.yaml parses as valid YAML with expected top-level keys."""
    cfg_path = project_root / "config" / "datasources.yaml"
    assert cfg_path.exists(), "config/datasources.yaml must exist after bootstrap"

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    required_keys = {
        "refseq",
        "card",
        "amrfinder",
        "bvbrc_ast",
        "enterobase",
        "phase1_drugs",
        "compute",
        "foundation_models",
    }
    missing = required_keys - cfg.keys()
    assert not missing, f"datasources.yaml missing keys: {missing}"


def test_phase1_drugs_have_known_loci(project_root: Path):
    """Each Phase 1 drug has at least one known resistance locus + a precision target."""
    cfg_path = project_root / "config" / "datasources.yaml"
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    drugs = cfg["phase1_drugs"]
    assert len(drugs) == 3, "Phase 1 specifies cipro + ceftriaxone + tetracycline"

    for drug in drugs:
        assert drug["known_loci"], f"{drug['name']} must have ≥1 known locus"
        assert "attribution_precision_target" in drug
        assert 0 < drug["attribution_precision_target"] <= 1


def test_foundation_models_complete(project_root: Path):
    """All 4 leaderboard foundation models are declared with required fields.

    Phase 2 added `mock` for plumbing-only smoke runs — it must be a subset of
    `foundation_models` but is NOT a leaderboard contender.
    """
    cfg_path = project_root / "config" / "datasources.yaml"
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    models = cfg["foundation_models"]
    leaderboard = {"evo", "dnabert2", "nucleotide_transformer", "gena_lm"}
    assert leaderboard.issubset(set(models.keys())), (
        f"Expected leaderboard models {leaderboard}, got {set(models.keys())}"
    )
    # Extra entries (e.g., `mock`) are permitted but must satisfy the same schema.
    for name, meta in models.items():
        assert "huggingface_id" in meta
        assert "embedding_dim" in meta and meta["embedding_dim"] > 0
        assert "max_context" in meta and meta["max_context"] > 0
