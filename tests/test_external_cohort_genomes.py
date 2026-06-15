"""Unit tests for external-cohort genome resolution (no network — fake resolver)."""
from __future__ import annotations

from dna_decode.data import external_cohort_genomes as ecg


class _FakeResolver:
    def __init__(self, mapping: dict[str, list[str]]):
        self.mapping = mapping
        self.calls: list[str] = []

    def biosample_to_assemblies(self, biosample: str) -> list[str]:
        self.calls.append(biosample)
        return list(self.mapping.get(biosample, []))


# --------------------------------------------------------------------------- #
# pick_assembly
# --------------------------------------------------------------------------- #
def test_pick_prefers_gcf():
    assert ecg.pick_assembly(["GCA_1.1", "GCF_1.1"]) == "GCF_1.1"


def test_pick_gca_when_no_gcf():
    assert ecg.pick_assembly(["GCA_2.1", "GCA_1.1"]) == "GCA_1.1"  # sorted-smallest GCA


def test_pick_empty_is_none():
    assert ecg.pick_assembly([]) is None
    assert ecg.pick_assembly(["", None]) is None


def test_pick_is_deterministic():
    a = ecg.pick_assembly(["GCF_9.1", "GCF_3.1", "GCA_1.1"])
    b = ecg.pick_assembly(["GCA_1.1", "GCF_3.1", "GCF_9.1"])
    assert a == b == "GCF_3.1"


# --------------------------------------------------------------------------- #
# resolve_cohort_genomes
# --------------------------------------------------------------------------- #
def test_resolve_splits_free_and_assembly_required():
    resolver = _FakeResolver({
        "SAMN_a": ["GCF_a.1"],
        "SAMN_b": ["GCA_b.1"],
        "SAMN_reads": [],            # reads-only -> ASSEMBLY_REQUIRED
    })
    out = ecg.resolve_cohort_genomes(["SAMN_a", "SAMN_b", "SAMN_reads"], resolver)
    assert out["free"] == {"SAMN_a": "GCF_a.1", "SAMN_b": "GCA_b.1"}
    assert out["assembly_required"] == ["SAMN_reads"]
    assert out["n_free"] == 2
    assert out["n_assembly_required"] == 1


def test_reads_only_not_dropped():
    resolver = _FakeResolver({"SAMN_reads": []})
    out = ecg.resolve_cohort_genomes(["SAMN_reads"], resolver)
    assert out["free"] == {}
    assert "SAMN_reads" in out["assembly_required"]  # counted, not silently dropped


def test_keys_are_gca_accessions_for_ensure_run():
    # The free values must be downloadable GCA/GCF accessions (the ensure_run contract).
    resolver = _FakeResolver({"SAMN_a": ["GCF_a.1", "GCA_a.1"]})
    out = ecg.resolve_cohort_genomes(["SAMN_a"], resolver)
    gca = out["free"]["SAMN_a"]
    assert gca.startswith(("GCA_", "GCF_"))


def test_dedup_and_sorted_input():
    resolver = _FakeResolver({"SAMN_a": ["GCF_a.1"], "SAMN_b": ["GCF_b.1"]})
    out = ecg.resolve_cohort_genomes(["SAMN_b", "SAMN_a", "SAMN_a"], resolver)
    assert out["n_free"] == 2  # deduped
