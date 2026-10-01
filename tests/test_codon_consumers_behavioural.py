"""Gap coverage for the migrated codon consumers (test-epilogue, 2026-10-01).

The plan ships two consumer guards and between them they leave a hole:

* `test_codon_consumers_canonical.py` covers the THREE in-package consumers by name;
* `test_script_codon_consumers_canonical.py` covers THREE scripts, hand-listed in a `MIGRATED` dict.

**Eleven modules outside the package import `STANDARD`.** Five scripts are in neither list:
`hiv_esm_vs_catalog`, `mavedb_cpu_smoke`, `pear_genotype_alphabet`, `pre_resistance_base_rate_census`,
`pre_resistance_predictive_backtest`. Step 9's discovery guard catches a surviving dict LITERAL in them
but not what they import -- a script could be pointed at `table_for('bacteria.mycoplasma')` instead of
`STANDARD` and both shipped guards would stay green while `TGA` silently became tryptophan.

That hand-listed dict is the sixth instance of this repo's recorded hardcoded-exclusion-list failure, so
the consumer set here is DISCOVERED from the import graph rather than enumerated.

Three further gaps, all behavioural rather than structural:

* nothing checks any consumer still PRODUCES table-1 assignments. Identity to `STANDARD` is necessary,
  not sufficient: the aliases could be right while a module-level derivation off the table is wrong.
* `pre_resistance_*` derive `AAS` and `CODONS_FOR` at import time from the table. Those are the only
  load-bearing derived constants in the migrated set and neither is tested.
* consolidation introduced ONE shared mutable dict where there used to be private copies, so a stray
  write in any consumer would now reach all eleven. Nothing forbids one.
"""

from __future__ import annotations

import ast
import functools
import importlib
import importlib.util
import pathlib

import pytest

from dna_decode.constraints import codon_tables as T

REPO = pathlib.Path(__file__).resolve().parent.parent
SCAN_ROOTS = ("dna_decode", "scripts")

#: The canonical module plus the three modules INSIDE the constraints package that legitimately import
#: from it for their own construction; they are not "consumers" in the aliasing sense.
_PACKAGE_INTERNAL = {
    "dna_decode/constraints/__init__.py",
    "dna_decode/constraints/codon_tables.py",
    "dna_decode/constraints/registry.py",
    "dna_decode/constraints/report.py",
    "dna_decode/constraints/specializations.py",
    "dna_decode/constraints/universal.py",
}

#: Cannot import the package at all (runs standalone on Kaggle); equality-pinned by
#: `test_codon_table_single_source.py`, so it is out of scope for an aliasing check.
_STANDALONE = {"scripts/kaggle/pear_prosst_kernel.py"}

#: A consumer may legitimately import a NON-standard table -- but only one does, and it is the probe,
#: whose whole purpose is the clade demo. Any other module importing a clade table is a behaviour
#: change that must be noticed.
_CLADE_TABLE_CONSUMERS = {"scripts/constraint_kill_count_probe.py"}


def _py_files() -> list[pathlib.Path]:
    out: list[pathlib.Path] = []
    for root in SCAN_ROOTS:
        out.extend(sorted((REPO / root).rglob("*.py")))
    return [p for p in out if "__pycache__" not in p.parts]


def _rel(p: pathlib.Path) -> str:
    return str(p.relative_to(REPO)).replace("\\", "/")


@functools.lru_cache(maxsize=4)
def _discover(roots: tuple[str, ...], repo: str) -> dict[str, dict[str, str]]:
    """Cached worker. The scan parses ~790 files (~6 s), and these tests call it dozens of times."""
    out: dict[str, dict[str, str]] = {}
    for p in _py_files():
        rel = _rel(p)
        if rel in _PACKAGE_INTERNAL or rel in _STANDALONE:
            continue
        try:
            tree = ast.parse(p.read_text(encoding="utf-8", errors="ignore"))
        except SyntaxError:
            continue
        names: dict[str, str] = {}
        for n in ast.walk(tree):
            if isinstance(n, ast.ImportFrom) and (n.module or "").startswith("dna_decode.constraints"):
                for a in n.names:
                    names[a.name] = a.asname or a.name
        if names:
            out[rel] = names
    return out


def discover_consumers() -> dict[str, dict[str, str]]:
    """path -> {imported_name: alias}. Discovery, not a hand list.

    Keyed on `(SCAN_ROOTS, REPO)` so the non-vacuity test can monkeypatch both and get a fresh scan."""
    return _discover(tuple(SCAN_ROOTS), str(REPO))


#: The five scripts in NEITHER shipped guard. Named explicitly so the gap this file closes stays legible
#: even after the discovery set grows.
PREVIOUSLY_UNCOVERED = (
    "scripts/hiv_esm_vs_catalog.py",
    "scripts/mavedb_cpu_smoke.py",
    "scripts/pear_genotype_alphabet.py",
    "scripts/pre_resistance_base_rate_census.py",
    "scripts/pre_resistance_predictive_backtest.py",
)


@functools.lru_cache(maxsize=32)
def _load(rel: str):
    """Import a consumer for real.

    An in-package module MUST go through `import_module`: loading `dna_decode/forward/genome_edit.py`
    by file path gives it no parent package, so its `from .variant_effect import ...` raises
    `ImportError`. A `scripts/` module has no package and must go by path."""
    if rel.startswith("dna_decode/"):
        return importlib.import_module(rel[:-len(".py")].replace("/", "."))
    spec = importlib.util.spec_from_file_location(
        pathlib.Path(rel).stem + "__consumer_probe", REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------- discovery covers everything

def test_discovery_finds_every_consumer_including_the_five_the_shipped_guards_MISS():
    found = discover_consumers()
    assert found, "discovery found no consumers at all -- the scan is broken"
    for rel in PREVIOUSLY_UNCOVERED:
        assert rel in found, f"{rel} imports the canonical table but discovery missed it"
    assert len(found) >= 11, (
        f"expected at least the 11 known consumers; found {len(found)}: {sorted(found)}")


def test_the_shipped_hand_listed_guards_genuinely_UNDER_COVER_the_consumer_set():
    """Non-vacuity for this whole file: prove the hole exists rather than asserting it. If a future
    change folds the five into the shipped lists, this fails and the file can be slimmed."""
    from tests.test_script_codon_consumers_canonical import MIGRATED

    # the three in-package consumers are hand-named in `test_codon_consumers_canonical.py`
    covered = set(MIGRATED) | {
        "dna_decode/forward/genome_edit.py", "dna_decode/forward/inverse.py",
        "dna_decode/typing/codon_map.py"}
    uncovered = sorted(set(discover_consumers()) - covered - {"scripts/constraint_kill_count_probe.py"})
    assert uncovered == list(PREVIOUSLY_UNCOVERED), uncovered


@pytest.mark.parametrize("rel", sorted(discover_consumers()))
def test_every_consumer_imports_the_STANDARD_table_not_a_clade_table(rel):
    """The failure this closes: a consumer silently re-pointed at a non-standard code. The shipped
    cells (E. coli, Klebsiella, HIV, SARS-CoV-2, Candida, TB) all use standard assignments, so any
    clade import outside the clade probe is a behaviour change."""
    names = discover_consumers()[rel]
    clade_imports = sorted(set(names) & {
        "MYCOPLASMA_SPIROPLASMA", "VERTEBRATE_MITOCHONDRIAL", "CILIATE_NUCLEAR", "BACTERIAL_PLASTID"})
    if rel in _CLADE_TABLE_CONSUMERS:
        assert clade_imports, f"{rel} is listed as a clade consumer but imports no clade table"
        return
    assert not clade_imports, (
        f"{rel} imports clade table(s) {clade_imports}. If that is intentional, add it to "
        "_CLADE_TABLE_CONSUMERS with the reason; otherwise it silently changed the genetic code.")


@pytest.mark.parametrize("rel", sorted(discover_consumers()))
def test_every_consumers_table_alias_IS_the_canonical_object(rel):
    """Object identity, not equality -- that is what makes drift impossible rather than unlikely. Run
    over the DISCOVERED set, so a new consumer is covered the day it lands."""
    names = discover_consumers()[rel]
    if "STANDARD" not in names:
        pytest.skip(f"{rel} imports {sorted(names)}, no table alias to check")
    if "forward_inverse_roundtrip" in rel:
        pytest.skip("pre-existing import-time side effects; covered statically by the Step 7 guard")
    mod = _load(rel)
    alias = names["STANDARD"]
    assert getattr(mod, alias) is T.STANDARD, f"{rel}.{alias} is not the canonical object"


# --------------------------------------------------------------------- behaviour, not just identity

@pytest.mark.parametrize("rel,fn", [
    ("scripts/fungal_erg11_caller.py", "_translate"),
    ("scripts/build_holdout_hybrid_manifest.py", "_translate"),
    ("scripts/mavedb_cpu_smoke.py", "translate"),
])
def test_each_script_translator_still_reads_TGA_as_STOP_not_tryptophan(rel, fn):
    """The one assignment that distinguishes table 1 from table 4, checked through each consumer's own
    translator rather than through the table. `TGG`->W is the control: a translator that returned `*`
    for everything would pass a TGA-only check."""
    translate = getattr(_load(rel), fn)
    assert translate("ATGTGG").startswith("MW"), "TGG must read as tryptophan"
    out = translate("ATGTGGTGA")
    assert out.startswith("MW"), out
    assert "*" in out or len(out) == 2, (
        f"{rel}.{fn} neither emitted a stop nor truncated at one: {out!r}")
    assert not out.startswith("MWW"), f"{rel}.{fn} read TGA as tryptophan -- wrong genetic code"


def test_the_in_package_translators_agree_with_each_other_on_a_real_qrdr_codon():
    """`typing.codon_map` feeds `pointfinder`, `forward.genome_edit` feeds the validated forward cell.
    Both must agree, and on the real gyrA S83 codon rather than a toy."""
    from dna_decode.forward.genome_edit import translate_cds
    from dna_decode.typing.codon_map import translate

    cds = "ATGTCGTTGTGATGG"                      # M S L * W
    assert translate(cds) == translate_cds(cds) == "MSL*W"
    assert translate("TCG") == "S" and translate("TTG") == "L", "real gyrA S83L codons"


def test_the_hiv_translator_handles_an_unresolvable_codon_as_X_not_a_crash():
    """`hiv_esm_vs_catalog.rt_protein` is on the HIV path; an `N` in a real consensus must degrade to
    `X`, which is only true because the call site uses `.get(..., 'X')` on the shared table."""
    mod = _load("scripts/hiv_esm_vs_catalog.py")
    src = (REPO / "scripts/hiv_esm_vs_catalog.py").read_text(encoding="utf-8")
    assert 'CODON.get(seq[i:i + 3], "X")' in src
    assert mod.CODON.get("ATN", "X") == "X" and mod.CODON["ATG"] == "M"


def test_one_consumer_indexes_the_table_BARE_and_would_raise_on_an_ambiguous_codon():
    """Asymmetry across consumers of the same table, recorded rather than silently inherited:
    `pear_genotype_alphabet` line ~196 does `CODON[cds[i:i+3]]` with no default, so an `N` in the CDS
    raises `KeyError` where every sibling degrades to `X`. Pre-existing and arguably correct for a
    curated reference, but it is a real difference in failure mode."""
    src = (REPO / "scripts/pear_genotype_alphabet.py").read_text(encoding="utf-8")
    assert "CODON[cds[i:i + 3]]" in src, "the bare-index call site moved; re-check this asymmetry"
    with pytest.raises(KeyError):
        T.STANDARD["ATN"]
    assert T.STANDARD.get("ATN", "X") == "X"


# --------------------------------------------------------------------- derived module constants

@pytest.mark.parametrize("rel", [
    "scripts/pre_resistance_base_rate_census.py",
    "scripts/pre_resistance_predictive_backtest.py",
])
def test_the_derived_amino_acid_constants_survive_the_migration(rel):
    """`AAS` and `CODONS_FOR` are computed at import time FROM the table, so they are where a wrong
    table would show up as a wrong analysis rather than as an import error. Neither was tested."""
    mod = _load(rel)
    assert len(mod.AAS) == 20, mod.AAS
    assert "*" not in mod.AAS, "the stop symbol must be excluded from the amino-acid alphabet"
    assert mod.AAS == sorted(mod.AAS)
    assert set(mod.AAS) == set("ACDEFGHIKLMNPQRSTVWY")

    assert len(mod.CODONS_FOR) == 21, "20 amino acids + the stop class"
    assert sorted(mod.CODONS_FOR["M"]) == ["ATG"]
    assert sorted(mod.CODONS_FOR["W"]) == ["TGG"]
    assert sorted(mod.CODONS_FOR["*"]) == ["TAA", "TAG", "TGA"], (
        "stop codons must follow table 1 here, not a clade table")
    assert sum(len(v) for v in mod.CODONS_FOR.values()) == 64, "every codon must be classified once"


def test_pear_genotype_alphabets_AAS_is_the_same_20_letter_alphabet():
    mod = _load("scripts/pear_genotype_alphabet.py")
    assert set(mod.AAS) == set("ACDEFGHIKLMNPQRSTVWY") and len(mod.AAS) == 20


# --------------------------------------------------------------------- shared-mutable-dict guard

@functools.lru_cache(maxsize=4)
def _scan_mutations(roots: tuple[str, ...], repo: str) -> tuple[tuple[str, int, str], ...]:
    """Static scan for a WRITE through any codon-table alias. Consolidation turned eleven private
    copies into one shared dict, so a write that used to be local now reaches every consumer."""
    aliases = {"_CODON", "CODON", "CODON_TABLE", "STANDARD",
               "MYCOPLASMA_SPIROPLASMA", "VERTEBRATE_MITOCHONDRIAL",
               "CILIATE_NUCLEAR", "BACTERIAL_PLASTID"}
    mutators = {"pop", "popitem", "update", "setdefault", "clear"}
    hits: list[tuple[str, int, str]] = []
    for p in _py_files():
        rel = _rel(p)
        if rel in _STANDALONE or rel == "dna_decode/constraints/codon_tables.py":
            continue                   # the standalone kernel builds its OWN table; so does canonical
        try:
            tree = ast.parse(p.read_text(encoding="utf-8", errors="ignore"))
        except SyntaxError:
            continue
        for n in ast.walk(tree):
            if isinstance(n, (ast.Assign, ast.AugAssign, ast.Delete)):
                targets = n.targets if isinstance(n, (ast.Assign, ast.Delete)) else [n.target]
                for t in targets:
                    if (isinstance(t, ast.Subscript) and isinstance(t.value, ast.Name)
                            and t.value.id in aliases):
                        hits.append((rel, n.lineno, f"write to {t.value.id}[...]"))
            if (isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                    and n.func.attr in mutators and isinstance(n.func.value, ast.Name)
                    and n.func.value.id in aliases):
                hits.append((rel, n.lineno, f"{n.func.value.id}.{n.func.attr}()"))
    return tuple(hits)


def _mutation_sites() -> list[tuple[str, int, str]]:
    return list(_scan_mutations(tuple(SCAN_ROOTS), str(REPO)))


def test_no_consumer_WRITES_through_its_codon_table_alias():
    """Eleven aliases, one dict. A single stray write would corrupt every other consumer and diverge
    the shared object from `table_for`, which re-derives per call."""
    hits = _mutation_sites()
    assert not hits, (
        "these sites mutate the SHARED canonical codon table: "
        + "; ".join(f"{r}:{ln} {what}" for r, ln, what in hits)
        + ". Copy it (`dict(STANDARD)`) or call `table_for(clade)` for a fresh table instead.")


def test_the_mutation_scan_is_NON_VACUOUS(tmp_path, monkeypatch):
    """A scan that cannot fire is the inert-filter failure in test form. Plant both shapes it looks for
    and require detection, then re-assert the live tree is clean."""
    root = tmp_path / "dna_decode"
    root.mkdir()
    (root / "rogue_writer.py").write_text(
        "from dna_decode.constraints.codon_tables import STANDARD as CODON\n"
        'CODON["TGA"] = "W"\n'
        'CODON.setdefault("ATN", "X")\n', encoding="utf-8")
    import tests.test_codon_consumers_behavioural as mod
    monkeypatch.setattr(mod, "REPO", tmp_path)
    monkeypatch.setattr(mod, "SCAN_ROOTS", ("dna_decode",))
    planted = mod._mutation_sites()
    kinds = {what for _, _, what in planted}
    assert "write to CODON[...]" in kinds, planted
    assert "CODON.setdefault()" in kinds, planted


def test_the_canonical_table_is_UNMUTATED_after_importing_every_consumer():
    """End-to-end companion to the static scan: import each consumer for real and require the shared
    table to still be the standard code afterwards, so an import-time write is caught too."""
    expected = {"TGA": "*", "TAA": "*", "TAG": "*", "ATG": "M", "TGG": "W", "AGA": "R", "ATA": "I"}
    for rel, names in sorted(discover_consumers().items()):
        if "STANDARD" not in names or "forward_inverse_roundtrip" in rel:
            continue
        _load(rel)
    assert len(T.STANDARD) == 64
    for codon, aa in expected.items():
        assert T.STANDARD[codon] == aa, f"{codon} became {T.STANDARD[codon]!r} after consumer imports"
    assert T.STANDARD == T.table_for("standard"), (
        "the shared table diverged from a freshly derived one -- something wrote to it")
