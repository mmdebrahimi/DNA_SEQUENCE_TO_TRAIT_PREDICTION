"""Reference-integrity gate for the committed C. auris FKS1 CDS reference.

The FKS1 echinocandin cell (micafungin/caspofungin) BLASTs `data/fungal_ref/Cauris_FKS1_cds.fna` against a
genome and checks the catalog hotspots (S639F/P/Y, F635del, R1354H). If the committed reference has the
wrong frame or coordinates, every echinocandin call is silently mis-numbered. This test translates the
committed reference and asserts it carries the catalog WT residues at the hotspot positions (F635 / S639 /
R1354) — a frame/coordinate error fails loudly here before any call. Mirrors the HIV RT / SARS-CoV-2 Mpro
reference self-checks.

Offline-safe: pure translation of the committed FASTA (no BLAST, no network). The reference was extracted
from the B8441 (GCF_002759435.1) annotation, where FKS1 is annotated GSC1 (locus B9J08_02922), and the
numbering was verified against the fungal_amr catalog on 2026-07-20.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.fungal_amr import FUNGAL_RESISTANCE_MUTATIONS  # noqa: E402
from scripts.fungal_erg11_caller import _CODON  # noqa: E402

_REF = Path(__file__).resolve().parent.parent / "data" / "fungal_ref" / "Cauris_FKS1_cds.fna"


def _read_cds(path: Path) -> str:
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln and not ln.startswith(">")]
    return "".join(lines).upper()


def _translate(cds: str) -> str:
    return "".join(_CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def test_fks1_reference_exists_and_in_frame():
    assert _REF.exists(), f"committed FKS1 reference missing at {_REF}"
    cds = _read_cds(_REF)
    assert cds.startswith("ATG"), f"FKS1 reference not in-frame (starts {cds[:3]})"
    assert len(cds) % 3 == 0, f"FKS1 reference length {len(cds)} not a multiple of 3"


def test_fks1_reference_wt_residues_match_catalog_numbering():
    # Catalog hotspots: S639F/P/Y (HS1) + F635del + R1354H (HS2). The reference must carry the WT
    # residue at each position, else the caller mis-numbers every echinocandin call.
    prot = _translate(_read_cds(_REF))
    assert len(prot) >= 1354, f"FKS1 protein only {len(prot)} aa; need >=1354 for R1354"
    assert prot[634] == "F", f"pos 635 = {prot[634]!r}, expected F (catalog F635del)"
    assert prot[638] == "S", f"pos 639 = {prot[638]!r}, expected S (catalog S639F/P/Y)"
    assert prot[1353] == "R", f"pos 1354 = {prot[1353]!r}, expected R (catalog R1354H)"


def test_catalog_hotspot_positions_are_covered_by_reference():
    # Every FKS1 substitution in the echinocandin catalog references a position within the reference length.
    prot_len = len(_translate(_read_cds(_REF)))
    for drug in ("micafungin", "caspofungin"):
        for sub in FUNGAL_RESISTANCE_MUTATIONS[drug]["FKS1"]:
            # sub forms: S639F / F635del / R1354H
            digits = "".join(c for c in sub if c.isdigit())
            pos = int(digits)
            assert pos <= prot_len, f"{drug} FKS1 hotspot {sub} at pos {pos} > reference len {prot_len}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
