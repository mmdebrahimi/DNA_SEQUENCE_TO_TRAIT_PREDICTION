"""`dna-pgx --all` — the one-command 14-gene decode + integrated interpretations from a GRCh38 VCF."""
import io
import json
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.pgx import cli  # noqa: E402
from dna_decode.pgx import (  # noqa: E402
    cyp2b6_catalog, cyp2c8_catalog, cyp2c19_catalog, cyp2c9_catalog, cyp2d6_catalog,
    cyp3a5_catalog, dpyd_catalog, nudt15_catalog, tpmt_catalog, ugt1a1_catalog,
    vkorc1, slco1b1, cyp4f2, abcg2,
)

_DIPLO_CATS = [cyp2b6_catalog, cyp2c8_catalog, cyp2c19_catalog, cyp2c9_catalog, cyp2d6_catalog,
               cyp3a5_catalog, dpyd_catalog, nudt15_catalog, tpmt_catalog, ugt1a1_catalog]
_SNP_MODS = [vkorc1, slco1b1, cyp4f2, abcg2]


def _full_reference_vcf() -> Path:
    """Build a VCF carrying a 0/0 (reference) record at EVERY defining position of all 14 genes,
    so the full decode resolves to *1/*1 / reference everywhere (no absent-position None)."""
    rows = {}
    for cat in _DIPLO_CATS:
        for d in cat.CORE_DEFINING:
            rows[(d.chrom, d.pos)] = (d.rsid, d.ref, d.alt)
        for s in getattr(cat, "SENTINELS", []):
            rows[(s.chrom, s.pos)] = (getattr(s, "rsid", "."), s.ref, s.alt)
    for m in _SNP_MODS:
        rows[(m.CHROM, m.POS)] = (m.RSID, m.REF, m.ALT)
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    for (chrom, pos), (rsid, ref, alt) in sorted(rows.items()):
        lines.append(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/0")
    p = Path(tempfile.mktemp(suffix=".vcf"))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _run_all_json(vcf: Path) -> dict:
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = cli.main([str(vcf), "--all", "--json-only", "--sample-id", "T"])
    assert rc == 0
    return json.loads(buf.getvalue())


def test_all_decodes_14_genes_reference_as_normal():
    out = _run_all_json(_full_reference_vcf())
    assert out["n_genes"] == 14
    # every diplotype gene reference -> *1/*1 Normal Metabolizer
    for g in ("cyp2c19", "cyp2c9", "tpmt", "nudt15", "dpyd", "ugt1a1"):
        assert out["results"][g]["diplotype"] == "*1/*1", g
        assert out["results"][g]["phenotype"] == "Normal Metabolizer", g
    # SNP genes reference -> Normal
    assert out["results"]["vkorc1"]["sensitivity"] == "Normal sensitivity"
    assert out["results"]["cyp4f2"]["function"] == "Normal Function"
    assert out["results"]["abcg2"]["function"] == "Normal Function"


def test_all_interpretations_resolve_not_indeterminate():
    out = _run_all_json(_full_reference_vcf())
    i = out["interpretations"]
    assert i["warfarin"]["status"] == "ok"
    assert i["warfarin"]["dose_direction"] == "standard_dose_requirement"
    assert i["statins"]["status"] == "ok"
    assert i["statins"]["simvastatin_myopathy_risk"] == "typical_risk"
    assert i["thiopurines"]["status"] == "ok"
    assert i["thiopurines"]["toxicity_risk"] == "normal_risk"


def test_all_variant_changes_interpretation():
    """Flip CYP4F2 to *3/*3 (homozygous ALT) -> warfarin dose direction shifts higher."""
    vcf = _full_reference_vcf()
    text = vcf.read_text(encoding="utf-8")
    text = text.replace(f"{cyp4f2.CHROM}\t{cyp4f2.POS}\t{cyp4f2.RSID}\t{cyp4f2.REF}\t{cyp4f2.ALT}\t.\tPASS\t.\tGT\t0/0",
                        f"{cyp4f2.CHROM}\t{cyp4f2.POS}\t{cyp4f2.RSID}\t{cyp4f2.REF}\t{cyp4f2.ALT}\t.\tPASS\t.\tGT\t1/1")
    vcf.write_text(text, encoding="utf-8")
    out = _run_all_json(vcf)
    assert out["results"]["cyp4f2"]["function"] == "Reduced Function"
    assert out["interpretations"]["warfarin"]["dose_direction"] == "higher_dose_requirement"


def test_all_caveat_not_clinical():
    out = _run_all_json(_full_reference_vcf())
    assert "NOT a clinical tool" in out["caveat"]
