"""Offline test for the MaveDB adapter's pure parsers (no torch/network)."""
import importlib.util, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
def _load():
    spec = importlib.util.spec_from_file_location("mv", ROOT / "scripts" / "mavedb_cpu_smoke.py")
    m = importlib.util.module_from_spec(spec); sys.modules["mv"]=m; spec.loader.exec_module(m); return m
def test_parse_hgvs_pro():
    M=_load()
    assert M.parse_hgvs_pro("p.Met1Ala")==("M",1,"A")
    assert M.parse_hgvs_pro("p.Val39Gly")==("V",39,"G")
    assert M.parse_hgvs_pro("p.Leu2=") is None          # synonymous
    assert M.parse_hgvs_pro("p.Lys4Ter") is None         # nonsense (Ter not an AA)
    assert M.parse_hgvs_pro("p.Met1Met") is None         # wt==mut
    assert M.parse_hgvs_pro("p.Gly10del") is None        # indel
def test_translate():
    M=_load()
    assert M.translate("ATGGTTTAA")=="MV"                # Met-Val-STOP
    assert M.translate("atgaaa")=="MK"
def test_build_dms():
    M=_load()
    csv_txt="hgvs_pro,score\np.Met1Ala,1.5\np.Leu2=,9.9\np.Val3Gly,-0.4\np.Lys4Ter,0.0\nbad,x\n"
    assert M.build_dms(csv_txt)=={"M1A":1.5,"V3G":-0.4}
