"""Per-organism essentiality validation: pure SGD parser + registry + real yeast smoke (slow)."""
from __future__ import annotations

import pytest

from dna_decode.fba.essentiality_labels import (
    ESSENTIALITY_LABEL_SOURCES,
    LABEL_WALLED,
    parse_essential,
    parse_sgd_essential,
)


def _row(orf, mutant_type, phenotype):
    # SGD phenotype_data.tab: col0=ORF, col6=mutant_type, col9=phenotype (14 cols)
    c = [orf, "ORF", "GENE", "S000000001", "ref", "classical genetics",
         mutant_type, "allele", "S288C", phenotype, "", "", "", ""]
    return "\t".join(c)


_TAB = "\n".join([
    _row("YAL001C", "null", "inviable"),                 # essential
    _row("YBR002W", "null", "viable"),                   # viable
    _row("YCL003C", "overexpression", "inviable"),       # NOT null -> not counted essential
    _row("YDL004W", "null", "inviable ascus"),           # 'inviable' substring -> essential
    "not\ta\tvalid\trow",                                # junk
])


def test_parse_sgd_essential():
    ess = parse_sgd_essential(_TAB)
    assert ess == {"YAL001C", "YDL004W"}                 # null + inviable only
    assert "YBR002W" not in ess and "YCL003C" not in ess


def test_parse_essential_dispatch():
    assert parse_essential("sgd", _TAB) == {"YAL001C", "YDL004W"}
    with pytest.raises(ValueError):
        parse_essential("nonsense-kind", _TAB)


def test_registry_shapes():
    # yeast is a scored source; S. aureus / P. aeruginosa are honest walls (disjoint sets)
    assert ESSENTIALITY_LABEL_SOURCES["yeast"][0] == "sgd"
    assert "saureus" in LABEL_WALLED and "paeruginosa" in LABEL_WALLED
    assert not (set(ESSENTIALITY_LABEL_SOURCES) & set(LABEL_WALLED))


@pytest.mark.slow
def test_yeast_essentiality_validation_real():
    pytest.importorskip("cobra")
    import urllib.request

    from dna_decode.fba.keio import confusion, metrics_from_confusion
    from dna_decode.fba.model import gene_essentiality, load_model

    url = ESSENTIALITY_LABEL_SOURCES["yeast"][1]
    txt = urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "test"}), timeout=90
    ).read().decode("utf-8", "replace")
    ess = parse_sgd_essential(txt)
    assert len(ess) > 800                                # yeast essentialome ~1100-1200

    m = load_model(organism="yeast")
    pred = {g: v[1] for g, v in gene_essentiality(m).items()}
    exp = {g: (g in ess) for g in pred}
    cm = confusion(exp, pred)
    met = metrics_from_confusion(cm)
    # the number is WEAK (a real finding) -- assert only that a real, non-degenerate metric came out
    assert cm["n"] > 800                                 # ~905 model genes scored
    assert 0.0 < met["accuracy"] <= 1.0
    assert met["mcc"] > 0.0                              # positive (better than chance), even if weak
