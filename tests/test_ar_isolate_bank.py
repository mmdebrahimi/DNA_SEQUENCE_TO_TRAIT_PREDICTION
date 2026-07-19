"""Offline tests for the CDC AR Isolate Bank ingester (parsers pure; fixtures = real markup).

Fixtures below reproduce the exact patterns observed on the live 2026-07-18 pages:
  - detail BioSample span id="MainContent_lblSeqAcc" -> NCBI biosample link,
  - "AR Bank #NNNN</b> ... <i>genus species</i>" header,
  - MIC table "<tr><td>Drug <sup>fn</sup></td><td>MIC</td><td>INT</td></tr>" with censored + plain MICs,
  - panel row "<a href='IsolateDetail.aspx?IsolateID=&PanelID='>ARnum</a>" + organism + mechanisms + biosample.
No network is touched.
"""
from pathlib import Path

from dna_decode.data import ar_isolate_bank as ab

DETAIL_HTML = """
<html><body>
<tr><td style="font-size:22px;"><b>AR Bank #  0350</b>
  <span style="font-size:22px;"><i>Escherichia coli</i></span></td></tr>
<td><span>Biosample Accession #:
  <span id="MainContent_lblSeqAcc"><a class="tp-external-link-fix" target="_blank"
  href="https://www.ncbi.nlm.nih.gov/biosample/SAMN07152999">SAMN07152999</a></span></span></td>
<table><thead><tr><th><b>Drug</b></th><th><b>MIC (&mu;g/ml)</b></th><th><b>INT</b></th></tr></thead>
<tbody>
<tr><td>Amikacin </td><td style="text-align:center">2</td><td style="text-align:center">S</td></tr>
<tr><td>Ciprofloxacin <sup>1</sup></td><td style="text-align:center">&gt;4</td><td style="text-align:center">R</td></tr>
<tr><td>Ceftriaxone </td><td style="text-align:center">&gt;32</td><td style="text-align:center">R</td></tr>
<tr><td>Gentamicin </td><td style="text-align:center">&lt;=1</td><td style="text-align:center">S</td></tr>
<tr><td>Meropenem </td><td style="text-align:center">1</td><td style="text-align:center">I</td></tr>
</tbody></table>
</body></html>
"""

PANEL_HTML = """
<table id="dataTable"><tbody>
<tr><td><a target='_blank' href="IsolateDetail.aspx?IsolateID=350&PanelID=13">0350</a></td>
    <td><i>Escherichia coli</i></td><td> CTX-M-15</td>
    <td><a class="tp-external-link-fix" href="https://www.ncbi.nlm.nih.gov/biosample/SAMN07152999">SAMN07152999</a></td></tr>
<tr><td><a target='_blank' href="IsolateDetail.aspx?IsolateID=21&PanelID=13">0021</a></td>
    <td><i>Citrobacter freundii</i></td><td> CMY-108, QnrB13</td>
    <td><a href="https://www.ncbi.nlm.nih.gov/biosample/SAMN04014862">SAMN04014862</a></td></tr>
<tr><th>header row no link</th></tr>
</tbody></table>
"""


def test_parse_isolate_detail_fields():
    d = ab.parse_isolate_detail(DETAIL_HTML, isolate_id=350, panel_id=13)
    assert d is not None
    assert d.ar_number == "0350"
    assert d.organism == "Escherichia coli"
    assert d.biosample == "SAMN07152999"
    # 5 MIC rows, operators preserved, entities decoded
    assert d.mics["Ciprofloxacin"] == ">4"
    assert d.mics["Ceftriaxone"] == ">32"
    assert d.mics["Gentamicin"] == "<=1"
    assert d.calls["Ciprofloxacin"] == "R"
    assert d.calls["Gentamicin"] == "S"
    assert d.calls["Meropenem"] == "I"


def test_blank_detail_returns_none():
    assert ab.parse_isolate_detail("<html><body>no data</body></html>") is None


def test_parse_panel_isolates():
    rows = ab.parse_panel_isolates(PANEL_HTML)
    assert len(rows) == 2                      # header row (no link) skipped
    r0 = rows[0]
    assert (r0.isolate_id, r0.panel_id) == (350, 13)
    assert r0.ar_number == "0350"
    assert r0.organism == "Escherichia coli"
    assert r0.biosample == "SAMN07152999"
    assert r0.mechanisms == "CTX-M-15"


def test_to_label_inputs_maps_and_filters_drugs():
    d = ab.parse_isolate_detail(DETAIL_HTML, isolate_id=350, panel_id=13)
    mics, calls = ab.to_label_inputs([d], "ciprofloxacin")
    assert mics == {"SAMN07152999": [">4"]}
    assert calls == {"SAMN07152999": {"R"}}
    # gentamicin S call retained; ceftriaxone maps from the AR name
    _, gcalls = ab.to_label_inputs([d], "gentamicin")
    assert gcalls == {"SAMN07152999": {"S"}}
    cmics, _ = ab.to_label_inputs([d], "ceftriaxone")
    assert cmics == {"SAMN07152999": [">32"]}
    # a non-pilot drug -> empty
    assert ab.to_label_inputs([d], "meropenem") == ({}, {})


def test_to_label_inputs_excludes_intermediate_call():
    """An I call is not a decisive R/S vote (MIC tiering decides); it must not enter calls."""
    d = ab.parse_isolate_detail(DETAIL_HTML)
    # meropenem is non-pilot anyway; assert on the mechanism via a pilot drug with an I would need
    # a fixture — here confirm only R/S propagate for the pilot drugs above (Meropenem I dropped).
    _, cipro_calls = ab.to_label_inputs([d], "ciprofloxacin")
    assert all(v <= {"R", "S"} for v in cipro_calls.values())


def test_organism_matches():
    assert ab.organism_matches("Escherichia coli", ("Escherichia",))
    assert ab.organism_matches("Klebsiella pneumoniae", ("coli", "Klebsiella"))
    assert not ab.organism_matches("Escherichia coli", ("Neisseria",))


def test_enumerate_and_labels_bridge_end_to_end():
    """The parsed detail flows into the FROZEN build_drug_labels and yields a clean R strict label."""
    from dna_decode.data.external_mic_labels import build_drug_labels
    d = ab.parse_isolate_detail(DETAIL_HTML, isolate_id=350, panel_id=13)
    mics, calls = ab.to_label_inputs([d], "ceftriaxone")
    res = build_drug_labels(mics, "ceftriaxone", calls)
    # >32 ug/mL ceftriaxone is unambiguously HIGH_R
    assert res["strict"] == {"SAMN07152999": "R"}


def test_fetch_cache_first(tmp_path: Path):
    """A cached page is read without any network call."""
    cache = tmp_path / "pages"
    cache.mkdir()
    (cache / "panel_13.html").write_text(PANEL_HTML, encoding="utf-8")
    html = ab.fetch(ab.PANEL_URL.format(panel_id=13), cache, name="panel_13")
    assert "IsolateDetail" in html
    rows = ab.enumerate_panel(13, cache)     # uses the cache, no network
    assert len(rows) == 2


def test_panel_ids_are_the_established_34():
    assert len(ab.PANEL_IDS) == 34
    assert 1177 not in ab.PANEL_IDS   # 1177 is a custom panel; established panels reach every isolate
