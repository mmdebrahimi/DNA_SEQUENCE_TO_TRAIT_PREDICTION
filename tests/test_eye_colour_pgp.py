"""Offline pins for the PGP replication harness parsers (no network).

Covers the profile download-link regex (on a real HTML snippet captured 2026-06-30), the Basic-Phenotypes
survey eye-colour parser, and the DTC-line SNP extractor. The real-surface run is scripts/eye_colour_pgp_validate.py --limit.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.eye_colour_pgp_validate import _DL_RE, _snps_from_lines, survey_eye_labels  # noqa: E402

# real fragment from https://my.pgp-hms.org/profile/hu826751 (2026-06-30)
_HTML = (
    'x<td data-summarize-as="name">genome_23andme_full_20141029132101</td> '
    '<td data-summarize-as="size"> <a href="/user_file/download/1285" rel="nofollow">Download</a>'
    '<td data-summarize-as="name">dna-data-2016-06-25.zip</td> '
    '<td data-summarize-as="size"> <a href="/user_file/download/2009" rel="nofollow">Download</a>'
)


def test_profile_link_regex():
    found = _DL_RE.findall(_HTML)
    assert ("genome_23andme_full_20141029132101", "1285") in found
    assert ("dna-data-2016-06-25.zip", "2009") in found


def test_survey_eye_labels():
    csv = (b"Participant,Timestamp,x,1.1 - Blood Type,2.3 - Left Eye Color - Text Description,"
           b"2.4 - Right Eye Color - Text Description\n"
           b"hu1,t,z,O,blue,blue\n"
           b"hu2,t,z,A,Dark Brown,brown\n"
           b"hu3,t,z,B,green,green\n"          # -> other
           b"hu4,t,z,O,,brown\n")             # left blank -> falls back to right
    labels = survey_eye_labels(csv)
    assert labels["hu1"] == "blue"
    assert labels["hu2"] == "brown"
    assert labels["hu3"] == "other"
    assert labels["hu4"] == "brown"


def test_snps_from_lines_23andme_and_ancestry():
    andme = ["# 23andMe", "rsid\tchromosome\tposition\tgenotype",
             "rs12913832\t15\t28365618\tGG", "rs1800407\t15\t1\tCC"]
    got = _snps_from_lines(andme, {"rs12913832", "rs1800407"})
    assert got["rs12913832"] == "GG" and got["rs1800407"] == "CC"
    ancestry = ["#AncestryDNA", "rsid\tchromosome\tposition\tallele1\tallele2",
                "rs12913832\t15\t28365618\tA\tG"]
    assert _snps_from_lines(ancestry, {"rs12913832"})["rs12913832"] == "AG"
