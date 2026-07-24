#!/usr/bin/env python3
"""Temp fetch script: BASEL E. coli phage collection genomes MZ501046.1..MZ501113.1.

Fetches FASTA + GenBank taxonomy, writes a manifest TSV. Polite to NCBI E-utilities.
"""
import os
import re
import sys
import time
import urllib.request
import urllib.error

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ROOT = os.path.join(os.path.dirname(__file__), "..", "data", "phage_ref")
BASEL_DIR = os.path.join(ROOT, "basel")
MANIFEST = os.path.join(ROOT, "basel_manifest.tsv")

ACCS = [f"MZ50{n}.1" for n in range(1046, 1114)]  # MZ501046.1 .. MZ501113.1
assert ACCS[0] == "MZ501046.1" and ACCS[-1] == "MZ501113.1" and len(ACCS) == 68, ACCS[:2]


def _get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent": "dna_decode-basel-fetch/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def fetch_text(acc, rettype):
    url = f"{BASE}?db=nucleotide&id={acc}&rettype={rettype}&retmode=text"
    for attempt in range(2):
        try:
            txt = _get(url)
            if txt and txt.strip():
                return txt
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
            sys.stderr.write(f"  {acc} {rettype} attempt {attempt} err: {e}\n")
        time.sleep(1.0)
    return ""


def parse_defline_name(fasta):
    """First line like '>MZ501046.1 Escherichia phage AdolfPortmann, complete genome'."""
    line = fasta.splitlines()[0] if fasta else ""
    line = line.lstrip(">").strip()
    # drop accession token
    parts = line.split(None, 1)
    desc = parts[1] if len(parts) > 1 else ""
    # trim trailing ', complete genome' etc.
    desc = re.split(r",\s*(?:complete|partial)\b", desc)[0].strip()
    return desc


def parse_gb_taxonomy(gb):
    """Extract ORGANISM line + taxonomy lineage from a GenBank flat file header.

    ORGANISM block:
      ORGANISM  Escherichia phage AdolfPortmann
                Viruses; ...; Caudoviricetes; ...; Straboviridae; ...; Tequatrovirus.
    """
    organism = ""
    lineage_tokens = []
    lines = gb.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("  ORGANISM"):
            organism = ln.split("ORGANISM", 1)[1].strip()
            j = i + 1
            lin = []
            while j < len(lines) and lines[j].startswith("            "):
                lin.append(lines[j].strip())
                j += 1
            lineage = " ".join(lin)
            lineage = lineage.rstrip(".")
            lineage_tokens = [t.strip() for t in lineage.split(";") if t.strip()]
            break
    return organism, lineage_tokens


def classify(lineage_tokens):
    """Return (genus, family). Virus ranks: family ends -viridae, genus is typically the
    last token in the lineage (ICTV lineage terminates at genus). subfamily ends -virinae."""
    family = "unclassified"
    genus = "unclassified"
    for t in lineage_tokens:
        if t.endswith("viridae"):
            family = t
    # ICTV lineage terminates at genus then a binomial species ("Genus epithet").
    # Genus is the last single-word token that isn't a family/higher rank; if only a
    # binomial species token is present, its first word is the genus.
    def _is_higher(t):
        return (t.endswith("viridae") or t.endswith("virinae") or
                t.endswith("viricetes") or t.endswith("viricota") or
                t.endswith("virales") or t.endswith("virae") or
                t.endswith("viria") or t == "Viruses")
    for t in lineage_tokens:
        if _is_higher(t):
            continue
        if " " in t:  # binomial species -> genus is first word
            genus = t.split()[0]
        else:
            genus = t
    if genus == "unclassified" and lineage_tokens:
        genus = lineage_tokens[-1]  # lowest-rank available
    return genus, family


def main():
    os.makedirs(BASEL_DIR, exist_ok=True)
    rows = []
    ok = 0
    for acc in ACCS:
        sys.stderr.write(f"{acc}...\n")
        fasta = fetch_text(acc, "fasta")
        fna_path = os.path.join(BASEL_DIR, f"{acc}.fna")
        if not fasta or not fasta.lstrip().startswith(">"):
            rows.append((acc, "FETCH_FAILED", "FETCH_FAILED", "FETCH_FAILED"))
            time.sleep(0.4)
            continue
        with open(fna_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(fasta if fasta.endswith("\n") else fasta + "\n")
        ok += 1
        name = parse_defline_name(fasta)
        time.sleep(0.4)
        gb = fetch_text(acc, "gb")
        organism, lineage = parse_gb_taxonomy(gb)
        if organism and not name:
            name = organism
        genus, family = classify(lineage)
        rows.append((acc, name or "unknown", genus, family))
        time.sleep(0.4)

    with open(MANIFEST, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tphage_name\tgenus\tfamily\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    sys.stderr.write(f"\nDONE ok={ok}/68 manifest={MANIFEST}\n")


if __name__ == "__main__":
    main()
