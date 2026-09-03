"""Test the gentamicin `rmt` rescue's specificity on BV-BRC — an archive independent of NCBI-PD.

WHY BV-BRC AND NOT SOMETHING ELSE. The 2026-09-02 hunt exhausted NCBI Pathogen Detection: 60 S-labelled
`rmt` carriers, all from one BioProject, killed as a LABEL_ARTIFACT by an aac(3) control. The obvious
alternatives fail on structure rather than effort (see wiki/rmt_independent_archive_search_2026-09-03.md):

  - `rmt` surveillance studies ascertain isolates ON high-level aminoglycoside resistance (e.g. "MIC >256
    to both amikacin and gentamicin, then PCR"), so by construction they hold ZERO susceptible carriers.
  - NARMS uploads its genomes to NCBI-PD weekly -- a feeder to the archive already exhausted.
  - ATLAS/Vivli has ~634k isolates with MICs but NO paired sequence, so the genotype side is unavailable.

BV-BRC is the one candidate that is genuinely different on BOTH axes:
  - PHENOTYPE ingestion: BioSample antibiograms PLUS ~300 publications curated by hand, carrying `pmid`,
    `laboratory_typing_method` and `testing_standard` per record. The publication-curated part is content
    NCBI-PD's `AST_phenotypes` field does not contain.
  - GENOTYPE calling: `sp_gene` is **CARD/BLAT**, not AMRFinder. So the carrier call is made by a
    different tool with different thresholds from both our pipeline and PD's -- a genuinely independent
    determination, not the same call re-served.

OVERLAP IS MEASURED, NOT ASSUMED. Both archives ultimately draw on public assemblies, so this resolves
genome_id -> assembly accession and reports how many carriers are NEW relative to the PD sweep. An
"independent archive" that returned the same isolates would be independent in name only.

Filters to `evidence = Laboratory Method` throughout: BV-BRC also ships ML-PREDICTED phenotypes, and
scoring a rule against a model's output is the circular-label gate (G1) firing.

Network-only. Writes wiki/gentamicin_rmt_bvbrc_hunt.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.eval.amr_rules import DRUG_RULE  # noqa: E402

API = "https://www.bv-brc.org/api"
RESCUE = DRUG_RULE["gentamicin"]["symbol_rescue"]
RESCUE_RE = re.compile(RESCUE)
# armA is EXCLUDED on purpose: AMRFinder files it under Subclass GENTAMICIN, so the frozen rule already
# counts it. Including it here would credit the rescue with carriers it never needed to rescue.
FAMILY_RE = re.compile(r"^(rmt[A-Z]\d*|npm[A-Z]\d*)$", re.IGNORECASE)


class BvbrcUnavailable(RuntimeError):
    """BV-BRC answers an outage with HTTP 200 wrapping a 503 envelope -- never read that as an empty
    result. This repo has been bitten by exactly that before."""


def fetch(collection: str, query: str, limit: int = 25000, offset: int = 0) -> list[dict]:
    url = f"{API}/{collection}/?{query}&limit({limit},{offset})"
    # BV-BRC 403s the default python-urllib User-Agent while serving curl fine -- an identifying UA is
    # required, and without it the failure looks like an access problem rather than a header problem.
    req = urllib.request.Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "dna_decode/0.13 (research; genotype-phenotype validation)",
    })
    for attempt in range(3):
        try:
            raw = urllib.request.urlopen(req, timeout=300).read().decode("utf8", "replace")
        except Exception as e:
            if attempt == 2:
                raise BvbrcUnavailable(f"{type(e).__name__}: {e}") from e
            time.sleep(5)
            continue
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as e:
            raise BvbrcUnavailable(f"non-JSON response: {raw[:200]}") from e
        if isinstance(data, dict):                    # the 200-wrapping-503 shape
            raise BvbrcUnavailable(f"error envelope: {str(data)[:200]}")
        return data
    raise BvbrcUnavailable("exhausted retries")


def page_all(collection: str, query: str, page: int = 25000, cap: int = 400000) -> list[dict]:
    out, offset = [], 0
    while offset < cap:
        batch = fetch(collection, query, limit=page, offset=offset)
        out.extend(batch)
        if len(batch) < page:
            break
        offset += page
        print(f"    ...{len(out)} rows", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json")
    a = ap.parse_args()

    print(f"deployed rescue (imported): {RESCUE}\n")

    # 1. Every genome BV-BRC calls as carrying a 16S-RMTase family member, via CARD.
    print("fetching RMTase carriers from sp_gene (CARD)...", flush=True)
    carriers: dict[str, set[str]] = {}
    for stem in ("rmtA", "rmtB", "rmtC", "rmtD", "rmtE", "rmtF", "rmtG", "rmtH", "rmtI",
                 "npmA", "npmB", "npmC"):
        try:
            rows = page_all("sp_gene", f"eq(gene,{stem})&select(genome_id,gene)", page=25000)
        except BvbrcUnavailable as e:
            print(f"  {stem}: UNAVAILABLE ({e})")
            return 2
        for r in rows:
            g = str(r.get("gene", ""))
            if FAMILY_RE.match(g):
                carriers.setdefault(str(r["genome_id"]), set()).add(g)
        if rows:
            print(f"  {stem:6} {len(rows):5d} rows", flush=True)
    print(f"  -> {len(carriers)} genomes carry an RMTase family member\n", flush=True)

    # 2. Every MEASURED gentamicin phenotype. Laboratory Method only -- BV-BRC also ships ML predictions,
    #    and scoring a rule against a model's output is the circular-label gate firing.
    print(f"fetching MEASURED {a.drug} AST (evidence=Laboratory Method)...", flush=True)
    q = (f"and(eq(antibiotic,{a.drug}),eq(evidence,Laboratory%20Method))"
         "&select(genome_id,genome_name,resistant_phenotype,measurement,measurement_value,"
         "measurement_sign,laboratory_typing_method,testing_standard,testing_standard_year,pmid,taxon_id)")
    try:
        ast = page_all("genome_amr", q)
    except BvbrcUnavailable as e:
        print(f"  UNAVAILABLE ({e})")
        return 2
    print(f"  -> {len(ast)} measured {a.drug} records\n", flush=True)

    by_genome: dict[str, dict] = {}
    for r in ast:
        gid = str(r.get("genome_id"))
        ph = str(r.get("resistant_phenotype") or "").strip()
        if ph in ("Resistant", "Susceptible", "Intermediate"):
            by_genome.setdefault(gid, r)

    hits = []
    for gid, genes in carriers.items():
        rec = by_genome.get(gid)
        if not rec:
            continue
        hits.append({"genome_id": gid, "genome_name": rec.get("genome_name"),
                     "rmt_genes": sorted(genes),
                     "rescued_by_deployed_rule": any(RESCUE_RE.match(g) for g in genes),
                     "phenotype": rec.get("resistant_phenotype"),
                     "measurement": rec.get("measurement"),
                     "measurement_value": rec.get("measurement_value"),
                     "laboratory_typing_method": rec.get("laboratory_typing_method"),
                     "testing_standard": rec.get("testing_standard"),
                     "testing_standard_year": rec.get("testing_standard_year"),
                     "pmid": rec.get("pmid"), "taxon_id": rec.get("taxon_id")})

    counts = {p: sum(1 for h in hits if h["phenotype"] == p)
              for p in ("Resistant", "Susceptible", "Intermediate")}
    print(f"RMTase carriers with a MEASURED {a.drug} phenotype: {len(hits)}")
    print(f"  R={counts['Resistant']}  S={counts['Susceptible']}  I={counts['Intermediate']}")
    n = counts["Resistant"] + counts["Susceptible"]
    if n:
        print(f"  PPV(rmt -> R) = {counts['Resistant']}/{n} = {counts['Resistant']/n:.4f}")

    s_hits = [h for h in hits if h["phenotype"] == "Susceptible"]
    if s_hits:
        print(f"\n  SUSCEPTIBLE carriers ({len(s_hits)}) -- the sought counter-examples:")
        for h in s_hits[:25]:
            print(f"    {h['genome_id']:16} {str(h['genome_name'])[:44]:46} {','.join(h['rmt_genes']):12} "
                  f"MIC={h['measurement'] or '-'} pmid={h['pmid']}")

    # 3. Overlap with the PD sweep -- an "independent archive" returning the same isolates is not one.
    overlap = {"checked": False}
    pd_art = ROOT / "wiki" / "gentamicin_rmt_specificity_hunt.json"
    if pd_art.is_file() and hits:
        pd = json.loads(pd_art.read_text(encoding="utf-8"))
        pd_acc = {r["acc"] for r in pd["rmt_R_records"] + pd["rmt_S_records"] if r["acc"].startswith("GC")}
        gids = [h["genome_id"] for h in hits]
        accs: dict[str, str] = {}
        for i in range(0, len(gids), 200):
            chunk = ",".join(gids[i:i + 200])
            try:
                for g in fetch("genome", f"in(genome_id,({chunk}))&select(genome_id,assembly_accession)"):
                    if g.get("assembly_accession"):
                        accs[str(g["genome_id"])] = g["assembly_accession"]
            except BvbrcUnavailable as e:
                overlap["error"] = str(e)
                break
        shared = {gid for gid, acc in accs.items() if acc in pd_acc}
        overlap = {"checked": True, "n_with_assembly_accession": len(accs),
                   "n_shared_with_pd_sweep": len(shared),
                   "n_new_relative_to_pd": len(hits) - len(shared),
                   "shared_genome_ids": sorted(shared)[:50]}
        print(f"\noverlap with the PD sweep: {len(shared)} of {len(accs)} resolvable carriers shared "
              f"-> {len(hits) - len(shared)} NEW")

    out = {"schema": "gentamicin-rmt-bvbrc-hunt-v1", "drug": a.drug, "archive": "BV-BRC",
           "deployed_rescue": RESCUE,
           "genotype_caller": "BV-BRC sp_gene (CARD/BLAT) -- a DIFFERENT tool from our AMRFinder and "
                              "from NCBI-PD's AMR_genotypes",
           "phenotype_filter": "evidence = Laboratory Method (excludes BV-BRC's ML-predicted phenotypes, "
                               "which would be a circular label)",
           "n_rmtase_carrier_genomes": len(carriers), "n_measured_ast_records": len(ast),
           "n_carriers_with_measured_phenotype": len(hits), "counts": counts,
           "ppv_rmt_to_R": (counts["Resistant"] / n) if n else None,
           "susceptible_carriers": s_hits, "all_hits": hits, "pd_overlap": overlap,
           "honest_limits": [
               "The CARRIER call is CARD/BLAT, not AMRFinder. That makes it independent of our pipeline "
               "but it is still a tool-derived feature, and CARD's thresholds differ -- a carrier here is "
               "not guaranteed to be a carrier under AMRFinder.",
               "BV-BRC and NCBI-PD both ultimately draw on public assemblies, so overlap is expected; "
               "the measured overlap count is what makes the independence claim checkable.",
               "Only `evidence = Laboratory Method` rows are used. BV-BRC also ships ML-predicted "
               "phenotypes; scoring a deterministic rule against those would be circular.",
               "A susceptible carrier found here still needs the same scrutiny the PD hunt applied: "
               "source concentration, and an aac(3) control on the submitting study's own labels.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
