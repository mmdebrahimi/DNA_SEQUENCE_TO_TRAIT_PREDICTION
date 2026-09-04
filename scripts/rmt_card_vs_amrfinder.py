"""Does AMRFinder agree with CARD that the 67 BV-BRC susceptible isolates carry `rmt`?

WHY THIS IS THE DECISIVE FOLLOW-UP. The BV-BRC finding (2026-09-03) rests on carrier calls from
**CARD/BLAT** via BV-BRC's `sp_gene`. The DEPLOYED rule consumes **AMRFinder**. If AMRFinder does not
call these isolates as `rmt` carriers, the rule never fires on them and they are not counter-examples to
it at all -- the discordance would be between two gene callers, not between our rule and reality.

THE FREE ROUTE, instead of re-running AMRFinder on 67 downloaded genomes (~1.7 h of Docker):
NCBI Pathogen Detection's `AMR_genotypes` field IS an AMRFinder call, made by NCBI on the same isolates.
So looking those isolates up in PD metadata reads AMRFinder's verdict directly, at zero compute cost.

THE JOIN KEY IS BIOSAMPLE, NOT ASSEMBLY ACCESSION. Only 3 of the 67 carry an `assembly_accession` in
BV-BRC; 64 carry a BioSample (SAMEA...) plus an SRA run. PD metadata carries `biosample_acc`, so that is
the key that actually connects them. Joining on assembly accession would have resolved 3 isolates and
looked like an answer.

Three outcomes, all informative and none of them a shrug:
  - AMRFinder AGREES on most      -> the counter-examples stand against the deployed rule.
  - AMRFinder DISAGREES on most   -> the finding is a CARD-vs-AMRFinder caller discordance; the deployed
                                     rule would never have fired, and its specificity is untouched.
  - The isolates are ABSENT from PD -> unresolvable by this route; the real AMRFinder run is required,
                                     and that is reported as unresolved rather than as agreement.

Network-only. Writes wiki/rmt_card_vs_amrfinder.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.data.pd_ast import ast_label_for  # noqa: E402
from gentamicin_rmt_specificity_hunt import (  # noqa: E402
    RESCUE_RE, latest_metadata_url, parse_amr_genotypes,
)

UA = {"Accept": "application/json",
      "User-Agent": "dna_decode/0.13 (research; genotype-phenotype validation)"}
FAMILY_RE = re.compile(r"^(rmt[A-Z]\d*|npm[A-Z]\d*)$", re.IGNORECASE)
ARMA_RE = re.compile(r"^armA\d*$", re.IGNORECASE)


def bvbrc(coll: str, query: str, limit: int = 25000) -> list[dict]:
    url = f"https://www.bv-brc.org/api/{coll}/?{query}&limit({limit},0)"
    raw = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300).read()
    data = json.loads(raw.decode("utf8", "replace"))
    if isinstance(data, dict):
        raise RuntimeError(f"BV-BRC error envelope: {str(data)[:200]}")
    return data


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups", default="Klebsiella,Escherichia_coli_Shigella,Salmonella")
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" / "rmt_card_vs_amrfinder.json")
    a = ap.parse_args()

    hunt = json.loads((ROOT / "wiki" / "gentamicin_rmt_bvbrc_hunt.json").read_text(encoding="utf-8"))
    S = [h for h in hunt["all_hits"] if h["phenotype"] == "Susceptible"]
    gids = [h["genome_id"] for h in S]
    print(f"{len(S)} BV-BRC susceptible rmt carriers (CARD calls)\n")

    # Resolve every identifier BV-BRC has -- BioSample is the one that reaches PD.
    ids: dict[str, dict] = {}
    for i in range(0, len(gids), 200):
        chunk = ",".join(gids[i:i + 200])
        for g in bvbrc("genome", f"in(genome_id,({chunk}))&select(genome_id,assembly_accession,"
                                 "biosample_accession,sra_accession,genome_name)"):
            ids[str(g["genome_id"])] = g
    biosamples = {str(g.get("biosample_accession")): gid for gid, g in ids.items()
                  if g.get("biosample_accession")}
    print(f"  resolvable by BioSample: {len(biosamples)}   "
          f"by assembly accession: {sum(1 for g in ids.values() if g.get('assembly_accession'))}")

    # Read AMRFinder's own verdict from PD metadata for those BioSamples.
    found: dict[str, dict] = {}
    for grp in [x.strip() for x in a.groups.split(",") if x.strip()]:
        try:
            r = urllib.request.urlopen(latest_metadata_url(grp), timeout=600)
        except Exception as e:
            print(f"  {grp}: fetch failed ({type(e).__name__}) -- skipped")
            continue
        cols = r.readline().decode("utf8", "replace").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        bi, gi, ai = idx.get("biosample_acc"), idx.get("AMR_genotypes"), idx.get("AST_phenotypes")
        if bi is None or gi is None:
            print(f"  {grp}: missing biosample_acc/AMR_genotypes -- skipped")
            continue
        n = 0
        for line in r:
            f = line.decode("utf8", "replace").rstrip("\n").split("\t")
            if len(f) <= max(bi, gi):
                continue
            bs = f[bi].strip()
            if bs not in biosamples:
                continue
            syms = parse_amr_genotypes(f[gi])
            found[bs] = {
                "group": grp,
                "amrfinder_rmt_family": sorted(s for s in syms if FAMILY_RE.match(s)),
                "amrfinder_would_be_rescued": any(RESCUE_RE.match(s) for s in syms),
                "amrfinder_has_armA": any(ARMA_RE.match(s) for s in syms),
                "amrfinder_has_aac3": any(s.startswith("aac(3)") for s in syms),
                "pd_gentamicin_label": ast_label_for(f[ai], "gentamicin") if ai is not None else None,
            }
            n += 1
        print(f"  {grp:26} matched {n}", flush=True)

    n_found = len(found)
    agree = [b for b, v in found.items() if v["amrfinder_rmt_family"]]
    rescued = [b for b, v in found.items() if v["amrfinder_would_be_rescued"]]
    print(f"\nresolved in PD (AMRFinder): {n_found} of {len(biosamples)} BioSamples")
    if n_found:
        print(f"  AMRFinder ALSO calls an RMTase : {len(agree)}/{n_found} ({len(agree)/n_found:.0%})")
        print(f"  would be RESCUED by the rule   : {len(rescued)}/{n_found}")
        lab = [v["pd_gentamicin_label"] for v in found.values() if v["pd_gentamicin_label"]]
        print(f"  of these, PD also has a gentamicin label for {len(lab)}: "
              f"{ {x: lab.count(x) for x in set(lab)} }")

    if n_found == 0:
        verdict = "UNRESOLVED_NOT_IN_PD"
        why = ("none of the susceptible carriers resolve in PD metadata, so AMRFinder's verdict cannot be "
               "read this way. The real AMRFinder run on the downloaded genomes is required; this is "
               "UNRESOLVED, not agreement.")
    elif len(rescued) / n_found >= 0.8:
        verdict = "CALLERS_AGREE"
        why = (f"AMRFinder independently calls an RMTase in {len(agree)}/{n_found} of them and the "
               f"deployed rule would fire on {len(rescued)}/{n_found} -- so these ARE counter-examples "
               "to the deployed rule, not a caller artefact.")
    elif len(rescued) / n_found <= 0.2:
        verdict = "CALLER_DISCORDANCE"
        why = (f"AMRFinder calls an RMTase in only {len(agree)}/{n_found} -- the deployed rule would "
               "mostly never fire on these isolates, so they are a CARD-vs-AMRFinder discordance rather "
               "than counter-examples to the rule.")
    else:
        verdict = "MIXED"
        why = (f"AMRFinder agrees on {len(agree)}/{n_found}; neither reading is clean and the subset "
               "where the callers agree is the only part that bears on the deployed rule.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "rmt-card-vs-amrfinder-v1",
           "question": "do CARD (BV-BRC sp_gene) and AMRFinder (NCBI-PD AMR_genotypes) agree that the 67 "
                       "gentamicin-susceptible BV-BRC isolates carry an RMTase?",
           "n_susceptible_carriers": len(S), "n_with_biosample": len(biosamples),
           "n_resolved_in_pd": n_found,
           "n_amrfinder_agrees": len(agree), "n_would_be_rescued": len(rescued),
           "per_biosample": found, "verdict": verdict, "why": why,
           "honest_limits": [
               "PD's AMR_genotypes is NCBI's AMRFinder run, not ours -- same tool, possibly a different "
               "version and database than the deployed pipeline uses.",
               "Isolates absent from PD are UNRESOLVED by this route, never counted as agreement.",
               "This tests CARRIER-CALL agreement only. It does not re-examine the phenotypes, which the "
               "aac(3) control already cleared (wiki/gentamicin_rmt_bvbrc_control.json).",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
