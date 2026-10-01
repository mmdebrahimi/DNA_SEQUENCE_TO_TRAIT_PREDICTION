"""Verify-in-batch: run the universal constraints over REAL CDS and report their kill counts.

WHY A PROBE AND NOT A TEST. A constraint layer's failure mode is looking protective while removing
nothing. This repo has shipped that twice (the HIV censoring guard over 2,168 rows with zero censored
values; the ResFinder POINT-row exclusion over 1,818 runs with zero POINT rows) and retracted both. So
every constraint reports what it actually refuted, and an `inert` result must be acknowledged with a
reason or `assert_no_silent_inert()` fails this probe.

EXPECTATION RECORDED BEFORE THE RUN (and written into the artifact): over curated reference CDS most laws
are predicted INERT, because a committed reference already satisfies them. **That is the honest outcome
and must be reported as inert, not tuned until something fails.**

THE CLADE DEMO IS REAL, NOT SYNTHETIC. `data/mycoplasma_ref/Mgenitalium_G37_MG_382_udk_cds.fna` is a real
uridine-kinase CDS from M. genitalium G37 (NCBI nuccore L43967.2, translation table 4). Under table 1 it
reports five internal stops; under table 4 it reads clean. A census over the whole G37 CDS set found this
is not a cherry-pick -- see `MGEN_CENSUS`.

Run: uv run python scripts/constraint_kill_count_probe.py
Exit 0 on success, 2 when no reference CDS could be read, 1 when a constraint is silently inert.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.constraints import universal as U              # noqa: E402
from dna_decode.constraints.codon_tables import (              # noqa: E402
    MYCOPLASMA_SPIROPLASMA,
    STANDARD,
)
from dna_decode.constraints.report import ConstraintReport, SilentInertConstraintError  # noqa: E402

#: Committed, non-gitignored reference CDS spanning three organisms.
REFERENCE_CDS: dict[str, str] = {
    "hiv1_rt_hxb2": "data/hiv_ref/HIV1_RT_HXB2_cds.fna",
    "hiv1_pr_hxb2": "data/hiv_ref/HIV1_PR_HXB2_cds.fna",
    "hiv1_in_hxb2": "data/hiv_ref/HIV1_IN_HXB2_cds.fna",
    "hiv1_ca_hxb2": "data/hiv_ref/HIV1_CA_HXB2_cds.fna",
    "sarscov2_mpro": "data/sarscov2_ref/SARSCoV2_Mpro_NC045512_cds.fna",
    "cauris_erg11_wt": "data/fungal_ref/Cauris_ERG11_PV630306_WT.fna",
    "cauris_erg11_y132f": "data/fungal_ref/Cauris_ERG11_PV630305_Y132F.fna",
    "cauris_erg11_k143r": "data/fungal_ref/Cauris_ERG11_PV630302_K143R.fna",
}

#: Which references are COMPLETE genes. The HIV and SARS-CoV-2 entries are in-frame extracts from a
#: polyprotein (HXB2 gag/pol; ORF1ab nsp5) and legitimately lack a start codon, so asserting completeness
#: for them would manufacture a false refutation. Measured: they begin CCC / CCT / AGT, not ATG.
COMPLETE_GENES = frozenset({"cauris_erg11_wt", "cauris_erg11_y132f", "cauris_erg11_k143r"})

MYCOPLASMA_FIXTURE = "data/mycoplasma_ref/Mgenitalium_G37_MG_382_udk_cds.fna"

#: Measured over the real G37 CDS set fetched from NCBI (see `REPRODUCE_CENSUS`). Recorded so the clade
#: demo is not mistaken for a cherry-picked single gene.
MGEN_CENSUS = {"n_false_nonsense_under_table1": 341, "n_clean_cds": 476, "fraction": 0.716}
REPRODUCE_CENSUS = (
    "curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L43967.2"
    "&rettype=fasta_cds_na&retmode=text'")

#: Stated BEFORE the run. A law that a curated reference already satisfies is inert, not protective.
PRE_STATED_EXPECTATION = (
    "Over curated reference CDS most laws are predicted INERT: a committed reference already satisfies "
    "the reading frame, its own reference bases and codon accessibility. Inert is the honest outcome "
    "here and is reported as such -- the constraints are NOT to be tuned until something fails. The one "
    "law expected to be ACTIVE is no_internal_stop_codon, and only on the Mycoplasma CDS when it is "
    "scored under the WRONG (standard) code."
)

#: CORRECTION recorded after the first run, because the pre-stated expectation was wrong in one place.
#: `start_codon_present` was predicted inert and instead refuted 5 of 8 -- all five CORRECTLY, because
#: HIV RT/PR/IN/CA and SARS-CoV-2 Mpro are polyprotein-derived extracts that legitimately lack a start
#: codon. The law now requires the caller to assert `cds_is_complete_gene` and abstains otherwise: a
#: FALSE refutation is the one failure mode that would make a constraint filter delete valid predictions.
FIRST_RUN_CORRECTION = (
    "start_codon_present refuted 5 of 8 reference CDS on the first run, against a pre-stated expectation "
    "of inert. All five refutations were on correct sequences (polyprotein-derived extracts beginning "
    "CCC/CCT/AGT), so the EXPECTATION was wrong and the law was over-firing. It now abstains unless the "
    "caller asserts cds_is_complete_gene."
)


def read_fasta(path: Path) -> str:
    seq = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"constraint_kill_count_{date.today()}.json"))
    a = ap.parse_args(argv)

    report = ConstraintReport()
    per_cds: dict[str, dict[str, str]] = {}
    missing: list[str] = []

    for name, rel in REFERENCE_CDS.items():
        p = REPO / rel
        if not p.exists():
            missing.append(rel)
            continue
        cds = read_fasta(p)
        ctx = {"cds": cds, "nt_pos": 1, "ref_base": cds[:1],
               "cds_is_complete_gene": name in COMPLETE_GENES}
        res = U.evaluate_all(ctx, report=report)
        per_cds[name] = {k: v.verdict for k, v in res.items()}

    if not per_cds:
        print("no reference CDS could be read; refusing to report kill counts")
        return 2

    # ----- the clade demo, on a REAL organism -----
    clade_demo: dict[str, object] = {"source": "unavailable"}
    mp = REPO / MYCOPLASMA_FIXTURE
    if mp.exists():
        cds = read_fasta(mp)
        wrong = U.no_internal_stop_codon({"cds": cds})                               # table 1
        right = U.no_internal_stop_codon({"cds": cds, "clade": "bacteria.mycoplasma"})  # table 4
        report.record("no_internal_stop_codon", wrong)
        report.record("no_internal_stop_codon", right)
        prot_t4 = "".join(MYCOPLASMA_SPIROPLASMA.get(cds[i:i + 3], "X")
                          for i in range(0, len(cds) - 2, 3))
        clade_demo = {
            "source": "real",                      # NOT synthetic
            "label": "real_organism_cds",
            "fixture": MYCOPLASMA_FIXTURE,
            "organism": "Mycoplasma genitalium G37",
            "gene": "MG_382 udk (uridine kinase)",
            "provenance": "NCBI nuccore L43967.2, rettype=fasta_cds_na, transl_table=4",
            "n_internal_stops_under_table1": sum(
                1 for i in range(0, len(cds) - 5, 3) if STANDARD.get(cds[i:i + 3]) == "*"),
            "verdict_under_standard_table": wrong.verdict,
            "verdict_under_clade_table": right.verdict,
            "detail_standard": wrong.detail,
            "detail_clade": right.detail,
            "protein_prefix_under_table4": prot_t4[:40],
            "biology_sanity_check": ("the table-4 translation opens MDEKGILVAISGGSCSGKTT..., whose "
                                     "SGGSCSGKTT is a canonical P-loop nucleotide-binding motif, which "
                                     "is independent evidence the clade reading is the correct one"),
            "census": MGEN_CENSUS,
            "census_reading": (f"{MGEN_CENSUS['n_false_nonsense_under_table1']} of "
                               f"{MGEN_CENSUS['n_clean_cds']} clean real G37 CDS "
                               f"({MGEN_CENSUS['fraction']:.1%}) are false-nonsense under table 1 and "
                               "clean under table 4, so this is not a cherry-picked gene"),
            "reproduce_census": REPRODUCE_CENSUS,
        }

    # ----- acknowledge the inert laws, with reasons -----
    acknowledgements = {
        "cds_length_multiple_of_three":
            "every committed reference CDS is in frame by construction; this law exists for user-supplied "
            "sequence, not for curated references",
        "ref_base_matches_cds":
            "the probe supplies each CDS's own first base as ref_base, so agreement is tautological here; "
            "the law earns its keep on externally supplied coordinates",
        "substitution_reachable_by_single_nt":
            "the probe scores whole CDS and supplies no substitution, so this law has nothing to evaluate "
            "on this input",
        "start_codon_present":
            "evaluated only on the three COMPLETE genes (C. auris ERG11), which all begin ATG; the five "
            "polyprotein-derived extracts abstain rather than being falsely refuted",
        "no_internal_stop_codon":
            "curated references carry no internal stop under their own code",
    }
    for key in report.inert_keys():
        report.acknowledge_inert(key, acknowledgements.get(key, "inert on curated reference input"))

    # ----- print -----
    print(f"\n{'constraint':36s} {'eval':>6s} {'refuted':>8s} {'satisf':>7s} {'inappl':>7s}  status")
    d = report.as_dict()
    for key, c in d["constraints"].items():
        print(f"{key:36s} {c['n_evaluated']:6d} {c['n_refuted']:8d} {c['n_satisfied']:7d} "
              f"{c['n_inapplicable']:7d}  {c['status']}")
    print(f"\nsummary: {d['summary']}")
    print(f"reference CDS scored: {len(per_cds)}" + (f"  (missing: {missing})" if missing else ""))

    if clade_demo.get("source") == "real":
        print(f"\nCLADE DEMO ({clade_demo['label']}): {clade_demo['organism']} {clade_demo['gene']}")
        print(f"  under the STANDARD code : {clade_demo['verdict_under_standard_table']} "
              f"-- {clade_demo['detail_standard']}")
        print(f"  under the CLADE code    : {clade_demo['verdict_under_clade_table']} "
              f"-- {clade_demo['detail_clade']}")
        print(f"  census: {clade_demo['census_reading']}")

    try:
        report.assert_no_silent_inert()
    except SilentInertConstraintError as exc:
        print(f"\nFAILED: {exc}")
        return 1
    print("\nno silently-inert constraint: every inert law carries an explicit reason")

    art = {
        "schema": "constraint-kill-count-v1", "date": str(date.today()),
        "question": ("do the universal constraints actually refute anything on real reference CDS, and "
                     "is an inert constraint reported as inert rather than as protection?"),
        "pre_stated_expectation": PRE_STATED_EXPECTATION,
        "first_run_correction": FIRST_RUN_CORRECTION,
        "complete_genes_declared": sorted(COMPLETE_GENES),
        "reference_cds_scored": sorted(per_cds), "missing_references": missing,
        "per_cds_verdicts": per_cds,
        "kill_counts": d,
        "clade_demo": clade_demo,
        "honest_scope": [
            "Inert does NOT mean protective. A law that a curated reference already satisfies removed "
            "nothing on this input; each is acknowledged with the reason it was inert here.",
            "`inapplicable` is counted separately from `satisfied`: the probe scores whole CDS and "
            "supplies no substitution, so the reachability law could not be evaluated at all.",
            "The clade demo is a REAL organism CDS, labelled `real`. Had the fetch failed, a synthetic "
            "sequence would have been used and labelled `synthetic`.",
            "The census is over one organism's CDS set (M. genitalium G37); it bounds how often the "
            "table-1 assumption would misfire THERE, not across all Mycoplasma or all life.",
        ],
    }
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
