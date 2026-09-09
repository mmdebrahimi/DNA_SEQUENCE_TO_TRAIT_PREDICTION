"""Salmonella enterica serovar caller — antigen allele DB via the shared blastn engine.

The Kauffmann-White-Le Minor scheme defines a serovar by an ANTIGENIC FORMULA: O-antigen group : H1
(phase-1 flagellin = fliC) : H2 (phase-2 flagellin = fljB). E.g. Typhimurium = `1,4,[5],12 : i : 1,2`.
SeqSero2 determines this from WGS by detecting the O-group + fliC + fljB allele sequences, then looks the
formula up in the White-Kauffmann table. This is the SAME shape as the E. coli serotype caller (best O
antigen + best H antigen) with a serovar-name lookup on top (the way ktype adds a wzi->KL map).

DB layout (a directory):
  * `salmonella_antigens.fasta` -- allele sequences, headers `<axis>__<antigen>__<id>` where
    axis in {O, H1, H2}, e.g. `O__9__01`, `H1__d__01`, `H2__1,2__01`. (The shared engine keys on the
    full header; this runner parses axis + antigen from it.)
  * `serovar_table.tsv` -- `O<tab>H1<tab>H2<tab>Serovar`, the White-Kauffmann-Le Minor formula table.

HONESTY (load-bearing): faithful to the SeqSero2 / Kauffmann-White method (blastn over the antigen DB +
formula lookup); NOT an independent baseline. Many formulas are shared by >1 serovar / are phase-incomplete
-> the call is the FORMULA + a serovar IFF the formula resolves uniquely (else serovar=None, like O?/H?).
Salmonella enterica only. Offline-safe via the shared engine (missing blastn/DB -> status 'unavailable').
"""
from __future__ import annotations

from pathlib import Path

from dna_decode.typing.blast_caller import call_alleles

# Salmonella antigen alleles are near-exact; SeqSero2 uses high-identity matching.
SEROVAR_IDENTITY_THRESHOLD = 90.0

# COVERAGE lowered 80.0 -> 40.0 on 2026-09-04, against a bar registered BEFORE the sweep ran
# (scripts/salmserovar_threshold_tradeoff.py::PREREGISTERED: adopt iff >=7 rescued AND <=2 newly wrong
# AND net >= +5). Measured on the 200-isolate wet-lab-labelled cohort: 36 abstentions became CORRECT
# calls for exactly 1 new error, net +35, plus 7 wrong->correct the rule did not even count.
#
# WHY COVERAGE AND NOT IDENTITY. The O-antigen alleles that were being discarded hit at near-perfect
# identity (median 99.8) but partial coverage (median 58.4, max 78.9) -- a partial-alignment/allele-
# length mismatch, concentrated on the O7 wzx/wzy reference. Identity stays at 90 because identity was
# never the failing axis; relaxing it would admit genuinely different alleles.
#
# THE TRADE WAS NOT ASSUMED FREE. Coverage rose 0.705 -> 0.900 AND accuracy-on-covered rose
# 0.702 -> 0.783 -- both improved, which is not the usual selective-classification trade and is the
# reason the change is safe rather than merely favourable.
SEROVAR_COVERAGE_THRESHOLD = 40.0


def parse_axis_antigen(allele_id: str) -> tuple[str, str] | None:
    """`O__9__01` -> ('O', '9'); `H1__d__01` -> ('H1', 'd'); `H2__1,2__01` -> ('H2', '1,2'). None if malformed.

    Deliberately returns None for the `SP__` axis (the two `tyr` genes the O procedure consults), so a
    specific gene can never be selected as an O/H1/H2 antigen by `_best_per_axis`."""
    parts = allele_id.split("__")
    if len(parts) < 2 or parts[0] not in ("O", "H1", "H2"):
        return None
    return parts[0], parts[1]


def parse_ss2_key(allele_id: str) -> str | None:
    """`O__9,46__239__O-9,46_wbaV__1002` -> `O-9,46_wbaV__1002`; a 3-field legacy header -> None.

    The 4th field is the ORIGINAL SeqSero2 header. It embeds `__` itself, so it is recovered by
    re-joining every field from index 3 on -- which round-trips exactly. Returning None on a legacy
    3-field header is load-bearing: the DB is gitignored, so a machine holding an older build must
    fall back to the pre-role behaviour visibly rather than silently mis-call every branch."""
    parts = allele_id.split("__")
    if len(parts) < 4:
        return None
    return "__".join(parts[3:])


def parse_source_index(allele_id: str) -> int:
    """The 3rd field is the record's index in SeqSero2's own FASTA. Used to reproduce upstream's
    file-order tiebreak (see `_o_argmax`). Unparseable -> a large sentinel so it sorts last."""
    parts = allele_id.split("__")
    try:
        return int(parts[2])
    except (IndexError, ValueError):
        return 10**9


# ---------------------------------------------------------------------------------------------------
# The O decision procedure, ported from SeqSero2 1.3.2 `SeqSero2_package.py::call_O_and_H_type`
# (lines 1070-1131), read from the pinned biocontainer. Source reading + the threshold-translation
# argument: `wiki/seqsero2_o_algorithm_reading_2026-09-09.md`.
#
# WHY A PORT AND NOT A HEURISTIC. Six of SeqSero2's O entries are NOT standalone positive alleles --
# two are DIFFERENTIAL MARKERS, two are the O9-vs-O9,46 discriminator `wbaV`, two are explicitly
# partial. Our previous caller read all six as ordinary positive alleles, which is one cause with two
# measured symptoms: `1,3,19` asserted whenever its marker hit (60% of a diverse cohort against a 2%
# true prevalence) and plain O9 structurally unemittable. There is no threshold that fixes that --
# the entries mean different things, so the procedure that consumes them has to know which is which.
#
# THRESHOLD TRANSLATION (argued, NOT measured -- see the memo's honest limits). SeqSero2's score is
# `len(ref_kmers & input_kmers) / len(ref_kmers) * 100`: coverage OF THE REFERENCE. Our
# `percent_coverage` is `alignment_length / qlen * 100` with the allele as query -- the same axis. So
# its `>70` / `>40` branch gates are compared against our coverage. Its membership bar is `>25`;
# ours is the caller's own coverage threshold (40 by default), i.e. STRICTER. That difference is left
# in place deliberately: the 40 was itself adopted against a pre-registered bar, and changing two
# things at once would make neither measurable.
WBAV_KEYS = ("O-9,46_wbaV__1002", "O-9,46_wbaV-from-II-9,12:z29:1,5-SRR1346254__1002")
WZY_946_KEYS = ("O-9,46_wzy__1191", "O-9,46_wzy_partial__216")
WZY_94627_KEY = "O-9,46,27_partial_wzy__1019"
WZX_310_KEY = "O-3,10_wzx__1539"
WZY_946_FULL_KEY = "O-9,46_wzy__1191"
NOT_IN_1319_KEY = "O-3,10_not_in_1,3,19__1519"   # used POSITIVELY: selects O-3,10
NOT_IN_310_KEY = "O-1,3,19_not_in_3,10__130"     # never a positive call; excluded from the argmax
WBAV_GATE = 70.0
WZY_GATE = 40.0


def _o_argmax(cands: list[tuple[str, str, int, float, float]],
              exclude: str | None = None) -> tuple[str, str] | None:
    """Best O candidate as (antigen, ss2_key), ranked IDENTITY-primary then coverage, or None.

    `cands` is (antigen, ss2_key, source_index, identity, coverage).

    WHY NOT UPSTREAM'S SCORE ORDER. SeqSero2 ranks on ONE number (k-mer recovery), which conflates
    how much of the reference is present with how well it matches. BLAST hands us those as TWO
    numbers, and which one leads was already settled here by measurement: coverage-primary picks the
    WRONG antigen when near-identical alleles both align full-length, so `_best_per_axis` ranks
    identity first. A first cut of this port ranked branch C on coverage alone and regressed six
    isolates -- five where the multi-gene `O:23-gene2` entry out-covered the correct `O-13_wzx`, and
    one where `9,46` out-covered `3,10` -- each of which our identity-primary ranking had called
    right. So the BRANCHES are ported and the SCORING is not: the branch preconditions come from
    upstream, the generic ranking stays the one this caller has evidence for.

    Ties are broken by source index (the entry's position in SeqSero2's own FASTA) rather than by
    whatever order blastn reported hits in, so the same genome cannot call differently between runs."""
    best: tuple[str, str] | None = None
    best_key = (-1.0, -1.0)
    for antigen, key, _idx, ident, cov in sorted(cands, key=lambda c: c[2]):
        if exclude is not None and key == exclude:
            continue
        if (ident, cov) >= best_key:
            best_key, best = (ident, cov), (antigen, key)
    return best


def call_o_antigen(per_allele: dict) -> dict:
    """SeqSero2's three-branch O procedure over called O alleles. Returns the antigen + which rule fired.

    `o_antigen_rule` is part of the output on purpose: a differential-marker call and an ordinary
    best-allele call are not the same kind of evidence, and collapsing them to a bare antigen string is
    what hid the defect in the first place."""
    o_cands: list[tuple[str, str, int, float, float]] = []
    scores: dict[str, float] = {}
    special: dict[str, tuple[float, float]] = {}
    legacy_headers = False
    for allele_id, hit in per_allele.items():
        if not hit["called"]:
            continue
        key = parse_ss2_key(allele_id)
        parts = allele_id.split("__")
        axis = parts[0] if parts else ""
        if axis not in ("O", "SP"):
            continue
        if key is None:
            legacy_headers = True
            continue
        ident, cov = float(hit["percent_identity"]), float(hit["percent_coverage"])
        if axis == "SP":
            special[key] = max(special.get(key, (0.0, 0.0)), (ident, cov))
        else:
            scores[key] = max(scores.get(key, 0.0), cov)
            o_cands.append((parts[1], key, parse_source_index(allele_id), ident, cov))

    if legacy_headers and not scores:
        # An older, role-less DB build. Say so rather than reporting a confident branch result off an
        # empty dict -- the caller then falls back to the legacy best-allele pick.
        return {"antigen": None, "o_antigen_rule": "legacy_db_no_role_field"}
    if not scores:
        return {"antigen": None, "o_antigen_rule": "no_o_allele_called"}

    # Branch A -- the O9 family, gated on wbaV coverage.
    if any(scores.get(k, 0.0) > WBAV_GATE for k in WBAV_KEYS):
        if any(scores.get(k, 0.0) > WZY_GATE for k in WZY_946_KEYS):
            return {"antigen": "9,46", "o_antigen_rule": "A_wbaV_with_wzy"}
        if WZY_94627_KEY in scores:
            return {"antigen": "9,46,27", "o_antigen_rule": "A_wbaV_with_partial_wzy_27"}
        # THE O9-vs-O2 REFINEMENT IS DELIBERATELY NOT PERFORMED, and this is a MEASURED refusal
        # rather than an omission. Upstream splits O-9 from O-2 by comparing k-mer scores for the two
        # `tyr` genes. Those two references are ~99% identical to each other, so under BLAST they BOTH
        # align essentially full length and the comparison collapses:
        #   real O9 genome GCA_009152955 -> tyr-O-9 id 99.30 cov 100.00 | tyr-O-2 id 99.20 cov 100.20
        #   real O9 genome GCA_009474265 -> tyr-O-9 id 100.0 cov 100.00 | tyr-O-2 id  99.80 cov 100.20
        # Coverage is degenerate (and exceeds 100 on the gapped alignment, so a 0.2-point artifact
        # decides an antigen); identity separates them by 0.1-0.2 points, which is noise. A first cut
        # of this port compared coverage and fired `O-2` on 16 isolates: 13 miss, 3 no-call, ZERO
        # hits -- every one a true O9 (Enteritidis, Typhi, Berta) that SeqSero2 called `9`.
        # A k-mer-recovery score genuinely carries this discrimination; a near-full-length homologous
        # alignment does not. So we take the branch's own DEFAULT (O-9) and name the gap.
        # KNOWN COST, stated rather than hidden: a true O-2 genome (Paratyphi A, Nitra, Kiel,
        # Koessen) will be called O-9. The cohort this was measured on contains NO O-2 serovar at all,
        # so the positive direction was never testable here -- which is itself a reason not to ship a
        # 0.1-point margin rule. The two `tyr` entries stay in the DB, and their scores are returned,
        # so a future k-mer-based refinement has what it needs.
        tyr = {k: v for k, v in special.items() if "tyr-O-" in k}
        out = {"antigen": "9", "o_antigen_rule": "A_wbaV_default_O9_tyr_not_discriminable_by_blast"}
        if tyr:
            out["tyr_scores"] = {k: {"percent_identity": v[0], "percent_coverage": v[1]}
                                 for k, v in sorted(tyr.items())}
        return out

    # Branch B -- the 3,10 / 1,3,19 family. The ONLY route to a positive `1,3,19`, and it is reached
    # by the ABSENCE of the reciprocal marker, never by a marker hit.
    if WZX_310_KEY in scores and WZY_946_FULL_KEY in scores:
        if NOT_IN_1319_KEY in scores:
            return {"antigen": "3,10", "o_antigen_rule": "B_not_in_1319_present"}
        return {"antigen": "1,3,19", "o_antigen_rule": "B_not_in_1319_absent"}

    # Branch C -- generic argmax, with upstream's two guards.
    win = _o_argmax(o_cands)
    if win is None:
        return {"antigen": None, "o_antigen_rule": "no_o_allele_called"}
    antigen, key = win
    if key in WBAV_KEYS and not any(k in scores for k in WZY_946_KEYS):
        return {"antigen": "9", "o_antigen_rule": "C_guard_wbaV_without_wzy"}
    if antigen == "1,3,19":
        # Guard 2: the differential marker cannot win on its own. Upstream re-runs the argmax with
        # that exact key excluded, and leaves the call EMPTY if nothing else is in the dict.
        win2 = _o_argmax(o_cands, exclude=NOT_IN_310_KEY)
        if win2 is None:
            return {"antigen": None, "o_antigen_rule": "C_guard_1319_marker_only"}
        antigen, key = win2
        if key in WBAV_KEYS and not any(k in scores for k in WZY_946_KEYS):
            return {"antigen": "9", "o_antigen_rule": "C_guard_wbaV_without_wzy"}
        return {"antigen": antigen, "o_antigen_rule": "C_argmax_after_1319_exclusion"}
    return {"antigen": antigen, "o_antigen_rule": "C_argmax"}


def load_serovar_table(tsv: str | Path) -> dict[tuple[str, str, str], str]:
    """Parse `O<tab>H1<tab>H2<tab>Serovar` -> {(O, H1, H2): serovar}. Header row (`O`/`Serovar`) skipped."""
    table: dict[tuple[str, str, str], str] = {}
    for line in Path(tsv).read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 4 or parts[0].lower() == "o" or parts[3].lower() == "serovar":
            continue
        o, h1, h2, serovar = parts[0], parts[1], parts[2], parts[3]
        if o and serovar:
            table[(o, h1, h2)] = serovar
    return table


def _best_per_axis(per_allele: dict) -> dict[str, dict]:
    """For each axis (O/H1/H2), the best CALLED antigen — ranked by IDENTITY then coverage.

    Identity-PRIMARY is load-bearing for the H antigens: flagellin (fliC/fljB) alleles cross-hybridize at
    near-full COVERAGE across different antigen types, so a coverage-only tiebreak picks the WRONG antigen
    (e.g. Typhimurium fliC=i blasts ~100% id while the cross-hybridizing fliC=r hits ~92% id at equal
    coverage). The true antigen is the highest-IDENTITY hit. (Verified on S. Typhimurium LT2: coverage-only
    gave 4:r:1,5,7; identity-primary gives the correct 4:i:1,2.)"""
    axis_best: dict[str, dict] = {}
    for allele_id, hit in per_allele.items():
        if not hit["called"]:
            continue
        pa = parse_axis_antigen(allele_id)
        if pa is None:
            continue
        axis, antigen = pa
        key = (hit["percent_identity"], hit["percent_coverage"])
        cur = axis_best.get(axis)
        if cur is None or key > (cur["percent_identity"], cur["percent_coverage"]):
            axis_best[axis] = {"axis": axis, "antigen": antigen, "best_allele": allele_id,
                               "percent_identity": hit["percent_identity"], "percent_coverage": hit["percent_coverage"]}
    return axis_best


def call_serovar(fasta: str | Path, db_dir: str | Path, *,
                 identity_threshold: float = SEROVAR_IDENTITY_THRESHOLD,
                 coverage_threshold: float = SEROVAR_COVERAGE_THRESHOLD,
                 blastn_bin: str | None = None, timeout: int = 600) -> dict:
    """blastn the Salmonella antigen DB vs `fasta`; assemble the O:H1:H2 formula + look up the serovar."""
    db_dir = Path(db_dir)
    antigens_fasta, serovar_tsv = db_dir / "salmonella_antigens.fasta", db_dir / "serovar_table.tsv"
    if not antigens_fasta.exists() or not serovar_tsv.exists():
        return {"status": "unavailable", "serovar": None, "antigenic_formula": None,
                "reason": f"Salmonella antigen DB not found in {db_dir} "
                          "(need salmonella_antigens.fasta + serovar_table.tsv)"}
    res = call_alleles(fasta, antigens_fasta, identity_threshold=identity_threshold,
                       coverage_threshold=coverage_threshold, blastn_bin=blastn_bin, timeout=timeout)
    if res["status"] != "ok":
        return {"status": "unavailable", "tool": res.get("tool"), "serovar": None,
                "antigenic_formula": None, "reason": res.get("reason")}

    axis_best = _best_per_axis(res["per_allele"])
    o_call = call_o_antigen(res["per_allele"])
    o_rule = o_call["o_antigen_rule"]
    if o_rule == "legacy_db_no_role_field":
        # Role-less DB build: fall back to the pre-2026-09-09 best-allele pick and SAY so, so a run on
        # an out-of-date DB is never mistaken for a run of the ported procedure.
        o = axis_best.get("O", {}).get("antigen")
    else:
        o = o_call["antigen"]
    h1 = axis_best.get("H1", {}).get("antigen")
    h2 = axis_best.get("H2", {}).get("antigen")
    formula = f"{o or 'O?'}:{h1 or 'H?'}:{h2 or '-'}" if (o or h1) else None

    serovar = None
    if o and h1:
        table = load_serovar_table(serovar_tsv)
        serovar = table.get((o, h1, h2 or "-")) or table.get((o, h1, h2 or ""))
        # phase-incomplete fallback: match O+H1 ignoring H2 IFF that resolves uniquely
        if serovar is None:
            cands = {sv for (to, th1, _), sv in table.items() if to == o and th1 == h1}
            serovar = next(iter(cands)) if len(cands) == 1 else None

    # The evidence rows are the BEST HIT per axis. Since the O procedure landed, the best O hit is
    # frequently NOT the O call -- a `wbaV` hit at 99.7% identity is the evidence FOR calling plain O9,
    # and rendering it beside a `9` call with no marking reads as a contradiction. Each row therefore
    # says whether it IS the call, and the O best hit is surfaced separately.
    for _row in axis_best.values():
        _row["is_call"] = (_row["antigen"] == {"O": o, "H1": h1, "H2": h2}.get(_row["axis"]))
    o_best = axis_best.get("O", {}).get("antigen")

    base = {"status": "ok", "tool": "blastn", "method": "seqsero2_o_procedure_port_v1",
            "parameters": {"identity_threshold": identity_threshold, "coverage_threshold": coverage_threshold},
            "o_antigen": o, "h1_antigen": h1, "h2_antigen": h2, "o_antigen_rule": o_rule,
            "o_antigen_best_hit": o_best,
            "antigenic_formula": formula, "serovar": serovar,
            "antigens": sorted(axis_best.values(), key=lambda v: v["axis"])}
    return base
