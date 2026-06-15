"""Alias -> BioSample crosswalk with an explicit conflict taxonomy.

A MIC table keys isolates by some study-local identifier (run accession, sample
alias, secondary accession, or the BioSample itself). The scorer requires BioSample
keys, so we build a crosswalk from ENA read_run records (Step 1's
`read_run_records_for_project`). The canonical BioSample is the `sample_accession`.

CONFLICT TAXONOMY (brainstorm C2):
  - many runs -> 1 BioSample          = OK (normal; not a conflict)
  - 1 native key -> >1 BioSample      = HARD CONFLICT (e.g. a sample_alias reused
                                         across two BioSamples)
  - a candidate value colliding across fields to different BioSamples = HARD CONFLICT

Conflicts carry field PROVENANCE (`mic_key`, `candidate_field`, `candidate_value`,
`resolved_biosample`, `source_row_id`) so a collision is debuggable, never a silent
collapse. `build_oxford_labels` (Step 6) ABORTS on any conflict touching its keys.
"""
from __future__ import annotations

import json
from pathlib import Path

# A candidate value in any of these fields can be a MIC-table native key.
# sample_accession is the canonical BioSample (the resolution TARGET).
CANDIDATE_FIELDS = ("sample_accession", "run_accession", "sample_alias", "secondary_sample_accession")


def build_crosswalk(records: list[dict]) -> dict:
    """Build {native_key: BioSample} + a conflict list from ENA read_run records.

    A record's BioSample is its `sample_accession`. Each candidate field's value
    becomes a native key pointing at that BioSample. A value that points at >1
    distinct BioSample (across rows / fields) is a HARD CONFLICT and is EXCLUDED
    from the clean crosswalk.
    """
    keymap: dict[str, dict[str, list[dict]]] = {}
    for i, rec in enumerate(records):
        bs = (rec.get("sample_accession") or "").strip()
        if not bs:
            continue
        for field in CANDIDATE_FIELDS:
            val = (rec.get(field) or "").strip()
            if not val:
                continue
            keymap.setdefault(val, {}).setdefault(bs, []).append(
                {"candidate_field": field, "source_row_id": i})
    crosswalk: dict[str, str] = {}
    conflicts: list[dict] = []
    for key, bsmap in keymap.items():
        if len(bsmap) == 1:
            crosswalk[key] = next(iter(bsmap))
        else:                       # one key -> multiple BioSamples = HARD CONFLICT
            for bs, prov in bsmap.items():
                for p in prov:
                    conflicts.append({
                        "mic_key": key,
                        "candidate_field": p["candidate_field"],
                        "candidate_value": key,
                        "resolved_biosample": bs,
                        "source_row_id": p["source_row_id"],
                    })
    return {"crosswalk": crosswalk, "conflicts": conflicts}


def resolve_keys(native_keys, records: list[dict]) -> dict:
    """Resolve MIC-table native keys to BioSamples.

    Returns {"resolved": {key: BioSample}, "unresolved": [key,...],
    "conflicts": [provenance-row,...]}. A key implicated in a conflict is NOT
    resolved (it appears only under conflicts). Unresolved = no ENA match at all.
    """
    built = build_crosswalk(records)
    cw = built["crosswalk"]
    conflict_keys = {c["mic_key"] for c in built["conflicts"]}
    wanted = {k.strip() for k in native_keys if k and k.strip()}
    resolved = {k: cw[k] for k in wanted if k in cw and k not in conflict_keys}
    unresolved = sorted(k for k in wanted if k not in cw and k not in conflict_keys)
    relevant_conflicts = [c for c in built["conflicts"] if c["mic_key"] in wanted]
    return {"resolved": resolved, "unresolved": unresolved, "conflicts": relevant_conflicts}


def write_crosswalk(path: str | Path, built: dict) -> None:
    """Persist the crosswalk + conflicts JSON (provenance for reuse + audit)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(built, indent=2), encoding="utf-8")
