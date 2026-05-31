# Two-Machine Capability Map + Transfer-Minimization Long-Term Plan — 2026-05-31

> Soraya planning deliverable. Goal: map what THIS laptop can do, and sequence laptop↔workhorse work so file transfers stop as soon as possible and stay stopped.
>
> **Core thesis:** transfers happen only because of *capability asymmetry* + *un-synced artifacts*. Erase both and the machines never need to swap files again — git carries code/manifests/decisions; public accessions re-derive the data on each machine independently; the CPU runtime is replicated locally. The GPU is the only genuine asymmetry, and the active track (EP-4 pathotype) doesn't use it.

---

## 1. Capability map

| Capability | This laptop (GTX 860M) | Workhorse (Precision 7780) | Notes |
|---|---|---|---|
| Git / commit / push (sync channel) | ✅ | ✅ | Both pull/push `main`. THE transfer replacement. |
| Python/uv venv (torch+cu118, xgboost, sklearn, biopython, h5py, pandas, openpyxl) | ✅ | ✅ | Laptop venv verified 2026-05-31 |
| CUDA GPU inference | ✅ GTX 860M (CC5.0, 4GB) — small/medium N | ✅ stronger GPU | NT v2 100M runs on 860M (slow); large N → Databricks cloud, NOT workhorse |
| Docker (Mash/AMRFinder/Bakta biocontainers) | ✅ Docker 29.4.3 | ✅ | `tools/docker_runner.py` already wraps these |
| VirulenceFinder + blastn (pathotype caller) | ⚠️ not installed — **replicable** (pip VF + BLAST+ Win binary OR Docker) | ✅ installed (Gate A) | **CPU-based** — no GPU needed. This is the unlock. |
| Web research (`/research`, WebFetch/WebSearch) | ✅ | ✅ (network-dependent) | Laptop did all the accession resolution this session |
| Disk for caches/DB/genomes | ✅ **D: 4.4 TB free** (C: 28 GB — avoid) | ✅ | Route ALL caches/DB/Docker volumes to D: on laptop |
| `pathotype_horesh_*` scripts + Gate A config | ❌ absent (workhorse-only, **not in git**) | ✅ | **The one real blocker** — see Sub-plan A2 |
| Path-gated skills (`/project-state`) | ⚠️ refuses from `rca_engine/articles` cwd | ✅ (from dna_decode cwd) | Launch sessions from `dna_decode/` to fix |

**Bottom line:** the laptop can do ~everything EP-4 needs once two gaps close (replicate VF runtime + sync the scripts). The workhorse's only irreplaceable asset is heavier GPU — used only by the parked AMR-classifier track, whose scale-out target is Databricks (cloud) anyway.

---

## 2. The transfer-minimization protocol (the rules that keep it transfer-free)

1. **Code / scripts / manifests / memos / ledger → git tracked dirs only.** Never zip. (`research_outputs/`, `plans/`, `scripts/`, `project_state/`, `soraya_runs/`, `wiki/` are tracked; `data/` and `reports/` are gitignored.)
2. **Data is NEVER transferred — it is re-derived from public accessions.** Manifests carry `wgs_accession` / `GCA_` / `LR` accessions; each machine fetches its own genomes from NCBI/ENA. Genomes, VF DB, caches stay machine-local (on D:), reproducible from public sources.
3. **Anything a machine needs from the other = commit it to a tracked dir.** If it's gitignored (data/DB), it must be re-derivable from a committed manifest + a documented fetch command, not shipped.
4. **Each machine: pull → run locally → commit results → push.** Weekly sync = `git pull && git push`. No bundles.
5. **Sessions launch from `dna_decode/`** so `/project-state` works (avoids the cwd path-gate that forced direct-Edit workarounds).

---

## 3. Long-term plan — tracks + phased sub-plans

### Track A — EP-4 Pathotype (ACTIVE) → make fully laptop-local

**Terminal:** the bounded 3-class slice (ExPEC+EPEC+ETEC) materialized + called + provenance-split-evaluated, entirely on the laptop, with Gate B sent. After A1+A2, **zero transfers** for this track.

| Sub-plan | What | Machine | Transfer? | Status |
|---|---|---|---|---|
| **A1** Replicate pathotype runtime on laptop | Install VirulenceFinder (pip 3.2.0) + BLAST+ Win binary (or use Docker VF), clone VF DB → **D:\dna_decode_stage2\virulencefinder_db**, pin versions | laptop | none (DB clones from CGE) | TODO — laptop-solo, dep-install gate |
| **A2** Sync workhorse pathotype scripts into git | Workhorse commits `scripts/pathotype_horesh_*.py` + `config/pathotype_gate_a_panel.json` + gate_a contract (NOT data/DB) → laptop `git pull` gets them | **workhorse (one-time)** | **last required workhorse→laptop sync** | EMIT to user |
| **A3** Materialize bounded slice by-accession | Fetch ExPEC 135 + EPEC 125 (Horesh WGS) + ETEC (8 refs/558 GCA von Mentzer) → D: cache, by accession | laptop | none (public accessions) | ready (manifests committed) |
| **A4** Run callers + provenance-split eval | VF/blastn over the slice; report resolver-conformance vs external-validity separately (ExPEC strong / EPEC medium / ETEC near-conformance) | laptop | none | after A1-A3 |
| **A5** Gate B outreach | Send 5 emails, score replies (≥2 yes = PASS) | **user-only** | none | EMIT — still unsent |

**After A1 + A2, Track A is permanently transfer-free.**

### Track B — AMR / NT classifier (PARKED, v0/v0.1 shipped)

**State:** cipro genome-input + cef cached shipped; tet falsified; v0 closed. Residual = Stage 2 N=150 (GPU-heavy).

| Sub-plan | What | Machine | Transfer? |
|---|---|---|---|
| B1 (parked) | Stage 2 N=150 NT-embedding populate | Databricks **cloud** (GTX 860M too slow; workhorse optional) | none (cloud) |
| B2 | Any re-analysis of existing results | laptop | none (results in git) |

This track has **no hard workhorse dependency** — its scale-out is cloud, not the workhorse.

### Track C — Roadmap Phases 3–6 (LONG-HORIZON, mostly research/planning)

Multi-organism (Klebsiella) → non-AMR phenotypes → multimodal → eukaryotic (`plans/Trait_Decoding_Roadmap.md`). Until execution, these are **research + substrate-survey + planning** = 100% laptop-solo (web research, feasibility, accession census). Execution phases re-evaluate compute then (likely cloud burst, not workhorse).

---

## 4. Critical path to "transfer-free"

```
A2 (workhorse commits pathotype scripts)  ──┐
                                             ├──► A3 → A4 (laptop runs EP-4 end-to-end, no transfers)
A1 (laptop installs VF+blastn+DB on D:)  ───┘
A5 (user sends Gate B) ── gates heavy build, runs in parallel
```

**Two moves unlock permanent transfer-freedom for the active track:** A1 (laptop-solo, I can do) + A2 (one last workhorse push). Everything after is pull/run/push.

---

## 5. What I (Soraya, laptop) can execute solo right now — the capability inventory

**Money-free, in-cwd, no-workhorse actions available to advance the plan:**
- A1: install VF + BLAST+ + clone VF DB to D: (dep-install gate; un-gated default). **Highest-VOI next.**
- A3: write the by-accession materialize driver (pure python; reuse `tools/docker_runner.py` patterns) — even before scripts sync, I can write a laptop-native materializer from the committed manifests.
- A4-prep: draft the provenance-split eval-harness spec.
- Track C: any roadmap-phase research / substrate census.
- Ongoing: ledger/decision/memo upkeep, manifest building, web research.

**Cannot do (emit / defer):**
- A2: workhorse must commit its scripts (emit).
- A5: user sends Gate B (emit).
- Heavy GPU N=150: Databricks (cloud, separate).

---

## 6. Emitted actions (user / workhorse — Soraya can't self-invoke)

1. **WORKHORSE (one-time, closes the last forced transfer):** from `dna_decode/` on the workhorse:
   `git add scripts/pathotype_horesh_*.py config/pathotype_gate_a_panel.json <gate_a contract md> && git commit && git push`
   (commit CODE + CONFIG only; never the VF DB or genomes — those re-derive from CGE + accessions.)
2. **USER:** send Gate B outreach (`research_outputs/pathotype_gate_b_send_kit_2026-05-29.md`).

---

## 7. Recommended immediate next (if you tell Soraya to advance)
A1 — replicate the VF+blastn runtime on the laptop (to D:), so EP-4 caller execution stops needing the workhorse. Single self-contained laptop action; the linchpin of the whole transfer-minimization plan.
