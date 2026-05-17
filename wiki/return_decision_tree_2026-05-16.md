# Return Decision Tree — 2026-05-16

> 30-second triage card for what to do when you return. Detail in appendix below.

## TL;DR (read this first)

**Two jobs ran while you were out.** Both write verdict packets at known paths; packets are atomic file-writes (independent of stdout buffering). **Packet presence/absence is the authoritative signal.**

```bash
# triage block — paste this; reports a 4-letter state code
cd C:/Users/Farshad/PythonProjects/dna_decode
S1=$(ls wiki/stage1_n40_cipro_mean-plus-max_*.md 2>/dev/null | wc -l)
PF=$(ls wiki/cipro_attribution_preflight_*.md 2>/dev/null | wc -l)
[ $S1 -gt 0 ] && S1S="S1+" || S1S="S1-"
[ $PF -gt 0 ] && PFS="PF+" || PFS="PF-"
echo "STATE: $S1S$PFS"
```

**Timestamp gate:** Stage 1b launched 2026-05-16 14:43 EDT. Expected verdict packet: ~19:43 EDT.
- Returning **>7 hr** after launch → packet present OR process dead (no third state).
- Returning **3-7 hr** after launch → packet may not exist yet. Check `tasklist | grep python.exe` — alive = wait; dead = packet should exist (if not, recovery).

| state | preflight | stage1b | next action |
|---|---|---|---|
| `S1+PF+` | landed | landed | Read both packets. Apply PIVOT TRIGGER (below). |
| `S1+PF-` | missing | landed | Stage 1b verdict alone is decisive. Apply trigger w/o attribution evidence. |
| `S1-PF+` | landed | missing | Stage 1b still running OR crashed. Check `tasklist`. Don't act on preflight alone — preflight already landed INCONCLUSIVE_MISS today; see appendix. |
| `S1-PF-` | missing | missing | Worst case. See appendix §D recovery. |

## 2×2 verdict matrix (when both packets exist)

Cross-reference Stage 1b `stage2_action` (in verdict packet) × preflight verdict (today's was **INCONCLUSIVE_MISS**):

| Stage 1b → ↓ Preflight | `BURST_STAGE_2` (CLEAN) | `HOLD_STAGE_2_CI_DEGENERATE` (NOISY) | `ALTERNATIVE_POOLING_RERUN` (FAIL) | `PIVOT_TO_BAKTA` |
|---|---|---|---|---|
| **STRONG_POSITIVE** | Stage 2 burst justified | Hold burst; tighten classical (AMRFinder POINT*) | apply PIVOT TRIGGER | Bakta pathway |
| **WEAK_POSITIVE** | Stage 2 burst w/ caveat | Hold burst; rerun w/ tighter comparator | apply PIVOT TRIGGER | Bakta pathway |
| **INCONCLUSIVE_MISS** (today) | Stage 2 burst BUT require attribution at N=150 verdict. **Sleep on it before firing burst — $$$ commit.** | Hold burst. Build curated AMR baseline FIRST per PIVOT TRIGGER. | apply PIVOT TRIGGER | Bakta pathway |

## PIVOT TRIGGER (post-sunk-cost-check 2026-05-16)

Declare frozen NT whole-genome pooling **falsified** for Phase 1 if ALL FOUR hold:

1. NT mean+max AUROC does NOT exceed k-mer by ≥0.05 absolute, AND CI includes ≤0 lift.
2. NT remains below 0.70 AUROC on cipro.
3. Attribution still fails to recover ANY cipro mechanism locus (already true for mean-pool; mean+max requires `gene_level_mutagenesis` refactor — defer).
4. Curated AMR baseline (AMRFinderPlus POINT + Bakta gene presence + k-mer + MLST features) reaches ≥0.80 AUROC at N=38 OR beats NT by ≥0.10.

If all 4 → stop running more pooling variants. Curated baseline + targeted-NT (per-gene window) are the next experiments.

If Stage 1b passes only under random/LOO BUT fails under MLST/phylogeny-aware splits → **lineage signal, not scientific validation.** Validation gate is phylogeny-aware.

## Cleanup ledger ops (after deciding next action)

```
/project-state dna-decode-2026-05-11 --append-action --class run-tests --description "Stage 1b N=38 cipro outcome (mean+max + scaled NT-logreg)" --outcome "<verdict text: stage2_action + gap_pp + CI + per-variant AUROCs>"

/project-state dna-decode-2026-05-11 --refresh-frame --current-state "<post-Stage-1b one-line state>"
```

Then ledger commit batch + push. **Expected untracked-but-OK in `git status`:** `heartbeat_stage1b.bat`, `stage1b_heartbeat.log`, `project_state/dna-decode-2026-05-11-scratch.md`, `wiki/smoke_gate_12strain_cipro_2026-05-15.md`. Anything else surprising → inspect before staging.

---

## Appendix — Recovery procedures (read only if needed)

### What was launched
| job | mechanism | output | ETA from 14:43 launch |
|---|---|---|---|
| Stage 1b | detached via `run_stage1b_detached.bat` (PID 4736) | `wiki/stage1_n40_cipro_mean-plus-max_<date>.md` + `stage1b_detached.log` | ~5 hr (verdict packet written at end only) |
| Preflight v2 | Bash bg job `bx13oua9k` | `wiki/cipro_attribution_preflight_<date>.md` + `.per_strain.json` sidecar | ~5-10 min (✓ COMPLETED 15:26 EDT — verdict INCONCLUSIVE_MISS) |
| Heartbeat sidecar | `heartbeat_stage1b.bat 4736` (broken sleep — loops fast, ~2 lines/sec; harmless but noisy) | `stage1b_heartbeat.log` | self-terminates on verdict-packet landing OR PID 4736 death |

**Heartbeat is NOT load-bearing for any decision branch** — kept as informational aliveness signal only.

### Preflight v2 verdict (already landed 2026-05-16 15:26)
INCONCLUSIVE_MISS. No QRDR + no expanded-set cipro loci recovered. Per packet's own rubric: NOT damning. Mean-pool dilution (1/N_genes per knockout, N≈5000) + RefSeq symbol coverage (~11% CDSs carry `gene=`) are unaddressed confounders. Next escalation per rubric: refactor `gene_level_mutagenesis` to accept aggregation kwarg + retest mean+max preflight matching Stage 1b's pooling. Damning-MISS verdict requires BOTH mean-pool AND mean+max preflights to miss.

### B. Stage 1b packet missing, preflight verdict exists
- Check `tasklist | grep python.exe` (or `MSYS_NO_PATHCONV=1 tasklist /FI "PID eq 4736"` for PID-specific filter — MSYS otherwise mangles the quoted arg per CLAUDE.md gotcha).
- If process alive past 7 hr from launch → suspect hang; consider killing + relaunching.
- If process dead but no verdict packet:
  - `taskkill /F` does **NOT** flush Python stdio buffers — it's `TerminateProcess` (Windows SIGKILL); partial stdout is DISCARDED, not recovered. If python.exe was launched WITHOUT `-u` / `PYTHONUNBUFFERED=1`, the `stage1b_detached.log` is unrecoverable. Verdict packet (atomic file write via `write_text`) is the only authoritative signal.
  - Relaunch via `run_stage1b_detached.bat` (now has `-u` flag after Step 3 patch); rewrite heartbeat with PowerShell `Start-Sleep` for reliable cadence.

### C. Preflight packet missing, Stage 1b packet exists
- Preflight likely crashed silently; check `C:/Users/Farshad/AppData/Local/Temp/preflight_v2.log` for traceback.
- Re-run: `uv run python -u scripts/cipro_attribution_preflight.py` (foreground, no `| tail` pipe).
- Apply Stage 1b verdict alone via 2×2 matrix (preflight column unknown → assume INCONCLUSIVE_MISS row pending re-run).

### D. Both packets missing
- Worst case. Inspect:
  - `tasklist | grep python.exe` — anything alive?
  - `stage1b_detached.log` + `preflight_v2.log` — partial output? traceback?
- Recovery: relaunch BOTH with observability:
  - Stage 1b: patched `.bat` (already has `-u`); add PowerShell-sleep heartbeat for monitoring.
  - Preflight: foreground, NOT bg, NOT piped through `| tail`.

### Open questions surfacing on return (from sunk-cost brainstorm)
1. **Are the 19 R strains' resistance mechanisms actually CHECKED?** Confirmed for ST131 (`1328433.3 / GCA_000522345.1`) via install smoke (gyrA-S83L/D87N + parC-E84V detected). The other 18 R strains: unknown. If many are R via non-QRDR mechanisms (qnr-plasmid / efflux / regulatory), the preflight's miss may partly reflect TRUE biology — NT shouldn't rank gyrA if gyrA isn't the mechanism in those strains. **Worth running AMRFinder against all 19 R strains as a cheap diagnostic before Stage 2 spend.**
2. **AST noise near MIC breakpoints** — how many cipro labels are right at the R/S cutoff?
3. **What specifically justifies Databricks burst spend?** NT > k-mer ≥5 pp? OR NT > curated AMR baseline? The curated baseline experiment has not been run; running it is the load-bearing PIVOT TRIGGER pre-req.
