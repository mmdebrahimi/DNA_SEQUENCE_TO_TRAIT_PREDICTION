<!-- intent-contract.md for 2026-06-03-0000-ep4-v0-cohort-calibrate -->
# Intent Contract — 2026-06-03-0000-ep4-v0-cohort-calibrate

**Persona:** Soraya `--advance` (money-only; UNATTENDED overnight — user asleep, full autonomous mandate "do any long runs you need now, all recommendations").
**Terminal goal:** advance the v0 pathotype resolver per the post-ship recommendation set.

## Batch (all 4 recommendations)
| # | recommendation | gate | decision |
|---|---|---|---|
| 1 | full-cohort H4 eval (24 genomes) with per-genome coverage cache | auto | run (the long run) |
| 2 | calibrate ExPEC sensitivity (lift recall vs abstention) | auto | run |
| 3 | real CGE VirulenceFinder side-by-side diff | dep-install (irreversible) | **DEFER** — needs KMA/BLAST+ external binaries; not unattended-installable on this Windows/uv host (no pip/conda, no aligner). Documented + gated to user. |
| 4 | cache per-genome k-mer coverage (speed) | auto | run (built into the eval) |

## Unattended-safety posture
- Money: none. Deletes: none. **Dep-installs: declined** (rec 3) — judgment widens the gate; an unattended BLAST/KMA install could break the env and the user can't approve. This is the correct fail-safe direction.
- Commit + push at milestones (user authorized push this session + "move forward").

## Stop condition
All executable recommendations done + committed/pushed; rec 3 documented as gated; handoff + retrospective written.
