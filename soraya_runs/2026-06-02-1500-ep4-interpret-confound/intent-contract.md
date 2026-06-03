<!-- intent-contract.md for 2026-06-02-1500-ep4-interpret-confound -->
# Intent Contract — 2026-06-02-1500-ep4-interpret-confound

**Persona:** Soraya `--advance` (money-only gate; attended dogfood flight).
**Terminal goal (north star):** read a DNA sequence → identify which sequence regions map to which trait. For this project (EP-4 / `ecoli-pathotype-prediction-cli-2026-05-26`): emit a pathotype call WITH provenance — which acquired virulence-gene clusters drove it.

## Why this run
EP-4 learned-model thread just closed: pooled whole-genome NT = 0.38 (≤chance) ≪ classical full-spectrum k-mer **presence/absence = 0.729**. Signal is accessory-genome GENE PRESENCE — which VINDICATES the ledger's LOCKED v0 design (deterministic virulence-gene-cluster resolver), not the learned pooled track.

But the **dominant unresolved risk** is the study confound: ExPEC=Salipante, EPEC=Hazen (no single study has both classes). The 0.729 could be study/assembler batch, not biology. **Highest-VOI move = de-risk before building: interpret the winning model and test whether its discriminative k-mers are biology or batch.** (expected_info_gain high, cost low — reuses the 24 cached genomes + sklearn.)

## Planned batch (≤ dynamic cap 8–10; this run uses 4)
| # | action | gate verdict | decision |
|---|---|---|---|
| 1 | git commit today's completed EP-4 work (local) | auto | run |
| 2 | build `scripts/pathotype_model_interpret_confound.py` | auto | run |
| 3 | run it: top discriminative presence/absence k-mers + per-study presence pattern → biology-vs-batch verdict | auto | run |
| 4 | record findings to ledger (`/project-state` emit — user-only) + write result/audit + dogfood notes | mixed | run local; EMIT `/project-state` |

**Held (judgment widens gate):** `git push origin main` — irreversible; Soraya money-only default would run it, but held per the user's ~weekly manual-sync preference + the unanswered commit-vs-push question. Emitted as a recommendation, not executed.

**Money actions:** none. **Deletes:** none.

## Stop condition
Batch complete, OR money-declined (n/a), OR a user-only command blocks the only path (emit + continue).
