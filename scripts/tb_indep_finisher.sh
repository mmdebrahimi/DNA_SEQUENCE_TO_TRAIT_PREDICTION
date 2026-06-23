#!/bin/bash
# Self-healing finisher for the independent-TB full run (unattended window). Loops the checkpointed runner
# until scoring stops progressing (remaining = transient fetch-fails), then runs the lineage-disclosure, then
# writes a DONE marker. Safe to leave running for hours: each pass resumes via checkpoint + skip-existing.
set +e
cd "C:/Users/Farshad/PythonProjects/dna_decode" || exit 1
W="D:/dna_decode_cache/tb_indep"
export TB_INDEP_WORK="$W"
LOG="$W/finish.log"
echo "=== finisher start $(date) ===" >> "$LOG"
prev=-1
for i in $(seq 1 40); do
  uv run python -m scripts.run_tb_independent_amr_portal --max 0 --work "$W" >> "$LOG" 2>&1
  scored=$(wc -l < "$W/results.jsonl" 2>/dev/null)
  vcf=$(ls "$W/vcf/"*.vcf 2>/dev/null | wc -l)
  echo "--- pass $i: scored=$scored vcf=$vcf $(date) ---" >> "$LOG"
  if [ "$scored" = "$prev" ]; then echo "no progress -> stop (remaining are fetch-fails)" >> "$LOG"; break; fi
  prev=$scored
done
echo "=== lineage disclosure $(date) ===" >> "$LOG"
uv run python -m scripts.tb_independent_lineage_collapse --work "$W" >> "$LOG" 2>&1
echo "=== DONE $(date) | scored=$(wc -l < "$W/results.jsonl") ===" | tee "$W/DONE" >> "$LOG"
