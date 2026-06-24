#!/bin/bash
# Self-healing finisher for the full 447-isolate ktype cohort validation (TB-finisher pattern).
# Loops the checkpointed runner until all 447 ERR isolates are scored (or stalls), then marks DONE.
cd /c/Users/Farshad/PythonProjects/dna_decode || exit 1
W="D:/dna_decode_cache/ktype_build"; R="$W/cohort_results.jsonl"; LOG="$W/cohort.log"
for pass in $(seq 1 60); do
  before=$(wc -l < "$R" 2>/dev/null | tr -d ' ' || echo 0); before=${before:-0}
  uv run --with openpyxl --no-sync python scripts/ktype_cohort_validate.py \
    --xlsx "$W/zenodo_suppl.xlsx" --limit 0 >> "$LOG" 2>&1
  after=$(wc -l < "$R" 2>/dev/null | tr -d ' ' || echo 0); after=${after:-0}
  echo "$(date +%H:%M) pass $pass: $before -> $after" >> "$LOG"
  if [ "$after" -ge 447 ]; then echo "DONE $after" > "$W/COHORT_DONE"; break; fi
  if [ "$after" -le "$before" ]; then
    stall=$((stall+1)); [ "$stall" -ge 3 ] && { echo "STALLED $after" > "$W/COHORT_DONE"; break; }
  else stall=0; fi
done
