#!/bin/bash
# Full-grid NICEOPS run + numpy blend-variant footprint for every track of a vehicle.
# usage: batch.sh VEHICLE DBNAME JSON
set -e
W=/Users/evg5332/Code/D/NICEOPS/.claude/worktrees/hemisphere-review
B=/Users/evg5332/Code/D/NICEOPS/build/niceops
P=/Users/evg5332/Code/Python/niceops-validation/.venv/bin/python
D=/Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend
OUT=$D/fp/$1
mkdir -p "$OUT"
cd "$D"
for f in "$W"/tracks/"$1"/*AC.csv; do
  n=$(basename "$f" AC.csv)
  "$B" -d "$W/databases/$2.nod" -v "$W/databases/$3" -t "$f" -o "$OUT/lae_$n.csv" \
       -S "$OUT/state_$n.csv" --no-spectral -l "$OUT/log_$n.log" -L warn >/dev/null
  echo "== $n" >> "$OUT/summary.txt"
  "$P" footprint.py "$OUT/state_$n.csv" "$OUT/lae_$n.csv" "$OUT/fp_$n.npz" "$2" >> "$OUT/summary.txt" 2>&1
done
grep -c validation "$OUT/summary.txt"
