#!/bin/bash
# Write NICEOPS trajectory state (-S) for every shipped track of a vehicle.
set -e
W=/Users/evg5332/Code/D/NICEOPS/.claude/worktrees/hemisphere-review
B=/Users/evg5332/Code/D/NICEOPS/build/niceops
OUT=/Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend/states
mkdir -p "$OUT/$1"
cd /Users/evg5332/.claude/jobs/a090ec0c/tmp/dbblend
for f in "$W"/tracks/"$1"/*AC.csv; do
  n=$(basename "$f" AC.csv)
  "$B" -d "$W/databases/$2" -v "$W/databases/$3" -t "$f" -o "$OUT/$1/lae_$n.csv" \
       -S "$OUT/$1/state_$n.csv" --gridstep 4000 --no-spectral -l "$OUT/$1/log_$n.log" -L warn >/dev/null
done
ls "$OUT/$1" | grep -c state_
