 rm -f ./integrated/B205/metrics.csv
 for file in ./tracks/B205/*.csv; do
     [ -f "$file" ] || continue
     casename=$(basename "$file" .csv)
     ./build/niceops -t "$file" -v databases/B205.json -d databases/B205_M1.nod  --qsam=true  -o./integrated/B205/"$casename"_LAE_grid.csv -f ./figures/B205/"$casename"_LAE_grid.eps -F ./figures/B205/"$casename".eps -S ./state/B205/"$casename"_state.csv -m ./integrated/B205/metrics.csv --smooth=5,1e6 &
 done
 wait
