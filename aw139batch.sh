  rm -f ./integrated/AW139/metrics.csv
 for file in ./tracks/AW139/*.csv; do
     [ -f "$file" ] || continue
     casename=$(basename "$file" .csv)
     ./niceops -t "$file" -v databases/AW139.json -d databases/AW139_M1.nod  --qsam=true  -o./integrated/AW139/"$casename"_LAE_grid.csv -f ./figures/AW139/"$casename"_LAE_grid.eps -F ./figures/AW139/"$casename".eps -S ./state/AW139/"$casename"_state.csv -m ./integrated/AW139/metrics.csv &
 done
 wait
