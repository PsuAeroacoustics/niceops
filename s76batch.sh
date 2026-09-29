 rm -f ./integrated/S-76D/metrics.csv
 for file in ./tracks/S-76D/*.csv; do
     [ -f "$file" ] || continue
     casename=$(basename "$file" .csv)
     build/niceops -t "$file" -v databases/S-76D.json -d databases/S-76D_M3.nod  --qsam=true  -o./integrated/S-76D/"$casename"_LAE_grid.csv -f ./figures/S-76D/"$casename"_LAE_grid.eps -F ./figures/S-76D/"$casename".eps -S ./state/S-76D/"$casename"_state.csv -m ./integrated/S-76D/metrics.csv &
 done
 wait
