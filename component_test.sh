
file=./tracks/B205/190353AC.csv
casename=$(basename "$file" .csv)
./build/niceops -t "$file" -v databases/tandem.json -d databases/AW139_M1.nod  --qsam=true  -o./integrated/tandem/"$casename"_LAE_grid.csv -f ./figures/tandem/"$casename"_LAE_grid.eps -F ./figures/tandem/"$casename".eps -S ./state/tandem/"$casename"_state.csv -m ./integrated/tandem/metrics.csv -L trace

