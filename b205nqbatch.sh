 for file in ./tracks/B205/*.csv; do
     [ -f "$file" ] || continue
     casename=$(basename "$file" .csv)
     build/niceops -t "$file" -v databases/B205.json -d databases/B205_M1.nod  --qsam=false  -o./integrated/B205/"$casename"_noq_LAE_grid.csv -f ./figures/B205/"$casename"_noq_LAE_grid.eps &
 done
 wait
