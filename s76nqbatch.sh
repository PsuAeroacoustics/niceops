 for file in ./tracks/S-76D/*.csv; do
     [ -f "$file" ] || continue
     casename=$(basename "$file" .csv)
     build/niceops -t "$file" -v databases/S-76D.json -d databases/S-76D_M3.nod  --qsam=false  -o./integrated/S-76D/"$casename"_noq_LAE_grid.csv -f ./figures/S-76D/"$casename"_noq_LAE_grid.eps &
 done
 wait
