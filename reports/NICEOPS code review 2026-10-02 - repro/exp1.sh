set -u
B=../tree/build
G=golden/197421_LAE.csv
echo "### E1 golden_compare with NaN rows"
python3 - <<'PY'
lines=open("golden/197421_LAE.csv").read().split("\n")
out=[]
for i,l in enumerate(lines):
    if 4<=i<=20 and l.strip():
        out.append(",".join("nan" for _ in l.split(",")))
    else: out.append(l)
open("integrated/nan_golden.csv","w").write("\n".join(out))
PY
$B/golden_compare integrated/nan_golden.csv $G; echo "exit=$?"
echo "### E2 two trajectories, -i / -S"
T1=tracks/S-76D/197421AC.csv; T2=$(ls tracks/S-76D/*.csv | sed -n 5p); echo "T2=$T2"
rm -rf hist hist1 hist2; mkdir hist hist1 hist2
$B/niceops -d databases/S-76D_M3.nod -v databases/S-76D.json -t $T1 -t $T2 --gridstep 1640 -i hist/h -L warn -o integrated/m1.csv -o integrated/m2.csv; echo "exit=$?"
$B/niceops -d databases/S-76D_M3.nod -v databases/S-76D.json -t $T1 --gridstep 1640 -i hist1/h -L warn -o integrated/s1.csv; 
$B/niceops -d databases/S-76D_M3.nod -v databases/S-76D.json -t $T2 --gridstep 1640 -i hist2/h -L warn -o integrated/s2.csv
echo "files in hist: $(ls hist | wc -l)"; 
same1=0; same2=0; for f in hist/*; do b=$(basename $f); cmp -s $f hist1/$b && same1=$((same1+1)); cmp -s $f hist2/$b && same2=$((same2+1)); done
echo "hist files equal to track1-alone: $same1 ; equal to track2-alone: $same2"
echo "### E3 exception in worker thread (unwritable -i dir)"
rm -rf ro; mkdir ro; chmod 555 ro
$B/niceops -S state/ci_state.csv --gridstep 1640 -i ro/h -L warn -o integrated/e3.csv -j 1 > e3a.out 2>&1; echo "-j 1 exit=$?"; tail -2 e3a.out
$B/niceops -S state/ci_state.csv --gridstep 1640 -i ro/h -L warn -o integrated/e3.csv -j 4 > e3b.out 2>&1; echo "-j 4 exit=$?"; tail -2 e3b.out
chmod 755 ro
