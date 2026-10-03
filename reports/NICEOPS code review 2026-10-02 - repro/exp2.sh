set -u
B=../tree/build
D="-d databases/S-76D_M3.nod -v databases/S-76D.json"
python3 - <<'PY'
import csv
rows=open("tracks/S-76D/197421AC.csv").read().split("\n")
hdr=rows[0]; body=[r for r in rows[1:] if r.strip()]
def write(name,f):
    out=[hdr]
    for k,r in enumerate(body):
        c=r.split(","); f(k,c); out.append(",".join(c))
    open(name,"w").write("\n".join(out)+"\n")
def flat(k,c): c[1]="0.0"
def desc(k,c): c[1]="%.6f"%(-0.02*k)
def nanalt(k,c):
    if k==3000: c[3]="nan"
write("integrated/t_flat.csv",flat); write("integrated/t_desc.csv",desc); write("integrated/t_nan.csv",nanalt)
PY
for t in flat desc nan; do
  echo "### E4/E7 track $t"; rm -f integrated/o_$t.csv
  $B/niceops $D -t integrated/t_$t.csv --gridstep 492 -L warn -o integrated/o_$t.csv > out_$t.txt 2>&1; echo "exit=$? stdout/stderr bytes=$(wc -c < out_$t.txt)"
  sed -n 4,6p integrated/o_$t.csv 2>/dev/null | cut -c1-80
done
echo "### E5 --trim past the end (ft)"
for s in 331 333 340; do $B/niceops $D -t tracks/S-76D/197421AC.csv --trim $s --gridstep 492 -L warn -o integrated/trim.csv > out_trim.txt 2>&1; echo "trim $s exit=$? out=$(head -c 200 out_trim.txt | tr '\n' ' ') log: $(grep -i 'skipping' log/niceops.log | tail -1)"; done
echo "### E6 grid counts: --grid -1500 1500 -1500 1500 --gridstep 100"
for u in ft m; do $B/niceops $D -t tracks/S-76D/197421AC.csv --units $u --grid -1500 1500 -1500 1500 --gridstep 100 -L warn -o integrated/g_$u.csv > /dev/null 2>&1; echo "units $u: rows=$(grep -vc '^[A-Za-z#]' integrated/g_$u.csv) $(head -c 0 integrated/g_$u.csv)"; done
echo "### E9 error visibility"
$B/niceops $D -t no_such_track.csv -o integrated/x.csv > e9.out 2>&1; echo "missing track: exit=$? bytes(stdout+stderr)=$(wc -c < e9.out)"; tail -2 log/niceops.log
$B/niceops --gridstep -5 > e9b.out 2>&1; echo "bad gridstep: exit=$? bytes=$(wc -c < e9b.out)"; head -c 200 e9b.out; tail -1 log/niceops.log
$B/niceops --no_such_option > e9c.out 2>&1; echo "bad option: exit=$? bytes=$(wc -c < e9c.out)"
