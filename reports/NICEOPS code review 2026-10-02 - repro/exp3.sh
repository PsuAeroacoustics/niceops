set -u
B=../tree/build
D="-d databases/S-76D_M3.nod -v databases/S-76D.json"
python3 - <<'PY'
rows=open("tracks/S-76D/197421AC.csv").read().split("\n")
hdr=rows[0]; body=[r for r in rows[1:] if r.strip()]
out=[hdr]
for k,r in enumerate(body):
    c=r.split(",")
    if k==3000: c[5]="nan"
    out.append(",".join(c))
open("integrated/t_nanx.csv","w").write("\n".join(out)+"\n")
PY
echo "### E7 one NaN x in a 6551-row track"
$B/niceops $D -t integrated/t_nanx.csv --gridstep 492 -L warn -o integrated/o_nanx.csv > out_nanx.txt 2>&1; echo "exit=$? bytes=$(wc -c < out_nanx.txt)"
python3 - <<'PY'
import re
vals=[]
for l in open("integrated/o_nanx.csv"):
    for t in l.strip().split(","):
        try: vals.append(float(t))
        except: pass
import math
n=len(vals); bad=sum(1 for v in vals if math.isnan(v) or math.isinf(v))
print(f"cells={n} non-finite={bad}")
PY
echo "### E6 columns"
for u in ft m; do head -4 integrated/g_$u.csv | tail -1 | tr ',' '\n' | grep -c .; done
echo "### E10 DNL on a WGS84 track with no clock"
python3 - <<'PY'
rows=open("tracks/S-76D/197421AC.csv").read().split("\n")
body=[r.split(",") for r in rows[1:] if r.strip()]
out=["time,lat,lon,alt"]
for k,c in enumerate(body):
    if k%10==0: out.append("%s,%s,%s,%s"%(c[1].strip(),c[2].strip(),c[3].strip(),c[4].strip()))
open("integrated/geo_noclock.csv","w").write("\n".join(out)+"\n")
PY
$B/niceops $D -t integrated/geo_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 492 -L info -o integrated/dnl_a.csv -D integrated/dnl_dnl.csv > out_dnl.txt 2>&1; echo "exit=$? bytes=$(wc -c < out_dnl.txt)"
grep -i 'night\|utc\|clock\|day' log/niceops.log | tail -4
python3 - <<'PY'
def read(p):
    v=[]
    for l in open(p):
        for t in l.strip().split(","):
            try: v.append(float(t))
            except: pass
    return v
a=read("integrated/dnl_a.csv"); d=read("integrated/dnl_dnl.csv")
import statistics
print("n",len(a),len(d)); 
pairs=[(x,y) for x,y in zip(a,d) if x==x and y==y and abs(x)<1e3 and abs(y)<1e3]
print("DNL-SEL median diff", statistics.median(y-x for x,y in pairs) if pairs else None)
PY
