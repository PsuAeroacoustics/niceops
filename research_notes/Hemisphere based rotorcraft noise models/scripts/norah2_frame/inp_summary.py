import numpy as np,sys
A="/Users/evg5332/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Outbox/NORAH2/NORAH2_V2.0.74_public/ARP/"
def load(path):
    rows=[];hdr=None
    for l in open(path):
        t=l.split()
        if t and t[0]=='TEM': hdr=t; continue
        if hdr and len(t)==len(hdr):
            try: rows.append([float(x.strip('"')) for x in t])
            except: pass
    return hdr,np.array(rows)
for p in sys.argv[1:]:
    hdr,R=load(A+p)
    c={h:i for i,h in enumerate(hdr)}
    t=R[:,0];x=R[:,1];y=R[:,2];H=R[:,c['HEADING']];ro=R[:,c['ROLL']]
    dx=np.gradient(x);dy=np.gradient(y)
    trk=np.degrees(np.arctan2(dx,dy))%360
    dtrk=np.gradient(np.unwrap(np.radians(trk)))
    print(p,hdr)
    # segment summary: group by (roll sign, rounded H-trk, turn sense)
    prev=None
    for i in range(len(t)):
        sense='R' if dtrk[i]>1e-5 else ('L' if dtrk[i]<-1e-5 else 'S')
        key=(sense,round(ro[i],2),round((H[i]-trk[i]+180)%360-180,1))
        if key!=prev:
            print(f"  t={t[i]:8.2f} x={x[i]:9.1f} y={y[i]:9.1f} z={R[i,3]:7.1f} track={trk[i]:6.1f} HEADING={H[i]:6.1f} H-trk={key[2]:6.1f} turn={sense} ROLL={ro[i]:6.2f} mode={R[i,4]:.0f} V={R[i,5]}")
            prev=key
