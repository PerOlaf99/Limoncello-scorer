"""Per-run offset of the 125 bp band, with run speed and data folder (helpers are reused by era.py).  ET550-R band-length check, docs/FRAGMENT_SIZING_QC.md."""
import pickle, numpy as np, collections, re
wells,clean=pickle.load(open("et550.pkl","rb"))
K=3;DEG=3
def pred_matrix(sc):
    n=len(sc); P=np.zeros((n,n))
    for i in range(n):
        nb=[j for j in range(max(0,i-K),min(n,i+K+1)) if j!=i]
        if len(nb)<DEG+3: nb=sorted([j for j in range(n) if j!=i],key=lambda j:abs(j-i))[:DEG+3]
        P[i,nb]=np.linalg.pinv(np.vander((sc[nb]-sc[i])/100.0,DEG+1))[-1]
    return P
run=collections.defaultdict(lambda: collections.defaultdict(list)); allw=collections.defaultdict(list)
for k,v in wells.items():
    v=sorted(v); bp=np.array([a[0] for a in v]); sc=np.array([a[1] for a in v]); h=np.array([a[2] for a in v])
    if len(bp)<20 or np.any(np.diff(sc)<=0): continue
    r=pred_matrix(sc)@bp-bp; d=dict(zip(bp,r))
    if 125 in d:
        allw[125].append(d[125]); run[k[0]]["r125"].append(d[125])
        run[k[0]]["sc125"].append(sc[list(bp).index(125)])
        run[k[0]]["sp"].append((sc[list(bp).index(150)]-sc[list(bp).index(100)])/50 if 150 in d and 100 in d else np.nan)
        run[k[0]]["path"]=clean[k]["path"]; run[k[0]]["rescued"].append(clean[k]["rescued"]=="1")
a=np.array(allw[125]); print("125 per-well hist:",np.histogram(a,bins=[-5,-1,0,0.5,1,1.5,2,2.5,3,5])[0])
rows=[]
for rn,d in run.items():
    m=re.search(r"MB\d|Paulo|Data_\d\d_\d\d_\d\d",d["path"]); 
    rows.append((np.median(d["r125"]),rn,len(d["r125"]),np.median(d["sc125"]),np.nanmedian(d["sp"]),"/".join(d["path"].split("/")[-3:-1])))
rows.sort()
for r in rows: print(f"{r[0]:+5.2f} n{r[2]:3d} scan125 {r[3]:6.0f} scans/bp {r[4]:5.2f}  {r[1][:40]:40s} {r[5][:45]}")
