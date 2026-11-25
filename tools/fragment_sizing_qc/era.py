"""Band offsets and heights for runs before/after the 2012 ladder change; writes era.pkl.  ET550-R band-length check, docs/FRAGMENT_SIZING_QC.md."""
import pickle, numpy as np, collections
exec(open("region.py").read().split("run=collections")[0])
L=np.array([60,75,100,125,150,175,200,225,250,275,300,310,325,350,375,400,425,450,475,500,525,550],float); idx={b:i for i,b in enumerate(L)}
# era per run from median 125 residual
runr=collections.defaultdict(list); cache={}
for k,v in wells.items():
    v=sorted(v); bp=np.array([a[0] for a in v]); sc=np.array([a[1] for a in v]); h=np.array([a[2] for a in v])
    if len(bp)<20 or np.any(np.diff(sc)<=0): continue
    P=pred_matrix(sc); r=P@bp-bp; cache[k]=(bp,sc,h,P,r)
    if 125 in bp: runr[k[0]].append(r[list(bp).index(125)])
era={rn:("late" if np.median(v)>1 else "early" if np.median(v)<0.5 else "mid") for rn,v in runr.items()}
print(collections.Counter(era.values()))
G={e:[np.zeros((22,22)),np.zeros(22),0] for e in ["early","late"]}
H={e:collections.defaultdict(list) for e in G}
for k,(bp,sc,h,P,r) in cache.items():
    e=era.get(k[0]); 
    if e not in G: continue
    M=np.zeros((len(bp),22)); M[:,[idx[b] for b in bp]]=np.eye(len(bp))-P
    G[e][0]+=M.T@M; G[e][1]+=M.T@r; G[e][2]+=1
    # band height relative to well median; spacing scans/bp locally
    med=np.median(h)
    for b,hh in zip(bp,h): H[e][b].append(hh/med)
print("band   offset_early  offset_late   rel.height early/late")
D={e:np.linalg.solve(g[0]+0.2*g[2]*np.eye(22),g[1]) for e,g in G.items()}
for i,b in enumerate(L):
    print(f"{b:4.0f}   {D['early'][i]:+6.2f}      {D['late'][i]:+6.2f}       {np.median(H['early'][b]):.2f} / {np.median(H['late'][b]):.2f}")
pickle.dump(era,open("era.pkl","wb"))
