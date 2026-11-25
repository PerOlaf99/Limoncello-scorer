"""Per-band length offsets (local cubic, ridge-solved), pooled and per run, from et550.pkl.  ET550-R band-length check, docs/FRAGMENT_SIZING_QC.md."""
import pickle, numpy as np, collections
wells,clean=pickle.load(open("et550.pkl","rb"))
L=np.array([60,75,100,125,150,175,200,225,250,275,300,310,325,350,375,400,425,450,475,500,525,550],float)
idx={b:i for i,b in enumerate(L)}
K=3; DEG=3
def pred_matrix(sc):
    n=len(sc); P=np.zeros((n,n))
    for i in range(n):
        nb=[j for j in range(max(0,i-K),min(n,i+K+1)) if j!=i]
        if len(nb)<DEG+3: nb=sorted([j for j in range(n) if j!=i],key=lambda j:abs(j-i))[:DEG+3]
        x=(sc[nb]-sc[i])/100.0
        P[i,nb]=np.linalg.pinv(np.vander(x,DEG+1))[-1]
    return P
AtA=np.zeros((22,22)); Aty=np.zeros(22); runAt={}
perrun=collections.defaultdict(lambda:[np.zeros((22,22)),np.zeros(22),0])
R=collections.defaultdict(list)
for k,v in wells.items():
    v=sorted(v); bp=np.array([a[0] for a in v]); sc=np.array([a[1] for a in v])
    if len(bp)<20 or np.any(np.diff(sc)<=0): continue   # near-complete wells only
    P=pred_matrix(sc); r=P@bp-bp
    # true = listed+d smooth => P(bp+d)=bp+d => (I-P)d = r
    M=np.zeros((len(bp),22)); M[:,[idx[b] for b in bp]]=np.eye(len(bp))-P
    AtA+=M.T@M; Aty+=M.T@r
    pr=perrun[k[0]]; pr[0]+=M.T@M; pr[1]+=M.T@r; pr[2]+=1
    for b,ri in zip(bp,r): R[b].append(ri)
nw=sum(p[2] for p in perrun.values()); print("wells used",nw,"runs",len(perrun))
for lam_per_well in [0.05,0.2,1.0]:
    d=np.linalg.solve(AtA+lam_per_well*nw*np.eye(22),Aty)
    print(f"ridge {lam_per_well}: "+" ".join(f"{b:.0f}:{x:+.2f}" for b,x in zip(L,d)))
# per-run ridge estimates -> spread across runs
lam=0.2
D=np.array([np.linalg.solve(p[0]+lam*p[2]*np.eye(22),p[1]) for p in perrun.values() if p[2]>=10])
print("runs with >=10 wells:",len(D))
print("band  pooled  run-median  run-IQR        runs>+0.5  runs<-0.5")
d=np.linalg.solve(AtA+lam*nw*np.eye(22),Aty)
for i,b in enumerate(L):
    c=D[:,i]; print(f"{b:4.0f}  {d[i]:+6.2f}   {np.median(c):+6.2f}   [{np.percentile(c,25):+5.2f},{np.percentile(c,75):+5.2f}]   {np.sum(c>0.5):3d}  {np.sum(c<-0.5):3d}")
# single-band hypothesis: which one band shift explains most of the pooled residual
y0=Aty; base=0
print("single-band shift fits (reduction of residual SS, per well):")
SS=[]
for i in range(22):
    a=AtA[i,i]; dk=Aty[i]/a; SS.append((Aty[i]**2/a/nw,L[i],dk))
for s,b,dk in sorted(SS,reverse=True)[:6]: print(f"  {b:.0f} bp shift {dk:+.2f} explains {s:.3f} bp^2/well")
