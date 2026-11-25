"""Collect clean ET550-R wells and their band scans from the QC CSVs into et550.pkl (run first, in a work folder).  ET550-R band-length check, docs/FRAGMENT_SIZING_QC.md."""
import csv, collections, pickle, numpy as np
D="/run/media/per/NVMe Disk/limoncello-data/fragment_sizing/qc_20261125/"
rd=lambda f: list(csv.DictReader(open(D+f)))
W=rd("auto_wells.csv"); B=rd("auto_bands.csv")
for f in ["rerun","rerun2"]:
    rw=rd(f+"_wells.csv"); rb=rd(f+"_bands.csv"); runs={r["run"] for r in rw}
    W=[r for r in W if r["run"] not in runs]+rw; B=[r for r in B if r["run"] not in runs]+rb
clean={(r["run"],r["well"]):r for r in W if r["ladder"]=="MegaBACE ET550-R" and r["clean"]=="1"}
print("clean ET550 wells",len(clean),"runs",len({k[0] for k in clean}),"rescued",sum(r["rescued"]=="1" for r in clean.values()))
print(collections.Counter(r["kind"] for r in clean.values()))
wells=collections.defaultdict(list)
for r in B:
    k=(r["run"],r["well"])
    if k in clean and r["ladder"]=="MegaBACE ET550-R":
        wells[k].append((float(r["band_bp"]),float(r["scan"]),float(r["height"]),float(r["residual_bp"]) if r["residual_bp"] else np.nan))
print("wells with bands",len(wells))
pickle.dump((dict(wells),clean),open("et550.pkl","wb"))
res=collections.defaultdict(list)
for v in wells.values():
    for bp,sc,h,rr in v: res[bp].append(rr)
print("band  n  medianLOO  IQR")
for bp in sorted(res):
    a=np.array(res[bp]); a=a[~np.isnan(a)]
    print(f"{bp:5.0f} {len(a):5d} {np.median(a):7.2f}  {np.percentile(a,25):6.2f} {np.percentile(a,75):6.2f}")
