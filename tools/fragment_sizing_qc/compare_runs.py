"""Per run: kit/channel and clean wells, earlier engine (auto + reruns) vs a new validate_auto.py output.  python compare_runs.py NEW_wells.csv (run in the QC data folder)."""
import csv, collections, sys
D = "/run/media/per/NVMe Disk/limoncello-data/fragment_sizing/qc_20261125/"
def load(files):
    R = {}
    for f in files:
        rows = collections.defaultdict(list)
        for r in csv.DictReader(open(D + f)): rows[r["run"]].append(r)
        R.update(rows)
    return R
def label(rows):
    d = rows[0]["detection"] or ""
    if not d or d.startswith(("no size", "run error")): return "none"
    return d.split(" in channel")[0].replace("ladder ", "").replace("MegaBACE ", "") + " ch" + rows[0]["channel"] + (" ?" if "uncertain" in d else "")
clean = lambda rows: sum(int(r["clean"] or 0) for r in rows)
old = load(["auto_wells.csv", "rerun_wells.csv", "rerun2_wells.csv"])
new = load([sys.argv[1] if len(sys.argv) > 1 else "v4_wells.csv"])
same = changed = 0; tot = collections.Counter()
for run in sorted(set(old) | set(new)):
    o, n = old.get(run), new.get(run)
    if not n: continue
    lo, ln = (label(o) if o else "?"), label(n)
    co, cn = (clean(o) if o else 0), clean(n)
    tot["old"] += co; tot["new"] += cn
    if lo == ln and abs(co - cn) <= max(2, 0.1 * max(co, 1)):
        same += 1; continue
    changed += 1
    print(f"{run[:46]:46s} {lo:28s} {co:3d} -> {ln:28s} {cn:3d}")
print(f"\nruns unchanged {same}, changed {changed}; clean wells {tot['old']} -> {tot['new']}")
