"""Merge the rebuild_*.csv scans into runs_current.csv and rebuild by_kind/.
Old runs_current.csv and by_kind/ are kept with a _20261119 suffix."""
import csv, os, sys, collections
from pathlib import Path

D = Path("/run/media/per/NVMe Disk")
L = D / "limoncello-data"
C = L / "catalog"
STAMP = "20261119"

cols = next(csv.reader(open(C / "runs_current.csv")))
rows = []
for part in ("rebuild_disk.csv", "rebuild_sids.csv", "rebuild_mybook.csv"):
    for r in csv.DictReader(open(C / part, encoding="utf-8", errors="replace")):
        rows.append({k: r.get(k, "") for k in cols})

old = list(csv.DictReader(open(C / "runs_current.csv", encoding="utf-8", errors="replace")))
oldp, newp = {r["path"] for r in old}, {r["path"] for r in rows}
added, gone = newp - oldp, oldp - newp

def source(p):
    rel = Path(p).relative_to(D).parts
    return rel[1] if rel[0] == "limoncello-data" else rel[0]

dup = collections.Counter((r["fingerprint"], r.get("machine", "")) for r in rows)
print(f"rows {len(rows)} (was {len(old)}): +{len(added)} new paths, -{len(gone)} gone")
print("by source of new paths:", collections.Counter(source(p) for p in added).most_common())
print("new rows whose fingerprint+machine matches another row:",
      sum(1 for r in rows if r["path"] in added and dup[(r['fingerprint'], r.get('machine', ''))] > 1))
for p in sorted(gone)[:10]:
    print("  gone:", p)

if "--write" not in sys.argv:
    sys.exit()

os.rename(C / "runs_current.csv", C / f"runs_current_{STAMP}.csv")
with open(C / "runs_current.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

bk = L / "by_kind"
os.rename(bk, L / f"by_kind_{STAMP}")
bk.mkdir()
for r in rows:
    d = bk / r["kind"]; d.mkdir(exist_ok=True)
    link = d / r["name"]
    if link.is_symlink() or link.exists():
        link = d / f"{r['name']}__{source(r['path'])}_{r.get('machine', '')}"
        i = 2
        while link.is_symlink():
            link = d / f"{r['name']}__{source(r['path'])}_{r.get('machine', '')}_{i}"; i += 1
    link.symlink_to(r["path"])
print("links:", sum(1 for _ in bk.rglob("*") if _.is_symlink()))
print("kinds:", collections.Counter(r["kind"] for r in rows).most_common())
