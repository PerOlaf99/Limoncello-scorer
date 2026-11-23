"""Testis 5, 3D LCM: answer key from a fractions workbook.  Usage: build_testis_3d.py P2F2|P2F1
Slice sheets S1..S24: Wells | M1..M5, per Per (2026-11-23): M1 = WT %, M2.. = mutant species %.
'pos' sheet: per spot, x/y/z in the same order (slice after slice, A1..H12)."""
import zipfile, re, csv, sys, xml.etree.ElementTree as ET
from collections import defaultdict
FRAG = sys.argv[1].upper()
SRC = {"P2F2": "/run/media/per/NVMe Disk/Excel/3D_Testis_P2F2_fractions.xlsx",
       "P2F1": "/run/media/per/NVMe Disk/Excel/P2F1_fractions3D.xlsx"}[FRAG]
CAT = "/run/media/per/NVMe Disk/limoncello-data/catalog/runs_current.csv"
OUT = f"/run/media/per/NVMe Disk/limoncello-data/answer_keys/testis5_3d_{FRAG.lower()}_scoring.csv"
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
z = zipfile.ZipFile(SRC)
ss = ["".join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}

def cells(target):
    x = ET.fromstring(z.read("xl/" + target.split("xl/")[-1]))
    for row in x.find("m:sheetData", NS):
        d = {}
        for c in row:
            v = c.find("m:v", NS)
            if v is None: continue
            val = ss[int(v.text)] if c.get("t") == "s" else v.text
            d[re.match(r"[A-Z]+", c.get("r")).group(0)] = val
        yield d

sheets = [(s.get("name"), rels[s.get(RID)]) for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS)]
slices = [(n, t) for n, t in sheets if re.fullmatch(r"S\d+", n)]
pos = [d for d in cells(dict(sheets)["pos"]) if d.get("A")]

runs = defaultdict(list)
for r in csv.DictReader(open(CAT)):
    m = re.search(r"Testis5?_3D_S(\d+)_" + FRAG + r"(?!\d)", r["name"], re.I)
    if m: runs[int(m.group(1))].append(r["path"])

rows, k = [], 0
for name, t in slices:
    data = list(cells(t))
    hdr, body = data[0], data[1:]
    mcols = [c for c in "BCDEF" if hdr.get(c, "").startswith("M")]
    for d in body:
        well = d.get("A"); p = pos[k] if k < len(pos) else {}; k += 1
        if p.get("A") != well: sys.exit(f"pos out of step at {name} {well} vs {p.get('A')}")
        vals = {hdr[c]: d.get(c, "").strip() for c in mcols}
        has = any(v != "" for v in vals.values())
        f = lambda v: round(float(v), 2) if v != "" else ""
        muts = [f(vals.get(f"M{i}", "")) for i in range(2, 6)]
        rows.append(dict(study="testis5_3d", fragment=FRAG, slice=name, well=well,
            x=p.get("B"), y=p.get("C"), z=p.get("D"),
            wt_percent=f(vals.get("M1", "")),
            mut1_percent=muts[0], mut2_percent=muts[1], mut3_percent=muts[2], mut4_percent=muts[3],
            mut_total_percent=round(sum(m for m in muts if m != ""), 2) if has else "",
            n_mutants=sum(1 for m in muts if m not in ("", 0, 0.0)) if has else "",
            scored="yes" if has else "no (blank in workbook)",
            runs=" | ".join(sorted(set(runs.get(int(name[1:]), [])))),
            source=f"Excel/{SRC.split('/')[-1]} sheet {name}"))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "rows ->", OUT)
