"""3D LCM mutant fractions for the breast case 170 (T170 tumour, N170 lymph node, L170 liver metastasis)
the bladder case 231 (T231 tumour, M231 = LN231 lymph node) and liver tumour #13 (Liver13), from Excel/ (2016-17).

Per slice sheet: Wells | Fractions.  Fractions = WT share (0-1): the workbooks' own formulas give
WT % = Fractions*100 and mutant % = 100 - WT %  (M1 = WT, as in the testis files; Per 2026-11-23).
x/y from the well grid (x by row, y by column, 180 apart).  z_workbook = the workbook's own z, which spaces the
slices evenly whatever the gap between them; z_by_section = -(section-1) * 12, from T170's consecutive sections."""
import zipfile, re, csv, collections, xml.etree.ElementTree as ET

EX = "/run/media/per/NVMe Disk/Excel/"
CAT = "/run/media/per/NVMe Disk/limoncello-data/catalog/runs_current.csv"
OUT = "/run/media/per/NVMe Disk/limoncello-data/answer_keys/tumor_3d_scoring.csv"
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"

BOOKS = [  # file, sample, tissue, run-name prefixes, 3D/stack sheet
    ("Breast_T170_3D.xlsx", "T170", "breast tumour", ["T170"], "3D_breast_P1F1"),
    ("Node_breast_N170_3D.xlsx", "N170", "breast case, lymph node", ["N170", "Node170"], "Breast_Lymphnode_3D"),
    ("Liver_Met_breast_3D.xlsx", "L170", "breast case, liver metastasis", ["L170"], "Breast_Liver_Met_3D"),
    ("T231_3d_P2F6_P5F6.xlsx", "T231", "bladder tumour", ["T231"], None),
    ("M231_3D_P2F6.xlsx", "LN231", "bladder case, lymph node (files say M231/Met231)", ["M231", "Met231"], "Stack"),
    ("M231_3D_P5F6.xlsx", "LN231", "bladder case, lymph node (files say M231/Met231)", ["M231", "Met231"], "Stack"),
    ("#13_Liver_P1F15_ALL.xlsx", "Liver13", "liver tumour #13", ["Liver13"], "Sheet1"),
]
# parts of these sheets hold exactly the same values as each other (~30 wells each): likely copy-paste
SUSPECT = {("LN231", "P2F6", 41), ("LN231", "P2F6", 46), ("LN231", "P5F6", 46)}


def book(path):
    z = zipfile.ZipFile(path)
    ss = ["".join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)] \
        if "xl/sharedStrings.xml" in z.namelist() else []
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for s in ET.fromstring(z.read("xl/workbook.xml")).find("m:sheets", NS):
        t = rels[s.get(RID)]
        if "chart" in t: continue
        sd = ET.fromstring(z.read("xl/" + t.split("xl/")[-1])).find("m:sheetData", NS)
        rows = []
        for row in (sd if sd is not None else []):
            d = {}
            for c in row:
                v = c.find("m:v", NS)
                if v is not None and v.text is not None:
                    d[re.sub(r"\d", "", c.get("r"))] = (ss[int(v.text)] if c.get("t") == "s" else v.text).strip()
            rows.append(d)
        out[s.get("name")] = rows
    return out


def num(v):
    try:
        f = float(v)
        return None if f != f else f
    except (TypeError, ValueError):
        return None


def wellkey(w):
    m = re.fullmatch(r"([A-H])0?(\d{1,2})", w or "")
    return (m.group(1), int(m.group(2))) if m else None


runs = defaultdict = collections.defaultdict(list)
for r in csv.DictReader(open(CAT)):
    runs[r["name"]].append(r["path"])


def find_runs(prefixes, sl, frag):
    hits = []
    for name, paths in runs.items():
        for p in prefixes:
            if re.match(rf"{p}_?(3D_)?S{sl}[_\-]", name, re.I) and re.search(rf"{frag}(?!\d)", name, re.I):
                hits += paths
    return sorted(set(hits))


rows = []
for fn, sample, tissue, prefixes, stackname in BOOKS:
    B = book(EX + fn)
    slices = [(n, r) for n, r in B.items() if re.search(r"_S\d+_|_S\d+$", n)]
    # z per slice: Stack sheets name the slice in col C; 3D sheets are blocks of 96 in slice-sheet order
    zof = {}
    st = B.get(stackname) if stackname else None
    if st and any(r.get("C", "").startswith("S") for r in st):
        for r in st[1:]:
            if r.get("C", "").startswith("S") and num(r.get("G")) is not None:
                zof.setdefault(int(r["C"][1:]), num(r["G"]))
    elif st:
        for i, (n, _) in enumerate(slices):
            blk = st[1 + 96 * i: 1 + 96 * (i + 1)]
            if blk: zof[int(re.search(r"_S(\d+)", n).group(1))] = num(blk[0].get("D"))
    for n, data in slices:
        sl = int(re.search(r"_S(\d+)", n).group(1))
        frag = re.search(r"P\dF\d+", n).group(0)
        rr = find_runs(prefixes, sl, frag)
        for d in data[1:]:
            wk = wellkey(d.get("A"))
            if not wk: continue
            f = num(d.get("B"))
            wt = round(f * 100, 2) if f is not None else ""
            rows.append(dict(
                sample=sample, tissue=tissue, fragment=frag, slice=f"S{sl}", well=f"{wk[0]}{wk[1]:02d}",
                x=1350 - 180 * "ABCDEFGH".index(wk[0]), y=2070 - 180 * (wk[1] - 1), z_workbook=zof.get(sl, ""), z_by_section=-(sl - 1) * 12,
                wt_percent=wt, mut_percent=round(100 - wt, 2) if wt != "" else "",
                scored="yes" if wt != "" else "no (NaN/blank in workbook)",
                flag="SUSPECT: ~30 wells hold exactly the same values as another sheet (copy-paste?)"
                     if (sample, frag, sl) in SUSPECT else "",
                runs=" | ".join(rr), source=f"Excel/{fn} sheet {n}"))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "rows ->", OUT)
c = collections.defaultdict(lambda: [0, 0, set(), set()])
for r in rows:
    k = (r["sample"], r["fragment"]); c[k][0] += 1; c[k][1] += r["scored"] == "yes"
    (c[k][2] if r["runs"] else c[k][3]).add(r["slice"])
for k, (n, s, wr, nr) in c.items():
    print(f"  {k[0]:6s} {k[1]:5s} spots {n:5d} scored {s:5d}  slices with run {len(wr):2d}  without: {sorted(nr, key=lambda v: int(v[1:]))}")
