"""Add Christian's electronic score sheet (Excel/Lavage colorectal_V2.xls, Apr 2015; converted copy
excel_sources_Lavage_colorectal_V2.xlsx) to the paper-based Lavage answer key.
Per (2026-11-23): the paper is the primary record. A blank Excel cell where the paper has a result = not entered.
Excel legend: 0 WT, 1 MUT, 2 MUTlmf, 3 neg (in practice: no PCR product), 4 -PCR."""
import zipfile, csv, xml.etree.ElementTree as ET
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
K = "arstad_lavage_scoring.csv"
z = zipfile.ZipFile("excel_sources_Lavage_colorectal_V2.xlsx")
ss = ["".join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS)]
g = {}
for row in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).find("m:sheetData", NS):
    for c in row:
        v = c.find("m:v", NS)
        if v is not None and v.text is not None:
            g[c.get("r")] = (ss[int(v.text)] if c.get("t") == "s" else v.text).strip()
blocks = {"E": "F G H I J K", "L": "M N", "O": "P Q R S", "T": "U V W X", "Y": "Z AA AB"}  # sample col -> fragment cols
MEAN = {"0": "WT", "1": "MUT", "2": "LOW", "3": "NOPCR", "4": "NOPCR", "": "not entered"}
xl = {}
for scol, fcols in blocks.items():
    for c in fcols.split():
        fr = g[c + "7"]
        for r in range(8, 200):
            plate = 1 if r < 104 else 2
            well = f"{g.get(f'C{r}', '')}{int(float(g.get(f'D{r}', '0'))):02d}"
            xl[(fr, plate, well)] = (g.get(f"{c}{r}", ""), g.get(f"{scol}{r}", ""))
rows = list(csv.DictReader(open(K)))
cols = list(rows[0]) + [c for c in ("sample_id", "excel_code", "excel_meaning", "excel_vs_paper") if c not in rows[0]]
SAME = {("WT", "not entered")}
for p in rows:
    k = (p["fragment"], int(p["plate"] or 1), p["well"])
    if k not in xl:
        p.update(sample_id="", excel_code="", excel_meaning="", excel_vs_paper="fragment not in Excel"); continue
    code, sid = xl[k]; m = MEAN.get(code, f"unknown code {code}")
    paper = {"VERYLOW": "MUT", "UNCLEAR": "?"}.get(p["code"], p["code"])
    if m == paper or (paper == "WT" and code == ""):
        flag = "agree"
    elif code == "":
        flag = "not entered in Excel"
    elif (p["fragment"], k[1]) == ("P1F5", 1) and code == "3":
        flag = "Excel put plate 2's failed PCR under plate 1; run 120315_Lavage_P1F5_Plate1 has good signal, paper correct"
    else:
        flag = f"DIFFERENT: paper {p['code']}, Excel {m}"
    p.update(sample_id=sid, excel_code=code, excel_meaning=m, excel_vs_paper=flag)
with open(K, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
import collections
print(collections.Counter(p["excel_vs_paper"][:40] for p in rows).most_common())
for p in rows:
    if p["excel_vs_paper"].startswith("DIFF"): print("  p", p["page"], p["fragment"], p["plate"], p["well"], p["sample_id"], "|", p["raw_mark"][:40], "|", p["excel_vs_paper"])
