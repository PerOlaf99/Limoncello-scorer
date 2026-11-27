"""Build one SNP genotype answer key from Per's hand-scored sheets.

Sources (see docs/GENOTYPING_ENGINE.md):
  * ~/Dokumenter/Genotype_*.ods      plate grids, 0/1/2/3 (+ notes AI, Flagg, no STD)
  * ~/CE/90_genotyping_tools/Genotyping.xlsx   training labels of the old GUI,
    well/code pairs, 0/1/2/4 (4 = heterozygote)
  * tests/data/rs1695_expected.csv   the T9 plate

Writes run,well,call,ai,flag,note,source with call in hom1/hom2/het/none.

    python build_snp_key.py OUT.csv
"""
import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

HOME = Path.home()
REPO = Path(__file__).resolve().parents[2]
CODES = {"0": "none", "1": "hom1", "2": "hom2", "3": "het", "4": "het"}

# Sheet title -> run folder.  "Run1" sheets name the first run of the day;
# T1 and ABCC2 N11 have two candidate runs, checked against the data
# (docs/GENOTYPING_ENGINE.md).
RUN_NAMES = {
    "OY_rs1695_T1_230910Run1": "OY_rs1695_T1_230910Run01",
    "OY_rs1695_T5_270910Run1": "OY_rs1695_T5_270910Run01",
    "OY_IL10_N2_080910Run1": "OY_IL10_N2_080910Run01",
    "OY_ABCC2_N1_160910Run01": "OY_ABCC2_N1_160910Run01",
    "OY_CYBA_N1_200910Run01": "OY_CYBA_N1_200910Run01",
    "OY_CYBA_N1_200910Run2": "OY_CYBA_N1_200910Run02",
    "OY_GSTA1_N1_060910Run1": "OY_GSTA1_N1_060910Run01",
    "OY_RS1695_N1210910Run01": "OY_rs1695_N1_210910Run01",
    "OY_CYBA_N2_60gr": "OY_CYBA_N2_60gr_210910Run01",
    "OY_CYBA_N2_59gr": "OY_CYBA_N2_59gr_200910Run01",
    "Genotype_ABCC2_N11": "OY_ABCC2_N11_010910Run01",
}


def parse_code(text):
    """'3, AI, no het, flagg' -> (call, ai, flag, note)."""
    t = (text or "").strip()
    m = re.match(r"\s*([0-4])\s*(\??)", t)
    if not m:
        return None
    low = t.lower()
    return (CODES[m.group(1)], "ai" in low.replace("baad", ""),
            "flag" in low or bool(m.group(2)), t)


# ---------------------------------------------------------------- .ods grids
_T = "urn:oasis:names:tc:opendocument:xmlns:table:1.0"
_X = "urn:oasis:names:tc:opendocument:xmlns:text:1.0"


def _ods_rows(path):
    root = ET.fromstring(zipfile.ZipFile(path).read("content.xml"))
    for row in root.iter(f"{{{_T}}}table-row"):
        cells = []
        for c in row.findall(f"{{{_T}}}table-cell"):
            n = int(c.get(f"{{{_T}}}number-columns-repeated", "1"))
            txt = " ".join("".join(p.itertext()) for p in c.findall(f"{{{_X}}}p"))
            cells += [txt] * min(n, 30)
        yield cells


def read_ods(path):
    """Plate grids: a title row (run), a header 'WELL Column 1..12', rows A-H."""
    out, title = [], Path(path).stem
    for cells in _ods_rows(path):
        first = cells[0].strip() if cells else ""
        if first.startswith("OY_"):
            title = first
            continue
        if len(cells) > 2 and cells[1].strip() in "ABCDEFGH" and cells[1].strip():
            row = cells[1].strip()
            for col in range(1, 13):
                txt = cells[col + 1] if col + 1 < len(cells) else ""
                p = parse_code(txt)
                if p:
                    out.append((title, f"{row}{col:02d}") + p)
    return out


# ---------------------------------------------------------------- .xlsx pairs
_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col_index(letters):
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx(path):
    """Runs as a title cell in column A, then rows of (well, code) pairs."""
    z = zipfile.ZipFile(path)
    ss = ["".join(t.itertext()) for t in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(f"{_M}si")]
    grid = {}
    for c in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(f"{_M}c"):
        v = c.find(f"{_M}v")
        if v is None:
            continue
        val = ss[int(v.text)] if c.get("t") == "s" else v.text
        col, row = re.match(r"([A-Z]+)(\d+)", c.get("r")).groups()
        grid.setdefault(int(row), {})[_col_index(col)] = val
    out, title = [], None
    for r in sorted(grid):
        cells = grid[r]
        if str(cells.get(0, "")).startswith("OY_"):
            title = cells[0]
            continue
        for k, v in cells.items():
            if re.fullmatch(r"[A-H]\d\d", str(v)):
                p = parse_code(str(cells.get(k + 1, "")))
                if p:
                    out.append((title, v) + p)
    return out


def read_t9():
    out = []
    with open(REPO / "tests/data/rs1695_expected.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            note = r.get("note") or ""
            out.append(("OY_rs1695_T9_270910Run01", r["well"], r["expected"],
                        "AI" in note, "Borderline" in note, note))
    return out


def main(out_path):
    rows = {}

    def add(items, source):
        for title, well, call, ai, flag, note in items:
            run = RUN_NAMES.get(title, title)
            rows[(run, well)] = (run, well, call, int(ai), int(flag), note, source)

    xlsx = HOME / "CE/90_genotyping_tools/Genotyping.xlsx"
    add(read_xlsx(xlsx), "Genotyping.xlsx")
    for ods in sorted((HOME / "Dokumenter").glob("Genotype_*.ods")):
        add(read_ods(ods), ods.name)
    add(read_t9(), "tests/data/rs1695_expected.csv")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "well", "call", "ai", "flag", "note", "source"])
        for k in sorted(rows):
            w.writerow(rows[k])
    from collections import Counter
    per = Counter(r[0] for r in rows.values())
    for run, n in sorted(per.items()):
        calls = Counter(r[2] for r in rows.values() if r[0] == run)
        print(f"{run:34s} {n:3d}  {dict(calls)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "snp_genotypes.csv")
