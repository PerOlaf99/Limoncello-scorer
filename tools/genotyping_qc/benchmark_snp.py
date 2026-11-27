"""Score SNP genotyping engines against Per's hand calls (snp_genotypes.csv).

    python benchmark_snp.py [--engine auto|seeded] [--runs SUBSTR ...]

Reads the raw OY runs, genotypes every well, and prints per plate the
agreement and the confusion of hand call -> engine call.  Writes the per-well
table to --out (default: scratch CSV).
"""
import argparse
import collections
import csv
import sys
import time
import warnings
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
warnings.filterwarnings("ignore")

from analyzer_core import load_trace  # noqa: E402
import genotyping as g  # noqa: E402

DATA = Path("/run/media/per/78B0C7DE1FA7081C/OY")
KEY = Path("/run/media/per/NVMe Disk/limoncello-data/answer_keys/snp_genotypes.csv")
MARKS = Path("/run/media/per/NVMe Disk/limoncello-data/answer_keys/snp_std_marks.json")
BASE_ORDER = "ACTG"
NORM = {"hom-1": "hom1", "hom-2": "hom2", "hom1": "hom1", "hom2": "hom2", "het": "het", "no-call": "none"}


def load_key():
    key = collections.defaultdict(dict)
    with open(KEY, newline="") as fh:
        for r in csv.DictReader(fh):
            key[r["run"]][r["well"]] = r
    return key


def load_plate(run):
    docs = {}
    for f in sorted((DATA / run).glob("*.rsd")):
        try:
            docs[f.stem] = load_trace(f, base_order=BASE_ORDER)
        except Exception as e:  # a broken file is one well, not the plate
            print(f"  {run} {f.stem}: {e}")
    return docs


def engine_unseeded(docs, run):
    return {w: g.auto_genotype(d, base_order=BASE_ORDER, run_name=run)
            for w, d in docs.items()}


def engine_plate(docs, run):
    import snp_genotyping as sg
    return sg.genotype_plate(docs, base_order=BASE_ORDER, run_name=run)


def engine_plate_matrix(docs, run):
    """The plate engine on traces corrected by the run's own matrix."""
    import snp_genotyping as sg
    return sg.genotype_plate(docs, base_order=BASE_ORDER, run_name=run, matrix="peaks")


ENGINES = {"unseeded": engine_unseeded, "plate": engine_plate,
           "plate_matrix": engine_plate_matrix}


def score_placement(engine, sel):
    """Standard placed where Per marked it: every mark within 6 scans of a band."""
    import json
    marks = json.load(open(MARKS))
    have = {p.name for p in DATA.iterdir()}
    tot = [0, 0]
    for run in sorted(marks):
        if run not in have or (sel and not any(x in run for x in sel)):
            continue
        rows = ENGINES[engine](load_plate(run), run)
        ok = n = 0
        for w, m in marks[run].items():
            if w not in rows:
                continue
            s = [int(x) for x in str(rows[w].get("std_scans") or "").split("/") if x]
            n += 1
            ok += bool(s) and all(min(abs(x - p) for p in s) <= 6 for x, _ in m)
        tot[0] += ok
        tot[1] += n
        print(f"  placement {run:34s} {ok:3d}/{n:3d}", flush=True)
    print(f"PLACEMENT {tot[0]}/{tot[1]} = {100*tot[0]/max(1,tot[1]):.1f}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="plate", choices=sorted(ENGINES))
    ap.add_argument("--runs", nargs="*", default=[])
    ap.add_argument("--placement", action="store_true",
                    help="also score standard placement on every plate with Per's marks")
    ap.add_argument("--out", default="/tmp/claude-1000/-home-per/ce00d501-b635-470a-b1e4-2b0f4fb1ced0/scratchpad/snp_bench.csv")
    a = ap.parse_args()
    key = load_key()
    rows, total = [], collections.Counter()
    for run in sorted(key):
        if a.runs and not any(s in run for s in a.runs):
            continue
        t = time.time()
        docs = load_plate(run)
        calls = ENGINES[a.engine](docs, run)
        conf = collections.Counter()
        for w, k in key[run].items():
            r = calls.get(w) or {}
            got = NORM.get(r.get("call"), r.get("call") or "missing")
            conf[(k["call"], got)] += 1
            rows.append({"run": run, "well": w, "hand": k["call"], "engine": got,
                         "ai": k["ai"], "flag": k["flag"], "note": k["note"],
                         "frac": r.get("frac"), "flags": r.get("flags"),
                         "std_source": r.get("std_source") or r.get("std_scans"),
                         "reason": r.get("reason"), "homo1": r.get("homo1"),
                         "homo2": r.get("homo2"), "hetero": r.get("hetero")})
        n = sum(conf.values())
        ok = sum(v for (h, e), v in conf.items() if h == e)
        total["ok"] += ok
        total["n"] += n
        wrong = {f"{h}->{e}": v for (h, e), v in sorted(conf.items()) if h != e}
        print(f"{run:32s} {ok:3d}/{n:3d}  {time.time()-t:5.0f}s  {wrong}", flush=True)
    print(f"TOTAL {total['ok']}/{total['n']} = {100*total['ok']/max(1,total['n']):.1f}%")
    if a.placement:
        score_placement(a.engine, a.runs)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
