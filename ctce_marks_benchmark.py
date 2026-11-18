#!/usr/bin/env python3
"""Score the CTCE species engine against hand-dragged areas and replicate plates.

    python3 ctce_marks_benchmark.py --data "/path/to/CTCE data" \
        --marks docs/ctce_marks/More_measure_#19_20.csv

1. Every marked well (the GUI's "Save areas" CSV: one row per drag, sorted
   along the trace as Homo1, Homo2, Hetero1, Hetero2 -- fewer rows are read
   as Homo1, Homo2 and one merged heteroduplex area) is compared run by run
   with the engine: minor fraction p = (Homo2 + heteroduplexes/2) / total.
2. Plates where every well holds the same DNA (021123: 96 replicates;
   010823 / 040723: a dilution curve across the columns) should give one
   fraction everywhere: their median, robust CV and column medians.

Run folders are looked up under --data by the run names in the marks file;
plates are given with --plate (repeatable), as the run-folder prefix before
"RunNN".
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import os
import sys
import warnings

import numpy as np

import ctce


def _num(s):
    return float(str(s).replace(",", "."))


def _p(a):
    return (a[1] + sum(a[2:]) / 2) / sum(a) if len(a) >= 2 and sum(a) > 0 else float("nan")


def _rcv(v):
    v = np.asarray(v, float)
    return 1.4826 * np.median(np.abs(v - np.median(v))) / np.median(v) if len(v) else float("nan")


def main(argv=None) -> int:
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", required=True, help="folder holding the run folders")
    ap.add_argument("--marks", help="the drags CSV (GUI: Peak area by drag, Save areas)")
    ap.add_argument("--plate", action="append", default=[],
                    help="run-folder prefix of a same-DNA plate, e.g. 021123_mt_nucl_#19_20")
    args = ap.parse_args(argv)
    base = args.data.rstrip("/") + "/"

    if args.marks:
        marks = collections.defaultdict(list)
        with open(args.marks, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                marks[(r["run"], r["well"])].append((int(r["peak_scan"]), _num(r["area_Vscan"])))
        print("Marked wells (engine per run vs the drags):")
        cache, diffs = {}, []
        for run, well in sorted(marks):
            if not os.path.isdir(base + run):
                print(f"  {run} not under --data; skipped")
                continue
            if run not in cache:
                cache[run] = {x.well: x for x in ctce.analyse_plate([base + run])}
            res = cache[run].get(well)
            ua = [a for _s, a in sorted(marks[(run, well)])]
            u = _p(ua)
            e = None if res is None else res.minor_fraction
            if e is not None and np.isfinite(u):
                diffs.append(e - u)
            print(f"  {run[-30:]:<30} {well} drags {u:.4f}{'' if len(ua) == 4 else f' ({len(ua)} marks)'}"
                  f"  engine {'below detection' if e is None else f'{e:.4f}'}"
                  f"  {[] if res is None else res.species}")
        if diffs:
            print(f"  engine - drags: mean {np.mean(diffs):+.4f}  mean |diff| {np.mean(np.abs(diffs)):.4f}"
                  f"  ({len(diffs)} runs)")

    for plate in args.plate:
        runs = sorted(glob.glob(base + plate + "Run0*"))
        if not runs:
            print(f"== {plate}: no runs under --data")
            continue
        res = ctce.analyse_plate(runs)
        mf = {r.well: r.minor_fraction for r in res if r.minor_fraction is not None}
        v = np.array(list(mf.values()))
        print(f"== {plate}: {len(v)}/{len(res)} wells with a fraction "
              f"({len(res) - len(v)} below detection)  median {np.median(v):.4f}  robust CV {_rcv(v):.2f}")
        cols = []
        for col in range(1, 13):
            c = [mf[f"{r}{col:02d}"] for r in "ABCDEFGH" if f"{r}{col:02d}" in mf]
            cols.append(f"{col}:{np.median(c):.3f}" if c else f"{col}:-")
        print("   column medians", " ".join(cols))
    return 0


if __name__ == "__main__":
    sys.exit(main())
