#!/usr/bin/env python3
"""Benchmark the CTCE species engine against hand-dragged areas.

The drags come from the GUI's Area mode (Save areas -> CSV or JSON): one row
per drag, with run, well, start_scan, end_scan, peak_scan, channel and
area_Vscan.  Each drag is matched to the engine block it overlaps most on the
same run and well, and both are printed side by side:

    python3 ctce_benchmark.py drags.csv /path/to/plate/Run01 Run02 Run03

The whole plate is analysed (the species pattern is learnt over every well),
then only the dragged wells are reported.  Drags on a run the engine did not
measure (another temperature) are listed but cannot be compared.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import ctce


def load_drags(path: Path) -> list[dict]:
    if path.suffix.lower() == ".json":
        rows = json.loads(path.read_text(encoding="utf-8"))
    else:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        try:
            out.append(dict(run=str(r["run"]), well=str(r["well"]),
                            start=int(float(r["start_scan"])), stop=int(float(r["end_scan"])),
                            peak=int(float(r.get("peak_scan") or 0)),
                            channel=int(float(r.get("channel") or 0)),
                            area=float(r["area_Vscan"])))
        except (KeyError, ValueError) as e:
            print(f"skipping row {r}: {e}", file=sys.stderr)
    return out


def engine_blocks(res: ctce.WellResult) -> dict:
    """run name -> (channel, quant trace, [(name, lo, hi, a, b)])"""
    out = {}
    for run in res.runs:
        ov = getattr(run, "_species_overlay", None)
        if ov:
            out[run.run] = (ov["channel"], ov["yq"], ov["blocks"])
    return out


def overlap(a0, a1, b0, b1) -> int:
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def minor_fraction(by_name: dict) -> float | None:
    if "Hom2" not in by_name or "Hom3" not in by_name:
        return None
    return ctce.species_minor_fraction([by_name.get(n, 0.0) for n in ctce.SPECIES_NAMES
                                        if n in by_name])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("drags", type=Path, help="areas saved from the GUI (CSV or JSON)")
    ap.add_argument("runs", nargs="+", type=Path, help="run folders of one plate")
    ap.add_argument("--channel", type=int, default=None, help="sample channel (default: auto)")
    args = ap.parse_args(argv)

    drags = load_drags(args.drags)
    if not drags:
        print("no drags found", file=sys.stderr)
        return 1
    results = {r.well: r for r in ctce.analyse_plate(args.runs, channel=args.channel)}

    by_well = defaultdict(list)
    for d in drags:
        by_well[(d["run"], d["well"])].append(d)

    ratios, bound_err, mf_pairs, lines = [], [], [], []
    hdr = (f"{'run':<44} {'well':<4} {'drag':>11} {'drag area':>10} | {'engine':<5} "
           f"{'block':>11} {'eng area':>10} {'ratio':>6} {'d_lo':>5} {'d_hi':>5} {'same span':>10}")
    print(hdr)
    print("-" * len(hdr))
    for (run, well), ds in sorted(by_well.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        res = results.get(well)
        blocks = engine_blocks(res).get(run) if res else None
        drag_named, eng_named = {}, {}
        for d in sorted(ds, key=lambda d: d["start"]):
            pre = f"{run[-44:]:<44} {well:<4} {d['start']:>5}-{d['stop']:<5} {d['area']:>10.0f} | "
            if not blocks:
                why = "well not analysed" if res is None else "engine did not measure this run"
                print(pre + why)
                continue
            ch, yq, bl = blocks
            if d["channel"] and d["channel"] != ch:
                pre += f"(drag ch{d['channel']}, engine ch{ch}) "
            best = max(bl, key=lambda b: overlap(d["start"], d["stop"], b[1], b[2]))
            if overlap(d["start"], d["stop"], best[1], best[2]) == 0:
                print(pre + "no engine block here (missed species or extra drag)")
                continue
            name, lo, hi, a, b = best
            ea = ctce.block_area(yq, lo, hi, a, b)
            same = ctce.block_area(yq, d["start"], d["stop"])   # engine trace, user's span
            r = ea / d["area"] if d["area"] else float("nan")
            ratios.append((name, r))
            bound_err.append((lo - d["start"], hi - d["stop"]))
            drag_named[name], eng_named[name] = d["area"], ea
            print(pre + f"{name:<5} {lo:>5}-{hi:<5} {ea:>10.0f} {r:>6.2f} {lo - d['start']:>+5} "
                        f"{hi - d['stop']:>+5} {same:>10.0f}")
        if blocks:
            dragged = {x for x in drag_named}
            for name, lo, hi, a, b in blocks[2]:
                if name not in dragged:
                    print(f"{'':<44} {well:<4} {'(no drag)':>11} {'':>10} | {name:<5} "
                          f"{lo:>5}-{hi:<5} {ctce.block_area(yq, lo, hi, a, b):>10.0f}")
            mdr, men = minor_fraction(drag_named), minor_fraction(eng_named)
            h4, h5 = drag_named.get("Het4"), drag_named.get("Het5")
            note = f"   minor fraction: drags {mdr:.4f}  engine {men:.4f}" if mdr and men else ""
            if h4 and h5:
                note += f"   drag Het5/Het4 {h5 / h4:.2f}"
            if note:
                print(note)
            if mdr and men:
                mf_pairs.append((mdr, men))
        print()

    if ratios:
        print("Summary (engine area / drag area):")
        for name in ctce.SPECIES_NAMES:
            rs = [r for n, r in ratios if n == name]
            if rs:
                print(f"  {name}: n={len(rs):<3} median {statistics.median(rs):.3f}  "
                      f"min {min(rs):.3f}  max {max(rs):.3f}")
        lo_e = [abs(a) for a, _ in bound_err]
        hi_e = [abs(b) for _, b in bound_err]
        print(f"  block ends vs drag: median |start| {statistics.median(lo_e):.0f} scans, "
              f"median |stop| {statistics.median(hi_e):.0f} scans")
        if mf_pairs:
            d = [e - m for m, e in mf_pairs]
            print(f"  minor fraction, engine - drags: median {statistics.median(d):+.4f} "
                  f"over {len(d)} wells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
