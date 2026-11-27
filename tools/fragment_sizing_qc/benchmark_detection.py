"""Benchmark ladder recognition on synthetic plates of every bundled kit.

Each plate is a run of WELLS wells sharing one kit, one ladder channel and one
instrument (4-channel MegaBACE-like or 5-channel LIZ-capable).  Wells differ:
capillary speed, band heights, a dropped band, fading long bands, dye bleed,
pull-up from the sample, busy or sparse sample channels.  Plates without any
ladder check that nothing is claimed.

    python tools/fragment_sizing_qc/benchmark_detection.py [plates_per_kit] [kit,kit,...]

Reports, per kit: kit right (or an identical-length twin), channel right,
"kit uncertain", a nested kit that gives the same sizes, wrong answers, and seconds per plate.
"""
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import fragment_sizing as fs  # noqa: E402

N = 9000
WELLS = 12


class Doc:
    def __init__(self, channels, well):
        self.channels = channels
        self.acgt = channels[:, :4]
        self.well = well
        self.path = f"{well}.rsd"
        self.base_order = "ACGT"


def peaks_at(scans, heights, width, rng, noise=2.0, base=20.0):
    x = np.arange(N)
    y = base + rng.normal(0, noise, N)
    for s, h in zip(scans, heights):
        if 0 < s < N:
            y += h * np.exp(-0.5 * ((x - s) / width) ** 2)
    return y


def plate(kit, rng, *, five=None, ladder=True):
    lad = fs.load_ladder(kit) if kit else None
    is_liz = bool(lad and lad.dye.upper() == "LIZ")
    five = is_liz if five is None else five
    nch = 5 if five else 4
    lch = 5 if is_liz else int(rng.choice([1, 2, 3, 4]))
    # run-level migration: scans for 1 bp, curvature, span of the window
    top = lad.lengths[-1] if lad else 600
    s_lo = rng.uniform(900, 1800)
    per_bp = rng.uniform(4.0, 11.0) * (600.0 / max(top, 300)) ** 0.5
    curv = rng.uniform(0.0, 0.35)
    run_end = rng.uniform(0.9, 1.3)       # <1: the run stops before the last bands

    def mig(bp, stretch):
        u = np.asarray(bp, float) / top
        return s_lo + stretch * per_bp * top * (u - curv * u * u / 2) / (1 - curv / 2)

    sample_kind = rng.choice(["sparse", "alleles", "busy"])
    sample_bp = {"sparse": lambda: rng.uniform(60, top, 3),
                 "alleles": lambda: np.repeat(rng.uniform(80, top * 0.8, 4), 2) + np.tile([0, 4], 4),
                 "busy": lambda: rng.uniform(40, top, 45)}[sample_kind]
    docs = []
    for w in range(WELLS):
        stretch = rng.uniform(0.9, 1.12)
        width = rng.uniform(3.0, 7.0)
        chans = np.zeros((N, nch))
        for c in range(nch):
            chans[:, c] = peaks_at([], [], width, rng)
        if lad is not None and ladder:
            bands = np.array(lad.lengths, float)
            h = rng.uniform(600, 2500) * rng.uniform(0.6, 1.4, len(bands))
            h *= np.linspace(1.0, rng.uniform(0.3, 1.0), len(bands))   # long bands fade
            if rng.random() < 0.4:
                h[rng.integers(len(bands))] = 0.0                      # a band missing
            sc = mig(bands, stretch)
            sc = np.where(sc < N * run_end, sc, -1)
            chans[:, lch - 1] = peaks_at(sc, h, width, rng)
            nb = lch % nch                                             # bleed next door
            chans[:, nb] += rng.uniform(0.04, 0.15) * (chans[:, lch - 1] - 20)
        # sample peaks in the other channels, pull-up into the ladder channel
        others = [c for c in range(1, nch + 1) if c != lch]
        for c in others[:2]:
            bp = sample_bp()
            hs = rng.uniform(200, 4000, len(bp))
            y = peaks_at(mig(bp, stretch), hs, width, rng, noise=0)
            chans[:, c - 1] += y - 20
            chans[:, lch - 1] += rng.uniform(0.0, 0.06) * (y - 20)
        if not ladder or lad is None:
            pass
        docs.append(Doc(chans, f"W{w:02d}"))
    return docs, lch


def twins():
    """Kits with identical length tables (cannot be told apart by peaks)."""
    by = defaultdict(set)
    for k in fs.list_ladders():
        by[tuple(fs.load_ladder(k).lengths)].add(k)
    return {k: grp for grp in by.values() for k in grp}


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    kits = sys.argv[2].split(",") if len(sys.argv) > 2 else [
        k for k in fs.BUILTIN_LADDERS]
    all_kits = tuple(fs.BUILTIN_LADDERS)
    tw = twins()
    rng = np.random.default_rng(12345)
    tot = defaultdict(lambda: defaultdict(int))
    for kit in kits + ["(none)"]:
        for i in range(per):
            docs, lch = plate(None if kit == "(none)" else kit, rng,
                              five=None if kit != "(none)" else bool(i % 2))
            t = time.time()
            det = fs.detect_ladder(iter(docs), kits=all_kits, enough=6, min_docs=8)
            dt = time.time() - t
            r = tot[kit]
            r["plates"] += 1
            r["sec"] += dt
            got = det["ladder"]
            if kit == "(none)":
                r["ok" if got is None else "wrong"] += 1
                if got:
                    print(f"  (none) -> {got} ch{det['channel']}")
                continue
            kit_ok = got is not None and got in tw.get(kit, {kit})
            ch_ok = det["channel"] == lch
            if not kit_ok and ch_ok and kit in det.get("same_sizes", []):
                r["same_sizes"] += 1      # a nested kit with the same lengths over the bands seen
                print(f"  {kit} ch{lch} -> {got} (same lengths over the bands seen)")
            elif kit_ok and ch_ok:
                r["ok"] += 1
                if det.get("rivals"):
                    r["uncertain"] += 1
                    print(f"  {kit}: uncertain, rivals {det['rivals']} same-sizes {det.get('same_sizes')}")
            else:
                r["wrong"] += 1
                print(f"  {kit} ch{lch} -> {got} ch{det['channel']} "
                      f"rivals={det.get('rivals')} wells_ok={det['wells_ok']}/{det['wells']}")
        r = tot[kit]
        print(f"{kit:18s} ok {r['ok']}/{r['plates']}  uncertain {r['uncertain']}  "
              f"same-sizes {r['same_sizes']}  wrong {r['wrong']}  {r['sec'] / r['plates']:.1f} s/plate", flush=True)


if __name__ == "__main__":
    main()
