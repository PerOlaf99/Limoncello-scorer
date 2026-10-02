"""Manual genotyping for capillary-electrophoresis reads (built from scratch).

Click on a peak (or just next to it) and the best available peak-recognition
method locates the peak and shades its area.  The recognition algorithm is
selectable; a polymerase-A-addition (+A) peak trailing the main peak by about
one repeat can be tagged automatically.  Nothing is tagged in front of the
main peak: for a single-base-extension product that leading shoulder is
another A-addition on the GC-clamp side, not stutter.

The four internal-standard peaks are captured through the CTC-CE duplex
pattern: peaks 1-2 are the two homoduplexes (they differ by the single SNP
base of the rs number) and peaks 3-4 are the two heteroduplexes made in the
PCR by Watson/Crick re-annealing (one mismatch base).  All four are the SAME
fragment (same bp), so a shared ``Fragment length (bp)`` applies to all four.

Peak position, channel/base, height, area and kind are collected in a table
that can be saved as CSV, Excel (.xlsx) or JSON — readable in Excel and usable
as ML training input.  For positions showing two peaks, the mutant (variant)
fraction is computed as small / (small + large).

The recognition logic lives in ``PeakPicker`` (headless, reusable from the main
window for batch picking); ``GenotypingEditor`` is the Tk widget that edits one
trace and ``GenotypingDialog`` the standalone Toplevel wrapper kept for
scripts/tests.
"""
import csv
import json
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from analyzer_core import acgt_index_for_channel, load_trace  # noqa: F401  (re-export)

# numpy dropped ``trapz`` in 2.0 (renamed ``trapezoid``); accept both so peak
# areas keep working on the whole supported numpy range (>=1.20).
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

CHANNEL_ORDER = "ACGT"
CLICK_RADIUS = 30          # scans searched around a click
HET_WINDOW = 8             # scans in which two mains count as heterozygote
SATELLITE_FRAC = 0.05      # min height of a +A satellite vs the main peak
# A +A product and a genuine minor allele one base downstream are close to
# indistinguishable on position and height alone: _shoulders() must pick one,
# and it currently always picks "+A", which drops the peak from clusters() and
# makes clust_frac() report 0.0.  A heterozygote then reads as a homozygote
# with no indication that anything was ambiguous -- the most damaging kind of
# wrong.
#
# The honest fix is not a better threshold (that needs a plate of known
# minor-allele hets to calibrate against, and a synthetic attempt regressed a
# genuine +A case) but to stop hiding the ambiguity.  Above this height
# fraction the tag is recorded as "could be a minor allele" so the caller can
# raise an uncertain instead of reporting a confident 0.0.
#
# 0.15 is NOT calibrated against known minor alleles -- it has never been
# tested against a plate of them, because none was available.  It is set from
# two facts: a genuine 50% minor allele one spacing from its parent measures
# only ~0.22 of that parent here, because the two peaks overlap badly and
# _sat_height reads the local valley, not the apex; and typical +A tails run
# 0.05-0.15.  So 0.15 sits above ordinary +A while still catching a minor allele
# of roughly 30% or more.  Expect to revise it once real minor-allele samples
# are available, and treat any flag it raises as "look at this", not as a call.
SATELLITE_AMBIGUOUS_FRAC = 0.15
MIN_PEAK_FRAC = 0.05       # peak candidates must stand off the segment floor
VALLEY_CAP = 1.5           # hard bound on how far a peak's area may reach, in x
                           # spacings, and the same reach as the heterozygote
                           # window: a peak can never measure into where a
                           # neighbouring allele could start.  The peak's own
                           # valley normally ends the search well before this
VALLEY_FLOOR = 0.2          # a valley counts only once the trace is down to
                           # this much of the peak's height above the baseline
                           # between peaks, so a step up the shoulder of a
                           # sharp peak is not mistaken for the valley beside it
VALLEY_FLAT = 0.01          # ...and a trace that has changed by less than this
                           # over a few scans has settled at the valley floor

PEAK_FINDERS = [
    ("best", "Best (prominence + area)"),
    ("max", "Simple local maxima"),
    ("gauss", "Gaussian fit"),
]


class _Record(dict):
    """One picked peak. Writable dict so rows feed straight into export."""


# --------------------------------------------------------------------------- #
# rs1695 internal-standard detection
# --------------------------------------------------------------------------- #
# On a CTC-CE run the internal standard is always present as four equimolar
# peaks, so it can be found without the click-pick dance that ``mark_std``
# otherwise needs.  The geometry is fixed by the chemistry: the two homoduplexes
# are one SNP base apart, the heteroduplexes are one mismatch base from
# Watson/Crick re-annealing, so d1 ~ d3 and d2 is roughly 2.5x either.
#
# The equimolar prior matters as much as the spacing.  Requiring the weakest of
# the four to be at least 40% of the strongest is what keeps an arbitrary run of
# four peaks from being read as a standard.
IS_MIN_SPACING = 25.0        # scans; d1 and d3 (one SNP base)
IS_OUTER_TOL = 0.75          # |d1 - d3| <= this * dm
IS_MID_LO, IS_MID_HI = 1.7, 3.3   # d2, in units of dm
IS_EQUIMOLAR_MIN = 0.40      # min(height) / max(height) across the four
IS_PROMINENCE_FRAC = 0.010
IS_HEIGHT_FRAC = 0.020
IS_SAT_SKIP = (10, 32)       # a +A tail sits this many scans behind its parent
IS_SAT_RATIO = 1.15
IS_MAX_CANDIDATES = 60       # keep the O(n^4) search bounded


# ``acgt_index_for_channel`` is re-exported from analyzer_core (imported at the
# top of this module) so tests and the GUI can keep calling
# ``genotyping.acgt_index_for_channel`` while the headless scorer uses the same
# implementation without pulling in Tk/matplotlib.


def _is_candidates(trace, cut=1900):
    """Despiked, baseline-corrected peak candidates of the standard channel."""
    from scipy.signal import find_peaks, medfilt, savgol_filter
    y = savgol_filter(medfilt(np.asarray(trace, float), 5), 9, 2, mode="interp")
    y = y - float(np.median(y[:400]))            # baseline before prominence
    mx = float(y.max())
    if mx <= 0:
        return []
    p, props = find_peaks(y, prominence=mx * IS_PROMINENCE_FRAC,
                          height=mx * IS_HEIGHT_FRAC)
    base = float(np.median(y[:400]))
    cands = [(int(i), float(y[i] - base),
              float(props["prominences"][k]))
             for k, i in enumerate(p) if i >= cut]
    # A +A tail is always shorter than and close behind its parent peak; drop it
    # so it cannot stand in for a real member of the quartet.
    kept = []
    for x, h, prom in sorted(cands, key=lambda c: -c[1]):
        if any(IS_SAT_SKIP[0] <= x - kx <= IS_SAT_SKIP[1] and kh > IS_SAT_RATIO * h
               for kx, kh, _ in kept):
            continue
        kept.append((x, h, prom))
    return sorted(kept, key=lambda c: c[0])[:IS_MAX_CANDIDATES]


def find_is_quartet(trace, cut=1900):
    """Locate the four internal-standard peaks, or ``None``.

    *trace* is one channel's samples.  Returns ``(peaks, heights)`` where
    *peaks* are the four scan positions in ascending order, chosen as the
    highest-prominence candidate quartet that satisfies the geometry.  Found in
    96/96 wells of the T9 rs1695 plate.
    """
    from itertools import combinations
    pk = _is_candidates(trace, cut)
    if len(pk) < 4:
        return None
    best = None
    for combo in combinations(range(len(pk)), 4):
        a, b, c, d = (pk[i] for i in combo)
        d1, d2, d3 = b[0] - a[0], c[0] - b[0], d[0] - c[0]
        if d1 < IS_MIN_SPACING or d3 < IS_MIN_SPACING:
            continue
        dm = (d1 + d3) / 2.0
        if abs(d1 - d3) > IS_OUTER_TOL * dm or not (IS_MID_LO * dm <= d2 <= IS_MID_HI * dm):
            continue
        heights = [a[1], b[1], c[1], d[1]]
        lo, hi = min(heights), max(heights)
        if lo <= 0 or lo / hi < IS_EQUIMOLAR_MIN:
            continue
        score = sum((a[2], b[2], c[2], d[2]))
        if best is None or score > best[0]:
            best = (score, [a[0], b[0], c[0], d[0]], heights)
    if best is None:
        return None
    return best[1], best[2]


# --------------------------------------------------------------------------- #
# batch auto-genotyping: one well -> one call, plus the reason when it cannot
# --------------------------------------------------------------------------- #
# The manual path is ``auto_mark_std`` for the standard and four hand-measured
# duplexes for the sample.  That is fine for one well and unusable for a plate,
# so this measures the same thing with no clicks: find the standard quartet,
# read the sample channel in the four windows it defines, and call the well.
#
# Which channel is which is a property of the KIT, not something a trace can be
# asked to work out.  Measured over the 96 rs1695 T9 wells with this module's
# own ``find_is_quartet``, per ``doc.acgt`` column (which is always A,C,G,T):
#
#     column 0 (A)  quartet in 62/96 wells, position spread 191 scans
#     column 1 (C)  quartet in 59/96 wells, position spread 190 scans
#     column 2 (G)  quartet in 96/96 wells, position spread 125 scans
#     column 3 (T)  quartet in 96/96 wells, position spread 125 scans
#
# Columns 2 and 3 carry the same standard -- a clean equimolar quartet present
# in every well at a consistent position is the signature, and it bleeds a
# little into G.  Column 3 is the dye itself, so that is the one to read.  The
# sample is column 1: reading it reproduces ``rs1695_measured.csv`` to within a
# few percent and its sample-to-standard area ratio is constant across the
# plate, which is what a real amplicon looks like and the standard never does.
#
# On this plate's "ACTG" dye order (Ch1=A, Ch2=C, Ch3=T, Ch4=G) that makes the
# standard Ch3 and the sample Ch2, i.e. ``auto_mark_std``'s ``channel=3``
# default.  Both constants stay explicit because getting this pair wrong does
# not fail loudly -- it scores the standard's own peaks as the sample and
# returns confident nonsense -- and every result row records the channels it
# was measured on.
#
# Beware when comparing against the analysis cache: ``t9raw.npz`` stores wells
# in the plate's physical channel order, not A,C,G,T, so a column index means
# something different there.  On this plate npz->acgt is [2, 3, 1, 0].
DEFAULT_IS_CHANNEL = 3
DEFAULT_SAMPLE_CHANNEL = 2
DEFAULT_IS_CUT = 1900        # scans; skip the injection front before looking

# A duplex is measured between the midpoints to its neighbours -- the standard's
# own spacing already says where one fragment stops and the next begins, which
# is tighter than any per-peak valley search and cannot wander into a
# neighbour's area.
NOISE_SGOLAY_WINDOW = 9      # noise = scatter left over by a 9-scan SG fit
SEGMENT_BASELINE_PAD = 80    # scans either side of a duplex to set its baseline
SEGMENT_APEX_RADIUS = 13     # scans either side of the standard peak to look for
                            # the sample apex: the two migrate close but not
                            # identically, and the residual shift is what the
                            # sample's own height must be read at
MIN_SEGMENT_SPAN = 3         # narrower than this and it is not a peak at all


def _noise_sigma(y):
    """Robust per-scan noise as 1.4826 x the MAD of a Savitzky-Golay residual.

    A plain standard deviation is wrong here: the four duplexes are large
    enough to dominate it, which would make the noise read far too high and
    every peak look insignificant.  The median absolute deviation of what the
    smooth fit fails to explain is not.
    """
    from scipy.signal import savgol_filter
    fit = savgol_filter(y, NOISE_SGOLAY_WINDOW, 2, mode="interp")
    d = y - fit
    return 1.4826 * float(np.median(np.abs(d - np.median(d))))


def _quartet_segments(scans):
    """Integration window per duplex, from the midpoints to its neighbours.

    Extrapolated half a spacing beyond the outer peaks so the first and last
    duplexes are measured over the same width as the inner two.
    """
    n = len(scans)
    if n < 2:
        return [(max(0, scans[0] - 1), min(len(scans), scans[0] + 1))] if n else []
    edge = (scans[1] - scans[0]) / 2.0
    bounds = [scans[0] - edge]
    bounds += [(scans[i] + scans[i + 1]) / 2.0 for i in range(n - 1)]
    bounds += [scans[-1] + edge]
    return [(int(round(bounds[i])), int(round(bounds[i + 1])))
            for i in range(n)]


def auto_genotype(doc, is_channel=DEFAULT_IS_CHANNEL,
                  sample_channel=DEFAULT_SAMPLE_CHANNEL,
                  base_order="ACTG", cut=DEFAULT_IS_CUT, run_name=""):
    """Genotype one well without any clicking -> a result row.

    Returns a plain dict (so it feeds straight into ``save_table``) with the
    call, the four duplex areas, their significances, and -- when there is no
    call -- a *reason* naming what was missing.  Never raises for a bad well:
    a plate is 96 chances to hit a bad well and one of them must not take the
    other 95 down with it.
    """
    import scorer

    row = {
        "run": run_name, "well": getattr(doc, "well", "") or "",
        "call": "no-call", "frac": 0.0, "flags": "",
        "hom1": 0.0, "hom2": 0.0, "het1": 0.0, "het2": 0.0,
        "snr1": 0.0, "snr2": 0.0, "snr3": 0.0, "snr4": 0.0,
        "is_channel": is_channel, "sample_channel": sample_channel,
        "std_scans": "", "reason": "",
    }

    try:
        is_col = acgt_index_for_channel(base_order, is_channel)
        samp_col = acgt_index_for_channel(base_order, sample_channel)
    except ValueError as e:
        row["reason"] = str(e)
        return row
    if is_col == samp_col:
        row["reason"] = (f"standard and sample are both on Ch{is_channel} "
                         "(same acgt column)")
        return row

    acgt = np.asarray(getattr(doc, "acgt", None), dtype=float)
    if acgt.ndim != 2 or acgt.shape[0] == 0 or acgt.shape[1] <= max(is_col, samp_col):
        row["reason"] = "trace has no usable channels"
        return row

    found = find_is_quartet(acgt[:, is_col], cut=cut)
    if found is None:
        row["reason"] = f"no standard quartet on Ch{is_channel}"
        return row

    scans, _heights = found
    row["std_scans"] = "/".join(str(x) for x in scans)
    segs = _quartet_segments(scans)
    if len(segs) != 4:
        row["reason"] = "standard quartet too short to measure"
        return row

    y = np.asarray(acgt[:, samp_col], dtype=float)
    n = y.size
    sigma = _noise_sigma(y)
    idx = np.arange(n)
    areas, snrs = [], []
    for k, (lo, hi) in enumerate(segs):
        lo, hi = max(0, lo), min(n, hi)
        if hi - lo < MIN_SEGMENT_SPAN:
            areas.append(0.0)
            snrs.append(0.0)
            continue
        # Baseline from the quiet trace either side of this duplex, not from its
        # own peak: a peak sitting on a raised baseline would otherwise measure
        # the step under it as signal.
        m = (idx >= max(0, lo - SEGMENT_BASELINE_PAD)) & \
            (idx <= min(n, hi + SEGMENT_BASELINE_PAD))
        base = float(np.median(y[m]))
        a = max(0, int(scans[k]) - SEGMENT_APEX_RADIUS)
        b = min(n, int(scans[k]) + SEGMENT_APEX_RADIUS + 1)
        apex = float(y[a:b].max()) if b > a else 0.0
        areas.append(max(0.0, float(_trapz(y[lo:hi] - base, dx=1.0))))
        snrs.append((apex - base) / sigma if sigma > 0 else 0.0)

    for i, name in enumerate(("hom1", "hom2", "het1", "het2")):
        row[name] = round(areas[i], 1)
        row["snr%d" % (i + 1)] = round(snrs[i], 1)

    call, frac, flags = scorer.t9_call(areas[0], areas[1], areas[2], areas[3], snrs)
    row["call"] = call
    row["frac"] = round(frac, 4)
    row["flags"] = ",".join(sorted(flags))
    if call == "no-call":
        if sigma <= 0:
            row["reason"] = f"no signal on Ch{sample_channel}"
        elif max(snrs) < scorer.T9_MIN_DOMINANT_SIGMA:
            row["reason"] = (f"weakest sample duplex is only "
                             f"{max(snrs):.0f}x the noise on Ch{sample_channel}")
        else:
            row["reason"] = "no standard quartet on Ch%d" % is_channel
    return row


DEFAULT_COLORS = {"A": "#00AA00", "C": "#0000DD", "G": "#111111", "T": "#DD0000"}


def channel_colors(base_order, colors=None, theme_mode="base"):
    """Map each ``doc.acgt`` column to its line colour.

    Shared by every trace tool here so a colour can never mean one thing in the
    picker and another in the area measure.  Keys are ACGT column indices; in
    "channel" mode the dye a channel carries is ignored and the colour follows
    the channel's position, which is what the genotyping themes want.
    """
    colors = dict(colors or DEFAULT_COLORS)
    out = {}
    for ci, base in enumerate((base_order or "ACTG").upper()[:4]):
        if base in CHANNEL_ORDER:
            key = CHANNEL_ORDER[ci] if theme_mode == "channel" else base
            out[CHANNEL_ORDER.index(base)] = colors.get(key, "#444444")
    return out


class PeakPicker:
    """Headless click-to-pick peak recognition for one trace.

    Holds the picked records, the CTC-CE internal-standard duplexes, the
    optional shared fragment length and all recognition/detection logic.  The
    UI reads the same attributes (records, std, length_bp, col_color) and
    calls the same methods, so one engine can drive a dialog, an embedded
    widget or the main window directly.
    """

    def __init__(self, doc, path, colors=None, base_order=None,
                 theme_mode="base", show=None, include_sh=True,
                 show_d2=True, finder="best"):
        self.doc = doc
        self.path = Path(path)
        self.show = show or (lambda col: True)
        self.colors = dict(colors or {"A": "#00AA00", "C": "#0000DD",
                                      "G": "#111111", "T": "#DD0000"})
        self.base_order = (base_order or "ACTG").upper()
        self.theme_mode = theme_mode
        self.include_sh = bool(include_sh)
        self.show_d2 = bool(show_d2)
        self.finder = finder or "best"

        self.records: list[_Record] = []
        self._gid = 0
        self._reject = None
        self.std = None
        self.length_bp = None
        self.col_color = channel_colors(self.base_order, self.colors,
                                        theme_mode)
        self._d2_cache: dict[int, tuple] = {}

    # ------------------------------------------------------------- detection
    def pick(self, scan, vol=None):
        """Pick the peak you clicked on: nearest scan wins, then nearest
        voltage (tells channels apart when two share a scan), then tallest.
        A peak whose area is already picked on the same channel is never
        picked again — undo it first if you need to (neighbouring peaks, e.g.
        the two alleles of a heterozygote, stay pickable).  Returns the new
        record (plus any +A record appended) or None; on a refused
        re-pick, ``_reject`` is set to ``"area"`` for the UI message."""
        self._reject = None
        radius = max(CLICK_RADIUS, int(self._spacing() * 2.0))
        cands = []
        for col, color in self.col_color.items():
            if not self.show(col):
                continue
            pk = self._detect(col, scan, radius)
            if pk is None:
                continue
            pk["col"] = col
            pk["color"] = color
            cands.append(pk)
        if not cands:
            self._reject = "none"
            return None
        if vol is None:
            cands.sort(key=lambda p: (abs(p["apex"] - scan), -p["height"]))
        else:
            cands.sort(key=lambda p: (abs(p["apex"] - scan),
                                      abs(p["height"] - float(vol)),
                                      -p["height"]))
        best = cands[0]
        tol = max(1, int(round(self._spacing() * 0.25)))
        for rec in self.records:
            if rec["col"] == best["col"] and abs(rec["scan"] - best["apex"]) <= tol:
                self._reject = "area"
                return None
        self._gid += 1
        rec = _Record(
            file=str(self.path),
            well=self.doc.well,
            scan=best["apex"],
            channel=best["col"] + 1,
            base=CHANNEL_ORDER[best["col"]],
            kind="main",
            height=best["height"],
            area=best["area"],
            left=best["left"],
            right=best["right"],
            onset=best.get("onset", best["left"]),
            end=best.get("end", best["right"]),
            color=best["color"],
            col=best["col"],
            gid=self._gid,
        )
        if self.include_sh:
            for sib in self._shoulders(best, radius):
                # Copy the detector's fields wholesale, then add the picker's
                # own.  Listing them one by one used to drop any field the
                # detector added later, which silently discarded the
                # satellite's ambiguity flag.
                rec2 = _Record(sib)
                rec2.update(
                    file=str(self.path), well=self.doc.well,
                    scan=sib["apex"], channel=sib["col"] + 1,
                    base=CHANNEL_ORDER[sib["col"]],
                    gid=self._gid,
                )
                self.records.append(rec2)
        self.records.append(rec)
        return rec

    def _detect(self, col, scan, radius):
        """Run the chosen peak-recognition algorithm near scan on one column."""
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        n = y.size
        if n < 3:
            return None
        c = int(np.clip(scan, 0, n - 1))
        w = np.arange(max(0, c - radius), min(n, c + radius + 1))
        if w.size < 3:
            return None
        alg = dict(PEAK_FINDERS).get(self.finder, "best")
        if alg == "gauss":
            v = self._numeric_peak(y, w, col)
            if v is None:
                return v
            g = self._gauss_fit(y, v["left"], v["right"], v["apex"], v["height"])
            if g is not None:
                v.update(g)
            return v
        return self._numeric_peak(y, w, col)

    def _d2(self, col):
        """Smoothed second derivative of one channel, plus a robust noise
        estimate (median-absolute-deviation * 1.4826).  Cached per channel."""
        if col in self._d2_cache:
            return self._d2_cache[col]
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        try:
            from scipy.signal import savgol_filter
            yy = savgol_filter(y, window_length=7, polyorder=2, mode="interp")
        except Exception:
            k = np.ones(5) / 5.0
            yy = np.convolve(y, k, mode="same")
            yy = np.convolve(yy, k, mode="same")
        d2 = np.gradient(np.gradient(np.asarray(yy, dtype=float)))
        s2 = 1.4826 * np.median(np.abs(d2 - np.median(d2)))
        self._d2_cache[col] = (d2, float(s2))
        return self._d2_cache[col]

    def _onset_end(self, col, apex, left, right):
        """Start/stop of the peak from the second derivative: the concave-up
        lift-off on the rising flank (onset) and the concave-up return on the
        falling flank (end).  Falls back to the valley boundaries (left/right)
        when no curvature threshold stands out."""
        d2, s2 = self._d2(col)
        th = 3.0 * s2
        apex = int(apex)
        left = max(0, int(left))
        right = min(self.doc.acgt.shape[0] - 1, int(right))
        onset = None
        for i in range(left, max(left, apex - 1)):
            if d2[i] > th and d2[i + 1] > th:
                onset = i
                break
        end = None
        for i in range(right, max(apex + 1, right - 1), -1):
            if d2[i] > th and d2[i - 1] > th:
                end = i
                break
        return (left if onset is None else onset,
                right if end is None else end)

    def _numeric_peak(self, y, w, col):
        seg = y[w]
        base = float(np.nanmin(seg))
        top = float(np.nanmax(seg))
        rng = top - base
        if not np.isfinite(rng) or rng <= 1e-12:
            return None
        m = (seg[1:-1] >= seg[:-2]) & (seg[1:-1] > seg[2:])
        idx = np.flatnonzero(m) + 1
        if idx.size == 0:
            return None
        on = idx[seg[idx] - base >= MIN_PEAK_FRAC * rng]
        if on.size == 0:
            on = idx
        # strongest candidate nearest to the window centre
        c = w.size // 2
        dist = np.abs(w[on] - w[c])
        order = np.lexsort((seg[on], dist))
        best = int(on[order[0]])
        apex = int(w[best])
        # A peak's own two valleys bound its area: the area is integrated out to
        # the valley on each side, so a variant fraction uses the peak's whole
        # hump rather than a fixed slice of it.  VALLEY_CAP is only a bound for
        # the two cases the valley itself cannot settle: a peak with no valley
        # beside it (a low minor allele sitting on a main's tail has no dip of
        # its own) must not measure into where its neighbour could start.
        sp = self._spacing()
        win_base = float(np.mean((y[w[0]], y[w[-1]])))
        left = self._own_valley(y, apex, -1, sp, win_base)
        right = self._own_valley(y, apex, +1, sp, win_base)
        if right <= left:
            right = min(y.size - 1, left + 2)
        xs = np.arange(left, right + 1)
        bl = np.linspace(float(y[left]), float(y[right]), right - left + 1)
        area = float(np.sum(np.clip(y[xs] - bl, 0.0, None)))
        height = float(y[apex] - max(y[left], y[right]))
        if height <= 0 or area <= 0:
            return None
        onset, end = self._onset_end(col, apex, left, right)
        return {"apex": apex, "left": left, "right": right,
                "onset": onset, "end": end,
                "area": area, "height": height}

    def _gauss_fit(self, y, left, right, apex, height):
        from scipy.optimize import curve_fit
        xs = np.arange(left, right + 1)
        if xs.size < 5:
            return None

        def gauss(x, a, mu, s, b):
            return a * np.exp(-(x - mu) ** 2 / (2.0 * s * s)) + b

        p0 = (float(height), float(apex), max(2.0, (right - left) / 3.0),
              float(y[apex] - height))
        try:
            popt, _ = curve_fit(gauss, xs, y[xs], p0=p0, maxfev=4000)
        except Exception:
            return None
        a, mu, s, _b = popt
        if not (np.isfinite(a) and np.isfinite(mu) and np.isfinite(s)) or s <= 0:
            return None
        return {"apex": int(round(mu)),
                "area": float(a * abs(s) * np.sqrt(2.0 * np.pi)),
                "height": float(a)}

    def _own_valley(self, y, apex, step, sp, base):
        """The scan where this peak's own trace turns back up (step -1 = left).

        A peak's area runs between its two valleys, so the peak has to keep every
        tail scan that is really its own -- a variant fraction is only as good as
        the two areas it is built from.  The dip that closes a peak is a broad
        valley, while the ripples riding its shoulder are a scan or two wide, so
        the search walks a slightly smoothed copy of the trace (the area itself
        is always integrated from the raw samples) and stops at the first real
        turn of that trace.  This keeps a peak's whole hump yet still ends at the
        dip between two dense alleles instead of running on into the next peak.

        A turn only counts once the trace is down near the baseline: one step up
        a scan or two from the top of a sharp peak is this same peak's shoulder,
        not a valley beside it.  A trace that has merely levelled off at that
        baseline counts as the valley floor too, which is what stops a peak from
        stretching back to a neighbour 40 scans away.  VALLEY_CAP bounds the
        search for a peak that has no valley of its own at all.
        """
        n = y.size
        w = max(3, int(round(sp / 3.0)) | 1)
        pad = w // 2
        ys = np.convolve(np.pad(y, pad, mode="edge"),
                         np.ones(w) / float(w), mode="valid")
        # start at the top of the hump, not one scan off it: a single-scan
        # spike puts the raw apex next to a higher sample, and a walk that
        # starts there reads that step up as the turn and stops at once
        a0 = max(0, apex - pad)
        b0 = min(n - 1, apex + pad)
        ap = a0 + int(np.argmax(ys[a0:b0 + 1]))
        top = float(ys[ap])
        span = max(top - base, 1e-9)
        floor = base + (1.0 - VALLEY_FLOOR) * span
        k = max(3, pad)
        cap = max(2, int(round(sp * VALLEY_CAP)))
        lo = max(0, ap + step)               # first sample outside the hump
        hi = min(n - 1, ap + step * cap)     # never reach a neighbour
        if (step > 0 and hi <= lo) or (step < 0 and hi >= lo):
            return max(0, min(n - 1, ap + step * 2))
        imin = lo
        vmin = float(ys[lo])
        armed = vmin <= floor
        i = lo + step
        while True:
            v = float(ys[i])
            if not armed and v <= floor:
                armed = True
            elif armed:
                back = i - k * step
                if 0 <= back < n and abs(float(ys[back]) - v) <= VALLEY_FLAT * span:
                    break                      # the trace settled at the floor
                inner = float(ys[i - step])          # toward the hump
                outer = float(ys[i + step]) if 0 <= i + step < n else v
                if v <= inner and v <= outer:        # the trace turns here
                    imin = i
                    j = i + step                     # ride the dip to its floor
                    while 0 <= j < n and float(ys[j]) <= v:
                        imin = j
                        j += step
                    break
            if v < vmin:
                vmin = v
                imin = i
            if i == hi:
                break
            i += step
        a = max(0, imin - pad)
        b = min(n, imin + pad + 1)
        return a + int(np.argmin(y[a:b]))    # snap onto the raw samples

    def _sat_height(self, y, x, r):
        """Height of the candidate satellite at scan x above its own two
        local minima, matching how ``_numeric_peak`` measures a main peak's
        height.  Comparing this against ``main["height"]`` puts both numbers on
        the same scale, unlike the raw window value."""
        left = int(np.argmin(y[max(0, x - r): x + 1])) + max(0, x - r)
        right = int(np.argmin(y[x: min(y.size, x + r + 1)])) + x
        if right <= left:
            return 0.0
        return float(y[x] - max(y[left], y[right]))

    def _shoulders(self, main, radius):
        """Tag the strongest TRAILING satellite: the +A polymerase A-addition
        product, found within ~(0.4–1.6) x one repeat (~one base) after the
        main apex.

        Nothing is tagged in front of the main peak.  The leading shoulder of a
        single-base-extension product is a polymerase A-addition too (it sits
        on the GC-clamp side and barely shifts migration), so it is a second
        +A view of the same fragment rather than stutter -- leaving it
        unmarked beats mislabelling it."""
        sp = self._spacing()
        y = np.asarray(self.doc.acgt[:, main["col"]], dtype=float)
        n = y.size
        out = []
        a, b = (main["apex"] + sp * 0.4, main["apex"] + sp * 1.6)
        lo, hi = (int(min(a, b)), int(max(a, b)))
        r = max(2, int(round(sp * 0.20)))
        cands = self._window_peaks(main["col"], lo, hi)
        # Compare like with like: a satellite's height has to be measured above
        # its own local baseline, the same way ``main["height"]`` is.  The raw
        # window value is an absolute voltage, so on a trace whose baseline sits
        # well above zero it reads far taller than it is -- the old filter then
        # rejected almost every real satellite (their ratio came out > 1) while
        # letting a genuine minor allele through as if it were a +A tail.
        cands = [(x, self._sat_height(y, x, r)) for x, _h in cands]
        cands = [(x, h) for x, h in cands
                 if SATELLITE_FRAC * main["height"] <= h
                 <= 0.9 * main["height"]]
        if not cands:
            return out
        x, h = max(cands, key=lambda c: c[1])
        if x <= main["apex"] + 2:                  # must trail the main
            return out
        left = int(np.argmin(y[max(0, x - r): x + 1])) + max(0, x - r)
        right = int(np.argmin(y[x: min(n, x + r + 1)])) + x
        if right <= left:
            right = min(n - 1, left + 2)
        xs = np.arange(left, right + 1)
        bl = np.linspace(y[left], y[right], right - left + 1)
        area = float(np.sum(np.clip(y[xs] - bl, 0.0, None)))
        hgt = float(y[x] - max(y[left], y[right]))
        if hgt <= 0 or area <= 0:
            return out
        onset, end = self._onset_end(main["col"], x, left, right)
        ambiguous = (main["height"] > 0
                     and hgt >= SATELLITE_AMBIGUOUS_FRAC * main["height"])
        out.append({"apex": x, "left": left, "right": right,
                    "onset": onset, "end": end,
                    "area": area, "height": hgt,
                    "col": main["col"], "color": main["color"],
                    "kind": "+A",
                    "ambiguous": ambiguous,
                    # Record which main this shadows.  The scorer needs it to
                    # attach the ambiguity flag to the right call, and the
                    # nearest-main search is no use here: a satellite sits up
                    # to ~1.6 spacings out, which is outside the allelic
                    # window, so it would find nothing and drop the flag.
                    "parent_scan": int(main["apex"]),
                    "height_frac": round(hgt / main["height"], 4)
                    if main["height"] > 0 else None})
        return out

    def _window_peaks(self, col, lo, hi):
        """All local maxima in a scan window of one channel, as (scan, value)."""
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        lo = max(0, int(lo))
        hi = min(y.size - 1, int(hi))
        if hi - lo < 3:
            return []
        w = np.arange(lo, hi + 1)
        seg = y[w]
        rng = float(seg.max()) - float(seg.min())
        if not np.isfinite(rng) or rng <= 1e-12:
            return []
        m = (seg[1:-1] >= seg[:-2]) & (seg[1:-1] > seg[2:])
        idx = [int(i) + 1 for i in np.flatnonzero(m)
               if seg[int(i) + 1] - seg.min() >= MIN_PEAK_FRAC * rng]
        return [(int(w[i]), float(seg[i])) for i in idx]

    def _spacing(self):
        p = np.asarray(getattr(self.doc, "peak_positions", []) or [], dtype=float)
        if p.size > 1:
            d = float(np.median(np.diff(p)))
            if np.isfinite(d) and 2.0 < d < 80.0:
                return d
        return float(max(6.0, self.doc.n_scans * 0.004))

    # ------------------------------------- internal-standard CTC-CE duplexes
    def auto_mark_std(self, channel=3, cut=1900):
        """Find the internal standard automatically and mark it as the duplex set.

        The CTC-CE standard is always present, so on a real run there is no need
        to click the four peaks by hand.  *channel* is the physical MegaBACE
        channel (3 = the T channel on a standard "ACTG" plate); it is translated
        to an ``acgt`` column here, because the two orderings differ and
        confusing them silently swaps sample and standard.

        Returns the same message :meth:`mark_std` returns, or raises ValueError
        when no quartet is present (in which case the caller should fall back to
        clicking).
        """
        col = acgt_index_for_channel(getattr(self.doc, "base_order", "ACTG"),
                                     channel)
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        found = find_is_quartet(y, cut=cut)
        if found is None:
            raise ValueError(
                f"No internal-standard quartet on channel {channel}; "
                "pick the standard peaks by hand.")
        scans, _heights = found
        self.std = [(x, n) for x, n in zip(scans, ["HOM1", "HOM2", "HET1", "HET2"])]
        self.length_bp = None
        return "Standard set: " + ", ".join(f"{n}@{x}" for x, n in self.std) \
            + "  (auto)"

    def mark_std(self, length_bp=None):
        """Tag the standard mains (scan order) as the four CTC-CE duplexes:
        HOM1, HOM2 (homoduplexes, one rs SNP base apart), HET1, HET2
        (heteroduplexes, one mismatch base from Watson/Crick re-annealing).
        All four are the same fragment, so a shared length (bp) is optional.
        The standard channel is the one of the last main you picked.

        Raises ValueError with a plain-language reason, or returns a message."""
        mains = [r for r in self.records if r["kind"] == "main"]
        if not mains:
            raise ValueError("Pick the standard main peaks first.")
        last_col = mains[-1]["col"]
        col_mains = sorted((m for m in mains if m["col"] == last_col),
                           key=lambda m: m["scan"])
        if len(col_mains) < 2:
            raise ValueError("Need ≥ 2 standard main peaks on the same channel.")
        names = ["HOM1", "HOM2", "HET1", "HET2"][:min(len(col_mains), 4)]
        self.std = [(m["scan"], n) for m, n in zip(col_mains[:4], names)]
        self.length_bp = length_bp
        msg = "Standard set: " + ", ".join(f"{n}@{x}" for x, n in self.std)
        if self.length_bp is not None:
            msg += f"  (len {self.length_bp:g} bp)"
        return msg

    def clear_std(self):
        self.std = None
        self.length_bp = None

    def duplex_of(self, rec):
        """Duplex label (HOM1/HOM2/HET1/HET2) for a main peak that is one of
        the marked standard peaks; '' otherwise."""
        if not self.std or rec["kind"] != "main":
            return ""
        for x, n in self.std:
            if abs(rec["scan"] - x) <= 2:
                return n
        return ""

    def _duplex_of(self, rec):
        """Back-compat alias for duplex_of (used by the editor table/tests)."""
        return self.duplex_of(rec)

    def _clust_frac(self, rec):
        """Back-compat alias for clust_frac (used by the editor table/tests)."""
        return self.clust_frac(rec)

    def _het_window(self):
        """How many scans count as one allelic position: a repeat (~ one base,
        from the channel spacing) with a little slack, never below HET_WINDOW.
        On a CTC-CE run the two alleles of a heterozygote sit one base apart,
        which is typically far more than the old fixed 8-scan window."""
        return max(HET_WINDOW, int(round(self._spacing() * 1.6)))

    def clusters(self):
        """Main peaks grouped into allelic positions: same channel and within
        one repeat (≈ one base) of another main in the group."""
        mains = [r for r in self.records if r["kind"] == "main"]
        win = self._het_window()
        out = []
        for r in sorted(mains, key=lambda m: m["scan"]):
            placed = False
            for cl in out:
                if any(r["col"] == m["col"]
                       and abs(r["scan"] - m["scan"]) <= win for m in cl):
                    cl.append(r)
                    placed = True
                    break
            if not placed:
                out.append([r])
        return out

    def clust_frac(self, rec):
        """Variant (mutant) fraction small/(small+large) for a main peak in a
        position with two peaks; 0.0 otherwise."""
        for cl in self.clusters():
            if rec in cl and len(cl) >= 2:
                areas = sorted(m["area"] for m in cl)
                if sum(areas) > 0:
                    return areas[0] / sum(areas)
        return 0.0

    # ---------------------------------------------------------------- mutate
    def undo_last(self):
        """Remove the most recently picked peak (plus its +A tag)."""
        if not self.records:
            return False
        gid = max(r["gid"] for r in self.records)
        self.records = [r for r in self.records if r["gid"] != gid]
        return True

    def clear_all(self):
        self.records = []
        self._d2_cache = {}

    # ---------------------------------------------------------------- export
    def export_rows(self):
        """One writable dict per picked peak, ready for CSV/JSON/XLSX.

        Rows come out in scan order, not the order they were clicked, so one
        sample's peaks read left to right down the migration axis."""
        rows = []
        for r in sorted(self.records,
                        key=lambda r: (int(r["scan"]), int(r["col"]),
                                       int(r.get("gid", 0)))):
            rows.append({
                "file": r["file"], "well": r["well"], "scan": r["scan"],
                "channel": r["channel"], "base": r["base"], "kind": r["kind"],
                "start_scan": r.get("onset"), "end_scan": r.get("end"),
                "height_V": round(r["height"], 4),
                "area_Vscan": round(r["area"], 3),
                "duplex": self.duplex_of(r),
                "length_bp": self.length_bp,
                "ambiguous": bool(r.get("ambiguous")),
                "parent_scan": r.get("parent_scan", ""),
                "height_frac": r.get("height_frac", ""),
                "fraction": round(self.clust_frac(r), 4)
                if r["kind"] == "main" else "",
            })
        return rows

    def ambiguous_satellites(self):
        """Satellites tall enough that they may really be a minor allele.

        These are the positions where a heterozygote can be silently reported
        as a homozygote: the peak is tagged ``+A`` so ``clusters()`` never sees
        it and ``clust_frac()`` returns 0.0.  Callers should treat any position
        near one of these as unresolved rather than confident.
        """
        return [r for r in self.records
                if r.get("kind") == "+A" and r.get("ambiguous")]

    # ---------------------------------------------------------------- overlay
    def plot_overlay(self, ax):
        """Draw the picked records onto a matplotlib axes: shaded area, apex
        ring, 2nd-derivative start/end marks and a small label."""
        n = self.doc.acgt.shape[0]
        tick = float(np.nanmax(self.doc.acgt)) * 0.03
        for r in self.records:
            if not (0 <= r["left"] <= r["right"] < n):
                continue
            xs = np.arange(r["left"], r["right"] + 1)
            y = self.doc.acgt[xs, r["col"]]
            ax.fill_between(xs, 0, y, color=r["color"], alpha=0.25, zorder=1)
            ax.plot([r["scan"]], [self.doc.acgt[r["scan"], r["col"]]],
                    marker="o", ms=5, mfc="none", mec=r["color"], zorder=4)
            if self.show_d2:
                for xx in (r["onset"], r["end"]):
                    if xx is None or not (0 <= int(xx) < n):
                        continue
                    xx = int(xx)
                    yy = self.doc.acgt[xx, r["col"]]
                    ax.vlines(xx, max(0, yy - tick), yy + tick,
                              color=r["color"], lw=1.0, zorder=4)
                    ax.plot([xx], [yy], marker="s", ms=2.5, mfc=r["color"],
                            mec=r["color"], zorder=4)
            txt = f"{r['scan']}·{r['base']}"
            if r["kind"] != "main":
                txt = f"{r['kind']} " + txt
            ax.annotate(txt, (r["scan"], y.max()), textcoords="offset points",
                        xytext=(0, 6), fontsize=6, color=r["color"],
                        ha="center", zorder=5)


def _write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        if not rows:
            return
        w.writerow(list(rows[0].keys()))
        for r in rows:
            w.writerow(list(r.values()))


def _write_xlsx(path, rows):
    from openpyxl import Workbook
    from openpyxl.utils import get_column_letter
    wb = Workbook()
    ws = wb.active
    ws.title = "peaks"
    if rows:
        headers = list(rows[0].keys())
        ws.append(headers)
        for r in rows:
            ws.append([r[h] for h in headers])
        for i, h in enumerate(headers, 1):
            ws.column_dimensions[get_column_letter(i)].width = \
                max(8, min(28, 6 + len(h)))
        ws.freeze_panes = "A2"
    wb.save(path)


def save_table(path, rows):
    """Write peak rows to CSV, JSON or XLSX depending on the file suffix.

    An empty *rows* list is written as an empty table rather than raising."""
    ext = Path(path).suffix.lower()
    if ext == ".json":
        Path(path).write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    elif ext == ".xlsx":
        _write_xlsx(path, rows)
    else:
        _write_csv(path, rows)


# --------------------------------------------------------------------------- #
# drag to measure: a deliberately plainer manual area tool
# --------------------------------------------------------------------------- #
# PeakPicker decides for itself where a peak begins and ends, which is what you
# want for a CTC-CE duplex run and exactly what you do not want when you are
# measuring one specific hump by hand.  This tool takes the operator's word for
# it: press, drag, release.  The two scans you drag between are the baseline
# endpoints, and the area is whatever the trace stands above that line.
DRAG_MIN_SPAN = 3           # scans; closer than this was a click, not a drag
DRAG_MIN_AREA = 1e-9        # V*scan


def region_area(y, start, stop):
    """Baseline-corrected area of *y* between two hand-placed scans.

    The two scans are the markers the operator dragged, and they are used
    exactly as placed: the baseline is the straight line joining the trace value
    at *start* to the value at *stop*, and the area is the part of the trace
    standing above that line.  Same recipe as the automatic pickers
    (:meth:`PeakPicker._numeric_peak`) -- sum the clipped differences, no
    trapezoid -- but nothing here goes hunting for valleys, because choosing
    where the hump begins is the whole point.

    Reversed spans are normalised, so dragging right-to-left works.  Returns a
    dict, or None when the span is degenerate or the trace never rises above
    the line at all.
    """
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 2:
        return None
    start, stop = int(start), int(stop)
    if start > stop:
        start, stop = stop, start
    start = max(0, min(n - 1, start))
    stop = max(0, min(n - 1, stop))
    if stop - start < DRAG_MIN_SPAN:
        return None
    xs = np.arange(start, stop + 1)
    bl = np.linspace(float(y[start]), float(y[stop]), stop - start + 1)
    above = np.clip(y[xs] - bl, 0.0, None)
    area = float(np.sum(above))
    if area <= DRAG_MIN_AREA:
        return None
    rel = int(np.argmax(above))
    return {"start": start, "stop": stop,
            "midpoint": (start + stop) // 2,
            "peak_scan": start + rel,
            "height": float(above[rel]),
            "area": area,
            "baseline_left": float(y[start]),
            "baseline_right": float(y[stop])}


class DragAreaPicker:
    """Headless store for hand-measured peak areas on one trace.

    Two scans and a channel in, one area out.  There is deliberately no peak
    finding, no valley walking, no ``+A`` handling and no genotype call here --
    a region measure should be exactly what the dragged span says, so the
    genotype logic in :mod:`scorer` has nothing to second-guess.
    """

    def __init__(self, doc, path, colors=None, base_order=None,
                 theme_mode="base", show=None, run_name=None):
        self.doc = doc
        self.path = Path(path)
        self.show = show or (lambda col: True)
        self.colors = dict(colors or DEFAULT_COLORS)
        self.base_order = (base_order or "ACTG").upper()
        self.col_color = channel_colors(self.base_order, self.colors,
                                        theme_mode)
        # A run is the folder a well came out of; wells in one run are what make
        # a table worth having, so keep them apart in every row.
        self.run_name = str(run_name or self.path.parent.name or "")
        self.well = str(getattr(doc, "well", None) or self.path.stem)
        self.records: list[_Record] = []

    # ------------------------------------------------------------- measuring
    def best_col(self, start, stop):
        """The channel carrying the biggest hump in the span.

        Defaults the measurement to the peak the operator was pointing at
        rather than to whatever channel happens to be first.  Ties go to the
        lowest column so the same drag always reports the same channel.
        """
        best, best_area = None, -1.0
        for col in sorted(self.col_color):
            if not self.show(col):
                continue
            res = region_area(self.doc.acgt[:, col], start, stop)
            if res and res["area"] > best_area:
                best, best_area = col, res["area"]
        return best

    def add(self, start, stop, col=None):
        """Measure the span *start*..*stop* and record it.

        *col* is an ``acgt`` column; None picks the tallest hump in the span.
        Returns the new record, or None when the span is unusable.
        """
        if col is None:
            col = self.best_col(start, stop)
        if col is None or not (0 <= col < self.doc.acgt.shape[1]):
            return None
        res = region_area(self.doc.acgt[:, col], start, stop)
        if res is None:
            return None
        rec = _Record(
            file=str(self.path), run=self.run_name, well=self.well,
            scan=res["midpoint"], start_scan=res["start"],
            end_scan=res["stop"], midpoint=res["midpoint"],
            peak_scan=res["peak_scan"], channel=col + 1,
            base=CHANNEL_ORDER[col], kind="region",
            height_V=round(res["height"], 4),
            area_Vscan=round(res["area"], 3),
            baseline_left_V=round(res["baseline_left"], 4),
            baseline_right_V=round(res["baseline_right"], 4),
        )
        self.records.append(rec)
        return rec

    def undo_last(self):
        if self.records:
            return self.records.pop()
        return None

    def clear_all(self):
        self.records = []

    def rows(self):
        """Recorded measurements, oldest first.

        Click order rather than scan order, matching the on-screen table: the
        sequence of measurements is itself the record of what was done.
        """
        return list(self.records)

    def export_rows(self):
        """Measurements left to right along the migration axis, for saving."""
        return [dict(r) for r in sorted(self.records,
                                        key=lambda r: (r["start_scan"],
                                                       r["channel"]))]

    # ---------------------------------------------------------------- drawing
    def plot_overlay(self, ax, pending=None):
        """Draw the baseline and markers for every recorded measurement.

        *pending* is the span currently being dragged as ``(start, stop, col)``
        and is drawn dashed, so what is about to be recorded is visible before
        the mouse is released.
        """
        for r in self.records:
            self._draw_span(ax, r["start_scan"], r["end_scan"],
                            r["channel"] - 1,
                            color=self.col_color.get(r["channel"] - 1, "#666666"),
                            label=f"{r['area_Vscan']:.0f}")
        if pending is not None and pending[0] is not None and pending[1] is not None:
            col = pending[2] if len(pending) > 2 else None
            self._draw_span(ax, pending[0], pending[1], col,
                            color="#999999", label="", dashed=True)

    def _draw_span(self, ax, start, stop, col, color="#666666", label="",
                   dashed=False):
        """Baseline line, shaded area and endpoint guides for one span."""
        a, b = sorted((int(start), int(stop)))
        if b - a < 1:
            return
        ls = "--" if dashed else "-"
        if col is None or not (0 <= col < self.doc.acgt.shape[1]):
            col = self.best_col(a, b) or 0
        y = np.asarray(self.doc.acgt[:, col], dtype=float)
        xs = np.arange(a, b + 1)
        bl = np.linspace(float(y[a]), float(y[b]), b - a + 1)
        ax.plot(xs, bl, color=color, lw=1.1, ls=ls, alpha=0.9, zorder=4)
        ax.fill_between(xs, bl, np.maximum(y[xs], bl), color=color,
                        alpha=0.16, linewidth=0, zorder=1)
        for x in (a, b):
            ax.axvline(x, color=color, lw=0.7, ls=":", alpha=0.7, zorder=3)
        if not dashed and label:
            k = int(np.argmax(y[xs] - bl))
            ax.annotate(label, (xs[k], y[xs[k]]), textcoords="offset points",
                        xytext=(0, 5), ha="center", fontsize=6,
                        color=color, zorder=5)


class GenotypingEditor(ttk.Frame):
    """Click-to-pick peak editor. Usable as a standalone Toplevel widget pack
    (GenotypingDialog wraps it); channels can be switched off one at a time.
    The recognition engine is a PeakPicker, so the editor and the main window
    share exactly the same picking behaviour."""

    def __init__(self, master, path, colors=None, base_order=None,
                 theme_mode="base", show=None, on_exit=None):
        super().__init__(master)
        self.show = show or (lambda col: True)
        self.path = Path(path)
        self.colors = dict(colors or {"A": "#00AA00", "C": "#0000DD",
                                      "G": "#111111", "T": "#DD0000"})
        self.base_order = (base_order or "ACTG").upper()
        self.theme_mode = theme_mode
        self.on_exit = on_exit

        try:
            self.doc = load_trace(self.path, base_order=self.base_order)
        except Exception as e:
            raise ValueError(f"Could not load trace: {e}")

        self.include_sh = tk.BooleanVar(value=True)
        self.show_d2 = tk.BooleanVar(value=True)
        self.finder = tk.StringVar(value="best")
        self.pk = PeakPicker(self.doc, self.path, colors=self.colors,
                             base_order=self.base_order,
                             theme_mode=self.theme_mode, show=self.show,
                             include_sh=self.include_sh.get(),
                             show_d2=self.show_d2.get())

        self._build()
        self.redraw()

    def __getattr__(self, name):
        pk = self.__dict__.get("pk")
        if pk is not None:
            return getattr(pk, name)
        raise AttributeError(name)

    @property
    def records(self):
        return self.pk.records

    # ------------------------------------------------------------------ UI
    def _build(self):
        pane = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        pane.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(pane)
        pane.add(left, weight=3)
        self.fig = Figure(figsize=(9, 5.5), dpi=100, facecolor="#FFFFFF")
        self.fig.patch.set_facecolor("#FFFFFF")
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, window=left)
        self.toolbar.update()
        self.canvas.mpl_connect("button_press_event", self._on_click)

        right = ttk.Frame(pane, width=360)
        pane.add(right, weight=1)
        info = ttk.LabelFrame(right, text="Run", padding=6)
        info.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(info, text=str(self.path), wraplength=330, justify=tk.LEFT).pack(anchor=tk.W)
        sub = (f"{self.doc.n_scans} scans, {self.doc.source.lower()} "
               f"· base order {self.doc.base_order}"
               + (f" · called {len(self.doc.sequence)} bp"
                  if self.doc.sequence else " · not called"))
        ttk.Label(info, text=sub, foreground="#555").pack(anchor=tk.W)

        det = ttk.LabelFrame(right, text="Peak recognition", padding=6)
        det.pack(fill=tk.X, padx=4, pady=4)
        self.finder.trace_add("write", lambda *a: self._sync_opts(redraw=False))
        ttk.Label(det, text="Algorithm").pack(anchor=tk.W)
        ttk.Combobox(det, textvariable=self.finder, state="readonly",
                     values=[v for _, v in PEAK_FINDERS], width=30).pack(anchor=tk.W)
        self.include_sh = tk.BooleanVar(value=True)
        ttk.Checkbutton(det, text="Add +A (A-addition) peak",
                        variable=self.include_sh,
                        command=self._sync_opts).pack(anchor=tk.W)
        self.show_d2 = tk.BooleanVar(value=True)
        ttk.Checkbutton(det, text="Mark start/end from the 2nd derivative",
                        variable=self.show_d2,
                        command=self._sync_opts).pack(anchor=tk.W)
        self.het_label = ttk.Label(det, text="Click a peak (or near it) to pick it.",
                                   foreground="#0F3A6E", wraplength=320, justify=tk.LEFT)
        self.het_label.pack(anchor=tk.W, pady=(4, 0))

        std = ttk.LabelFrame(right, text="Internal standard (CTC-CE duplex pattern)",
                             padding=6)
        std.pack(fill=tk.X, padx=4, pady=4)
        fw = 330
        ttk.Label(std, text="The four standard peaks are the SAME fragment (same bp),\n"
                            "separated by cycling-temperature CE according to sequence:\n"
                            "peaks 1-2 = the two homoduplexes (they differ by the single\n"
                            "SNP base of the rs number); peaks 3-4 = the two heteroduplexes\n"
                            "made in the PCR by Watson/Crick re-annealing (one mismatch base).",
                  justify=tk.LEFT, foreground="#555", wraplength=fw).pack(anchor=tk.W, pady=(0, 4))
        sz = ttk.Frame(std)
        sz.pack(fill=tk.X)
        ttk.Label(sz, text="Fragment length (bp)").pack(side=tk.LEFT)
        self.len_entry = ttk.Entry(sz, width=8, justify=tk.RIGHT)
        self.len_entry.pack(side=tk.LEFT, padx=4)
        ttk.Label(sz, text="optional — same for all four",
                  foreground="#888").pack(side=tk.LEFT)
        sbtn = ttk.Frame(std)
        sbtn.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(sbtn, text="Mark picked peaks as standard",
                   command=self._mark_std).pack(side=tk.LEFT)
        ttk.Button(sbtn, text="Clear", command=self._clear_std).pack(side=tk.LEFT, padx=4)
        ttk.Label(std, text="The four first main peaks (scan order) are tagged\n"
                            "HOM1, HOM2 (homoduplexes) then HET1, HET2 (heteroduplexes).",
                  foreground="#777", wraplength=fw, justify=tk.LEFT).pack(anchor=tk.W, pady=(4, 0))
        self.std_lbl = ttk.Label(std, text="No standard set.",
                                 foreground="#555", wraplength=fw, justify=tk.LEFT)
        self.std_lbl.pack(anchor=tk.W, pady=(4, 0))

        bars = ttk.Frame(right)
        bars.pack(fill=tk.X, padx=4, pady=2)
        ttk.Button(bars, text="Undo last",
                   command=self._undo_last).pack(side=tk.LEFT)
        ttk.Button(bars, text="Clear all",
                   command=self._clear_all).pack(side=tk.LEFT, padx=4)

        tblf = ttk.LabelFrame(right, text="Picked peaks",
                              padding=4)
        tblf.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        cols = ("#", "scan", "duplex", "ch", "kind", "height V", "area V·sc", "frac")
        self.tree = ttk.Treeview(tblf, columns=cols, show="headings", height=12)
        widths = {"#": 34, "scan": 54, "duplex": 58, "ch": 40, "kind": 66,
                  "height V": 70, "area V·sc": 78, "frac": 50}
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=widths[c], anchor="e" if c not in ("#", "kind") else "w",
                             stretch=(c in ("scan", "kind")))
        vs = ttk.Scrollbar(tblf, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=vs.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vs.pack(side=tk.RIGHT, fill=tk.Y)

        http_bar = ttk.Frame(right)
        http_bar.pack(fill=tk.X, padx=4, pady=4)
        ttk.Button(http_bar, text="Save table…", command=self._save).pack(side=tk.LEFT)
        ttk.Button(http_bar, text="Close / back to viewer",
                   command=self._request_close).pack(side=tk.RIGHT)

        self.frac_var = tk.StringVar(value="No positions with two peaks yet.")
        fout = ttk.Label(right, textvariable=self.frac_var,
                         foreground="#555", wraplength=340, justify=tk.LEFT)
        fout.pack(fill=tk.X, padx=4, pady=(0, 6))

    # ------------------------------------------------------------- picking
    def _on_click(self, event):
        if event.xdata is None or event.inaxes is None:
            return
        if self.toolbar.mode != "":
            return
        self._pick_peak(int(round(event.xdata)), vol=event.ydata)

    def _pick_peak(self, scan, vol=None):
        rec = self.pk.pick(scan, vol=vol)
        if rec is None:
            if self.pk._reject == "area":
                self._status("That area is already picked — undo it first to "
                             "pick it again.")
            else:
                self._status("No peak found near that scan — try again closer "
                             "to a hump.")
            return
        self._status(f"Peak at scan {rec['scan']} · {rec['base']}"
                     f" height {rec['height']:.3f} V, area {rec['area']:.2f} V·scan")
        self.redraw()
        self._sync_table()

    def _sync_opts(self, redraw=True):
        self.pk.include_sh = bool(self.include_sh.get())
        self.pk.show_d2 = bool(self.show_d2.get())
        self.pk.finder = self.finder.get()
        if redraw:
            self.redraw()
            self._sync_table()

    def _request_close(self):
        """Leave peak picking: back to the main viewer when embedded."""
        if self.on_exit is not None:
            self.on_exit()
        else:
            self.destroy()

    def _close(self):
        self._request_close()

    def _status(self, msg):
        self.het_label.config(text=msg)

    # ------------------------------------- internal-standard CTC-CE duplexes
    def _mark_std(self):
        try:
            val = float(self.len_entry.get())
        except ValueError:
            val = None
        try:
            msg = self.pk.mark_std(val)
        except ValueError as e:
            self.std_lbl.config(text=str(e), foreground="#A33")
            return
        self._sync_table()

    def _clear_std(self):
        self.pk.clear_std()
        self.std_lbl.config(text="No standard set.", foreground="#555")
        self._sync_table()

    # ---------------------------------------------------------------- view
    def redraw(self):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n = self.doc.acgt.shape[0]
        x = np.arange(n)
        for col, color in self.pk.col_color.items():
            if not self.show(col):
                continue
            ax.plot(x, self.doc.acgt[:, col], color=color, lw=0.7,
                    label=f"Ch{col + 1} {CHANNEL_ORDER[col]}")
        # base letters of an existing call, faint, for orientation
        if self.doc.sequence and self.doc.peak_positions:
            seq = self.doc.sequence
            pos = np.asarray(self.doc.peak_positions, dtype=int)
            top = float(np.nanmax(self.doc.acgt))
            for pi, b in enumerate(seq):
                if pi >= len(pos) or b not in CHANNEL_ORDER:
                    continue
                if not self.show(CHANNEL_ORDER.index(b)):
                    continue
                ax.text(pos[pi], top * 1.01, b, ha="center", va="bottom",
                        fontsize=5, color=self.colors.get(b, "#444"), alpha=0.85, zorder=2)
        self.pk.plot_overlay(ax)
        ax.set_xlim(0, n)
        ax.set_ylim(0, float(np.nanmax(self.doc.acgt)) * 1.08 or 1.0)
        ax.set_xlabel("scan")
        ax.set_ylabel("V")
        ax.set_title(f"{self.path.name} — click to pick a peak", fontsize=9)
        ax.grid(True, alpha=0.15)
        ax.legend(loc="upper right", fontsize=7, ncol=2, framealpha=0.6)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._sync_fractions()

    def _sync_fractions(self):
        lines = []
        for cl in self.pk.clusters():
            if len(cl) >= 2 and sum(m["area"] for m in cl) > 0:
                areas = sorted(m["area"] for m in cl)
                frac = areas[0] / sum(areas)
                chans = "+".join(sorted(m["base"] for m in cl))
                lines.append(f"scan {cl[0]['scan']}  {chans}:  "
                             f"variant {frac:.3f}")
        self.frac_var.set("\n".join(lines) if lines
                          else "No position with two peaks yet (mutant "
                               "fraction appears here).")

    # ---------------------------------------------------------------- table
    def _sync_table(self):
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(self.records, 1):
            fr = self._clust_frac(r) if r["kind"] == "main" else ""
            self.tree.insert("", tk.END, values=(
                i, r["scan"], self._duplex_of(r), r["base"], r["kind"],
                f"{r['height']:.3f}", f"{r['area']:.1f}",
                f"{fr:.3f}" if fr else ""))

    def _undo_last(self):
        if not self.pk.undo_last():
            return
        self.redraw()
        self._sync_table()
        self._status("Removed last picked peak.")

    def _clear_all(self):
        self.pk.clear_all()
        self.redraw()
        self._sync_table()
        self._status("Table cleared.")

    # -------------------------------------------------------------- export
    def _save(self):
        if not self.records:
            messagebox.showinfo("Save table", "Pick some peaks first.", parent=self)
            return
        types = [("CSV (Excel-compatible)", "*.csv"), ("JSON (ML)", "*.json")]
        try:
            import openpyxl  # noqa: F401
            types.insert(1, ("Excel workbook (.xlsx)", "*.xlsx"))
        except ImportError:
            pass
        path = filedialog.asksaveasfilename(parent=self, defaultextension=".csv",
                                            filetypes=types)
        if not path:
            return
        rows = self.pk.export_rows()
        try:
            save_table(path, rows)
        except Exception as e:
            messagebox.showerror("Save table", f"Could not write file:\n{e}", parent=self)
            return
        self._status(f"Saved {len(rows)} peak rows to {path}")


class GenotypingDialog(tk.Toplevel):
    """Standalone window wrapping a GenotypingEditor (kept for scripts/tests;
    the app picks peaks directly in the main window instead, batching several
    samples at a time)."""

    def __init__(self, parent, path, colors=None, base_order=None,
                 theme_mode="base"):
        super().__init__(parent)
        self.title("Manual genotyping — " + Path(path).name)
        self.geometry("1080x680")
        try:
            ed = GenotypingEditor(self, path, colors=colors, base_order=base_order,
                                  theme_mode=theme_mode)
        except Exception as e:
            messagebox.showerror("Genotyping", f"Could not load trace:\n{e}",
                                 parent=self)
            self.destroy()
            return
        self._editor = ed
        ed.pack(fill=tk.BOTH, expand=True)
        self.bind("<Control-z>", lambda e: self._editor._undo_last())
        self.bind("<Control-Z>", lambda e: self._editor._undo_last())

    def __getattr__(self, name):
        ed = self.__dict__.get("_editor")
        if ed is None:
            raise AttributeError(name)
        return getattr(ed, name)
