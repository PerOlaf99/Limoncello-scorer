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
from itertools import combinations
from tkinter import filedialog, messagebox, ttk
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from analyzer_core import load_trace

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

# --- plate-learned geometry -------------------------------------------------
#
# The four constants above are NOT assay-general.  They were fitted to the T9
# rs1695 plate, and on the ABCC2 N10 plate the operator's own hand marks break
# all three of them at once:
#
#   IS_MIN_SPACING  25.0   but ABCC2's d3 (the one-SNP gap HET1->HET2) is 12
#   IS_MID_HI        3.3   but ABCC2's d2/d1 is 3.36
#   IS_SAT_SKIP  (10, 32)  and ABCC2's HET2 sits 12 scans behind HET1, slightly
#                          shorter -- so it was being deleted as a +A satellite
#                          before the quartet search ever ran
#
# The measured spacings are very stable (d1 IQR 7, d3 IQR 2 scans across 22
# wells) while the whole quartet slides 418 scans across the plate.  So the
# *shape* is a property of the assay and belongs in a model fitted from the
# operator's marks; only the *position* drifts per capillary.  T9 stays the
# fallback so an unmarked plate still behaves exactly as it did.


class ISGeometry:
    """The expected shape of one internal standard.

    ``d1`` is H1->H2 and ``d2`` is H2->the heteroduplex region.  Those two gaps
    are what genotyping needs: the standard is read to establish the nature of
    the two homoduplexes, and everything downstream of H2 is only evidence that
    the sample carried both alleles.  ``d3`` (HET1->HET2) is therefore *not*
    load-bearing.  It is kept only as a refinement, because the heteroduplex
    pair is where allelic imbalance shows up, and that analysis needs both bands
    measured separately.

    So ``d3`` is ``None`` whenever the heteroduplexes co-migrate, which is the
    normal state rather than an edge case: the run temperature is optimised to
    separate the two homoduplexes and the mismatched heteroduplexes are the
    least stable species, so they merge whenever the setpoint is tuned for H1 vs
    H2.  Measured d3 spans 5 to 116 scans across fragments, so it is set by the
    mismatch and the temperature together and cannot be a constant either.
    """

    __slots__ = ("d1", "d2", "d3", "tol1", "tol2", "tol3", "center_tol")

    def __init__(self, d1, d2, d3, tol1=None, tol2=None, tol3=None,
                 center_tol=None):
        self.d1, self.d2 = float(d1), float(d2)
        self.d3 = None if d3 is None else float(d3)
        # A quarter of the gap, floored so a near-zero d3 (the SNP gap) still
        # admits an integer separation rather than collapsing to +-1.
        self.tol1 = tol1 if tol1 is not None else max(3.0, 0.25 * self.d1)
        self.tol2 = tol2 if tol2 is not None else max(5.0, 0.25 * self.d2)
        self.tol3 = (None if self.d3 is None
                     else (tol3 if tol3 is not None else max(3.0, 0.35 * self.d3)))
        self.center_tol = center_tol

    @property
    def merged(self):
        """True when the heteroduplexes co-migrate (three bands, not four).

        This says nothing about whether the well can be genotyped -- it always
        can; see :meth:`matches_hom`.
        """
        return self.d3 is None

    @classmethod
    def from_quartets(cls, quartets, center_tol=None):
        """Fit a template from the operator's marked quartets.

        Each entry is four ascending scan positions.  Robust spreads (median
        absolute deviation) rather than min/max, so one mis-marked well widens
        the window slightly instead of blowing it open.
        """
        qs = [sorted(int(x) for x in q) for q in quartets
              if q is not None and len(q) == 4]
        if len(qs) < 2:
            raise ValueError("need at least 2 marked quartets to fit a geometry")
        import statistics as _st

        d1 = _st.median([q[1] - q[0] for q in qs])
        d2 = _st.median([q[2] - q[1] for q in qs])
        d3 = _st.median([q[3] - q[2] for q in qs])
        ctr = [sum(q) / 4.0 for q in qs]
        med = _st.median(ctr)
        # cover every marked well, then a margin for an unmarked capillary
        span = max(abs(c - med) for c in ctr) + 25.0

        def _tol(i, floor):
            # Widest deviation actually observed, not a robust spread estimate.
            # MAD*2 came out narrower than the operator's own marks (d1 needed
            # 12 scans and the fit gave 10, d3 needed 5 and the fit gave 3), and
            # a tolerance that rejects a hand-marked well is worse than useless.
            return max(floor, max(abs(q[i + 1] - q[i] - _st.median(
                [x[i + 1] - x[i] for x in qs])) for q in qs) * 1.25)

        return cls(d1, d2, d3,
                   tol1=_tol(0, 4.0), tol2=_tol(1, 6.0), tol3=_tol(2, 4.0),
                   center_tol=span if center_tol is None else center_tol)

    @classmethod
    def from_triplets(cls, triplets, d1=None, d2=None):
        """Fit a template from three-band wells.

        *d1*/*d2* are normally taken from the same run's split wells, since the
        two homoduplexes separate identically whether or not the heteroduplexes
        do.  When a run has only three-band wells there is nothing to borrow
        from, so the gaps come from the triplets themselves -- which is the
        normal case for the fully merged runs (GSTA1 N11 and ABCC2 T1 are 37/38
        and 13/14 three-band).  Only the third gap is unmeasurable.
        """
        ts = [sorted(int(x) for x in t) for t in triplets
              if t is not None and len(t) == 3]
        if not ts:
            return None
        import statistics as _st
        if d1 is None:
            d1 = _st.median([t[1] - t[0] for t in ts])
        if d2 is None:
            d2 = _st.median([t[2] - t[1] for t in ts])
        g = cls(d1, d2, None)
        g.tol1 = max(4.0, max(abs(t[1] - t[0] - d1) for t in ts) * 1.25)
        g.tol2 = max(6.0, max(abs(t[2] - t[1] - d2) for t in ts) * 1.25)
        return g

    @classmethod
    def t9_default(cls):
        """The shape implied by the original T9 constants."""
        return cls(IS_MIN_SPACING, 2.5 * IS_MIN_SPACING, IS_MIN_SPACING)

    def describe(self):
        d3 = "het merged" if self.merged else f"d3={self.d3:.0f}±{self.tol3:.0f}"
        return (f"d1={self.d1:.0f}±{self.tol1:.0f} "
                f"d2={self.d2:.0f}±{self.tol2:.0f} {d3}")

    def matches_hom(self, a, b):
        """H1->H2 only.  The minimal test a genotyping call depends on."""
        return b - a > 0 and abs((b - a) - self.d1) <= self.tol1

    def matches(self, a, b, c, d=None):
        """*c* is the start of the heteroduplex region; *d* is HET2, if resolved.

        With no *d* the heteroduplexes are merged and the shape is fully
        specified, so a three-band standard is a complete match rather than a
        degraded one.
        """
        if d is None:
            return (b - a > 0 and c - b > 0
                    and abs((b - a) - self.d1) <= self.tol1
                    and abs((c - b) - self.d2) <= self.tol2)
        if self.merged:
            return False
        d1, d2, d3 = b - a, c - b, d - c
        if d1 <= 0 or d2 <= 0 or d3 <= 0:
            return False
        return (abs(d1 - self.d1) <= self.tol1
                and abs(d2 - self.d2) <= self.tol2
                and abs(d3 - self.d3) <= self.tol3)


def weakest_std_band_snr(acgt, is_col, scans):
    """Apex of the standard's weakest band, in units of the channel noise."""
    y = np.asarray(acgt[:, is_col], dtype=float)
    sigma = _noise_sigma(y)
    if sigma <= 0 or not len(scans):
        return 0.0
    base = float(np.median(y[:400]))
    return min(max(0.0, (float(y[max(0, int(s) - 6):int(s) + 7].max()) - base)
                      / sigma) for s in scans)


def prime_plate_std_snr(model, wells):
    """Fill a model's weak-standard reference from the wells that were marked.

    *wells* is an iterable of ``(acgt, is_col, scans)``.  This is the automated
    form of what the operator does by eye: judge a peak against the wells above
    and below it.  Call it once over the marked wells before calling any
    unmarked one, so a weak standard is flagged relative to this plate rather
    than against a fixed cutoff.
    """
    if model is None:
        return
    model.std_snr = [weakest_std_band_snr(a, c, s) for a, c, s in wells if s]
    model.std_snr = [v for v in model.std_snr if v > 0]


def prime_plate_co_migration(model, wells):
    """Fill a model's sample/standard migration offsets from the marked wells.

    *wells* is an iterable of ``(acgt, samp_col, scans)``, the same arguments
    :func:`prime_plate_std_snr` takes.  For each duplex the offset of the
    strongest peak near its band centre is recorded.

    Only wells where that duplex is the *clearest* band present contribute a
    reference.  This matters more than it sounds: on a plate that is mostly
    hom-1 there is no H2 product anywhere, yet the H2 neighbourhood is full of
    ordinary noise that clears any absolute threshold.  Letting those wells
    teach the model where H2 lives produced a gate 300 scans wide -- wider than
    the segment it was meant to constrain, and so wide that the very off-target
    peak it was built to exclude fell back inside.  Requiring the duplex to be
    the strongest band means H2 collects references only from the wells where
    H2 really is present, which on that plate is too few to gate it.  A duplex
    with too little reference keeps the old fixed-radius search, which is the
    status quo and cannot regress a plate that was already working.
    """
    if model is None:
        return
    import scorer
    collected = [[] for _ in range(4)]
    for acgt, samp_col, scans in wells:
        if not scans:
            continue
        y = np.asarray(acgt[:, samp_col], dtype=float)
        sigma = _noise_sigma(y)
        n = len(y)
        if sigma <= 0:
            continue
        idx = np.arange(n)
        segs = _quartet_segments(scans)
        snrs = []
        for k, (lo, hi) in enumerate(segs):
            lo, hi = max(0, lo), min(n, hi)
            if hi - lo < MIN_SEGMENT_SPAN:
                snrs.append(0.0)
                continue
            m = (idx >= max(0, lo - SEGMENT_BASELINE_PAD)) & \
                (idx <= min(n, hi + SEGMENT_BASELINE_PAD))
            base = float(np.median(y[m]))
            ctr = int(scans[k])
            a, b = max(0, ctr - CO_MIGRATION_SEARCH), min(n, ctr + CO_MIGRATION_SEARCH + 1)
            win = y[a:b]
            snrs.append((float(win.max()) - base) / sigma if b > a else 0.0)
        if not snrs:
            continue
        loudest = max(range(len(snrs)), key=lambda i: snrs[i])
        if snrs[loudest] < scorer.T9_MIN_DOMINANT_SIGMA:
            continue
        for k in range(len(model.co_mig)):
            if k == loudest and len(snrs) > k:
                ctr = int(scans[k])
                a, b = max(0, ctr - CO_MIGRATION_SEARCH), min(n, ctr + CO_MIGRATION_SEARCH + 1)
                lo, hi = max(0, segs[k][0]), min(n, segs[k][1])
                if b > a and hi - lo >= MIN_SEGMENT_SPAN:
                    m = (idx >= max(0, lo - SEGMENT_BASELINE_PAD)) & \
                        (idx <= min(n, hi + SEGMENT_BASELINE_PAD))
                    base = float(np.median(y[m]))
                    collected[k].append(
                        float(int(np.argmax(y[a:b])) + a - ctr))
    for k in range(len(model.co_mig)):
        model.co_mig[k] = collected[k]


def acgt_index_for_channel(base_order: str, channel: int) -> int:
    """Column of ``doc.acgt`` holding MegaBACE *channel* (1-based).

    ``doc.acgt`` is always in ACGT order, but the physical channel order is the
    plate's dye order -- "ACTG" on a MegaBACE, so Ch1=A, Ch2=C, Ch3=T, Ch4=G.
    Getting this backwards is easy and silently swaps the sample and standard
    channels, so route every channel lookup through here.
    """
    order = (base_order or "ACTG").upper()
    if len(order) != 4 or set(order) != set("ACGT"):
        raise ValueError(f"Unexpected base_order {base_order!r}")
    if not 1 <= channel <= 4:
        raise ValueError(f"Channel must be 1..4, got {channel!r}")
    return "ACGT".index(order[channel - 1])


def _is_candidates(trace, cut=1900, geometry=None):
    """Despiked, baseline-corrected peak candidates of the standard channel.

    *geometry* changes the +A suppression.  The fixed ``IS_SAT_SKIP`` window of
    10-32 scans is a T9 measurement; on ABCC2 the real HET1/HET2 pair is 12
    scans apart, so that filter was deleting HET2 before the search began.  When
    a fitted geometry is supplied, only offsets that do *not* line up with one
    of the fitted gaps are treated as satellites.
    """
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
    if geometry is None:
        gaps = None
    else:
        gaps = [(geometry.d1, geometry.tol1), (geometry.d2, geometry.tol2)]
        if not geometry.merged:
            gaps.append((geometry.d3, geometry.tol3))
    # A +A tail is always shorter than and close behind its parent peak; drop it
    # so it cannot stand in for a real member of the quartet.
    kept = []
    for x, h, prom in sorted(cands, key=lambda c: -c[1]):
        sat = False
        for kx, kh, _ in kept:
            dx = x - kx
            if not (IS_SAT_SKIP[0] <= dx <= IS_SAT_SKIP[1]):
                continue
            if kh <= IS_SAT_RATIO * h:
                continue
            if gaps and any(abs(dx - g) <= t for g, t in gaps):
                continue          # this offset is a real duplex gap, keep it
            sat = True
            break
        if sat:
            continue
        kept.append((x, h, prom))
    return sorted(kept, key=lambda c: c[0])[:IS_MAX_CANDIDATES]


def _legacy_matches(a, b, c, d):
    """The original hard-coded T9 geometry test, kept for unmarked plates."""
    d1, d2, d3 = b - a, c - b, d - c
    if d1 < IS_MIN_SPACING or d3 < IS_MIN_SPACING:
        return False
    dm = (d1 + d3) / 2.0
    return (abs(d1 - d3) <= IS_OUTER_TOL * dm
            and IS_MID_LO * dm <= d2 <= IS_MID_HI * dm)


def _equimolar(heights):
    lo, hi = min(heights), max(heights)
    return lo > 0 and lo / hi >= IS_EQUIMOLAR_MIN


def _search_legacy_quartet(pk):
    """Best quartet under the original hard-coded T9 geometry."""
    best = None
    for combo in combinations(range(len(pk)), 4):
        picks = [pk[i] for i in combo]
        scans = [p[0] for p in picks]
        if not _legacy_matches(*scans):
            continue
        h = [p[1] for p in picks]
        if not _equimolar(h):
            continue
        score = sum(p[2] for p in picks)
        if best is None or score > best[0]:
            best = (score, scans, h)
    return None if best is None else (best[1], best[2])


def _refine_het(pk, scans, heights, geometry):
    """Attach HET2 when the heteroduplex pair really is resolved.

    Purely additive.  A missing or unresolvable HET2 leaves the three-band core
    untouched, and the caller still has everything genotyping requires; what it
    loses is the separate HET2 area, which is the quantity an allelic-imbalance
    read needs.
    """
    if geometry is None or geometry.merged:
        return scans, heights
    c = scans[-1]
    best = None
    for x, h, prom in pk:
        if x <= c or abs((x - c) - geometry.d3) > geometry.tol3:
            continue
        if not _equimolar(heights + [h]):
            continue
        if best is None or prom > best[0]:
            best = (prom, x, h)
    if best is None:
        return scans, heights
    return scans + [best[1]], heights + [best[2]]


def _search_core(pk, geometry):
    """Best H1 + H2 + heteroduplex-region match, by summed prominence.

    This is the shape a genotyping call actually needs, and it is a complete
    shape: whether the heteroduplexes inside that region resolve into two bands
    is not part of the test.  Anything that only works when HET1/HET2 separate
    is rejecting a well the assay can still genotype, which is what made the
    IL10 N7 split template score 9/20 on wells whose H1/H2 were all correct.

    Candidates are scored *including* any HET2 they admit.  That matters
    because the third band of the core is otherwise ambiguous: on a resolved
    standard both HET1 and HET2 sit close enough to the H2->HET1 gap to satisfy
    the core on their own, and scoring the bare core picked whichever was
    taller.  Picking HET2 as the core left the search hunting for a fifth band
    and silently dropped the resolved pair back to three.
    """
    best = None
    for combo in combinations(range(len(pk)), 3):
        picks = [pk[i] for i in combo]
        scans = [p[0] for p in picks]
        if not geometry.matches(scans[0], scans[1], scans[2], None):
            continue
        heights = [p[1] for p in picks]
        if not _equimolar(heights):
            continue
        found_scans, found_heights = _refine_het(pk, scans, heights, geometry)
        score = sum(p[2] for p in picks)
        if len(found_scans) > len(scans):
            # A resolved pair is strictly better evidence than one band.
            score = sum(h for h in found_heights) * 2.0
        if best is None or score > best[0]:
            best = (score, found_scans, found_heights)
    return None if best is None else (best[1], best[2])


def _search_quartets(pk, geometry, want=4):
    """Backwards-compatible entry point used by the unmarked-plate path."""
    if geometry is None:
        return _search_legacy_quartet(pk)
    return _search_core(pk, geometry)


def find_is_quartet(trace, cut=1900, geometry=None, center=None, center_tol=None):
    """Locate the four internal-standard peaks, or ``None``.

    *trace* is one channel's samples.  Returns ``(peaks, heights)`` where
    *peaks* are the four scan positions in ascending order, chosen as the
    highest-prominence candidate quartet that satisfies the geometry.  Found in
    96/96 wells of the T9 rs1695 plate.

    With *geometry* supplied, the quartet is matched against that fitted shape
    instead of the hard-coded T9 ratios.  That matters because the T9 ratios
    reject every genuine quartet on other assays (see :class:`ISGeometry`).

    *center*/*center_tol* restrict the search to a window around where the
    quartet is expected.  This filters the candidate list rather than slicing
    the trace, because the baseline is estimated from the first 400 samples and
    a slice starting mid-trace would measure that region as "baseline" and
    corrupt every height in the window.
    """
    pk = _is_candidates(trace, cut, geometry=geometry)
    if center is not None and center_tol:
        span = (geometry.d1 + (geometry.d3 or 0.0) + geometry.d2) / 2.0 \
            if geometry else 0.0
        lo, hi = center - center_tol - span, center + center_tol + span
        pk = [c for c in pk if lo <= c[0] <= hi]
    if geometry is None:
        # Unmodelled plate: the legacy T9 test still requires a full quartet.
        return _search_legacy_quartet(pk) if len(pk) >= 4 else None
    if len(pk) < 3:
        return None
    core = _search_core(pk, geometry)
    if core is None:
        return None
    return core


class PlateISModel:
    """Where the standard quartet sits on one plate, learned from marks.

    Detection is only reliable once the plate's own shape and position range
    are known: across the operator-marked ABCC2 wells the spacings barely
    move (d1 IQR 14, d3 IQR 6 scans) while the quartet centre slides 418 scans
    capillary to capillary.  So the geometry is fitted from the marks, and each
    unmarked well is then searched only inside that window.

    A model can carry two geometries -- one for wells whose heteroduplexes
    split and one for wells whose heteroduplexes co-migrate -- because that
    difference is per capillary, not per plate.
    """

    def __init__(self, geometry, centers, n_marked=0, merged_geometry=None,
                 merged_seen=False):
        import statistics as _st
        self.geometry = geometry
        self.centers = sorted(float(c) for c in centers)
        self.n_marked = int(n_marked)
        self.merged_geometry = merged_geometry
        # Whether the operator's own marks on this plate ever show a merged
        # heteroduplex.  Decides what a three-band detection means here: routine
        # on a plate that runs merged (GSTA1 N11 is 37/38), and off-model on one
        # that never does (all 56 ABCC2_N10 marks are four-band, so its two
        # three-band reads are a standard that failed to form).
        self.merged_seen = bool(merged_seen)
        # Weakest standard band's apex over the channel noise, per marked well.
        # Judged against the plate rather than an absolute floor, because a weak
        # peak is read by comparison with its neighbours: the eye aligns the
        # pattern above and below and still recognises it.  An absolute cutoff
        # would throw away a well that is plainly readable in context.
        self.std_snr = []
        # Per-duplex sample/standard migration offsets, in scans, relative to
        # each fitted band centre.  The sample and its standard are in the same
        # tube, so a real product lands a short, repeatable distance from the
        # band it belongs to -- but the distance is not zero, and it is not the
        # same for every duplex (a heteroduplex has its own mobility).  Filled
        # by prime_plate_co_migration() from the marked wells.
        self.co_mig = [[] for _ in range(4)]

    @classmethod
    def from_marks(cls, marks):
        """Fit from every marked well, splitting by 4-band and 3-band wells.

        *marks* is an iterable of ascending scan positions per well.  Four-band
        wells give the full template; three-band wells describe wells whose
        heteroduplexes co-migrate and produce a second, merged template.
        """
        ms = [sorted(int(x) for x in m) for m in marks if m]
        qs = [m for m in ms if len(m) == 4]
        ts = [m for m in ms if len(m) == 3]
        if len(qs) < 2:
            if len(ts) < 2:
                return None
            # Only three-band wells available: d1/d2 come straight from them.
            merged = ISGeometry.from_triplets(ts)
            return cls(merged, [sum(t) / 3.0 for t in ts], len(ms), merged,
                       merged_seen=True)
        geom = ISGeometry.from_quartets(qs)
        centers = [sum(m) / 4.0 for m in ms]
        merged = None
        if ts:
            merged = ISGeometry.from_triplets(ts, d1=geom.d1, d2=geom.d2)
        return cls(geom, centers, len(ms), merged, merged_seen=bool(ts))

    @classmethod
    def from_quartets(cls, quartets):
        return cls.from_marks(quartets)

    @property
    def center(self):
        import statistics as _st
        return _st.median(self.centers)

    def weak_std_threshold(self):
        """Signal level below which a standard is worth a second look.

        Derived from this plate's own marked wells, not fixed: a quarter of the
        plate's typical weakest-band signal.  Nothing is rejected on it -- it
        only raises a flag, because a weak standard is still readable in the
        company of its neighbours.
        """
        if len(self.std_snr) < 4:
            return None
        import statistics as _st
        return max(5.0, 0.25 * _st.median(self.std_snr))

    def co_migration_window(self, duplex):
        """Scans either side of a fitted band centre where its product sits.

        Returns ``(lo, hi)`` as offsets from the band centre, or ``None`` when
        this plate has too little reference material to gate that duplex.  The
        gate exists because a duplex's segment is wide: with ``d1`` at 81 scans
        the H2 window spans ~170, and a strong peak from somewhere else entirely
        falls inside it and gets integrated as if it were the second allele.
        Searching only where the product belongs removes that failure at the
        source instead of trying to recognise the wrong peak afterwards.
        """
        refs = self.co_mig[duplex] if 0 <= duplex < len(self.co_mig) else []
        if len(refs) < CO_MIGRATION_MIN_REFS:
            return None
        import statistics as _st
        med = _st.median(refs)
        mad = _st.median([abs(v - med) for v in refs])
        span = max(CO_MIGRATION_MAD_K * mad, CO_MIGRATION_MIN_SPAN)
        # A gate wider than any plausible product spread is not a measurement,
        # it is a contaminated reference set.  Refuse it: the ungated fallback
        # is strictly the old behaviour, so declining can only miss the
        # improvement, never introduce a regression.
        if span > CO_MIGRATION_MAX_HALF_SPAN:
            return None
        return med - span, med + span

    def window(self):
        """Half-width of the search window around the predicted centre."""
        c = self.center
        lo = min(abs(c - self.centers[0]), abs(self.centers[-1] - c))
        return max(self.geometry.center_tol or 0.0, 40.0) + max(lo, 0.0)

    def find(self, trace, cut=1900):
        """Anchored search, then widened.  Heteroduplexes are never load-bearing.

        A single template covers both three- and four-band wells, so there is no
        split/merged branch to get wrong: the search matches H1 + H2 + the
        heteroduplex region, then attaches HET2 only if the trace actually
        resolves it.  ``merged_geometry`` is still consulted afterwards for a
        plate whose homoduplex gaps genuinely differ between capillaries.
        """
        found = find_is_quartet(trace, cut=cut, geometry=self.geometry,
                                center=self.center, center_tol=self.window())
        if found is not None:
            return found[0], found[1], "anchored"
        found = find_is_quartet(trace, cut=cut, geometry=self.geometry)
        if found is not None:
            return found[0], found[1], "anchored-wide"
        if self.merged_geometry is not None:
            found = find_is_quartet(trace, cut=cut, geometry=self.merged_geometry,
                                    center=self.center, center_tol=self.window())
            if found is not None:
                return found[0], found[1], "anchored-merged"
            found = find_is_quartet(trace, cut=cut, geometry=self.merged_geometry)
            if found is not None:
                return found[0], found[1], "anchored-merged-wide"
        return None


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

# Sample/standard co-migration.  The sample product does not sit on its
# standard band: on ABCC2_N10 it lands a median 11 scans *before* the fitted
# centre, and every genuine hom-1 there falls between -20 and -7.  These
# describe how wide to allow that band, measured per plate and per duplex.
CO_MIGRATION_MAD_K = 5.0     # robust spread multiplier for the gate
CO_MIGRATION_MIN_SPAN = 8.0  # scans; floor so a tight plate cannot go knife-edge
CO_MIGRATION_MIN_REFS = 8    # marked wells needed before a duplex is gated
# Eight is deliberately above the five that would have let a noisy duplex
# through.  At five, CYBA_N1 installed an H3 gate at +21..+41 from five
# reference wells and thereby moved two real heteroduplexes below the
# significance floor -- its own references contradicted the wells it changed.
# A gate moves calls, so it has to rest on more than a handful of wells, and
# too few references must leave the duplex ungated rather than guess.
CO_MIGRATION_SEARCH = 50     # scans either side to search when a duplex has no gate
CO_MIGRATION_MAX_HALF_SPAN = 40.0  # wider than this and the references are noise


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


def _sample_quadrature(trace, samp_col, sigma, n):
    """Detect the four-peak sample pattern directly, with no standard.

    Used only when a well has no internal standard.  hom1/hom2 are not
    scorable without it -- there is nothing to say which allele the observed
    species is -- but a four-peak pattern is a heterozygote on its own, so the
    no-IS case is not automatically a no-call.
    """
    from scipy.signal import find_peaks, savgol_filter
    y = np.asarray(trace, dtype=float)
    ys = savgol_filter(y, 9, 2, mode="interp")
    base = float(np.median(ys[:400]))
    rng = float(ys.max() - base)
    if rng <= 0 or sigma <= 0:
        return []
    # Primary peaks only: this assay carries a Taq-A satellite ~20 scans behind
    # every real peak, so counting all of them double-counts each species.
    p, _ = find_peaks(ys, height=base + 0.35 * rng,
                      prominence=max(2.0 * sigma, 0.01 * rng))
    return [int(x) for x in p if x >= 400]


def _no_is_het_chance(doc, row, acgt, is_col, samp_col, base_order,
                      is_channel, sample_channel, geometry=None, center=None,
                      center_tol=None):
    """A well with no internal standard: het only, never hom1/hom2.

    Without the standard there is nothing to anchor migration, so the observed
    species cannot be identified as allele 1 or 2 -- which is why a missing IS
    used to be an automatic no-call that threw away wells a four-peak pattern
    can still call.  Returns *row*.

    With a fitted *geometry* the four-peak test is the assay's own template
    rather than a bare peak count: sample and standard are the same fragment,
    so a het must show the same d1/d2/d3 spacing the standard does.  Verified
    on ABCC2 well B10, whose sample bands at 2419/2491/2728/2739 measure
    d1=72, d2=237, d3=11 against the standard's 76/257/12.
    """
    y = np.asarray(acgt[:, samp_col], dtype=float)
    sigma = _noise_sigma(y)
    hits = []
    if geometry is not None:
        from itertools import combinations
        pk = _is_candidates(y, DEFAULT_IS_CUT, geometry=geometry)
        if center is not None and center_tol:
            span = (geometry.d1 + geometry.d2 + geometry.d3) / 2.0
            lo, hi = center - center_tol - span, center + center_tol + span
            pk = [c for c in pk if lo <= c[0] <= hi]
        for combo in combinations(range(len(pk)), 4):
            a, b, c, d = (pk[i] for i in combo)
            if geometry.matches(a[0], b[0], c[0], d[0]):
                h = [a[1], b[1], c[1], d[1]]
                if min(h) > 0 and min(h) / max(h) >= IS_EQUIMOLAR_MIN:
                    hits.append((a[0], b[0], c[0], d[0]))
        if hits:
            best = min(hits, key=lambda q: abs(sum(q) / 4.0 - (center or q[0])))
            row["sample_peaks"] = "/".join(str(x) for x in best)
            row["call"] = "het"
            row["flags"] = ",".join(filter(None, [row["flags"], "no-is"]))
            row["reason"] = (f"no IS on Ch{is_channel}; het from a four-peak "
                             "pattern matching the fragment geometry")
            for key, val in zip(("hom1", "hom2", "het1", "het2"), best):
                row[key] = float(y[val])
            return row
    peaks = _sample_quadrature(y, samp_col, sigma, y.size)
    row["sample_peaks"] = "/".join(str(x) for x in peaks)
    if len(peaks) >= 4:
        row["call"] = "het"
        row["flags"] = ",".join(filter(None, [row["flags"], "no-is"]))
        row["reason"] = (f"no IS on Ch{is_channel}; het from {len(peaks)}-peak "
                         "sample pattern")
        if peaks:
            row["het1"] = float(y[peaks[0]])
            row["het2"] = float(y[peaks[1]])
        if len(peaks) > 2:
            row["hom1"] = float(y[peaks[2]])
        if len(peaks) > 3:
            row["hom2"] = float(y[peaks[3]])
    return row


def auto_genotype(doc, is_channel=DEFAULT_IS_CHANNEL,
                  sample_channel=DEFAULT_SAMPLE_CHANNEL,
                  base_order="ACTG", cut=DEFAULT_IS_CUT, run_name="",
                  std_scans_manual=None, is_model=None):
    """Genotype one well without any clicking -> a result row.

    *std_scans_manual* is the operator's own four internal-standard picks for
    this well.  When given it wins over ``find_is_quartet``: that detector
    accepts any equimolar set, and on a weak-IS plate it will happily lock
    onto four equal-magnitude noise ripples in preference to the real standard,
    which on ABCC2 N10 F12 meant reading ~680-sigma noise as the standard while
    the genuine 3000-5700-sigma peaks went unused.

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
        "std_scans": "", "std_scans_manual": "", "std_source": "",
        "het_resolved": True, "std_snr": 0.0, "reason": "",
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

    manual = [int(x) for x in (std_scans_manual or [])]
    if manual:
        # Operator picks are authoritative. Validate the shape, not the values:
        # four ascending, distinct scans is all the window logic needs.
        manual = sorted(set(manual))
        if len(manual) < 2:
            row["reason"] = ("manual standard picks need at least two scans; "
                             "falling back to detection")
            row["flags"] = "manual-std-rejected"
            manual = []
    if manual:
        scans = manual
        row["std_scans_manual"] = "/".join(str(x) for x in scans)
        row["std_source"] = "manual"
    else:
        # Fallback 1: the plate model fitted to the operator's own marks.
        hit = (is_model.find(acgt[:, is_col], cut=cut)
               if is_model is not None else None)
        if hit is not None:
            scans, _heights, how = hit
            row["std_source"] = how
        else:
            # Fallback 2: the original global T9 search -- but only on a plate
            # nobody marked.  Once a fitted model exists its geometry is known
            # to disagree with the T9 ratios, and falling back anyway is what
            # produced every remaining false quartet on ABCC2 (13 hom-1 wells
            # read as hom-2 off invented peaks near scan 2000).
            found = (None if is_model is not None
                     else find_is_quartet(acgt[:, is_col], cut=cut))
            if found is None:
                # No IS. hom1/hom2 are not scorable without one, but a four-peak
                # sample pattern on its own is a heterozygote.
                row["std_source"] = "none"
                row["reason"] = f"no standard quartet on Ch{is_channel}"
                return _no_is_het_chance(
                    doc, row, acgt, is_col, samp_col, base_order, is_channel,
                    sample_channel,
                    geometry=(is_model.geometry if is_model is not None else None),
                    center=(is_model.center if is_model is not None else None),
                    center_tol=(is_model.window() if is_model is not None else None))
            scans, _heights = found
            row["std_source"] = "detected"
    row["std_scans"] = "/".join(str(x) for x in scans)
    # A standard showing three bands carries a merged heteroduplex, which is a
    # complete and equally valid standard.  Only a standard too short to fix
    # even H1 + H2 + a heteroduplex region is unusable.
    if len(scans) < 3:
        row["reason"] = "standard too short to measure"
        return row

    # Geometry alone says where the bands should be, not that they are there.
    # The plate's own marks say whether a merged standard is routine or off-model
    # here: on the four plates that never merge it is not, and a three-band read
    # is a standard that failed to form rather than one that merged.  Those are
    # not genotypeable, because the sample windows read off a half-formed
    # standard are misplaced.  Manual marks bypass this -- the operator marking
    # the well is the explicit statement that the standard is usable there.
    merged_het = len(scans) < 4
    row["het_resolved"] = not merged_het
    if merged_het and is_model is not None and not is_model.merged_seen \
            and not manual:
        row["reason"] = ("standard shows an unresolved heteroduplex on a plate "
                         "whose marked wells always resolve it")
        return row

    y = np.asarray(acgt[:, samp_col], dtype=float)
    n = y.size
    sigma = _noise_sigma(y)
    idx = np.arange(n)
    segs = _quartet_segments(scans)
    # Three segments means the heteroduplexes merged; measure the single band
    # and leave the fourth slot empty rather than inventing a second one.  The
    # call then rests on that band alone, which is conclusive for the genotype
    # but reports no separate HET2 area.
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
        # Look for this duplex's apex only where its product belongs.  The
        # segment itself is far wider than the product (H2 spans ~170 scans
        # when d1 is 81), so taking the segment maximum hands the measurement to
        # whatever else happens to sit inside it -- which is how ABCC2 H01 grew
        # a second allele it does not have, and how F12's off-target became its
        # only "product".  With a calibrated gate the peak is sought inside the
        # co-migration band; without one the old fixed radius stands.
        gate = is_model.co_migration_window(k) if is_model is not None else None
        if gate is not None:
            a = max(0, int(round(scans[k] + gate[0])))
            b = min(n, int(round(scans[k] + gate[1])) + 1)
        else:
            a = max(0, int(scans[k]) - SEGMENT_APEX_RADIUS)
            b = min(n, int(scans[k]) + SEGMENT_APEX_RADIUS + 1)
        apex = float(y[a:b].max()) if b > a else 0.0
        areas.append(max(0.0, float(np.trapz(y[lo:hi] - base, dx=1.0))))
        snrs.append((apex - base) / sigma if sigma > 0 else 0.0)

    while len(areas) < 4:
        areas.append(0.0)
        snrs.append(0.0)
    for i, name in enumerate(("hom1", "hom2", "het1", "het2")):
        row[name] = round(areas[i], 1)
        row["snr%d" % (i + 1)] = round(snrs[i], 1)

    call, frac, flags = scorer.t9_call(
        areas[0], areas[1], areas[2], None if merged_het else areas[3], snrs)
    row["call"] = call
    row["frac"] = round(frac, 4)
    row["flags"] = ",".join(sorted(flags))
    if merged_het:
        # Legitimate on a plate that runs merged, and never score-affecting --
        # but it is exactly the well where an allelic-imbalance read cannot be
        # checked, because the two heteroduplex areas that would show the
        # imbalance are not separable here.
        row["flags"] = ",".join(sorted(flags | {"het-merged"}))

    # Weak-standard flag, relative to this plate.  The call stands: a weak
    # pattern still reads in context, and ABCC2 D07 is a hom-1 the operator
    # would call from the wells above and below it.  Deliberately applied to
    # manual marks too -- where the operator's positions are trusted is a
    # separate question from how strong the standard turned out to be.
    if is_model is not None:
        thr = is_model.weak_std_threshold()
        if thr is not None:
            weakest = weakest_std_band_snr(acgt, is_col, scans)
            row["std_snr"] = round(weakest, 1)
            if weakest < thr:
                row["flags"] = ",".join(sorted(
                    set(f for f in row["flags"].split(",") if f)
                    | {"std-weak"}))
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
    def pick(self, scan, vol=None, only_col=None):
        """Pick the peak you clicked on: nearest scan wins, then nearest
        voltage (tells channels apart when two share a scan), then tallest.

        *only_col* restricts the pick to one acgt column.  The channels are
        overlaid in a single axes, so a click carries no reliable channel
        information -- with the internal standard appearing on two channels a
        few scans apart, nearest-apex silently moved picks onto the wrong
        trace (it is what wrote ``channel=4`` for peaks the operator believed
        they had clicked on Ch3).  Locking the column removes the guess.

        A peak whose area is already picked on the same channel is never
        picked again — undo it first if you need to (neighbouring peaks, e.g.
        the two alleles of a heterozygote, stay pickable).  Returns the new
        record (plus any +A record appended) or None; on a refused
        re-pick, ``_reject`` is set to ``"area"`` for the UI message."""
        self._reject = None
        radius = max(CLICK_RADIUS, int(self._spacing() * 2.0))
        cands = []
        for col, color in self.col_color.items():
            if only_col is not None and col != only_col:
                continue
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

    # ------------------------------------- internal-standard CTC-CE duplexes
    def _duplex_names(self, n):
        """CTC-CE duplex names for `n` picked peaks, in migration order.

        In migration order the homoduplexes come first — the two alleles of a
        heterozygote one rs SNP base apart, or a single one for a homozygote —
        then the heteroduplexes the re-annealed strands form.  So the number of
        homoduplexes tells the two apart: 4 peaks = 2+2 (the internal standard,
        and a heterozygote whose heteroduplexes are baseline-resolved), 3 peaks
        = 1+2 (a homozygote, whose mutant strands have all re-annealed), 2
        peaks = 2+0.  Naming a 3-peak position as a truncated standard would
        call one of its heteroduplexes a homoduplex and halve a low MF.
        """
        n = min(n, 4)
        n_hom = 1 if n == 3 else min(n, 2)
        return (["HOM1", "HOM2"][:n_hom] + ["HET1", "HET2"][:n - n_hom])

    def _standard_mains(self):
        """The picked mains of the last-picked channel, in scan order.

        The internal standard's four duplexes are one fragment separated by
        cycling temperature, so they span many repeats of the migration axis —
        the whole channel is the set, and that is the point of tagging them.
        """
        mains = [r for r in self.records if r["kind"] == "main"]
        if not mains:
            raise ValueError("Pick the standard main peaks first.")
        last_col = mains[-1]["col"]
        col_mains = sorted((m for m in mains if m["col"] == last_col),
                           key=lambda m: m["scan"])
        if len(col_mains) < 2:
            raise ValueError("Need ≥ 2 standard main peaks on the same channel.")
        return col_mains[:4]

    def _position_mains(self):
        """The picked mains of ONE allelic position on the last-picked channel.

        A sample's duplex species sit within a few repeats of each other, so
        the set is the last-picked main's own cluster, grown outward one repeat
        at a time.  A separate position tens of repeats away stops the growth,
        which keeps two heterozygous positions in one well apart.
        """
        mains = [r for r in self.records if r["kind"] == "main"]
        if not mains:
            raise ValueError("Pick the main peaks first.")
        last = mains[-1]
        col_mains = sorted((m for m in mains if m["col"] == last["col"]),
                           key=lambda m: m["scan"])
        if len(col_mains) < 2:
            raise ValueError("Need ≥ 2 main peaks on the same channel — a "
                             "single peak cannot be told apart from a missed "
                             "pick.")
        win = self._het_window()
        i = col_mains.index(last)
        lo = hi = i
        while lo > 0 and col_mains[lo]["scan"] - col_mains[lo - 1]["scan"] <= win:
            lo -= 1
        while hi < len(col_mains) - 1 and \
                col_mains[hi + 1]["scan"] - col_mains[hi]["scan"] <= win:
            hi += 1
        pos = col_mains[lo:hi + 1]
        if len(pos) < 2:
            raise ValueError("Need ≥ 2 main peaks within one repeat of each "
                             "other to read a duplex position.")
        return pos[:4]

    def mark_std(self, length_bp=None):
        """Tag the standard mains (scan order) as the four CTC-CE duplexes:
        HOM1, HOM2 (homoduplexes, one rs SNP base apart), HET1, HET2
        (heteroduplexes, one mismatch base from Watson/Crick re-annealing).
        All four are the same fragment, so a shared length (bp) is optional.
        The standard channel is the one of the last main you picked.

        Raises ValueError with a plain-language reason, or returns a message."""
        col_mains = self._standard_mains()
        names = self._duplex_names(len(col_mains))
        self.std = [(m["scan"], n) for m, n in zip(col_mains, names)]
        self.length_bp = length_bp
        msg = "Standard set: " + ", ".join(f"{n}@{x}" for x, n in self.std)
        if self.length_bp is not None:
            msg += f"  (len {self.length_bp:g} bp)"
        return msg

    def mark_duplex(self):
        """Tag ONE allelic position's picked mains as its duplex species, in
        migration order — the same HOM1/HOM2/HET1/HET2 reading the internal
        standard gets, so a sample position can carry the mass-action MF.

        Raises ValueError with a plain-language reason, or returns a message."""
        pos = self._position_mains()
        names = self._duplex_names(len(pos))
        self.std = [(m["scan"], n) for m, n in zip(pos, names)]
        return "Duplex species: " + ", ".join(f"{n}@{x}" for x, n in self.std)

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

    def labelled_species(self):
        """The picked mains that carry a duplex label, in migration order.
        Tagging a well's position puts them in one list — two peaks for a
        homozygote, four when the heteroduplexes are baseline-resolved."""
        return sorted((m for m in self.records
                       if m["kind"] == "main" and self.duplex_of(m)),
                      key=lambda m: m["scan"])

    def mass_action(self, rec):
        """PCR mass-action mutant fraction of the duplex position `rec` is in:

            MF = (A_MUT + ½ × A_HET) / (A_WT + A_MUT + A_HET)

        A_WT and A_MUT are the homoduplex areas — with both alleles present the
        larger is taken as wild type, the right reading for a rare mutation,
        and a true heterozygote is symmetric either way; with one homoduplex
        (a homozygous WT allele) it is all wild type, and the mutant fraction
        lives entirely in the heteroduplex.  A_HET is the combined area of the
        labelled heteroduplex peaks.

        The half-heteroduplex term is the point of the formula: it makes a clean
        heterozygote read 0.5 instead of the 0.25 a plain area ratio of the two
        homoduplexes gives, and it is what carries a low mutant fraction, where
        the mutant strands are essentially all in heteroduplex and no mutant
        homoduplex is visible at all.

        Also returns the allelic imbalance of the two homoduplexes,
        AI = A_HOMO1 / (A_HOMO1 + A_HOMO2), which needs no wild-type choice
        and is None unless both are present.

        Returns None unless the position carries duplex labels (see
        mark_duplex / mark_std), so an unlabelled two-peak position keeps
        reporting the plain small/(small+large) fraction."""
        if rec is None or rec.get("kind") != "main":
            return None
        species = self.labelled_species()
        if not any(m is rec or m.get("gid") == rec.get("gid") for m in species):
            return None
        hom, het = {}, 0.0
        for m in species:
            name = self.duplex_of(m)
            if name in ("HOM1", "HOM2"):
                hom[name] = float(m["area"])
            elif name in ("HET1", "HET2"):
                het += float(m["area"])
        if not hom:
            return None
        # A single homoduplex is a homozygous WT allele: its mutant strands are
        # all in heteroduplex, which is exactly the low-MF case where no mutant
        # homoduplex is visible at all.
        a_wt = float(hom.get("HOM1", 0.0))
        a_mut = float(hom.get("HOM2", 0.0))
        if a_wt < a_mut:
            a_wt, a_mut = a_mut, a_wt
        total = a_wt + a_mut + het
        if total <= 0:
            return None
        h1, h2 = hom.get("HOM1"), hom.get("HOM2")
        ai = (h1 / (h1 + h2)
              if (h1 is not None and h2 is not None and (h1 + h2) > 0) else None)
        return {"mf": (a_mut + 0.5 * het) / total,
                "a_wt": a_wt, "a_mut": a_mut, "a_het": het, "ai": ai,
                "n_homoduplex": len(hom)}

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
            mf = self.mass_action(r) if r["kind"] == "main" else None
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
                "mf": round(mf["mf"], 4) if mf else "",
                "ai": round(mf["ai"], 4) if mf and mf["ai"] is not None else "",
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
    """Write rows keyed by column name, not by position.

    Rows here do not all carry the same keys -- a load-error row has no
    measured areas, and a no-IS row has no standard scans -- so writing
    ``list(r.values())`` against ``rows[0]`` headers silently shifts every
    value after the first gap into the wrong column.
    """
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        if not rows:
            return
        headers = []
        for r in rows:
            for k in r:
                if k not in headers:
                    headers.append(k)
        w.writerow(headers)
        for r in rows:
            w.writerow([r.get(h, "") for h in headers])


def _write_xlsx(path, rows):
    from openpyxl import Workbook
    from openpyxl.utils import get_column_letter
    wb = Workbook()
    ws = wb.active
    ws.title = "peaks"
    if rows:
        headers = []
        for r in rows:
            for k in r:
                if k not in headers:
                    headers.append(k)
        ws.append(headers)
        for r in rows:
            ws.append([r.get(h, "") for h in headers])
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
