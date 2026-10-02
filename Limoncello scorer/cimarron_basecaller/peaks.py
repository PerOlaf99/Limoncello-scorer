"""
Putative peak detection, and the OmitOkN / GapCheck fuzzy-logic blocks.

Follows the described method (EP0944739A1, Figs. 5-13): peaks are found
on the single envelope (max across the 4 aligned, deconvolved channels)
rather than per-channel, each peak's spacing/height/cross-banding/width
is measured relative to a quadratic fit of its neighbors' values, and two
fuzzy-logic blocks decide (a) whether a putative peak looks like a
spurious insertion (OmitOkN) and (b) whether a gap between two peaks is
large enough that one or more bases were likely missed (GapCheck).

This is an independent implementation built from the patent's prose
description and stated formulas/thresholds; exact membership-function
shapes for OmitOkN/GapCheck are not fully specified in the patent (unlike
BaseQual, whose breakpoints are given explicitly -- see quality.py), so
reasonable trapezoidal shapes consistent with the described behavior are
used here and are natural tuning targets against real trace data.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from .fuzzy import FuzzySet, AND, OR, NOT


@dataclass
class Band:
    """A called/putative band (peak) with its measured attributes."""
    position: int          # sample index of the peak, in the aligned segment
    channel: int            # which of the 4 channels dominates here (base id)
    height: float = 0.0
    spacing_left: float = 0.0    # distance to previous band
    spacing_right: float = 0.0   # distance to next band
    cross_banding: float = 1.0   # dominant / next-dominant amplitude ratio
    width: float = 0.0
    quality: float | None = None
    inserted: bool = False        # True if added by GapCheck splitting


def detect_putative_peaks(aligned_segment: np.ndarray) -> list[Band]:
    """Find local maxima of the 4-channel envelope; a peak is any sample
    taller than both neighbors (liberal definition, per the patent --
    over-detects on purpose, relying on OmitOkN to prune false positives)."""
    envelope = aligned_segment.max(axis=1)
    channel_id = aligned_segment.argmax(axis=1)
    bands = []
    for i in range(1, len(envelope) - 1):
        if envelope[i] > envelope[i - 1] and envelope[i] >= envelope[i + 1]:
            bands.append(Band(position=i, channel=int(channel_id[i]), height=float(envelope[i])))
    return bands


def measure_spacing_and_crossbanding(aligned_segment: np.ndarray, bands: list[Band]) -> None:
    """Fill in spacing_left/right, cross_banding and width for each band,
    in place."""
    envelope = aligned_segment.max(axis=1)
    positions = [b.position for b in bands]
    for i, b in enumerate(bands):
        b.spacing_left = float(positions[i] - positions[i - 1]) if i > 0 else float(positions[i])
        b.spacing_right = float(positions[i + 1] - positions[i]) if i < len(bands) - 1 else b.spacing_left

        # Cross banding: ratio of this channel's amplitude at the peak to
        # the next-largest channel's amplitude at the same sample.
        vals = np.sort(aligned_segment[b.position])[::-1]
        b.cross_banding = float(vals[0] / (vals[1] + 1e-9))

        # Width: samples on either side until amplitude falls to half max.
        half = b.height / 2.0
        left = b.position
        while left > 0 and envelope[left] > half:
            left -= 1
        right = b.position
        while right < len(envelope) - 1 and envelope[right] > half:
            right += 1
        b.width = float(right - left)


def expected_spacing_curve(bands: list[Band]) -> np.ndarray:
    """Quadratic fit of observed spacings vs. band index, used as the
    'expected spacing' at each position along the segment."""
    if len(bands) < 3:
        avg = np.mean([b.spacing_right for b in bands]) if bands else 20.0
        return np.full(len(bands), avg)
    idx = np.arange(len(bands))
    spacings = np.array([b.spacing_right for b in bands])
    coeffs = np.polyfit(idx, spacings, deg=2)
    return np.polyval(coeffs, idx)


def expected_width_curve(bands: list[Band]) -> np.ndarray:
    """Quadratic fit of observed band widths vs. band index."""
    if len(bands) < 3:
        avg = np.mean([b.width for b in bands]) if bands else 5.0
        return np.full(len(bands), avg)
    idx = np.arange(len(bands))
    widths = np.array([b.width for b in bands])
    coeffs = np.polyfit(idx, widths, deg=2)
    return np.polyval(coeffs, idx)


# ---------------------------------------------------------------------
# OmitOkN fuzzy logic block (Figs. 6-10): classify each putative peak as
# OK, AMBIGUOUS, or OMIT (spurious insertion to be discarded).
# ---------------------------------------------------------------------

def _height_membership_sets(med_intersect_pt: float) -> tuple[FuzzySet, FuzzySet]:
    """tinyHt / okHt membership functions, whose breakpoints scale with the
    median band-intersection amplitude, per the patent's stated relation:
    tinyHt breaks at 0.4x and reaches 0 by 1.1x; okHt leaves 0 at 0.5x and
    flattens to 1.0 by 1.5x (x = med_intersect_pt)."""
    m = med_intersect_pt
    tiny = FuzzySet(np.array([0.0, 0.4 * m, 1.1 * m]), np.array([1.0, 1.0, 0.0]), name="tinyHt")
    ok = FuzzySet(np.array([0.5 * m, 1.5 * m, 1.5 * m + 1e-6]), np.array([0.0, 1.0, 1.0]), name="okHt")
    return tiny, ok


def _spacing_membership_sets() -> tuple[FuzzySet, FuzzySet]:
    """okSp / badSp membership functions over normalized spacing
    (observed / expected, so 1.0 = perfect match to an integer multiple
    of expected spacing). "Accept" near integer multiples; "deprecate"
    away from them."""
    ok = FuzzySet(np.array([0.0, 0.85, 1.0, 1.15, 2.0]), np.array([0.0, 0.6, 1.0, 0.6, 0.0]), name="okSp")
    bad = FuzzySet(np.array([0.0, 0.5, 1.0, 1.5, 2.0]), np.array([1.0, 0.5, 0.0, 0.5, 1.0]), name="badSp")
    return ok, bad


def _crossbanding_membership_sets() -> tuple[FuzzySet, FuzzySet]:
    """badXb / negligibleXb over the cross-banding ratio. Ratios >= 1.5 are
    treated as negligible cross-banding (per the patent's stated example);
    ratios approaching 1.0 (two channels of near-equal height) are bad."""
    bad = FuzzySet(np.array([1.0, 1.5, 2.5]), np.array([1.0, 0.25, 0.0]), name="badXb")
    negligible = FuzzySet(np.array([1.0, 1.5, 2.5]), np.array([0.0, 1.0, 1.0]), name="negligibleXb")
    return bad, negligible


def omit_ok_n(bands: list[Band]) -> list[str]:
    """Classify every band as 'OK', 'AMBIGUOUS', or 'OMIT'. Returns a
    label list parallel to `bands`. Does not mutate `bands`."""
    if not bands:
        return []

    intersect_heights = [b.height for b in bands]
    med_intersect = float(np.median(intersect_heights)) if intersect_heights else 1.0
    tinyHt, okHt = _height_membership_sets(med_intersect)
    okSp, badSp = _spacing_membership_sets()
    badXb, negligibleXb = _crossbanding_membership_sets()

    expected_spacing = expected_spacing_curve(bands)
    labels = []

    for b, exp_sp in zip(bands, expected_spacing):
        norm_spacing = b.spacing_left / exp_sp if exp_sp > 0 else 1.0

        ok_sp = okSp.membership(norm_spacing)
        bad_sp = badSp.membership(norm_spacing)
        bad_xb = badXb.membership(b.cross_banding)
        negl_xb = negligibleXb.membership(b.cross_banding)
        tiny_h = tinyHt.membership(b.height)
        ok_h = okHt.membership(b.height)

        ok_conclusion = AND(negl_xb, OR(ok_h, ok_sp))
        ambiguous_conclusion = AND(bad_xb, OR(ok_h, ok_sp))
        omit_conclusion = AND(tiny_h, NOT(ok_sp))

        scores = {"OK": ok_conclusion, "AMBIGUOUS": ambiguous_conclusion, "OMIT": omit_conclusion}
        labels.append(max(scores, key=scores.get))

    return labels


# ---------------------------------------------------------------------
# GapCheck fuzzy logic block (Figs. 11-13): decide whether a gap between
# two neighboring bands should be split (i.e. one or more bases inserted).
# ---------------------------------------------------------------------

def gc_richness(bases: list[int], window: int = 5) -> float:
    """Fraction of G/C calls (channels 1 and 2, by convention) in the last
    `window` upstream bases -- used to temper GapCheck's aggressiveness in
    GC-rich stretches, per the patent."""
    if not bases:
        return 0.0
    recent = bases[-window:]
    gc = sum(1 for c in recent if c in (1, 2))
    return gc / len(recent)


def gap_check(bands: list[Band]) -> list[bool]:
    """Return a parallel list of booleans: True where the gap to the LEFT
    of that band should be split (i.e. a base was likely missed there)."""
    if len(bands) < 2:
        return [False] * len(bands)

    expected_spacing = expected_spacing_curve(bands)
    expected_width = expected_width_curve(bands)

    big_gap = FuzzySet(np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.0, 0.0, 0.7, 1.0]), name="bigGap")
    small_gap = FuzzySet(np.array([-1.0, -0.3, 0.0, 0.3]), np.array([1.0, 1.0, 0.3, 0.0]), name="smallGap")
    big_width = FuzzySet(np.array([-1.0, 0.0, 0.5, 1.0]), np.array([0.0, 0.0, 0.5, 1.0]), name="bigWidth")

    called_channels: list[int] = []
    split_flags = [False] * len(bands)

    for i in range(1, len(bands)):
        b_prev, b_cur = bands[i - 1], bands[i]
        exp_sp = expected_spacing[i]
        exp_wd = expected_width[i]

        norm_gap = (b_cur.spacing_left / exp_sp - 1.0) if exp_sp > 0 else 0.0
        norm_wid_cur = (b_cur.width / exp_wd - 1.0) if exp_wd > 0 else 0.0
        norm_wid_prev = (b_prev.width / exp_wd - 1.0) if exp_wd > 0 else 0.0

        big_gap_m = big_gap.membership(norm_gap)
        small_gap_m = small_gap.membership(norm_gap)
        big_wid_cur = big_width.membership(norm_wid_cur)
        big_wid_prev = big_width.membership(norm_wid_prev)
        gc_rich = gc_richness(called_channels)

        rule_norm = OR(
            AND(big_gap_m, gc_rich),
            AND(small_gap_m, NOT(OR(big_wid_cur, big_wid_prev))),
            NOT(big_gap_m),
        )
        rule_split = OR(
            AND(big_gap_m, OR(big_wid_prev, big_wid_cur)),
            AND(big_gap_m, AND(NOT(small_gap_m), NOT(gc_rich))),
        )

        split_flags[i] = rule_split > rule_norm
        called_channels.append(b_cur.channel)

    return split_flags


def apply_gap_splits(bands: list[Band], split_flags: list[bool]) -> list[Band]:
    """Insert a synthetic band at the centroid of any gap flagged by
    GapCheck, placing it on the 'shoulder' of the poorly-resolved region
    rather than in the trough, per the patent's stated placement rule."""
    if not bands:
        return bands
    out = [bands[0]]
    for i in range(1, len(bands)):
        if split_flags[i]:
            prev, cur = bands[i - 1], bands[i]
            mid_pos = (prev.position + cur.position) // 2
            inserted = Band(
                position=mid_pos,
                channel=prev.channel,  # best guess without re-scanning raw data
                height=min(prev.height, cur.height) * 0.5,
                inserted=True,
            )
            out.append(inserted)
        out.append(bands[i])
    return out
