"""
BaseQual fuzzy-logic quality scoring (EP0944739A1, Figs. 15-21).

Assigns each called band a quality value on [0, 100] from six measured
attributes -- height, cross-banding, width, shape, small (upstream) gap,
and large (downstream) gap -- combined through a set of nine graded
rules whose output "conclusion" fuzzy sets are centered at 0, 13, 25, 38,
50, 63, 75, 88 and 100.

The band-height breakpoints and the rule structure below reproduce the
*numeric parameters* disclosed in the patent (thresholds/breakpoints are
functional facts, not creative expression); the control flow and code
itself is an independent implementation, not a transcription of the
original C++.
"""

from __future__ import annotations
import numpy as np

from .fuzzy import FuzzySet, AND, OR, NOT
from .peaks import Band

# Height membership sets (as disclosed): tiny / small / moderate / normal / tall
_TINY_H = FuzzySet(np.array([0.00, 0.02, 0.03]), np.array([1.00, 1.00, 0.00]), "tinyH")
_SMALL_H = FuzzySet(np.array([0.015, 0.025, 0.04, 0.05]), np.array([0.0, 1.0, 1.0, 0.0]), "smallH")
_MODERATE_H = FuzzySet(np.array([0.035, 0.045, 0.14, 0.16]), np.array([0.0, 1.0, 1.0, 0.0]), "moderateH")
_NORMAL_H = FuzzySet(np.array([0.135, 0.145, 0.48, 0.55]), np.array([0.0, 1.0, 1.0, 0.0]), "normalH")
_TALL_H = FuzzySet(np.array([0.475, 0.485, 1.00]), np.array([0.0, 1.0, 1.0]), "tallH")

# Cross-banding: present (bad) vs. negligible
_PRESENT_XB = FuzzySet(np.array([1.0, 1.35]), np.array([1.0, 0.0]), "presentXb")
_NEGLIGIBLE_XB = FuzzySet(np.array([1.35, 1.5]), np.array([0.0, 1.0]), "negligibleXb")

# Normalized width: acceptable band around 0
_NORMAL_W = FuzzySet(np.array([-1.0, -0.4, 0.4, 1.0]), np.array([0.0, 1.0, 1.0, 0.0]), "normalW")

# Shape (linear correlation coefficient vs. an ideal band): good if high
_GOOD_SHAPE = FuzzySet(np.array([0.4, 0.5, 1.0]), np.array([0.0, 1.0, 1.0]), "goodShape")

# Gaps (normalized): good/OK near 0.3-0.7 in the disclosed example
_GOOD_GAP = FuzzySet(np.array([0.3, 0.5, 0.7]), np.array([1.0, 0.0, 1.0]), "goodGap")

# Baseline buzz: OK below ~0.2-0.4
_OK_BUZZ = FuzzySet(np.array([0.0, 0.2, 0.4]), np.array([1.0, 1.0, 0.0]), "okBuzz")

# Conclusion (output) fuzzy sets, nine graded quality levels, centroids
# at 0, 13, 25, 38, 50, 63, 75, 88, 100 as disclosed.
_CONCLUSION_CENTERS = [0, 13, 25, 38, 50, 63, 75, 88, 100]
_CONCLUSION_HALF_WIDTH = 12.5


def _conclusion_set(center: float) -> FuzzySet:
    lo, hi = center - _CONCLUSION_HALF_WIDTH, center + _CONCLUSION_HALF_WIDTH
    return FuzzySet(np.array([lo, center, hi]), np.array([0.0, 1.0, 0.0]), name=f"C{center}")


def score_band(
    height: float,
    cross_banding: float,
    norm_width: float,
    shape: float,
    small_gap: float,
    large_gap: float,
    buzz: float,
) -> float:
    """Compute the BaseQual quality value (0-100) for one band from its
    six measured attributes, using the nine-rule scheme disclosed in the
    patent."""
    tiny_h = _TINY_H.membership(height)
    small_h = _SMALL_H.membership(height)
    moderate_h = _MODERATE_H.membership(height)
    normal_h = _NORMAL_H.membership(height)
    tall_h = _TALL_H.membership(height)
    ok_height = OR(moderate_h, normal_h, tall_h)

    present_xb = _PRESENT_XB.membership(cross_banding)
    negligible_xb = _NEGLIGIBLE_XB.membership(cross_banding)

    normal_w = _NORMAL_W.membership(norm_width)
    good_shape = _GOOD_SHAPE.membership(shape)
    good_gap_small = _GOOD_GAP.membership(small_gap)
    good_gap_large = _GOOD_GAP.membership(large_gap)
    ok_gap = AND(good_gap_small, good_gap_large)
    ok_buzz = _OK_BUZZ.membership(buzz)

    # Count how many of the four "quality" measures (buzz, width, shape,
    # gap) are out of tolerance, per the disclosed 1Bad..4Bad combinatorics.
    measures_ok = [ok_buzz, normal_w, good_shape, ok_gap]
    n_bad = sum(1 for m in measures_ok if m < 0.5)

    run_channel = AND(NOT(good_shape), ok_height)
    run_fill = AND(run_channel, AND(negligible_xb, AND(normal_w, AND(ok_buzz, ok_gap))))
    four_ok = AND(ok_buzz, AND(normal_w, AND(good_shape, ok_gap)))

    rules = {
        0: AND(tiny_h, NOT(ok_gap)),
        13: AND(small_h, AND(present_xb, NOT(ok_gap))),
        25: AND(tiny_h, ok_gap),
        38: AND(OR(small_h, moderate_h), OR(present_xb, 1.0 if n_bad >= 3 else 0.0)),
        50: AND(OR(small_h, moderate_h), OR(NOT(negligible_xb), 1.0 if n_bad >= 2 else 0.0)),
        63: AND(ok_height, OR(OR(NOT(ok_buzz), NOT(negligible_xb)), 1.0 if n_bad >= 2 else 0.0)),
        75: OR(
            AND(ok_height, AND(ok_buzz, AND(negligible_xb, 1.0 if n_bad >= 2 else 0.0))),
            OR(run_fill, AND(small_h, AND(NOT(present_xb), four_ok))),
        ),
        88: AND(ok_buzz, AND(ok_height, AND(negligible_xb, 1.0 if n_bad == 1 else 0.0))),
        100: AND(ok_height, AND(negligible_xb, four_ok)),
    }

    conclusion_xs = np.linspace(-15, 115, 261)
    conclusion_ys = np.zeros_like(conclusion_xs)
    for center, strength in rules.items():
        cset = _conclusion_set(center)
        membership = np.array([cset.membership(x) for x in conclusion_xs])
        conclusion_ys = np.maximum(conclusion_ys, np.minimum(membership, strength))

    combined = FuzzySet(conclusion_xs, conclusion_ys)
    quality = combined.centroid()
    return float(np.clip(quality, 0.0, 100.0))


def score_bands(bands: list[Band]) -> None:
    """Score every band in place, filling in `.quality`. Bands must already
    have height/cross_banding/width populated (see peaks.py); gap/shape/
    buzz are derived here from neighboring bands."""
    from .peaks import expected_width_curve

    if not bands:
        return
    expected_width = expected_width_curve(bands)
    heights = np.array([b.height for b in bands])
    max_h = heights.max() if heights.max() > 0 else 1.0
    norm_heights = heights / max_h  # normalize into the [0,1] range the
    # disclosed breakpoints assume (they're expressed as fractions of a
    # normalized band amplitude).

    for i, b in enumerate(bands):
        norm_width = (b.width / expected_width[i] - 1.0) if expected_width[i] > 0 else 0.0
        # Shape: correlation between this band's local height profile and
        # an idealized triangular/gaussian bump of the expected width.
        shape = _estimate_shape(bands, i)
        small_gap = b.spacing_left / (b.spacing_left + b.spacing_right + 1e-9)
        large_gap = b.spacing_right / (b.spacing_left + b.spacing_right + 1e-9)
        buzz = _estimate_buzz(bands, i)

        b.quality = score_band(
            height=float(norm_heights[i]),
            cross_banding=b.cross_banding,
            norm_width=norm_width,
            shape=shape,
            small_gap=small_gap,
            large_gap=large_gap,
            buzz=buzz,
        )


def _estimate_shape(bands: list[Band], i: int) -> float:
    """Cheap shape proxy: how consistent this band's height is with its
    immediate neighbors (a real implementation would fit the underlying
    trace samples, not just peak heights -- left as a refinement point)."""
    if len(bands) < 3:
        return 1.0
    lo, hi = max(0, i - 1), min(len(bands), i + 2)
    local = [b.height for b in bands[lo:hi]]
    if len(local) < 2:
        return 1.0
    spread = np.std(local) / (np.mean(local) + 1e-9)
    return float(np.clip(1.0 - spread, 0.0, 1.0))


def _estimate_buzz(bands: list[Band], i: int) -> float:
    """Cheap baseline-buzz proxy from local cross-banding variability."""
    if len(bands) < 3:
        return 0.0
    lo, hi = max(0, i - 1), min(len(bands), i + 2)
    local_xb = [1.0 / (b.cross_banding + 1e-9) for b in bands[lo:hi]]
    return float(np.clip(np.mean(local_xb), 0.0, 1.0))
