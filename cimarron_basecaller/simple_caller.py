"""
Shared signal helpers for the spacing tracker (``spacing_caller.track_bases``):
rolling-minimum baseline subtraction, rolling local normalisation and
per-channel mobility correction.  The module once also held a simple
per-channel peak-picking caller; it scored far below track_bases and was
removed (see git history).
"""

from __future__ import annotations
import bisect
import numpy as np
from scipy.ndimage import minimum_filter1d, uniform_filter1d


def robust_baseline_subtract(trace: np.ndarray, window: int = 151) -> np.ndarray:
    """Rolling-minimum baseline subtraction per channel: a simple, robust
    alternative to the patent's rising-exponential-threshold algorithm.
    A rolling minimum tracks the noise floor; smoothing it avoids chasing
    individual low-noise samples."""
    out = np.empty_like(trace)
    for c in range(trace.shape[1]):
        col = trace[:, c]
        baseline = minimum_filter1d(col, size=window, mode="nearest")
        baseline = uniform_filter1d(baseline, size=window // 3 + 1, mode="nearest")
        out[:, c] = np.clip(col - baseline, 0, None)
    return out


def rolling_percentile(x: np.ndarray, size: int, percentile: float) -> np.ndarray:
    """``scipy.ndimage.percentile_filter(x, percentile, size, mode="nearest")``
    for a 1-D signal, computed with a sorted sliding window.

    Same rank rule, window alignment and edge padding as ndimage, so the
    result is identical value for value; the cost is O(n log size) instead of
    ndimage's O(n * size), which made this the bulk of a track_bases call at
    the 1800-scan windows the plate configs use (~150 ms -> ~11 ms per channel).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n == 0:
        return x.copy()
    size = max(1, int(size))
    if not np.isfinite(x).all():
        # NaN does not sort; leave that case to ndimage.
        from scipy.ndimage import percentile_filter
        return percentile_filter(x, percentile=percentile, size=size, mode="nearest")
    rank = int(float(size) * percentile / 100.0)
    if rank == size:
        rank -= 1
    left = size // 2
    xp = np.concatenate([np.full(left, x[0]), x,
                         np.full(size - left - 1, x[-1])]).tolist()
    win = sorted(xp[:size])
    out = np.empty(n)
    out[0] = win[rank]
    for i in range(1, n):
        del win[bisect.bisect_left(win, xp[i - 1])]
        bisect.insort(win, xp[i - 1 + size])
        out[i] = win[rank]
    return out


def normalize_channels_local(trace: np.ndarray, window: int = 400, percentile: float = 95.0) -> np.ndarray:
    """Like normalize_channels, but the scale factor is a ROLLING local
    high-percentile rather than one global value -- compensates for signal
    amplitude naturally rising/falling over the length of a read (e.g.
    weaker signal near the very start/end), which a single global
    normalization factor doesn't account for.
    """
    out = np.empty_like(trace)
    for c in range(trace.shape[1]):
        local_scale = rolling_percentile(trace[:, c], window, percentile)
        local_scale = np.maximum(local_scale, np.percentile(trace[:, c], 90) * 0.1 + 1e-6)
        out[:, c] = trace[:, c] / local_scale
    return out


def apply_mobility_correction(trace: np.ndarray, shifts: list[int] | tuple[int, ...]) -> np.ndarray:
    """Shift each channel's samples by a per-channel integer offset to
    correct for dye mobility differences (each dye migrates through the
    capillary at a slightly different rate, so the same physical base
    shows up at a slightly different scan index per channel). `shifts[c]`
    is the amount channel c's peaks are observed to be delayed by --
    positive means the channel's real peak appears `shifts[c]` samples
    LATER than it should, so we shift it earlier to correct.
    """
    out = np.zeros_like(trace)
    n = trace.shape[0]
    for c, s in enumerate(shifts):
        s = int(round(s))
        if s == 0:
            out[:, c] = trace[:, c]
        elif s > 0:
            out[:n - s, c] = trace[s:, c]
        else:
            out[-s:, c] = trace[:n + s, c]
    return out
