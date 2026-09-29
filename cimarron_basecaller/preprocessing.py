"""
Preprocessing stage: begin/end trimming, baseline subtraction, and
spectral/leakage separation, following the general approach described
in EP0944739A1 section IV.A ("Preprocessing").

Input convention: a trace is a (n_samples, 4) float array, one column
per fluorescent channel / lane (A, C, G, T order is up to the caller;
this module is channel-order agnostic).
"""

from __future__ import annotations
import numpy as np


def detect_begin_end(trace: np.ndarray, n_zones: int = 8) -> tuple[int, int]:
    """Find usable Begin/End sample indices for a trace.

    Per the patent: split the trace into zones, find each zone's max
    amplitude to locate the putative primer peak (first large-amplitude
    region), set Begin just before it, and set End near the trace's
    natural end (trimmed back to exclude trailing noise).

    This is a reasonable-fidelity reconstruction of the described method;
    the patent's exact zone geometry is not fully specified, so the zone
    count is exposed as a parameter for tuning against real data.
    """
    n = trace.shape[0]
    envelope = trace.max(axis=1)

    zone_bounds = np.linspace(0, n, n_zones + 1, dtype=int)
    zone_maxima = [
        envelope[zone_bounds[i]:zone_bounds[i + 1]].max()
        if zone_bounds[i + 1] > zone_bounds[i] else 0.0
        for i in range(n_zones)
    ]

    # First zone whose max amplitude is a large fraction of the global max
    # is treated as containing the primer peak.
    global_max = max(zone_maxima) if zone_maxima else 0.0
    threshold = 0.5 * global_max
    begin_zone = next(
        (i for i, m in enumerate(zone_maxima) if m >= threshold), 0
    )
    begin = zone_bounds[begin_zone]

    # Refine: begin point = first sample below the mean of the first half
    # of the signal (as described).
    half = max(1, n // 2)
    first_half_mean = envelope[:half].mean()
    below = np.where(envelope[begin:half] < first_half_mean)[0]
    if len(below):
        begin = begin + int(below[0])

    # End point: back off 350 samples from the very end, per the patent's
    # stated default -- clamp for short traces.
    end = max(begin + 1, n - 350)

    return int(begin), int(end)


def rising_exponential_baseline(
    channel: np.ndarray, rise_span: int = 75, rise_fraction: float = 1.0 / 3.0
) -> np.ndarray:
    """One-directional baseline estimate using a rising-exponential threshold.

    Follows the described algorithm: start the threshold at the minimum of
    the first 10 samples; walk forward, ratcheting the threshold up slowly;
    whenever a sample falls below the current threshold, anchor a new
    baseline segment there (piecewise-linear between anchors) and reset the
    threshold. If 100 samples pass with no new anchor, force a checkpoint
    and speed up the ramp.
    """
    n = len(channel)
    baseline = np.zeros(n)
    anchor_idx = 0
    anchor_val = float(channel[:10].min()) if n >= 10 else float(channel.min())
    threshold = anchor_val
    samples_since_anchor = 0
    rise_rate = rise_fraction / rise_span

    for i in range(n):
        if channel[i] <= threshold:
            # New anchor point: fill baseline linearly since last anchor.
            if i > anchor_idx:
                baseline[anchor_idx:i + 1] = np.linspace(
                    anchor_val, channel[i], i - anchor_idx + 1
                )
            anchor_idx = i
            anchor_val = float(channel[i])
            threshold = anchor_val
            samples_since_anchor = 0
        else:
            samples_since_anchor += 1
            threshold += rise_rate * max(channel[i] - anchor_val, 0.0)
            if samples_since_anchor >= 100:
                # Force a checkpoint and steepen the ramp, as described.
                baseline[anchor_idx:i + 1] = np.linspace(
                    anchor_val, channel[i], i - anchor_idx + 1
                )
                anchor_idx = i
                anchor_val = float(channel[i])
                threshold = anchor_val
                samples_since_anchor = 0
                rise_rate *= 2.0

    if anchor_idx < n - 1:
        baseline[anchor_idx:] = anchor_val
    return baseline


def subtract_baseline(trace: np.ndarray) -> np.ndarray:
    """Baseline-subtract every channel using the geometric mean of a
    left-to-right and a right-to-left rising-exponential baseline pass,
    as described in the patent (Fig. 22)."""
    out = np.empty_like(trace, dtype=float)
    for c in range(trace.shape[1]):
        col = trace[:, c].astype(float)
        fwd = rising_exponential_baseline(col)
        bwd = rising_exponential_baseline(col[::-1])[::-1]
        # Geometric mean of two non-negative baseline estimates.
        combined = np.sqrt(np.clip(fwd, 0, None) * np.clip(bwd, 0, None))
        out[:, c] = col - combined
    return out


def build_characteristic_matrix(
    calibration_traces: np.ndarray, peak_indices: Sequence_ = None
) -> np.ndarray:
    """Build the characteristic matrix (CHM) capturing spectral cross-talk
    between channels, from calibration data where each channel's peaks are
    known/dominant. Columns are normalized so the largest element is 1.0,
    as described.

    calibration_traces: (n_samples, 4) array from a run/region where each
    channel has clean, well-separated peaks.
    """
    import numpy as _np
    n_channels = calibration_traces.shape[1]
    chm = _np.zeros((n_channels, n_channels))
    for lane in range(n_channels):
        col = calibration_traces[:, lane]
        peak_idx = _np.argmax(col)
        ratios = calibration_traces[peak_idx, :] / (col[peak_idx] + 1e-12)
        chm[:, lane] = ratios
    # Normalize each column so its largest element is 1.0
    for lane in range(n_channels):
        col_max = chm[:, lane].max()
        if col_max > 0:
            chm[:, lane] /= col_max
    return chm


# typing helper (avoids importing Sequence at module top redundantly)
from typing import Sequence as Sequence_  # noqa: E402


def spectral_separate(trace: np.ndarray, chm: np.ndarray) -> np.ndarray:
    """Apply spectral/leakage separation given a characteristic matrix.

    Per the patent: SST = inv(CHM * CHM^T)^-1 ... practically this is
    solved as a linear unmixing problem: observed = CHM @ true, so
    true = inv(CHM) @ observed (per-sample).
    """
    chm_inv = np.linalg.pinv(chm)
    return trace @ chm_inv.T


def cubic_spline_upsample(trace: np.ndarray, factor: int = 2) -> np.ndarray:
    """Increase sample density via cubic spline interpolation, used when a
    trace has fewer than ~8 scanlines per band (Section IV.A)."""
    from scipy.interpolate import CubicSpline

    n = trace.shape[0]
    x = np.arange(n)
    x_new = np.linspace(0, n - 1, n * factor)
    out = np.empty((len(x_new), trace.shape[1]))
    for c in range(trace.shape[1]):
        cs = CubicSpline(x, trace[:, c])
        out[:, c] = cs(x_new)
    return out
